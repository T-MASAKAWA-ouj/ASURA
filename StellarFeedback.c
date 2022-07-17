#include "config.h"
#include "NeighborSearch.h"
#include "StellarFeedback.h"
#include "SizeDetermination.h"

/*! \file StellarFeedback.c
 * \brief Feedback routines for type II/Ia SNe.
 */


static int StellarFeedbackExportFlagsMaxAllocated = 0;
struct StructActiveSNParticle *ActiveSNParticle; 

static int MaxIterationTimes = 20;
static bool OverMaxIterationTimes = false;
static double LocalKernelMax = 0.e0;

#define USE_STELLARFEEDBACK_RADIUS_LOCAL_UPDATE
#define __CHECK_SUM__
#define __CHECK_WEIGHT__
// static int Test;

static double SNIINumberPerMass;
static double SNIIEnergy;
static double ECSNEnergy;
static double HNEnergy;
static double EnergyConvertToSimulationUnit;
static double InternalEnergyConvertToSimulationUnit;
static double MassConvertToSimulationUnit;
static double MassConvertToMsun;
static double KernelPeak;


static inline double __attribute__((always_inline)) KernelStellarFeedback(const double r, const double InvKerneli){

    double u = r*InvKerneli;
    double coef = M_1_PI*CUBE(InvKerneli);
    if(u<1.e0){
        return (coef*(1.e0 - 1.5*SQ(u) + 0.75*CUBE(u)));
    } else if (u<2.e0){
        return (coef*(0.25*CUBE(2.e0-u)));
    } else {
        return 0.e0;
    }
}

#ifdef USE_CELIB //{
/*!
 * This function initializes the whole functions of CELib.
 */
void InitializeStellarFeedback(void){

    CELibSetRunParameterIMFType(CHEMICALEVOLUTION_IMFTYPE);
    CELibSetRunParameterSNIIRange(CHEMICALEVOLUTION_SNII_UPPERMASS,CHEMICALEVOLUTION_SNII_LOWERMASS);
    CELibSetRunParameterSNIIYieldsTableID(CHEMICALEVOLUTION_SNII_YIELD_TYPE);
    CELibSetRunParameterSNIIHyperNovaFraction(CHEMICALEVOLUTION_SNII_HYPERNOVA_FRACTION);

#if CHEMICALEVOLUTION_SNIa_TYPE == -1
    CELibSetRunParameterSNIaType(1);
#else
    CELibSetRunParameterSNIaType(CHEMICALEVOLUTION_SNIa_TYPE);
#endif

    CELibSetRunParameterSNIaYieldsTableID(CHEMICALEVOLUTION_SNIa_YIELD_TYPE);
    CELibSetRunParameterSNIaYieldsTableModelID(CHEMICALEVOLUTION_SNIa_YIELD_MODEL);

#if CHEMICALEVOLUTION_SNIa_TYPE == 0
    CELibSetRunParameterSNIaRange(CHEMICALEVOLUTION_SNIa_UPPERMASS,CHEMICALEVOLUTION_SNIa_UPPERMASS);
#endif // CHEMICALEVOLUTION_SNIa_TYPE

    CELibSetRunParameterSNIaNassociation(CHEMICALEVOLUTION_SNIa_EVENTNUMBER);

#if CHEMICALEVOLUTION_POPIII_IMF == 0
    CELibSetRunParameterPopIIIIMF(0);
    CELibSetRunParameterPopIIISNe(1);
    CELibSetRunParameterPopIIIAGB(0);
    CELibSetRunParameterPopIIILifeTime(1);
#else
    CELibSetRunParameterPopIIIIMF(1);
    CELibSetRunParameterPopIIISNe(1);
    CELibSetRunParameterPopIIIAGB(1);
    CELibSetRunParameterPopIIILifeTime(1);
#endif


#ifdef USE_CELIB_AGB //{
    CELibSetRunParameterAGBBinNumber(CHEMICALEVOLUTION_AGB_NBIN);
    CELibSetRunParameterAGBBinType(CHEMICALEVOLUTION_AGB_BIN_TYPE);
    CELibSetRunParameterAGBBinTimeInterval(CHEMICALEVOLUTION_AGB_INTERVAL);
#endif // USE_CELIB_AGB //}

#ifdef USE_CELIB_NSM //{
    CELibSetRunParameterNSMDTDPowerLawIndex(CHEMICALEVOLUTION_NSM_DTD_INDEX);
    CELibSetRunParameterNSMDTDOffsetForPower(CHEMICALEVOLUTION_NSM_DTD_OFFSET);
#endif //USE_CELIB_NSM //}


#ifdef USE_CELIB_ECSN //{
#if USE_CELIB_ECSN_MASS_RANGE==0
    CELibSetRunParameterECSNMassRange(0);
#elif USE_CELIB_ECSN_MASS_RANGE==1
    CELibSetRunParameterECSNMassRange(1);
#else
    CELibSetRunParameterECSNMassRange(2);
#endif
#endif //USE_CELIB_ECSN //}


#ifdef USE_CELIB_HN //{
    CELibSetRunParameterHNYieldsTableID(CHEMICALEVOLUTION_SNII_YIELD_TYPE);
#endif //USE_CELIB_HN //

    if(MPIGetMyID() == MPI_ROOT_RANK){
        CELibSetRunParameterTestMode(true);
    } else {
        CELibSetRunParameterTestMode(false);
    }

    // CELibSetRunParameterTestMode(true);
    CELibInit();

    if(MPIGetMyID() == MPI_ROOT_RANK){
        CELibShowVersion();
        CELibShowCurrentStatus();
    }

    EnergyConvertToSimulationUnit = GetUnitEnergy();
    InternalEnergyConvertToSimulationUnit = GetUnitSpecificEnergy();

    MassConvertToSimulationUnit = MSUN_CGS/Pall.UnitMass;
    MassConvertToMsun = Pall.UnitMass/MSUN_CGS;

    // SNIINumberPerMass = CELibGetSNIINumberPerMass();
    //SNIIEnergy = CELibGetSNIIEnergy();
    SNIIEnergy = 1.e51;
    ECSNEnergy = 9.e49;
    HNEnergy = 1.e52;

#ifdef SNII_PEAK_TEMPERATURE //{
    KernelPeak = KernelStellarFeedback(0,1.0);
#endif // SNII_PEAK_TEMPERATURE //}

    return ;
}

int StellarFeedbackGetIMFType(void){
    return CELibGetRunParameterIMFType();
}

int StellarFeedbackGetSNIaType(void){
    return CELibGetRunParameterSNIaType();
}

#if 0
/*
 * This function returns the next explosion time in year.
 * It is necessary to be called just after a star particle forms and a Type II
 * or Ia SN explosion takes place in order to compute the next explosion time.
 */
double StellarFeedbackGetNextEventTime(const int mode, const double Metallicity, const double InitialMass_in_Msun, const int Count){

    double R = gsl_rng_uniform(RandomGenerator);
    if(mode == CELibFeedbackType_SNII){
        double Time = CELibGetSNIIEventTime(R,Metallicity);
        if(Time < 0.0){
            return 10*Pall.TEnd;
        } else { 
            return Time;
        }
    } else if(mode == CELibFeedbackType_SNII){
#if (CHEMICALEVOLUTION_SNIa_TYPE==-1) //{
        return 10*Pall.TEnd;
#else // CHEMICALEVOLUTION_SNIa_TYPE  //}//{
        // double Time = CELibGetSNIaEventTime(R,Metallicity,InitialMass_in_Msun,Count);
        // fprintf(stderr,"SNIa Next Exp Time = %g, %g, Tcurrent = %g, Count %d \n",Time,
                // Time*YEAR_CGS/Pall.UnitTime,Pall.TCurrent,Count);

        return CELibGetSNIaEventTime(R,Metallicity,InitialMass_in_Msun,Count);
#endif // CHEMICALEVOLUTION_SNIa_TYPE  //}
    } else if(mode == CELibFeedbackType_AGB){
#ifdef USE_CELIB_AGB //{
        double Time = CELibGetAGBFeedbackTime(R,Metallicity,Count);

        // fprintf(stderr,"AGB Next Exp Time = %g, %g, Tcurrent = %g, Count %d \n",Time,
                // Time*YEAR_CGS/Pall.UnitTime,
                // Pall.TCurrent,Count);

        if(Time < 0.0){
            return 10*Pall.TEnd;
        } else { 
            return Time;
        }
#endif // USE_CELIB_AGB  //}
    } else if(mode == CELibFeedbackType_NSM){
#ifdef USE_CELIB_NSM //{
        double Time = CELibGetNSMFeedbackTime(R,Metallicity,InitialMass_in_Msun,Count);
        if(Time < 0.0){
            return 10*Pall.TEnd;
        } else {
            return Time;
        }
#endif // USE_CELIB_NSM  //}
    }
    return NONE;
}
#endif

#if 0
static inline bool __attribute__((always_inline)) CheckEventTime(const int Index){

#ifdef USE_CELIB_AGB // {
    if(Pstar[Index]->EventTime < 0.e0) 
        return false;

    if(Pall.TCurrent > Pstar[Index]->EventTime){
        return true;
    }
    if(Pall.TCurrent > Pstar[Index]->EventTimeAGB){
        return true;
    }

    return false;
#else // USE_CELIB_AGB //}//{
    if(Pstar[Index]->EventTime < 0.e0) 
        return false;

    double Time = Pstar[Index]->EventTime;

    if(Pall.TCurrent > Time){
        return true;
    } else {
        return false;
    }
#endif // USE_CELIB_AGB //}

}
#else
static inline int __attribute__((always_inline)) CheckEventTime(const int Index){

#ifdef USE_CELIB_AGB // {
    if(Pstar[Index]->EventTime < 0.e0) 
        return NONE;

    if(Pall.TCurrent > Pstar[Index]->EventTime){
        if(Pstar[Index]->SNIaCount == -1){ // II 
            return CELibFeedbackType_SNII;
        } else {
            return CELibFeedbackType_SNIa;
        }
    }
    if(Pall.TCurrent > Pstar[Index]->EventTimeAGB){
        return CELibFeedbackType_AGB;
    }

#ifdef USE_CELIB_NSM //{
    if(Pall.TCurrent > Pstar[Index]->EventTimeNSM){
        return CELibFeedbackType_NSM;
    }
#endif //USE_CELIB_NSM //}

#ifdef USE_CELIB_ECSN //{
    if(Pall.TCurrent > Pstar[Index]->EventTimeECSN && Pstar[Index]->ECSNCount == 0){
        return CELibFeedbackType_ECSN;
    }
#endif //USE_CELIB_ECSN //}

#ifdef USE_CELIB_HN //{
    if(Pall.TCurrent > Pstar[Index]->EventTimeHN && Pstar[Index]->HNCount == 0){
        return CELibFeedbackType_HN;
    }
#endif //USE_CELIB_HN //}



    return NONE;
#else // USE_CELIB_AGB //}//{
    if(Pstar[Index]->EventTime < 0.e0) 
        return NONE;

    if(Pall.TCurrent > Pstar[Index]->EventTime){
        if(Pstar[Index]->SNIaCount == -1){ // II 
            return CELibFeedbackType_SNII;
        } else {
            return CELibFeedbackType_SNIa;
        }
    }
    return NONE;
#endif // USE_CELIB_AGB //}

}
#endif

int CountFeedbackNumber(void){
    int counter = 0;
    for(int i=0;i<Pall.Nstars;i++){
        if(PstarBody(i)->Active){
            int CurrentFeedBackType = CheckEventTime(i);
            if(CurrentFeedBackType != NONE){
                // List[counter] = i;
                // FeedBackType[counter] = CurrentFeedBackType;
                counter ++;
            }
        }
    }
    return counter;
}


static int StellarFeedbackNContactedDomains;
static int *StellarFeedbackContactedDomainID;

//#define StellarFeedbackRadiusFactInc   (1.14) // 1.5 ^ (1.3)
#define StellarFeedbackRadiusFactInc   (1.2596) //  2^(0.333)
#define StellarFeedbackRadiusFactDec   (0.79) // 0.75 ^ (1.3)
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER
static double SmoothedMassConversionFactor;
#endif //USE_SMOOTHED_NEIGHBOR_NUMBER

struct StructStellarFeedbackExport{
    double    Radius;  // Feedback radius (2xRadius is the actual feedback radius).
    double    Pos[3];  // Position.
    int       Leaf;
#ifdef __CHECK_WEIGHT__
    int       GlobalID;
#endif // __CHECK_WEIGHT__
};

struct StructStellarFeedbackImport{
    int    Leaf;
    int    Nlist; // Nlist.
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
    double SmoothedNumber; // Smoothed mas..
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
    double Density;     // Mass weighted normalization factor.
#ifdef SET_SNII_TEMPERATURE
    double GasMass;     // Local gas mass
#endif //SET_SNII_TEMPERATURE
#ifdef MAXIMUM_ENERGY_INPUT
    double DistanceMin;
    unsigned long int  DistanceMinGlobalID;
#endif
#ifdef __CHECK_SUM__ //{
    int CheckSum;
#endif // __CHECK_SUM__ //}
};

#if 0
struct StructActiveSNParticle{
    int Index; // 
    int Nlist; // Number of neighbors.
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
    double SmoothedNumber; // Smoothed mas..
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
    double Pos[3];
    double Radius;  // Feedback radius is 2*Radius.
    double Density; // The mass weighted nomalization factor.
#ifdef SET_SNII_TEMPERATURE
    double GasMass; // Local gas mass.
#endif //SET_SNII_TEMPERATURE
#ifdef MAXIMUM_ENERGY_INPUT
    double DistanceMin;
    unsigned long int DistanceMinGlobalID;
#endif
#ifdef __CHECK_SUM__ //{
    int CheckSum;
#endif // __CHECK_SUM__ //}
    bool LocalUpdateFlags; // Flag for the local update.
    // bool TypeII;
    // bool TypeIa;
    int Type; // 1=TypeII, 2=TypeIa, 3=AGB
    int Count;
    double InitialMass;
    double Metallicity;
    double Rvalue;
    double Lvalue;
    int IterationCount;
} *ActiveSNParticle; 
#endif


struct StructStellarFeedbackLocalInfo{
    int Nlist;
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER
    double SmoothedNumber;
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER
    double Density;
#ifdef SET_SNII_TEMPERATURE
    double GasMass;// Local gas mass.
#endif //SET_SNII_TEMPERATURE
#ifdef MAXIMUM_ENERGY_INPUT
    double DistanceMin;
    unsigned long int DistanceMinGlobalID;
#endif
#ifdef __CHECK_SUM__ //{
    int CheckSum;
#endif // __CHECK_SUM__ //}
}; 

static int HydroUpdateFlagSize = 0;
static bool *HydroUpdateFlag;

static inline void __attribute__((always_inline)) StellarFeedbackAllocateContanctedDomainID(void){
    StellarFeedbackContactedDomainID = malloc(sizeof(int)*MPIGetNumProcs());
    return ;
}


static inline bool __attribute__((always_inline)) CheckLocalExternalDomainsContacted(const int MyID, const int ExtID){

    for(int k=0;k<3;k++){
        if((EdgesForHydro[MyID].PosMax[k] < EdgesForHydro[ExtID].PosMin[k])||
           (EdgesForHydro[MyID].PosMin[k] > EdgesForHydro[ExtID].PosMax[k]))  return false;
    }
    return true;

}


/*
 * This function checkes the number of contacted domains by comparing the local
 * domain edge to the external domains. 
 */
static inline void __attribute__((always_inline)) CheckContactedDomain(void){
    int NProcs = MPIGetNumProcs();
    int MyID = MPIGetMyID();
    StellarFeedbackNContactedDomains = 0;
    for(int i=0;i<NProcs-1;i++){
        int NodeID = CommunicationTable[i].SendRank;
        assert(MPIGetMyID() != NodeID);
        if(CheckLocalExternalDomainsContacted(MyID,NodeID)){
            StellarFeedbackContactedDomainID[StellarFeedbackNContactedDomains] = i;
            StellarFeedbackNContactedDomains ++;
        }
    }
    return ;
}


static inline bool __attribute__((always_inline)) OverlapDomainStellarFeedback(double Pos[restrict], const double h, const int NodeID){ 

    double Dist2 = 0.e0;
    for(int k=0;k<3;k++){
        if(Pos[k] < EdgesForHydro[NodeID].PosMin[k]) 
            Dist2 += SQ(EdgesForHydro[NodeID].PosMin[k]-Pos[k]);
        if(Pos[k] > EdgesForHydro[NodeID].PosMax[k])
            Dist2 += SQ(EdgesForHydro[NodeID].PosMax[k]-Pos[k]);
    }
    return (Dist2 < SQ(h));
}

static inline bool __attribute__((always_inline)) CheckInLocalDomain(double Pos[], double Kernel, const int MyID){
    for(int k=0;k<3;k++){
        if(Pos[k]+2.e0*Kernel > EdgesForHydro[MyID].PosMax[k]) return false;
        if(Pos[k]-2.e0*Kernel < EdgesForHydro[MyID].PosMin[k]) return false;
    }
    return true;
}


/*
 * This function return number of particles which should export to the node ID of [NodeIndex].
 */
//#define __EXPORT_ALL__
#if 0
static inline int __attribute__((always_inline)) CheckStellarFeedbackExportFlags(const int NodeIndex, const int NProcs, bool StellarFeedbackExportFlags[][NProcs], const int NActives, struct StructActiveSNParticle ActiveSNParticle[restrict]){

    if(NActives == 0) 
        return 0;

    int ExportNodeID = CommunicationTable[NodeIndex].SendRank;
    int NExport = 0;
    for(int i=0;i<NActives;i++){
        if(StellarFeedbackExportFlags[i][NProcs-1]){
#ifndef __EXPORT_ALL__
            if(OverlapDomainStellarFeedback(ActiveSNParticle[i].Pos,2.0*ActiveSNParticle[i].Radius,ExportNodeID)){
#endif // __EXPORT_ALL__
                StellarFeedbackExportFlags[i][NodeIndex] = true;
                NExport ++;
#ifndef __EXPORT_ALL__
            }
#endif // __EXPORT_ALL__
        }
    }

    return NExport;
}
#endif

// #define __EXPORT_ALL__
static inline int __attribute__((always_inline)) CheckStellarFeedbackExportFlagsModified(const int NodeIndex, const int NProcs, bool StellarFeedbackExportFlags[restrict], const int NActives, struct StructActiveSNParticle ActiveSNParticle[restrict]){

    if(NActives == 0) 
        return 0;

    int ExportNodeID = CommunicationTable[NodeIndex].SendRank;
    // Node By Node Comparison
    double BoxCenter[] = {GravityNode[0].Pos[0],GravityNode[0].Pos[1],GravityNode[0].Pos[2]};
    if(!OverlapDomainStellarFeedback(BoxCenter,LocalKernelMax,ExportNodeID)){
        return 0;
    }

    int NExport = 0;
    for(int i=0;i<NActives;i++){
        int Offset = i*NProcs;
        //if(StellarFeedbackExportFlags[i][NProcs-1]){
        if(StellarFeedbackExportFlags[Offset+NProcs-1]){
#ifndef __EXPORT_ALL__
            if(OverlapDomainStellarFeedback(ActiveSNParticle[i].Pos,2.0*ActiveSNParticle[i].Radius,ExportNodeID)){
#endif // __EXPORT_ALL__
                //ExportFlags[i][NodeIndex] = true;
                StellarFeedbackExportFlags[Offset+NodeIndex] = true;
                NExport ++;
#ifndef __EXPORT_ALL__
            }
#endif // __EXPORT_ALL__
        }
    }

    return NExport;
}

#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
static int GetSmoothedNumberofNeighbors(double Pos[restrict], const double Kernel, int Neighbors[restrict], double *SmoothedNumber){

    int nlist = 0;
    int RootNodeID = 0;
    int CurrentNodeID = HydroNode[RootNodeID].Children;
    do {
        CurrentNodeID = GetNeighborsSmoothedNumberIterativeApproach(CurrentNodeID,
                Pos,2.e0*Kernel,&nlist,Neighbors,SmoothedNumber);
    }while(CurrentNodeID != RootNodeID);

    return nlist;
}
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}


/*
 * This function returns a structure which involves various values (the neighbor
 * number, the smootehd neighbor number, the density, the local gas mass, the
 * minimum distance, and the globalID of the closest particle).
 */
struct StructStellarFeedbackLocalInfo RetrunStellarFeedbackLocalInfo(double Pos[restrict], const double Radius){

    static int Neighbors[MaxNeighborSize];
    struct StructStellarFeedbackLocalInfo TempStellarFeedbackLocalInfo = {
        .Nlist = 0,
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
        .SmoothedNumber = 0.0,
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
        .Density = 0.e0,
#ifdef SET_SNII_TEMPERATURE //{
        .GasMass = 0.e0,
#endif //SET_SNII_TEMPERATURE //}
#ifdef __CHECK_SUM__ //{
        .CheckSum = 0,
#endif //__CHECK_SUM__ //}
    };

    double InvRadiusi = 1.e0/Radius;
    TempStellarFeedbackLocalInfo.Nlist = 
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
            GetSmoothedNumberofNeighbors(Pos,Radius,Neighbors,&(TempStellarFeedbackLocalInfo.SmoothedNumber));
#else // USE_SMOOTHED_NEIGHBOR_NUMBER 
            GetNeighborsLimited(Pos,2.0*Radius,Neighbors);
            //GetNeighborsDirect(Pos,2*Radius,Neighbors);
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}

    for(int i=0;i<TempStellarFeedbackLocalInfo.Nlist;i++){
        int leaf = Neighbors[i];
        double xij[3];
#ifdef PERIODIC_RUN
        xij[0] = PeriodicDistance(Pos[0],PhydroPosP(leaf)[0],0);
        xij[1] = PeriodicDistance(Pos[1],PhydroPosP(leaf)[1],1);
        xij[2] = PeriodicDistance(Pos[2],PhydroPosP(leaf)[2],2);
#else // PERIODIC_RUN
        xij[0] = Pos[0]-PhydroPosP(leaf)[0];
        xij[1] = Pos[1]-PhydroPosP(leaf)[1];
        xij[2] = Pos[2]-PhydroPosP(leaf)[2];
#endif // PERIODIC_RUN

        double r = NORM(xij);
        double w = KernelStellarFeedback(r,InvRadiusi);

        TempStellarFeedbackLocalInfo.Density += PhydroMass(leaf)*w;
        //assert(PhydroMass(leaf)*w > 0.e0);
#ifdef SET_SNII_TEMPERATURE
        TempStellarFeedbackLocalInfo.GasMass += PhydroMass(leaf);
#endif //SET_SNII_TEMPERATURE
#ifdef MAXIMUM_ENERGY_INPUT
        if(i==0){
            TempStellarFeedbackLocalInfo.DistanceMinID = PhydroBody(leaf)->GlobalID;
            TempStellarFeedbackLocalInfo.DistanceMin = r;
        } else {
            if (TempStellarFeedbackLocalInfo.DistanceMin > r){
                TempStellarFeedbackLocalInfo.DistanceMinID = PhydroBody(leaf)->GlobalID;
                TempStellarFeedbackLocalInfo.DistanceMin = r;
            }
        }
#endif
#ifdef __CHECK_SUM__ //{
        TempStellarFeedbackLocalInfo.CheckSum += PhydroBody(leaf)->GlobalID;
#endif // __CHECK_SUM__ //}
    }

    return TempStellarFeedbackLocalInfo;
}


/*
 * This function checks the convergence of the feedback radius, using bisetion
 * method. If the radius is convergered, this function returns a "true" flag.
 * If not, it returns a "false" flag.
 */
static inline bool __attribute__((always_inline)) CheckNeighborNumberAndUpdateFeedbackRadius_i(const int NProcs, bool StellarFeedbackExportFlags_i[restrict], struct StructActiveSNParticle *ActiveSNParticle_i){ 

#ifdef USE_SN_INPUT_PARTICLE_NUMBER //{
    int NBmin = SN_INPUT_PARTICLE_NUMBER-SN_INPUT_PARTICLE_NUMBER_MARGIN;
    int NBmax = SN_INPUT_PARTICLE_NUMBER+SN_INPUT_PARTICLE_NUMBER_MARGIN;
#else // USE_SN_INPUT_PARTICLE_NUMBER
    int NBmin = Pall.Ns-Pall.Npm;
    int NBmax = Pall.Ns+Pall.Npm;
#endif // USE_SN_INPUT_PARTICLE_NUMBER //}

#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
    double Nlist = ActiveSNParticle_i->SmoothedNumber*
        SmoothedMassConversionFactor*CUBE(ActiveSNParticle_i->Radius);
#else // USE_SMOOTHED_NEIGHBOR_NUMBER
    int Nlist = ActiveSNParticle_i->Nlist;
#endif  // USE_SMOOTHED_NEIGHBOR_NUMBER //}

    // Here, the convergence condition is satisfied. 
    if(((NBmin)<=Nlist)&&(Nlist<=(NBmax))){
        StellarFeedbackExportFlags_i[NProcs-1] = false;
        return true;
    //}else if((Nlist<=(NBmax))&&(ActiveSNParticle_i->Rvalue>0.e0)&&(ActiveSNParticle_i->Lvalue>0.e0)){
    }else if(((0.8*NBmin)<=Nlist)&&(Nlist<=(2.0*NBmax))&&(ActiveSNParticle_i->Rvalue>0.e0)&&(ActiveSNParticle_i->Lvalue>0.e0)){
        if(ActiveSNParticle_i->Rvalue-ActiveSNParticle_i->Lvalue < 1.e-6*ActiveSNParticle_i->Lvalue){
            StellarFeedbackExportFlags_i[NProcs-1] = false;
            return true;
        }
    }
    if(StellarFeedbackExportFlags_i[NProcs-1]){
        if(Nlist<NBmin){
            ActiveSNParticle_i->Lvalue = fmax(ActiveSNParticle_i->Lvalue,ActiveSNParticle_i->Radius);
        } else if(Nlist>NBmax){
            if(ActiveSNParticle_i->Rvalue > 0.e0){
                ActiveSNParticle_i->Rvalue = fmin(ActiveSNParticle_i->Rvalue,ActiveSNParticle_i->Radius);
            }else{
                ActiveSNParticle_i->Rvalue = ActiveSNParticle_i->Radius;
            }
        }

        if((ActiveSNParticle_i->Lvalue>0.e0)&&(ActiveSNParticle_i->Rvalue>0.e0)){
            ActiveSNParticle_i->Radius = cbrt(0.5*(CUBE(ActiveSNParticle_i->Lvalue)+CUBE(ActiveSNParticle_i->Rvalue)));
        }else{
            if((ActiveSNParticle_i->Rvalue == 0.e0)&&(ActiveSNParticle_i->Lvalue > 0.e0)){
                ActiveSNParticle_i->Radius *= StellarFeedbackRadiusFactInc;
            }else if((ActiveSNParticle_i->Rvalue > 0.e0)&&(ActiveSNParticle_i->Lvalue == 0.e0)){
                ActiveSNParticle_i->Radius *= StellarFeedbackRadiusFactDec;
            }
        }
    }
    return false;
}


/*
 * This function calls when the feedback radius of a SN particle is enclosed in
 * the local domain.
 */
#if 0
static inline void __attribute__((always_inline)) UpdateStellarFeedbackRadiusLocal(const int Index, struct StructActiveSNParticle ActiveSNParticle[restrict], int Neighbors[restrict], const int MyID, const int NProcs, bool StellarFeedbackExportFlags[restrict][NProcs]){

    if(CheckNeighborNumberAndUpdateFeedbackRadius_i(NProcs,
                StellarFeedbackExportFlags[Index],ActiveSNParticle+Index) == true)
        return;

    do{
        if(!CheckInLocalDomain(ActiveSNParticle[Index].Pos,ActiveSNParticle[Index].Radius,MyID)) return;
        for(int i=0;i<StellarFeedbackNContactedDomains;i++){
            int NodeID = StellarFeedbackContactedDomainID[i];
            if(OverlapDomainStellarFeedback(ActiveSNParticle[Index].Pos,2.0*ActiveSNParticle[Index].Radius,NodeID)) return;
        }
        struct StructStellarFeedbackLocalInfo TemporalData 
            = RetrunStellarFeedbackLocalInfo(ActiveSNParticle[Index].Pos,ActiveSNParticle[Index].Radius);
        ActiveSNParticle[Index].Nlist = TemporalData.Nlist;
        ActiveSNParticle[Index].Density = TemporalData.Density;
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
        ActiveSNParticle[Index].SmoothedNumber = TemporalData.SmoothedNumber;
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
#ifdef SET_SNII_TEMPERATURE //{
        ActiveSNParticle[Index].GasMass = TemporalData.GasMass;
#endif //SET_SNII_TEMPERATURE //}
#ifdef MAXIMUM_ENERGY_INPUT //{
        if(TemporalData.Nlist > 0){
            ActiveSNParticle[Index].DistanceMin = TemporalData.DistanceMin;
            ActiveSNParticle[Index].DistanceMinGlobalID = TemporalData.DistanceMinGlobalID;
        }
#endif // MAXIMUM_ENERGY_INPUT //}
#ifdef __CHECK_SUM__ //{
        ActiveSNParticle[Index].CheckSum = TemporalData.CheckSum;
#endif // __CHECK_SUM__ //}
        ActiveSNParticle[Index].IterationCount ++;
    }while(CheckNeighborNumberAndUpdateFeedbackRadius_i(NProcs,
                StellarFeedbackExportFlags[Index],ActiveSNParticle+Index) == false);

    return;
}
#endif

static inline void __attribute__((always_inline)) UpdateStellarFeedbackRadiusLocalModified(const int Index, const int NActives, 
        struct StructActiveSNParticle ActiveSNParticle[restrict], int Neighbors[restrict], const int MyID, const int NProcs, 
            bool StellarFeedbackExportFlags[restrict]){

    //if(CheckNeighborNumberAndUpdateFeedbackRadius_i(NProcs,
                //StellarFeedbackExportFlags[Index],ActiveSNParticle+Index) == true)
    if(CheckNeighborNumberAndUpdateFeedbackRadius_i(NProcs,
                StellarFeedbackExportFlags+Index*NProcs,ActiveSNParticle+Index) == true)
        return;

    do{
        if(!CheckInLocalDomain(ActiveSNParticle[Index].Pos,ActiveSNParticle[Index].Radius,MyID)) return;
        for(int i=0;i<StellarFeedbackNContactedDomains;i++){
            int NodeID = StellarFeedbackContactedDomainID[i];
            if(OverlapDomainStellarFeedback(ActiveSNParticle[Index].Pos,2.0*ActiveSNParticle[Index].Radius,NodeID)) return;
        }
        struct StructStellarFeedbackLocalInfo TemporalData 
            = RetrunStellarFeedbackLocalInfo(ActiveSNParticle[Index].Pos,ActiveSNParticle[Index].Radius);
        ActiveSNParticle[Index].Nlist = TemporalData.Nlist;
        ActiveSNParticle[Index].Density = TemporalData.Density;
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
        ActiveSNParticle[Index].SmoothedNumber = TemporalData.SmoothedNumber;
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
#ifdef SET_SNII_TEMPERATURE //{
        ActiveSNParticle[Index].GasMass = TemporalData.GasMass;
#endif //SET_SNII_TEMPERATURE //}
#ifdef MAXIMUM_ENERGY_INPUT //{
        if(TemporalData.Nlist > 0){
            ActiveSNParticle[Index].DistanceMin = TemporalData.DistanceMin;
            ActiveSNParticle[Index].DistanceMinGlobalID = TemporalData.DistanceMinGlobalID;
        }
#endif // MAXIMUM_ENERGY_INPUT //}
#ifdef __CHECK_SUM__ //{
        ActiveSNParticle[Index].CheckSum = TemporalData.CheckSum;
#endif // __CHECK_SUM__ //}
        ActiveSNParticle[Index].IterationCount ++;
    }while(CheckNeighborNumberAndUpdateFeedbackRadius_i(NProcs,
                StellarFeedbackExportFlags+Index*NProcs,ActiveSNParticle+Index) == false);

    return;
}

#if 0
static inline int __attribute__((always_inline)) CheckNeighborNumberAndUpdateFeedbackRadius(const int NActives, const int NProcs, bool StellarFeedbackExportFlags[restrict][NProcs], struct StructActiveSNParticle ActiveSNParticle[restrict]){ 

    int NBmin = Pall.Ns-Pall.Npm;
    int NBmax = Pall.Ns+Pall.Npm;

    int NLocalActiveLeaves = 0;
    for(int i=0;i<NActives;i++){
        if(StellarFeedbackExportFlags[i][NProcs-1]){ 
#ifdef USE_STELLARFEEDBACK_RADIUS_LOCAL_UPDATE //{
        if(ActiveSNParticle[i].LocalUpdateFlags == false){ 
#endif //USE_STELLARFEEDBACK_RADIUS_LOCAL_UPDATE //}
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
            double Nlist = ActiveSNParticle[i].SmoothedNumber*
                SmoothedMassConversionFactor*CUBE(ActiveSNParticle[i].Radius);
#else
            int Nlist = ActiveSNParticle[i].Nlist;
#endif //USE_MOOTHED_NEIGHBOR_NUMBER //}

            if(((NBmin)<=Nlist)&&(Nlist<=(NBmax))){
                StellarFeedbackExportFlags[i][NProcs-1] = false;
            }else if((Nlist<=(NBmax))&&(ActiveSNParticle[i].Rvalue>0.e0)&&(ActiveSNParticle[i].Lvalue>0.e0)){
                if(ActiveSNParticle[i].Rvalue-ActiveSNParticle[i].Lvalue < 1.e-3*ActiveSNParticle[i].Lvalue)
                    StellarFeedbackExportFlags[i][NProcs-1] = false;
            }
            if(StellarFeedbackExportFlags[i][NProcs-1]){
                if(Nlist<NBmin){
                    ActiveSNParticle[i].Lvalue = fmax(ActiveSNParticle[i].Lvalue,ActiveSNParticle[i].Radius);
                } else if(Nlist>NBmax){
                    if(ActiveSNParticle[i].Rvalue > 0.e0){
                        ActiveSNParticle[i].Rvalue = fmin(ActiveSNParticle[i].Rvalue,ActiveSNParticle[i].Radius);
                    }else{
                        ActiveSNParticle[i].Rvalue = ActiveSNParticle[i].Radius;
                    }
                }

                if((ActiveSNParticle[i].Lvalue>0.e0)&&(ActiveSNParticle[i].Rvalue>0.e0)){
                    ActiveSNParticle[i].Radius = cbrt(0.5*(CUBE(ActiveSNParticle[i].Lvalue)+CUBE(ActiveSNParticle[i].Rvalue)));
                }else{
                    if((ActiveSNParticle[i].Rvalue == 0.e0)&&(ActiveSNParticle[i].Lvalue > 0.e0)){
                        ActiveSNParticle[i].Radius *= StellarFeedbackRadiusFactInc;
                    }else if((ActiveSNParticle[i].Rvalue > 0.e0)&&(ActiveSNParticle[i].Lvalue == 0.e0)){
                        ActiveSNParticle[i].Radius *= StellarFeedbackRadiusFactDec;
                    }
                }
                NLocalActiveLeaves ++;
            }
            ActiveSNParticle[i].IterationCount ++;
#ifdef USE_STELLARFEEDBACK_RADIUS_LOCAL_UPDATE //{
        } else {
            NLocalActiveLeaves ++;
        }
#endif //USE_STELLARFEEDBACK_RADIUS_LOCAL_UPDATE //}
        }
    }

    return NLocalActiveLeaves;
}
#endif

static inline int __attribute__((always_inline)) CheckNeighborNumberAndUpdateFeedbackRadiusModified(const int NActives, const int NProcs, bool StellarFeedbackExportFlags[restrict], struct StructActiveSNParticle ActiveSNParticle[restrict]){ 

#ifdef USE_SN_INPUT_PARTICLE_NUMBER //{
    int NBmin = SN_INPUT_PARTICLE_NUMBER-SN_INPUT_PARTICLE_NUMBER_MARGIN;
    int NBmax;
    if(OverMaxIterationTimes){
        NBmax = SN_INPUT_PARTICLE_NUMBER+2*SN_INPUT_PARTICLE_NUMBER_MARGIN;
    } else {
        NBmax = SN_INPUT_PARTICLE_NUMBER+SN_INPUT_PARTICLE_NUMBER_MARGIN;
    }
#else // USE_SN_INPUT_PARTICLE_NUMBER
    int NBmin = Pall.Ns-Pall.Npm;
    int NBmax = Pall.Ns+Pall.Npm;
#endif // USE_SN_INPUT_PARTICLE_NUMBER //}

    int NLocalActiveLeaves = 0;
    for(int i=0;i<NActives;i++){
        int Offset = i*NProcs;
        //if(StellarFeedbackExportFlags[i][NProcs-1]){ 
        if(StellarFeedbackExportFlags[Offset+NProcs-1]){ 
#ifdef USE_STELLARFEEDBACK_RADIUS_LOCAL_UPDATE //{
        if(ActiveSNParticle[i].LocalUpdateFlags == false){ 
#endif //USE_STELLARFEEDBACK_RADIUS_LOCAL_UPDATE //}
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
            double Nlist = ActiveSNParticle[i].SmoothedNumber*
                SmoothedMassConversionFactor*CUBE(ActiveSNParticle[i].Radius);
#else
            int Nlist = ActiveSNParticle[i].Nlist;
#endif //USE_MOOTHED_NEIGHBOR_NUMBER //}

            if(((NBmin)<=Nlist)&&(Nlist<=(NBmax))){
                StellarFeedbackExportFlags[Offset+NProcs-1] = false;
            }else if(((0.8*NBmin)<=Nlist)&&(Nlist<=(2.0*NBmax))&&(ActiveSNParticle[i].Rvalue>0.e0)&&(ActiveSNParticle[i].Lvalue>0.e0)){
                if(ActiveSNParticle[i].Rvalue-ActiveSNParticle[i].Lvalue < 1.e-6*ActiveSNParticle[i].Lvalue){
                    StellarFeedbackExportFlags[Offset+NProcs-1] = false;
                }
            }
            if(StellarFeedbackExportFlags[Offset+NProcs-1]){
                if(Nlist<NBmin){
                    ActiveSNParticle[i].Lvalue = fmax(ActiveSNParticle[i].Lvalue,ActiveSNParticle[i].Radius);
                } else if(Nlist>NBmax){
                    if(ActiveSNParticle[i].Rvalue > 0.e0){
                        ActiveSNParticle[i].Rvalue = fmin(ActiveSNParticle[i].Rvalue,ActiveSNParticle[i].Radius);
                    }else{
                        ActiveSNParticle[i].Rvalue = ActiveSNParticle[i].Radius;
                    }
                }

                if((ActiveSNParticle[i].Lvalue>0.e0)&&(ActiveSNParticle[i].Rvalue>0.e0)){
                    ActiveSNParticle[i].Radius = cbrt(0.5*(CUBE(ActiveSNParticle[i].Lvalue)+CUBE(ActiveSNParticle[i].Rvalue)));
                }else{
                    if((ActiveSNParticle[i].Rvalue == 0.e0)&&(ActiveSNParticle[i].Lvalue > 0.e0)){
                        ActiveSNParticle[i].Radius *= StellarFeedbackRadiusFactInc;
                    }else if((ActiveSNParticle[i].Rvalue > 0.e0)&&(ActiveSNParticle[i].Lvalue == 0.e0)){
                        ActiveSNParticle[i].Radius *= StellarFeedbackRadiusFactDec;
                    }
                }
                NLocalActiveLeaves ++;
            }
            ActiveSNParticle[i].IterationCount ++;
#ifdef USE_STELLARFEEDBACK_RADIUS_LOCAL_UPDATE //{
        } else {
            NLocalActiveLeaves ++;
        }
#endif //USE_STELLARFEEDBACK_RADIUS_LOCAL_UPDATE //}
        }
    }

    return NLocalActiveLeaves;

}

static int MaxIterationTimesForInitialGuess = 10;
double StellarFeedbackRadiusInitialGuess(const int Index){

    int Neighbors[MaxNeighborSize];
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
    double SmoothedNumber;
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
    //double Radius = 5.0*PstarBody(Index)->Eps; 
    double Radius = 2.0*PstarBody(Index)->Eps; 

    if(Pall.Nhydro == 0)
        return Radius;

    double Pos[3] = {PstarBody(Index)->Pos[0],PstarBody(Index)->Pos[1],PstarBody(Index)->Pos[2]};
    int Iteration = 0;
    do{ 
        int Nlist = 
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
                GetSmoothedNumberofNeighbors(Pos,Radius,Neighbors,&SmoothedNumber);
#else // USE_SMOOTHED_NEIGHBOR_NUMBER 
                GetNeighborsLimited(Pos,2.0*Radius,Neighbors);
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
        if(Nlist > 0){
            double Dist2 = DISTANCE2(Pos,NBCache[Neighbors[0]].Pos);
            double RadiusMin = NBCache[Neighbors[0]].Kernel;
            for(int i=1;i<Nlist;i++){
                int leaf = Neighbors[i];
                double CurrentDist2 = DISTANCE2(Pos,NBCache[leaf].Pos);
                if(Dist2 > CurrentDist2){
                    Dist2 = CurrentDist2;
                    RadiusMin = NBCache[leaf].Kernel;
                }
            }
            // dprintlmpi(Nlist);
            // dprintlmpi(Iteration);
            return RadiusMin;
        }
        Radius *= StellarFeedbackRadiusFactInc; 
#if 0
        if(Nlist > 2*Pall.Ns){
            double Dist2 = DISTANCE2(Pos,NBCache[Neighbors[0]].Pos);
            double RadiusMin = NBCache[Neighbors[0]].Kernel;
            for(int i=1;i<Nlist;i++){
                int leaf = Neighbors[i];
                double CurrentDist2 = DISTANCE2(Pos,NBCache[Neighbors[i]].Pos);
                if(Dist2 > CurrentDist2){
                    Dist2 = CurrentDist2;
                    RadiusMin = Radius;
                }
            }
            return RadiusMin;
        } else if(Nlist == 0){
            return 2.0*Radius;
        } else {
            return NBCache[Neighbors[0]].Kernel;
        }
#endif
        Iteration ++;
    } while(Iteration < MaxIterationTimesForInitialGuess);
    return 2.0*Radius;
}

#define _ResetTiming_ 2
static void ResetKernelSize(const int NActives, const int Niteration, const int NProcs, 
        const int IndexList[restrict], bool StellarFeedbackExportFlags[restrict]){

    if((Niteration+1)%(_ResetTiming_*MaxIterationTimes) == 0){
        for(int i=0;i<NActives;i++){
            int Offset = i*NProcs;
            if(StellarFeedbackExportFlags[Offset+NProcs-1]){ 
                ActiveSNParticle[i].Rvalue = ActiveSNParticle[i].Lvalue = 0.e0;
                ActiveSNParticle[i].Radius = 2.0*PstarBody(IndexList[i])->Eps
                                            *(gsl_rng_uniform(RandomGenerator)+0.5);
            }
        }
    }

}


//static bool **StellarFeedbackExportFlagsLog;
//static bool *StellarFeedbackExportFlagsLog; //! Export data log for stellar feedback.
static bool *StellarFeedbackExportFlags;

static bool first = true;
static void CalcFeedbackRadius(int NActives, const int IndexList[restrict], const int TypeList[restrict]){

    if(first){
        StellarFeedbackAllocateContanctedDomainID();
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER
        SmoothedMassConversionFactor = (4.0*M_PI/3.0)*8.0;
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER
        first = false;
    }

    OverMaxIterationTimes = false;

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();
    int Neighbors[MaxNeighborSize];

    //static int StellarFeedbackExportFlagsMaxAllocated = 0;
    //static bool (*StellarFeedbackExportFlags)[NProcs];

    if(StellarFeedbackExportFlagsMaxAllocated < MAX(NActives,NAdditionUnit)){
        if(StellarFeedbackExportFlagsMaxAllocated > 0){
            free(StellarFeedbackExportFlags);
            free(ActiveSNParticle);
        }
        StellarFeedbackExportFlagsMaxAllocated = (int)(MAX(ForAngelsShare*NActives,NAdditionUnit));
        StellarFeedbackExportFlags = malloc(sizeof(bool)*StellarFeedbackExportFlagsMaxAllocated*NProcs);
        ActiveSNParticle = malloc(sizeof(struct StructActiveSNParticle)*StellarFeedbackExportFlagsMaxAllocated);
        //StellarFeedbackExportFlagsLog = (bool **)StellarFeedbackExportFlags;
        // StellarFeedbackExportFlagsLog = StellarFeedbackExportFlags;
    }

    for(int k=0;k<NActives;k++){
        int Offset = k*NProcs;
        StellarFeedbackExportFlags[Offset+NProcs-1] = ON;
        ActiveSNParticle[k].Index = IndexList[k];
        ActiveSNParticle[k].Pos[0] = PstarPos(IndexList[k])[0];
        ActiveSNParticle[k].Pos[1] = PstarPos(IndexList[k])[1]; 
        ActiveSNParticle[k].Pos[2] = PstarPos(IndexList[k])[2]; 
        //ActiveSNParticle[k].Radius = 2.0*PstarBody(IndexList[k])->Eps;
        ActiveSNParticle[k].Radius = StellarFeedbackRadiusInitialGuess(IndexList[k]);
        // ActiveSNParticle[k].TypeII = ActiveSNParticle[k].TypeIa = false;
        /*
        if(Pstar[IndexList[k]]->SNIaCount == -1){
            ActiveSNParticle[k].TypeII = true;
        } else {
            ActiveSNParticle[k].TypeIa = true;
        }
        */
        ActiveSNParticle[k].Type = TypeList[k];
        if(TypeList[k] == CELibFeedbackType_AGB){
            ActiveSNParticle[k].Count = Pstar[IndexList[k]]->AGBCount;
#ifdef USE_CELIB_NSM //{
        } else if(TypeList[k] == CELibFeedbackType_NSM){
            ActiveSNParticle[k].Count = Pstar[IndexList[k]]->NSMCount;
#endif // USE_CELIB_NSM //}
#ifdef USE_CELIB_ECSN //{
        } else if(TypeList[k] == CELibFeedbackType_ECSN){
            ActiveSNParticle[k].Count = Pstar[IndexList[k]]->ECSNCount;
#endif // USE_CELIB_ECSN //}
#ifdef USE_CELIB_HN //{
        } else if(TypeList[k] == CELibFeedbackType_HN){
            ActiveSNParticle[k].Count = Pstar[IndexList[k]]->HNCount;
#endif // USE_CELIB_HN //}

        } else {
            ActiveSNParticle[k].Count = Pstar[IndexList[k]]->SNIaCount;
        }
        ActiveSNParticle[k].InitialMass = Pstar[IndexList[k]]->InitialMass;
        ActiveSNParticle[k].Metallicity = Pstar[IndexList[k]]->Z;
        ActiveSNParticle[k].Rvalue = ActiveSNParticle[k].Lvalue = 0.e0;
        ActiveSNParticle[k].IterationCount = 0;
    }

    int BitMask = 0x01; 
    int NExportThisTime[NProcs-1];
    int NImportThisTime[NProcs-1];
    int NExportThisTimeNew[NProcs];
    int NImportThisTimeNew[NProcs];

    struct StructStellarFeedbackExport *StellarFeedbackExportSend[NProcs-1];
    struct StructStellarFeedbackExport *StellarFeedbackExportRecv = NULL;
    struct StructStellarFeedbackImport *StellarFeedbackImportSend = NULL;
    struct StructStellarFeedbackImport *StellarFeedbackImportRecv[NProcs-1];
    MPI_Status mpi_status_Export_Send[NProcs-1];
    MPI_Request mpi_request_Export_Send[NProcs-1];
    MPI_Status mpi_status_Export_Recv[NProcs-1];
    MPI_Request mpi_request_Export_Recv[NProcs-1];

#ifdef PRINT_LOG_KERNEL_ITERATION
    if(MPIGetMyID()==MPI_ROOT_RANK)
        fprintf(stderr,"Feedback Iteration");
#endif // PRINT_LOG_KERNEL_ITERATION

    int NActiveLeaves;
    MPI_Allreduce(&NActives,&NActiveLeaves,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

    CheckContactedDomain();
    // Node center
    double BoxCenter[] = {GravityNode[0].Pos[0],GravityNode[0].Pos[1],GravityNode[0].Pos[2]};

    int Niteration = 0;
    do{
#ifdef PRINT_LOG_KERNEL_ITERATION
        if(MPIGetMyID()==MPI_ROOT_RANK)
            fprintf(stderr,":[%d] = %d ",Niteration,NActiveLeaves);

#endif // PRINT_LOG_KERNEL_ITERATION
        for(int i=0;i<NActives;i++){ 
            //int Offset = i*NActives;
            int Offset = i*NProcs;
            //if(StellarFeedbackExportFlags[i][NProcs-1]){
            if(StellarFeedbackExportFlags[Offset+NProcs-1]){
                if(Niteration > MaxIterationTimes){
// #ifdef SET_SNII_TEMPERATURE
                    // fprintf(stderr,"%ld %d %g %g %g\n",PstarBody(ActiveSNParticle[i].Index)->GlobalID,
                        // ActiveSNParticle[i].Nlist,ActiveSNParticle[i].GasMass,
                        // ActiveSNParticle[i].Lvalue,ActiveSNParticle[i].Rvalue);
// #else 
                    // fprintf(stderr,"%ld %d | %g %g\n",PstarBody(ActiveSNParticle[i].Index)->GlobalID,
                        // ActiveSNParticle[i].Nlist,
                        // ActiveSNParticle[i].Lvalue,ActiveSNParticle[i].Rvalue);
// #endif 
                    // fflush(NULL);
                }
                ActiveSNParticle[i].Nlist = 0;
                ActiveSNParticle[i].Density = 0.0;
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
                ActiveSNParticle[i].SmoothedNumber = 0.e0;
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
#ifdef SET_SNII_TEMPERATURE //{
                ActiveSNParticle[i].GasMass = 0.e0;
#endif //SET_SNII_TEMPERATURE //}
#ifdef MAXIMUM_ENERGY_INPUT //{
                ActiveSNParticle[i].DistanceMin = 2.0*ActiveSNParticle[i].Radius;
#endif // MAXIMUM_ENERGY_INPUT //}
                for(int k=0;k<NProcs-1;k++)
                    StellarFeedbackExportFlags[Offset+k] = false;
                LocalKernelMax = fmax(LocalKernelMax,
                        DISTANCE(BoxCenter,ActiveSNParticle[i].Pos)+2.0*ActiveSNParticle[i].Radius);
#ifdef __CHECK_SUM__ //{
                ActiveSNParticle[i].CheckSum = 0;
#endif // __CHECK_SUM__ //}
            }
        }

        bool ExportFlags[NActives];

        for(int i=0;i<NProcs-1;i++){
            memset(ExportFlags,0,sizeof(bool)*NActives);

            NExportThisTime[i] = CheckStellarFeedbackExportFlagsModified(i,
                    NProcs,StellarFeedbackExportFlags,NActives,ActiveSNParticle);

            CheckSizeofBufferExportSendIndex(NExportThisTime[i],
                    sizeof(struct StructStellarFeedbackExport),i);
            CheckSizeofBufferImportRecvIndex(NExportThisTime[i],
                    sizeof(struct StructStellarFeedbackImport),i);
            StellarFeedbackExportSend[i] = BufferExportSend[i];
            StellarFeedbackImportRecv[i] = BufferImportRecv[i];

            int NExport = 0;
            if(NExportThisTime[i] > 0){
            for(int k=0;k<NActives;k++){
                int Offset = k*NProcs;
                if(StellarFeedbackExportFlags[Offset+NProcs-1]){
                    if(StellarFeedbackExportFlags[Offset+i]&BitMask){ 
                        StellarFeedbackExportSend[i][NExport].Pos[0] = ActiveSNParticle[k].Pos[0];
                        StellarFeedbackExportSend[i][NExport].Pos[1] = ActiveSNParticle[k].Pos[1];
                        StellarFeedbackExportSend[i][NExport].Pos[2] = ActiveSNParticle[k].Pos[2];
                        StellarFeedbackExportSend[i][NExport].Radius = ActiveSNParticle[k].Radius;
                        StellarFeedbackExportSend[i][NExport].Leaf = k;
#ifdef __CHECK_WEIGHT__
                        StellarFeedbackExportSend[i][NExport].GlobalID =
                            PstarBody(ActiveSNParticle[k].Index)->GlobalID;
#endif //__CHECK_WEIGHT__
                        NExport ++;
                    }
                }
            }
            }
            NExportThisTime[i] = NExport;
        }

        int NImportThisTime2[NProcs];
        int NExportThisTime2[NProcs];
        NImportThisTime2[MPIGetMyID()] = 0;
        for(int i=0;i<NProcs-1;i++){
            NExportThisTime2[CommunicationTable[i].SendRank] = NExportThisTime[i];
        }
        MPI_Alltoall(NExportThisTime2,1,MPI_INT,NImportThisTime2,1,MPI_INT,MPI_COMM_WORLD);
        int NImport = 0;
        for(int i=0;i<NProcs-1;i++){
            NImportThisTime[i] = NImportThisTime2[CommunicationTable[i].RecvRank];
            NImport += NImportThisTime[i];
        }
        int NImportAll = NImport;

        CheckSizeofBufferExportRecv(NImportAll,sizeof(struct StructStellarFeedbackExport));
        CheckSizeofBufferImportSend(NImportAll,sizeof(struct StructStellarFeedbackImport));
        StellarFeedbackExportRecv = BufferExportRecv;
        StellarFeedbackImportSend = BufferImportSend; 

        NImport = 0;
        int counter_send = 0;
        int counter_recv = 0;

        int SendFlag,RecvFlag;
        for(int i=0;i<NProcs-1;i++){
            if(NExportThisTime[i]>0){
                MPI_Isend(StellarFeedbackExportSend[i],
                    NExportThisTime[i]*sizeof(struct StructStellarFeedbackExport),
                        MPI_BYTE,CommunicationTable[i].SendRank,TAG_STELLARFEEDBACK_EXPORT+i,
                            MPI_COMM_WORLD,mpi_request_Export_Send+counter_send);
                MPI_Test(mpi_request_Export_Send+counter_send,&SendFlag,MPI_STATUS_IGNORE);
                counter_send ++;
            }
            if(NImportThisTime[i]>0){
                MPI_Irecv(StellarFeedbackExportRecv+NImport,
                    NImportThisTime[i]*sizeof(struct StructStellarFeedbackExport),
                        MPI_BYTE,CommunicationTable[i].RecvRank,TAG_STELLARFEEDBACK_EXPORT+i,
                            MPI_COMM_WORLD,mpi_request_Export_Recv+counter_recv);
                MPI_Test(mpi_request_Export_Recv+counter_recv,&RecvFlag,MPI_STATUS_IGNORE);
                counter_recv ++;
                NImport += NImportThisTime[i];
            }
        }

        for(int i=0;i<NActives;i++){  // Check local
            int Offset = i*NProcs;
            //if(StellarFeedbackExportFlags[i][NProcs-1]){
            if(StellarFeedbackExportFlags[Offset+NProcs-1]){
                struct StructStellarFeedbackLocalInfo TemporalData = 
                    RetrunStellarFeedbackLocalInfo(ActiveSNParticle[i].Pos,ActiveSNParticle[i].Radius);
                ActiveSNParticle[i].Nlist = TemporalData.Nlist;
                ActiveSNParticle[i].Density = TemporalData.Density;

#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
                ActiveSNParticle[i].SmoothedNumber = TemporalData.SmoothedNumber;
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
#ifdef SET_SNII_TEMPERATURE //{
                ActiveSNParticle[i].GasMass = TemporalData.GasMass;
#endif //SET_SNII_TEMPERATURE //}
#ifdef MAXIMUM_ENERGY_INPUT //{
                if(TemporalData.Nlist > 0){
                    ActiveSNParticle[i].DistanceMin = TemporalData.DistanceMin;
                    ActiveSNParticle[i].DistanceMinGlobalID = TemporalData.DistanceMinGlobalID;
                }
#endif // MAXIMUM_ENERGY_INPUT //}
#ifdef __CHECK_SUM__ //{
                ActiveSNParticle[i].CheckSum = TemporalData.CheckSum;
                // dprintlmpi(TemporalData.CheckSum);
#endif //__CHECK_SUM__ //}

#ifdef USE_STELLARFEEDBACK_RADIUS_LOCAL_UPDATE //{
                /// Insert Local Update Routine here.
                ActiveSNParticle[i].LocalUpdateFlags = false;
                int IsLocal = 0;
                for(int k=0;k<NProcs-1;k++){
                    //if(StellarFeedbackExportFlags[i][k]&BitMask){
                    if(StellarFeedbackExportFlags[Offset+k]&BitMask){
                        IsLocal ++;
                    }
                }
                if(IsLocal == 0){
                    ActiveSNParticle[i].LocalUpdateFlags = true;
                    //UpdateStellarFeedbackRadiusLocal(i,ActiveSNParticle,Neighbors,
                            //MyID,NProcs,StellarFeedbackExportFlags);
                    UpdateStellarFeedbackRadiusLocalModified(i,NActives,ActiveSNParticle,Neighbors,
                            MyID,NProcs,StellarFeedbackExportFlags);
                }
#endif // USE_KERNEL_LOCAL_UPDATE //}
            }
        }

        MPI_Waitall(counter_send,mpi_request_Export_Send,mpi_status_Export_Send);
        MPI_Waitall(counter_recv,mpi_request_Export_Recv,mpi_status_Export_Recv);

        for(int i=0;i<NImportAll;i++){ // For imported data.
            struct StructStellarFeedbackLocalInfo TemporalData = 
                RetrunStellarFeedbackLocalInfo(StellarFeedbackExportRecv[i].Pos,StellarFeedbackExportRecv[i].Radius);

            StellarFeedbackImportSend[i].Nlist = TemporalData.Nlist;
            StellarFeedbackImportSend[i].Density = TemporalData.Density;
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
            StellarFeedbackImportSend[i].SmoothedNumber = TemporalData.SmoothedNumber;
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
#ifdef SET_SNII_TEMPERATURE //{
            StellarFeedbackImportSend[i].GasMass = TemporalData.GasMass;
#endif //SET_SNII_TEMPERATURE //}
#ifdef MAXIMUM_ENERGY_INPUT //{
            StellarFeedbackImportSend[i].DistanceMin = TemporalData.DistanceMin;
            StellarFeedbackImportSend[i].DistanceMinGlobalID = TemporalData.DistanceMinGlobalID;
#endif // MAXIMUM_ENERGY_INPUT //}
#ifdef __CHECK_SUM__ //{
            StellarFeedbackImportSend[i].CheckSum = TemporalData.CheckSum;
#endif //__CHECK_SUM__ //}
            StellarFeedbackImportSend[i].Leaf = StellarFeedbackExportRecv[i].Leaf;
        }

        NImportAll = 0;
        int NImportAllNew = 0;
        for(int i=0;i<NProcs-1;i++){
            NImportThisTimeNew[i] = 0;
            for(int k=0;k<NImportThisTime[i];k++){
                if(StellarFeedbackImportSend[NImportAll].Nlist > 0){
                    StellarFeedbackImportSend[NImportAllNew] = StellarFeedbackImportSend[NImportAll];
                    NImportThisTimeNew[i] ++;
                    NImportAllNew ++;
                }
                NImportAll ++;
            }
        }

        int NImportThisTimeNew2[NProcs];
        int NExportThisTimeNew2[NProcs];
        NImportThisTimeNew2[MPIGetMyID()] = 0;
        for(int i=0;i<NProcs-1;i++){
            NImportThisTimeNew2[CommunicationTable[i].SendRank] = NImportThisTimeNew[i];
        }
        MPI_Alltoall(NImportThisTimeNew2,1,MPI_INT,NExportThisTimeNew2,1,MPI_INT,MPI_COMM_WORLD);
        for(int i=0;i<NProcs-1;i++){
            NExportThisTimeNew[i] = NExportThisTimeNew2[CommunicationTable[i].RecvRank];
        }

        NImport = 0;
        counter_send = counter_recv = 0;
        for(int i=0;i<NProcs-1;i++){
            if(NImportThisTimeNew[i]>0){
                MPI_Isend(StellarFeedbackImportSend+NImport,
                    NImportThisTimeNew[i]*sizeof(struct StructStellarFeedbackImport),
                        MPI_BYTE,CommunicationTable[i].SendRank,TAG_STELLARFEEDBACK_IMPORT+i,
                            MPI_COMM_WORLD,mpi_request_Export_Send+counter_send);
                MPI_Test(mpi_request_Export_Send+counter_send,&SendFlag,MPI_STATUS_IGNORE);
                counter_send ++;
            }
            if(NExportThisTimeNew[i]>0){
                MPI_Irecv(StellarFeedbackImportRecv[i],
                    NExportThisTimeNew[i]*sizeof(struct StructStellarFeedbackImport),
                        MPI_BYTE,CommunicationTable[i].RecvRank,TAG_STELLARFEEDBACK_IMPORT+i,
                            MPI_COMM_WORLD,mpi_request_Export_Recv+counter_recv);
                MPI_Test(mpi_request_Export_Recv+counter_recv,&RecvFlag,MPI_STATUS_IGNORE);
                counter_recv ++;
            }
            NImport += NImportThisTimeNew[i];
        }
        MPI_Waitall(counter_send,mpi_request_Export_Send,mpi_status_Export_Send);
        MPI_Waitall(counter_recv,mpi_request_Export_Recv,mpi_status_Export_Recv);

        for(int i=0;i<NProcs-1;i++){
            for(int k=0;k<NExportThisTimeNew[i];k++){ 
                int leaf = StellarFeedbackImportRecv[i][k].Leaf;
                ActiveSNParticle[leaf].Nlist += StellarFeedbackImportRecv[i][k].Nlist;
                ActiveSNParticle[leaf].Density += StellarFeedbackImportRecv[i][k].Density;
                
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
                ActiveSNParticle[leaf].SmoothedNumber += StellarFeedbackImportRecv[i][k].SmoothedNumber;
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
#ifdef SET_SNII_TEMPERATURE //{
                ActiveSNParticle[leaf].GasMass += StellarFeedbackImportRecv[i][k].GasMass;
#endif //SET_SNII_TEMPERATURE //}
#ifdef MAXIMUM_ENERGY_INPUT //{
                if(ActiveSNParticle[leaf].DistanceMin > StellarFeedbackImportRecv[i][k].DistanceMin){
                    ActiveSNParticle[leaf].DistanceMin = StellarFeedbackImportRecv[i][k].DistanceMin;
                    ActiveSNParticle[leaf].DistanceMinGlobalID = StellarFeedbackImportRecv[i][k].DistanceMinGlobalID;
                }
#endif // MAXIMUM_ENERGY_INPUT //}
#ifdef __CHECK_SUM__ //{
                ActiveSNParticle[leaf].CheckSum += StellarFeedbackImportRecv[i][k].CheckSum;
#endif //__CHECK_SUM__ //}
            }
        }

        // if(Niteration > MaxIterationTimes)
            // break;
        assert(Niteration < 1000);
        if(Niteration > MaxIterationTimes)
            OverMaxIterationTimes = true;
        ResetKernelSize(NActives,Niteration,NProcs,IndexList,StellarFeedbackExportFlags);

        int NLocalActiveLeaves = CheckNeighborNumberAndUpdateFeedbackRadiusModified(NActives,
                NProcs,StellarFeedbackExportFlags,ActiveSNParticle);
        MPI_Allreduce(&NLocalActiveLeaves,&NActiveLeaves,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

        Niteration ++;
    } while (0<NActiveLeaves);

#ifdef EVALUATE_KERNEL_BY_ITERATION //{
#ifdef PRINT_LOG_KERNEL_ITERATION //{
    if(MPIGetMyID()==MPI_ROOT_RANK)
        fprintf(stderr,"\n");
#else // PRINT_LOG_KERNEL_ITERATION 
    if(MPIGetMyID()==MPI_ROOT_RANK)
        fprintf(stderr,"%d iterations for kernel determination.\n",Niteration);
#endif // PRINT_LOG_KERNEL_ITERATION //}
#endif // EVALUATE_KERNEL_BY_ITERATION //}

    return;
}


struct StructEnergyHeavyElements{
    double Pos[3];
    double Radius;
    int    Type;
    // bool TypeII;
    // bool TypeIa;
    double Density;    // The mass weighted nomalization factor.
#ifdef SET_SNII_TEMPERATURE
    double GasMass;// Local gas mass.
#endif //SET_SNII_TEMPERATURE
#ifdef MAXIMUM_ENERGY_INPUT
    double DistanceMin;
    unsigned long int DistanceMinGlobalID;
#endif
    /*
    union {
        struct CELibStructSNIIFeedbackOutput SNII; // Energy, Ejecta, Remnant masses, and Elements.
        struct CELibStructSNIaFeedbackOutput SNIa; // Energy, Ejecta, Remnant masses, and Elements.
        struct CELibStructAGBFeedbackOutput  AGB;  // Reseased elements during AGB phase.
        struct CELibStructNSMFeedbackOutput  NSM;  // Reseased elements during NSM.
        struct CELibStructECSNFeedbackOutput  ECSN;  // Reseased elements during ECSN.
    } CELibData;
    */
    struct CELibStructFeedbackOutput  CELibData;
#ifdef  __CHECK_WEIGHT__
    int Leaf;
    int GlobalID;
#endif // __CHECK_WEIGHT__
};// *EnergyHeavyElements;

static inline void __attribute__((always_inline)) IncrementElements(double DestElements[restrict], double SrcElements[restrict], const double Weight){

    DestElements[CELibYield_H]  += SrcElements[CELibYield_H] *Weight;
    DestElements[CELibYield_He] += SrcElements[CELibYield_He]*Weight;
    DestElements[CELibYield_C]  += SrcElements[CELibYield_C] *Weight;
    DestElements[CELibYield_N]  += SrcElements[CELibYield_N] *Weight;
    DestElements[CELibYield_O]  += SrcElements[CELibYield_O] *Weight;
    DestElements[CELibYield_Ne] += SrcElements[CELibYield_Ne]*Weight;
    DestElements[CELibYield_Na] += SrcElements[CELibYield_Na]*Weight;
    DestElements[CELibYield_Mg] += SrcElements[CELibYield_Mg]*Weight;
    DestElements[CELibYield_Si] += SrcElements[CELibYield_Si]*Weight;
    DestElements[CELibYield_S]  += SrcElements[CELibYield_S] *Weight;
    DestElements[CELibYield_Ca] += SrcElements[CELibYield_Ca]*Weight;
    DestElements[CELibYield_Fe] += SrcElements[CELibYield_Fe]*Weight;
    DestElements[CELibYield_Ni] += SrcElements[CELibYield_Ni]*Weight;
    DestElements[CELibYield_Zn] += SrcElements[CELibYield_Zn]*Weight;
    DestElements[CELibYield_Sr] += SrcElements[CELibYield_Sr]*Weight;
    DestElements[CELibYield_Eu] += SrcElements[CELibYield_Eu]*Weight;

    return;
}


/*
 * This function distributes the released energy and heavy elements to the
 * surrouding ISM.
 * The definition of the weight is as follows:
 *     Weight = Mass*w/Density.
 * This weight function guarantees that \sum Weight = 1, except a rounding
 * error.
 */
static void DistributeEnergyHeavyElements(const int NExplosion, struct StructEnergyHeavyElements EnergyHeavyElements[restrict]){

    int Neighbors[MaxNeighborSize];
    for(int i=0;i<NExplosion;i++){
        double InvRadiusi = 1.e0/EnergyHeavyElements[i].Radius;
        double InvDensityi = 1.e0/EnergyHeavyElements[i].Density;
        int Nlist = GetNeighborsLimited(EnergyHeavyElements[i].Pos,2.e0*EnergyHeavyElements[i].Radius,Neighbors);

        for(int k=0;k<Nlist;k++){
            int leaf = Neighbors[k];
            double xij[3],Posj[3];

            Posj[0] = PhydroPosP(leaf)[0];
            Posj[1] = PhydroPosP(leaf)[1];
            Posj[2] = PhydroPosP(leaf)[2];
#ifdef PERIODIC_RUN 
            xij[0] = PeriodicDistance(EnergyHeavyElements[i].Pos[0],Posj[0],0);
            xij[1] = PeriodicDistance(EnergyHeavyElements[i].Pos[1],Posj[1],1);
            xij[2] = PeriodicDistance(EnergyHeavyElements[i].Pos[2],Posj[2],2);
#else // PERIODIC_RUN 
            xij[0] = EnergyHeavyElements[i].Pos[0]-Posj[0];
            xij[1] = EnergyHeavyElements[i].Pos[1]-Posj[1];
            xij[2] = EnergyHeavyElements[i].Pos[2]-Posj[2];
#endif // PERIODIC_RUN

            double r = NORM(xij); 
            double w = KernelStellarFeedback(r,InvRadiusi);

            double Weight = PhydroMass(leaf)*w*InvDensityi;
            // double m = PhydroMass(leaf);
#if 0
            if(Weight > 1.0){
#ifdef __CHECK_WEIGHT__
                fprintf(stderr,"Node[%d] Leaf = %d, GID = %ld |%d| %g %s:%d\n",MPIGetMyID(),EnergyHeavyElements[i].Leaf,
                        EnergyHeavyElements[i].GlobalID,Test,EnergyHeavyElements[i].Density,
                        __FUNCTION__,__LINE__);
                fflush(NULL);
#endif // __CHECK_WEIGHT__
                fprintf(stderr,"--- %ld %d %g * %g * %g = %g\n",PhydroBody(leaf)->GlobalID,leaf,
                        PhydroMass(leaf),w,InvDensityi,PhydroMass(leaf)*w*InvDensityi);
                fflush(NULL);
                assert(Weight < 1.0);
            }
#endif

            if(EnergyHeavyElements[i].Type == CELibFeedbackType_SNII){
#ifdef MAXIMUM_ENERGY_INPUT //{
                if(PhydroBody(leaf)->GlobalID == EnergyHeavyElements[i].DistanceMinGlobalID)
                    Phydro[leaf]->DQheat += EnergyHeavyElements[i].CELibData.Energy;
#else 
                Phydro[leaf]->DQheat += EnergyHeavyElements[i].CELibData.Energy*Weight;
#endif // MAXIMUM_ENERGY_INPUT //}
                IncrementElements(Phydro[leaf]->Elements,EnergyHeavyElements[i].CELibData.Elements,Weight);
            } else if(EnergyHeavyElements[i].Type == CELibFeedbackType_SNIa){
#ifdef MAXIMUM_ENERGY_INPUT //{
                if(PhydroBody(leaf)->GlobalID == EnergyHeavyElements[i].DistanceMinGlobalID)
                    Phydro[leaf]->DQheat += EnergyHeavyElements[i].CELibData.Energy;
#else
                Phydro[leaf]->DQheat += EnergyHeavyElements[i].CELibData.Energy*Weight;
#endif // MAXIMUM_ENERGY_INPUT  //}
                IncrementElements(Phydro[leaf]->Elements,EnergyHeavyElements[i].CELibData.Elements,Weight);
            } else if(EnergyHeavyElements[i].Type == CELibFeedbackType_AGB){
                IncrementElements(Phydro[leaf]->Elements,EnergyHeavyElements[i].CELibData.Elements,Weight);
#ifdef USE_CELIB_NSM //{
            } else if(EnergyHeavyElements[i].Type == CELibFeedbackType_NSM){
                IncrementElements(Phydro[leaf]->Elements,EnergyHeavyElements[i].CELibData.Elements,Weight);
#endif // USE_CELIB_NSM //}
 
#ifdef USE_CELIB_ECSN //{
            } else if(EnergyHeavyElements[i].Type == CELibFeedbackType_ECSN){
#ifdef MAXIMUM_ENERGY_INPUT //{
                if(PhydroBody(leaf)->GlobalID == EnergyHeavyElements[i].DistanceMinGlobalID)
                    Phydro[leaf]->DQheat += EnergyHeavyElements[i].CELibData.Energy;
#else 
                Phydro[leaf]->DQheat += EnergyHeavyElements[i].CELibData.Energy*Weight;
#endif // MAXIMUM_ENERGY_INPUT //}

                IncrementElements(Phydro[leaf]->Elements,EnergyHeavyElements[i].CELibData.Elements,Weight);
#endif // USE_CELIB_ECSN //}
#ifdef USE_CELIB_HN //{
            } else if(EnergyHeavyElements[i].Type == CELibFeedbackType_HN){
#ifdef MAXIMUM_ENERGY_INPUT //{
                if(PhydroBody(leaf)->GlobalID == EnergyHeavyElements[i].DistanceMinGlobalID)
                    Phydro[leaf]->DQheat += EnergyHeavyElements[i].CELibData.Energy;
#else 
                Phydro[leaf]->DQheat += EnergyHeavyElements[i].CELibData.Energy*Weight;
#endif // MAXIMUM_ENERGY_INPUT //}

                IncrementElements(Phydro[leaf]->Elements,EnergyHeavyElements[i].CELibData.Elements,Weight);
#endif // USE_CELIB_HN //}

            }
            // turn on Hydro Update Flag
            HydroUpdateFlag[leaf] = true;
        }
    }

    return ;
}


static inline void __attribute__((always_inline)) ConvertSNIIEjectaintoUnitMass(double Elements[restrict]){

    for(int i=0;i<CELibYield_Number;i++){
        Elements[i] *= MSUN_CGS/Pall.UnitMass;
    }
    return ;
}


static inline void __attribute__((always_inline)) ConvertSNIaEjectaintoUnitMass(double Elements[restrict]){

    for(int i=0;i<CELibYield_Number;i++){
        Elements[i] *= MSUN_CGS/Pall.UnitMass;
    }
    return ;
}

#if 0 
static inline void __attribute__((always_inline)) ConvertSNIIEjectaUnits(struct CELibStructSNIIFeedbackOutput *SNIIFeedback){

    SNIIFeedback->Energy *= EnergyConvertToSimulationUnit;
    for(int i=0;i<CELibYield_Number;i++){
        SNIIFeedback->Elements[i] *= MassConvertToSimulationUnit;
    }
    SNIIFeedback->EjectaMass *= MassConvertToSimulationUnit;
    SNIIFeedback->RemnantMass *= MassConvertToSimulationUnit;

    return ;
}


static inline void __attribute__((always_inline)) ConvertSNIaEjectaUnits(struct CELibStructSNIaFeedbackOutput *SNIaFeedback){

    SNIaFeedback->Energy *= EnergyConvertToSimulationUnit;
    for(int i=0;i<CELibYield_Number;i++){
        SNIaFeedback->Elements[i] *= MassConvertToSimulationUnit;
    }
    SNIaFeedback->EjectaMass *= MassConvertToSimulationUnit;
    SNIaFeedback->RemnantMass *= MassConvertToSimulationUnit;

    return ;
}

static inline void __attribute__((always_inline)) ConvertAGBEjectaUnits(struct CELibStructAGBFeedbackOutput *AGBFeedback){

    AGBFeedback->Energy = 0.e0;
    for(int i=0;i<CELibYield_Number;i++){
        AGBFeedback->Elements[i] *= MassConvertToSimulationUnit;
    }
    AGBFeedback->EjectaMass *= MassConvertToSimulationUnit;
    AGBFeedback->RemnantMass *= MassConvertToSimulationUnit;

    return ;
}

static inline void __attribute__((always_inline)) ConvertNSMEjectaUnits(struct CELibStructNSMFeedbackOutput *NSMFeedback){

    NSMFeedback->Energy = 0.e0;
    for(int i=0;i<CELibYield_Number;i++){
        NSMFeedback->Elements[i] *= MassConvertToSimulationUnit;
    }
    NSMFeedback->EjectaMass *= MassConvertToSimulationUnit;
    NSMFeedback->RemnantMass *= MassConvertToSimulationUnit;

    return ;
}
#endif

static inline void __attribute__((always_inline)) ConvertEjectaUnits(struct CELibStructFeedbackOutput *CELibData, const int type){

    if((type == CELibFeedbackType_SNII)||(type == CELibFeedbackType_SNIa)||(type == CELibFeedbackType_ECSN)||(type == CELibFeedbackType_HN)){
        CELibData->Energy *= EnergyConvertToSimulationUnit;
    } else {
        CELibData->Energy = 0.e0;
    }
    for(int i=0;i<CELibYield_Number;i++){
        CELibData->Elements[i] *= MassConvertToSimulationUnit;
    }
    CELibData->EjectaMass *= MassConvertToSimulationUnit;
    CELibData->RemnantMass *= MassConvertToSimulationUnit;

    return ;
}



/*
 * This function calculates the amounts of the release energy and heavy elements
 * of star particles, and then it scatters them to the surrounding ISM.
 */
static void ReleaseEnergyHeavyElements(const int NExplosion, const int IndexList[restrict]){
    
    int NProcs = MPIGetNumProcs();
    MPI_Status  mpi_status;

    struct StructEnergyHeavyElements EnergyHeavyElements[NExplosion];

    // Allocate HydroUpdateFlag
    if(Pall.Nhydro > HydroUpdateFlagSize){
        HydroUpdateFlagSize = ForAngelsShare*MAX(Pall.Nhydro,NAdditionUnit);
        HydroUpdateFlag = realloc(HydroUpdateFlag,sizeof(bool)*HydroUpdateFlagSize); 
    }
    memset(HydroUpdateFlag,0,HydroUpdateFlagSize*sizeof(bool));

    // Calc Energy and heavy elements.
    int Counter[CELibFeedbackType_Number];
    for(int i=0;i<CELibFeedbackType_Number;i++)
        Counter[i] = 0;

    for(int i=0;i<NExplosion;i++){
        Counter[ActiveSNParticle[i].Type] ++;
    }

    int WCounter[CELibFeedbackType_Number];
    MPI_Allreduce(Counter,WCounter,CELibFeedbackType_Number,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

    for(int i=0;i<CELibFeedbackType_Number;i++)
        Counter[i] = 0;

    // Pack release energy/mass/elements.
    for(int i=0;i<NExplosion;i++){
        EnergyHeavyElements[i].Pos[0] = ActiveSNParticle[i].Pos[0];
        EnergyHeavyElements[i].Pos[1] = ActiveSNParticle[i].Pos[1];
        EnergyHeavyElements[i].Pos[2] = ActiveSNParticle[i].Pos[2];
        EnergyHeavyElements[i].Radius = ActiveSNParticle[i].Radius;
        EnergyHeavyElements[i].Density = ActiveSNParticle[i].Density;
#ifdef __CHECK_WEIGHT__
        EnergyHeavyElements[i].Leaf = ActiveSNParticle[i].Index;
        EnergyHeavyElements[i].GlobalID = PstarBody(ActiveSNParticle[i].Index)->GlobalID;
#endif // __CHECK_WEIGHT__
#ifdef SET_SNII_TEMPERATURE
        EnergyHeavyElements[i].GasMass = ActiveSNParticle[i].GasMass;
#endif //SET_SNII_TEMPERATURE
#ifdef MAXIMUM_ENERGY_INPUT
        EnergyHeavyElements[i].DistanceMin = ActiveSNParticle[i].DistanceMin;
        EnergyHeavyElements[i].DistanceMinGlobalID = ActiveSNParticle[i].DistanceMinGlobalID;
#endif // MAXIMUM_ENERGY_INPUT

        int leaf = IndexList[i];
        if(ActiveSNParticle[i].Type == CELibFeedbackType_SNII){ // type II
            EnergyHeavyElements[i].CELibData 
                = CELibGetFeedback((struct CELibStructFeedbackInput){
                        .Mass = ActiveSNParticle[i].InitialMass,
                        .MassConversionFactor = MassConvertToMsun,
                        .Metallicity = ActiveSNParticle[i].Metallicity,
                        .Elements = Pstar[leaf]->Elements,
                        },CELibFeedbackType_SNII);

            // SNIIFeedback->Energy *= EnergyConvertToSimulationUnit;
            // 10^51 erg = g*cm^2/s^2 -> convert simulation unit
            // g=M/UnitMass; cm=L/UnitLength, s=T/UnitTime^2
            // ->
            // 10^51 erg = 10^51 * M/UnitMass (L/unitLength)^2 (UnitTime/T)^2
            // = 10^51 * UnitT^2/UnitLength^2/UnitMass M*L^2/T^2

            ConvertEjectaUnits(&(EnergyHeavyElements[i].CELibData),CELibFeedbackType_SNII);
            // Count ejecta mass.
            PstarBody(leaf)->Mass -= EnergyHeavyElements[i].CELibData.EjectaMass;
            Pstar[leaf]->Mass -= EnergyHeavyElements[i].CELibData.EjectaMass;

	    if(PstarBody(leaf)->Mass < 0)
		PstarBody(leaf)->Mass = 0;
	    if(Pstar[leaf]->Mass < 0)
		Pstar[leaf]->Mass = 0;

#if defined(PRESERVE_SNII_EVENTRATE)
            double MassInMsun = Pstar[leaf]->InitialMass*Pall.UnitMass/MSUN_CGS;
            double Energy_in_erg = EnergyHeavyElements[i].CELibData.Energy/EnergyConvertToSimulationUnit;
//#error should be checked.
            double p = gsl_rng_uniform(RandomGenerator);
            fprintf(stderr," !!!!!!!!!!!!!!!!!! %g %g %g\n",Energy_in_erg,SNIIEnergy,
                    EnergyHeavyElements[i].CELibData.Energy);
            if(Energy_in_erg < SNIIEnergy){
                fprintf(stderr," !!!!!!!!!!!!!!!!!! %g %g p=%g, %g\n",Energy_in_erg,SNIIEnergy,p,
                        EnergyHeavyElements[i].CELibData.Energy);
                if(Energy_in_erg > p*SNIIEnergy){ // accept.
                    EnergyHeavyElements[i].CELibData.Energy = SNIIEnergy*EnergyConvertToSimulationUnit;
                } else {
                    EnergyHeavyElements[i].CELibData.Energy = 0.e0;
                }
                fprintf(stderr," !!!!!!!!!!!!!!!!!! %g %g p=%g, %g\n",Energy_in_erg,SNIIEnergy,p,
                        EnergyHeavyElements[i].CELibData.Energy);
            }

#if 0 // old ver
            if(Pstar[leaf]->TypeIIProb){
#error SNIINumberPerMass is not evaluated.
                double prob = SNIINumberPerMass*MassInMsun;
                if(prob >= 1.0){ // Do nothing!
                    EnergyHeavyElements[i].CELibData.Energy = SNIIEnergy*MassInMsun;
                } else { // Put energy of one SN.
                    EnergyHeavyElements[i].CELibData.Energy = SNIIEnergy*EnergyConvertToSimulationUnit;
                }
            } else {
                EnergyHeavyElements[i].CELibData.Energy = 0.e0;
            }
#endif // old ver
#elif defined(SET_SNII_TEMPERATURE) //} PRESERVE_SNII_EVENTRATE //{
            double MassInMsun = ActiveSNParticle[i].InitialMass*MassConvertToMsun;
            double Usn_in_sim_unit = EnergyConvertToSimulationUnit
                        *CELibGetSNIIIntegratedEnergy(ActiveSNParticle[i].Metallicity)*MassInMsun
                        /(EnergyHeavyElements[i].GasMass);
            double Tsn = Pall.ConvertUtoT*Usn_in_sim_unit;

#if SNII_PEAK_TEMPERATURE //{
            double prob = Tsn/(KernelPeak*SNII_TEMPERATURE);
#else //SNII_PEAK_TEMPERATURE //} //{
            double prob = Tsn/SNII_TEMPERATURE;
#endif //SNII_PEAK_TEMPERATURE //}
            /*
            gprintlmpi(Tsn);
            gprintlmpi(prob);
            if(isinf(prob)){
                fprintf(stderr,"--- %g\n",prob);
                fflush(NULL);
                assert(!isinf(prob));
            }
            */
#if 1
            if(prob >= 1.0){ // Do nothing.
                //EnergyHeavyElements[i].CELibData.Energy = EnergyConvertToSimulationUnit*SNIINumberPerMass*SNIIEnergy*MassInMsun;
            }else{
                if(prob > gsl_rng_uniform(RandomGenerator)){
#if SNII_PEAK_TEMPERATURE //{
                    EnergyHeavyElements[i].CELibData.Energy = 
                        Pall.ConvertTtoU*KernelPeak*SNII_TEMPERATURE*EnergyHeavyElements[i].GasMass;
#else //SNII_PEAK_TEMPERATURE //} //{
                    EnergyHeavyElements[i].CELibData.Energy = 
                        Pall.ConvertTtoU*SNII_TEMPERATURE*EnergyHeavyElements[i].GasMass;
#endif //SNII_PEAK_TEMPERATURE //}
                }else{
                    EnergyHeavyElements[i].CELibData.Energy = 0.e0;
                }
            }
#endif
            /*
            fprintf(stderr,"p = %g, Tsn = %g [k], E = %g\n",
                    prob,Tsn,EnergyHeavyElements[i].CELibData.Energy);
            fflush(NULL);
            */
#endif // SET_SNII_TEMPERATURE //}

            EnergyHeavyElements[i].Type = CELibFeedbackType_SNII;
            Counter[CELibFeedbackType_SNII] ++;
        } else if (ActiveSNParticle[i].Type == CELibFeedbackType_SNIa){ // type Ia
            EnergyHeavyElements[i].CELibData 
                = CELibGetFeedback((struct CELibStructFeedbackInput){
                        .Mass = ActiveSNParticle[i].InitialMass,
                        .MassConversionFactor = MassConvertToMsun,
                        .Metallicity = ActiveSNParticle[i].Metallicity,
                        },CELibFeedbackType_SNIa);
                    
            ConvertEjectaUnits(&(EnergyHeavyElements[i].CELibData),CELibFeedbackType_SNIa);
            PstarBody(leaf)->Mass -= EnergyHeavyElements[i].CELibData.EjectaMass;
            Pstar[leaf]->Mass -= EnergyHeavyElements[i].CELibData.EjectaMass;
	    if(PstarBody(leaf)->Mass < 0)
		PstarBody(leaf)->Mass = 0;
	    if(Pstar[leaf]->Mass < 0)
		Pstar[leaf]->Mass = 0;


            EnergyHeavyElements[i].Type = CELibFeedbackType_SNIa;
            Counter[CELibFeedbackType_SNIa] ++;
        } 
#ifdef USE_CELIB_AGB
        else if (ActiveSNParticle[i].Type == CELibFeedbackType_AGB){ // AGB
            EnergyHeavyElements[i].CELibData
                = CELibGetFeedback((struct CELibStructFeedbackInput){
                        .Mass = ActiveSNParticle[i].InitialMass,
                        .MassConversionFactor = MassConvertToMsun,
                        .Metallicity = ActiveSNParticle[i].Metallicity,
                        .Count = ActiveSNParticle[i].Count, 
                        .Elements = Pstar[leaf]->Elements,
                        },CELibFeedbackType_AGB);

            // Convert into the unit mass.
            ConvertEjectaUnits(&(EnergyHeavyElements[i].CELibData),CELibFeedbackType_AGB);

            // Count ejecta mass.
            PstarBody(leaf)->Mass -= EnergyHeavyElements[i].CELibData.EjectaMass;
            Pstar[leaf]->Mass -= EnergyHeavyElements[i].CELibData.EjectaMass;
	    if(PstarBody(leaf)->Mass < 0)
		PstarBody(leaf)->Mass = 0;
	    if(Pstar[leaf]->Mass < 0)
		Pstar[leaf]->Mass = 0;


            EnergyHeavyElements[i].Type = CELibFeedbackType_AGB;
            Counter[CELibFeedbackType_AGB] ++;
        }
#endif // USE_CELIB_AGB
#ifdef USE_CELIB_NSM
    else if (ActiveSNParticle[i].Type == CELibFeedbackType_NSM){ // NSM
        EnergyHeavyElements[i].CELibData
            = CELibGetFeedback((struct CELibStructFeedbackInput){
                    .Mass = ActiveSNParticle[i].InitialMass,
                    .MassConversionFactor = MassConvertToMsun,
                    },CELibFeedbackType_NSM);

        // Convert into the unit mass.
        ConvertEjectaUnits(&(EnergyHeavyElements[i].CELibData),CELibFeedbackType_NSM);

        // Count ejecta mass.
        PstarBody(leaf)->Mass -= EnergyHeavyElements[i].CELibData.EjectaMass;
        Pstar[leaf]->Mass -= EnergyHeavyElements[i].CELibData.EjectaMass;
	    if(PstarBody(leaf)->Mass < 0)
		PstarBody(leaf)->Mass = 0;
	    if(Pstar[leaf]->Mass < 0)
		Pstar[leaf]->Mass = 0;


        EnergyHeavyElements[i].Type = CELibFeedbackType_NSM;
        Counter[CELibFeedbackType_NSM] ++;
    }
#endif // USE_CELIB_NSM

#ifdef USE_CELIB_ECSN
    else if (ActiveSNParticle[i].Type == CELibFeedbackType_ECSN){ // ECSN
        EnergyHeavyElements[i].CELibData
            = CELibGetFeedback((struct CELibStructFeedbackInput){
                    .Mass = ActiveSNParticle[i].InitialMass,
                    .MassConversionFactor = MassConvertToMsun,
                    .Metallicity = ActiveSNParticle[i].Metallicity,
                    .Elements = Pstar[leaf]->Elements,
                    },CELibFeedbackType_ECSN);

        // Convert into the unit mass.
        ConvertEjectaUnits(&(EnergyHeavyElements[i].CELibData),CELibFeedbackType_ECSN);

        // Count ejecta mass.
        PstarBody(leaf)->Mass -= EnergyHeavyElements[i].CELibData.EjectaMass;
        Pstar[leaf]->Mass -= EnergyHeavyElements[i].CELibData.EjectaMass;
	    if(PstarBody(leaf)->Mass < 0)
		PstarBody(leaf)->Mass = 0;
	    if(Pstar[leaf]->Mass < 0)
		Pstar[leaf]->Mass = 0;


        EnergyHeavyElements[i].Type = CELibFeedbackType_ECSN;
        Counter[CELibFeedbackType_ECSN] ++;
    }
#endif // USE_CELIB_ECSN
#ifdef USE_CELIB_HN
    else if (ActiveSNParticle[i].Type == CELibFeedbackType_HN){ // HN
        EnergyHeavyElements[i].CELibData
            = CELibGetFeedback((struct CELibStructFeedbackInput){
                    .Mass = ActiveSNParticle[i].InitialMass,
                    .MassConversionFactor = MassConvertToMsun,
                    .Metallicity = ActiveSNParticle[i].Metallicity,
                    .Elements = Pstar[leaf]->Elements,
                    },CELibFeedbackType_HN);

        // Convert into the unit mass.
        ConvertEjectaUnits(&(EnergyHeavyElements[i].CELibData),CELibFeedbackType_HN);

        // Count ejecta mass.
        PstarBody(leaf)->Mass -= EnergyHeavyElements[i].CELibData.EjectaMass;
        Pstar[leaf]->Mass -= EnergyHeavyElements[i].CELibData.EjectaMass;
	    if(PstarBody(leaf)->Mass < 0)
		PstarBody(leaf)->Mass = 0;
	    if(Pstar[leaf]->Mass < 0)
		Pstar[leaf]->Mass = 0;


        EnergyHeavyElements[i].Type = CELibFeedbackType_HN;
        Counter[CELibFeedbackType_HN] ++;
    }
#endif // USE_CELIB_HN


    }


    // Export
    struct StructEnergyHeavyElements *EnergyHeavyElementsExportSend[NProcs-1];
    struct StructEnergyHeavyElements *EnergyHeavyElementsExportRecv = NULL;
    MPI_Status mpi_status_Export_Send[NProcs-1];
    MPI_Request mpi_request_Export_Send[NProcs-1];
    MPI_Status mpi_status_Export_Recv[NProcs-1];
    MPI_Request mpi_request_Export_Recv[NProcs-1];


    int NExportThisTime[NProcs-1];
    for(int i=0;i<NProcs-1;i++){
        NExportThisTime[i] = 0;
        for(int k=0;k<NExplosion;k++){
            int Offset = k*NProcs;
            if(StellarFeedbackExportFlags[Offset+i]){
                NExportThisTime[i] ++;
            }
        }

        CheckSizeofBufferExportSendIndex(NExportThisTime[i],sizeof(struct StructEnergyHeavyElements),i);
        EnergyHeavyElementsExportSend[i] = BufferExportSend[i];

        NExportThisTime[i] = 0; 
        for(int k=0;k<NExplosion;k++){
            int Offset = k*NProcs;
            if(StellarFeedbackExportFlags[Offset+i]){
                EnergyHeavyElementsExportSend[i][NExportThisTime[i]] = EnergyHeavyElements[k];
                NExportThisTime[i] ++;
            }
        }
    }

    int NImport = 0;
    int NImportThisTime[NProcs-1];
    int NImportThisTime2[NProcs];
    int NExportThisTime2[NProcs];
    NImportThisTime2[MPIGetMyID()] = 0;
    for(int i=0;i<NProcs-1;i++){
        NExportThisTime2[CommunicationTable[i].SendRank] = NExportThisTime[i];
    }
    MPI_Alltoall(NExportThisTime2,1,MPI_INT,NImportThisTime2,1,MPI_INT,MPI_COMM_WORLD);
    for(int i=0;i<NProcs-1;i++){
        NImportThisTime[i] = NImportThisTime2[CommunicationTable[i].RecvRank];
        NImport += NImportThisTime[i];
    }

    CheckSizeofBufferExportRecv(NImport,sizeof(struct StructEnergyHeavyElements));
    EnergyHeavyElementsExportRecv = BufferExportRecv;

    NImport = 0;
    int counter_send = 0;
    int counter_recv = 0;
    int SendFlag,RecvFlag;
    for(int i=0;i<NProcs-1;i++){
        if(NExportThisTime[i]>0){
            MPI_Isend(EnergyHeavyElementsExportSend[i],
                NExportThisTime[i]*sizeof(struct StructEnergyHeavyElements),
                MPI_BYTE,CommunicationTable[i].SendRank,TAG_STELLARFEEDBACK_EXPORT+i,
                    MPI_COMM_WORLD,mpi_request_Export_Send+counter_send);
            MPI_Test(mpi_request_Export_Send+counter_send,&SendFlag,MPI_STATUS_IGNORE);
            counter_send ++;
        }
        if(NImportThisTime[i]>0){
            MPI_Irecv(EnergyHeavyElementsExportRecv+NImport,
                NImportThisTime[i]*sizeof(struct StructEnergyHeavyElements),
                MPI_BYTE,CommunicationTable[i].RecvRank,TAG_STELLARFEEDBACK_EXPORT+i,
                    MPI_COMM_WORLD,mpi_request_Export_Recv+counter_recv);
            MPI_Test(mpi_request_Export_Recv+counter_recv,&RecvFlag,MPI_STATUS_IGNORE);
            counter_recv ++;
            NImport += NImportThisTime[i];
        }
    }

    // Test = 0;
    // Distribution for local ISM particles
    DistributeEnergyHeavyElements(NExplosion,EnergyHeavyElements);

    MPI_Waitall(counter_send,mpi_request_Export_Send,mpi_status_Export_Send);
    MPI_Waitall(counter_recv,mpi_request_Export_Recv,mpi_status_Export_Recv);

    // Test = 1;
    // Distribution for external ISM particles.
    DistributeEnergyHeavyElements(NImport,EnergyHeavyElementsExportRecv);

    return ;
}

static void StellarFeedbackEndProcedure(const int NExplosion, const int IndexList[restrict], const int TypeList[restrict]){

    // Turn on the TypeII flags.
    for(int i=0;i<NExplosion;i++){
        double NextEventTime;
        int leaf = IndexList[i];
        int Type; 
        double Metallicity = Pstar[leaf]->Z;
        double Mass_in_Msun = Pstar[leaf]->InitialMass*MassConvertToMsun;
        int Count;

        if((TypeList[i] == CELibFeedbackType_SNII)||(TypeList[i] == CELibFeedbackType_SNIa)){ // II/Ia
            if(TypeList[i] == CELibFeedbackType_SNII){
                Pstar[leaf]->TypeII = true;
            } else if (TypeList[i] == CELibFeedbackType_SNIa){
                Pstar[leaf]->TypeIa = false;
            }
            Pstar[leaf]->SNIaCount ++;
            Pstar[leaf]->EventTime = Pstar[leaf]->FormationTime
                    +CELibGetNextEventTime((struct CELibStructNextEventTimeInput){
                        .R = gsl_rng_uniform(RandomGenerator),
                        .InitialMass_in_Msun = Pstar[leaf]->InitialMass*Pall.UnitMass/MSUN_CGS,
                        .Metallicity = Pstar[leaf]->Z,
                        .Count = Pstar[leaf]->SNIaCount,
                        },CELibFeedbackType_SNIa)
                        *YEAR_CGS/Pall.UnitTime;
        } else if(TypeList[i] == CELibFeedbackType_AGB){ // AGB
            Pstar[leaf]->AGBCount ++;
            Pstar[leaf]->EventTimeAGB = Pstar[leaf]->FormationTime
                    +CELibGetNextEventTime((struct CELibStructNextEventTimeInput){
                        .R = gsl_rng_uniform(RandomGenerator),
                        .InitialMass_in_Msun = Pstar[leaf]->InitialMass*Pall.UnitMass/MSUN_CGS,
                        .Metallicity = Pstar[leaf]->Z,
                        .Count = Pstar[leaf]->AGBCount,
                        },CELibFeedbackType_AGB)
                        *YEAR_CGS/Pall.UnitTime;
#ifdef USE_CELIB_NSM //{
        } else if(TypeList[i] == CELibFeedbackType_NSM){ // NSM
            Pstar[leaf]->NSMCount ++;
            Pstar[leaf]->EventTimeNSM = Pstar[leaf]->FormationTime
                    +CELibGetNextEventTime((struct CELibStructNextEventTimeInput){
                        .R = gsl_rng_uniform(RandomGenerator),
                        .InitialMass_in_Msun = Pstar[leaf]->InitialMass*Pall.UnitMass/MSUN_CGS,
                        .Metallicity = Pstar[leaf]->Z,
                        .Count = Pstar[leaf]->NSMCount,
                        },CELibFeedbackType_NSM)
                        *YEAR_CGS/Pall.UnitTime;
#endif // USE_CELIB_NSM //}

#ifdef USE_CELIB_ECSN //{
        } else if(TypeList[i] == CELibFeedbackType_ECSN){ // ECSN
            Pstar[leaf]->ECSNCount ++;
            Pstar[leaf]->EventTimeECSN = 10.0*Pall.TEnd;
#endif // USE_CELIB_ECSN //}
#ifdef USE_CELIB_HN //{
        } else if(TypeList[i] == CELibFeedbackType_HN){ // HN
            Pstar[leaf]->HNCount ++;
            Pstar[leaf]->EventTimeHN = 10.0*Pall.TEnd;
#endif // USE_CELIB_HN //}

        }


#if 0
    //\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

        //if(TypeList[i] != CELib_Feedback_AGB){ // II/Ia
        if((TypeList[i] == CELibFeedbackType_SNII)||(TypeList[i] == CELibFeedbackType_SNIa)){ // II/Ia
            if(TypeList[i] == CELibFeedbackType_SNII){
                Pstar[leaf]->TypeII = true;
            } else if (TypeList[i] == CELibFeedbackType_SNIa){
                Pstar[leaf]->TypeIa = false;
            }
            Pstar[leaf]->SNIaCount ++;
            Type = CELibFeedbackType_SNIa;
            Count = Pstar[leaf]->SNIaCount;
            // fprintf(stderr," +SNIa[ID=%ld] Count = %d\n",PstarBody(leaf)->GlobalID,Count);
            // dprintlmpi(Count);
        } else if(TypeList[i] == CELibFeedbackType_AGB){ // AGB
            Pstar[leaf]->AGBCount ++;
            Type = CELibFeedbackType_AGB;
            Count = Pstar[leaf]->AGBCount;
            //dprintlmpi(Count);
            // fprintf(stderr," +AGB[ID=%ld] Count = %d\n",PstarBody(leaf)->GlobalID,Count);
#ifdef USE_CELIB_NSM //{
        } else if(TypeList[i] == CELibFeedbackType_NSM){ // NSM
            Pstar[leaf]->NSMCount ++;
            Type = CELibFeedbackType_NSM;
            Count = Pstar[leaf]->NSMCount;
#endif // USE_CELIB_NSM //}
        }

        // Next explosion time.
        NextEventTime = StellarFeedbackGetNextEventTime(Type,Metallicity,Mass_in_Msun,Count);

        if(NextEventTime > 0.e0){
            //if(TypeList[i] != CELib_Feedback_AGB){ // II/Ia
            if((TypeList[i] == CELib_Feedback_TypeII)||(TypeList[i] == CELib_Feedback_TypeIa)){ // II/Ia
                Pstar[leaf]->EventTime = Pstar[leaf]->FormationTime+NextEventTime*YEAR_CGS/Pall.UnitTime;
                // fprintf(stderr,"SNIa end procedure [%d: %ld] T = %g %g /%g/, Count %d\n",leaf,PstarBody(leaf)->GlobalID,
                        // Pstar[leaf]->EventTime,Pall.TCurrent,NextEventTime,Pstar[leaf]->SNIaCount);
            //} else { // AGB
            } else if(TypeList[i] == CELib_Feedback_AGB){ // AGB
                Pstar[leaf]->EventTimeAGB = Pstar[leaf]->FormationTime+NextEventTime*YEAR_CGS/Pall.UnitTime;
                // fprintf(stderr,"AGB end procedure [%d: %ld] T = %g %g /%g/, Count %d\n",leaf,PstarBody(leaf)->GlobalID,
                        // Pstar[leaf]->EventTimeAGB,Pall.TCurrent,NextEventTime,Pstar[leaf]->AGBCount);
#ifdef USE_CELIB_NSM //{
            } else if(TypeList[i] == CELib_Feedback_NSM){ // NSM
                Pstar[leaf]->EventTimeNSM = Pstar[leaf]->FormationTime+NextEventTime*YEAR_CGS/Pall.UnitTime;
                // fprintf(stderr,"NSM end procedure [%d: %ld] T = %g %g /%g/, Count %d\n",leaf,PstarBody(leaf)->GlobalID,
                // Pstar[leaf]->EventTimeNSM,Pall.TCurrent,NextEventTime,Pstar[leaf]->AGBCount);
#endif // USE_CELIB_NSM //}
            }
        } else {
            Pstar[leaf]->EventTime = -1;
        }
#endif
    }

    // Increment SPH particle mass.
    for(int i=0;i<Pall.Nhydro;i++){
        if(HydroUpdateFlag[i]){
            double Mass = 0.e0;
            for(int k=0;k<CELibYield_Number;k++){
                Mass += Phydro[i]->Elements[k];
            }
            Phydro[i]->Mass = PhydroBody(i)->Mass = Mass;

            double MassLightElements = 
                Phydro[i]->Elements[CELibYield_H]+Phydro[i]->Elements[CELibYield_He];
            Phydro[i]->Z = (Mass-MassLightElements)/Mass;
           // Phydro[i]->Z = Phydro[i]->Elements[CELibYield_Fe]/(0.09767910863914096*Mass);
        }
    }

    return ;
}

static void AllocateActiveSNParticle(const int NActives){

    int NProcs = MPIGetNumProcs();
    if(StellarFeedbackExportFlagsMaxAllocated < MAX(NActives,NAdditionUnit)){
        if(StellarFeedbackExportFlagsMaxAllocated > 0){
            free(StellarFeedbackExportFlags);
            free(ActiveSNParticle);
        }
        StellarFeedbackExportFlagsMaxAllocated = (int)(MAX(ForAngelsShare*NActives,NAdditionUnit));
        StellarFeedbackExportFlags = malloc(sizeof(bool)*StellarFeedbackExportFlagsMaxAllocated*NProcs);
        ActiveSNParticle = malloc(sizeof(struct StructActiveSNParticle)*StellarFeedbackExportFlagsMaxAllocated);
    }

    return ;
}


void StellarFeedback(void){

#ifdef EVALUATE_SIZES_ALL_TOGETHER //{

    int NExplosion = ReturnCalcSizeElementNumber(CS_TypeSN,false);
    if(ReturnCalcSizeElementNumber(CS_TypeSN,true) == 0){
        return ;
    }
    int List[NExplosion];
    int FeedBackType[NExplosion];

    AllocateActiveSNParticle(NExplosion);
    CalcSizeSetSNInfo(ActiveSNParticle);
    for(int i=0;i<NExplosion;i++){
        List[i] = ActiveSNParticle[i].Index;
        FeedBackType[i] = ActiveSNParticle[i].Type;
    }

#else // EVALUATE_SIZES_ALL_TOGETHER //}//{

    // pick up star particles which explode during this time-step.
    int NExplosion = 0;
    for(int i=0;i<Pall.Nstars;i++){
        if(PstarBody(i)->Active){
            if(CheckEventTime(i) != NONE){
                NExplosion ++;
            }
        }
    }
    // fflush(NULL);

    int GlobalNExplosion = 0;
    MPI_Allreduce(&NExplosion,&GlobalNExplosion,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

    if(GlobalNExplosion == 0){
        return ;
    }
    if(MPIGetMyID()==MPI_ROOT_RANK)
        dprintlmpi(GlobalNExplosion);

    int List[NExplosion];
    int FeedBackType[NExplosion];
    int counter = 0;
    for(int i=0;i<Pall.Nstars;i++){
        if(PstarBody(i)->Active){
            int CurrentFeedBackType = CheckEventTime(i);
            if(CurrentFeedBackType != NONE){
                List[counter] = i;
                FeedBackType[counter] = CurrentFeedBackType;
                counter ++;
            }
        }
    }

    // fflush(NULL);

    // Feedback radius.
    CalcFeedbackRadius(NExplosion,List,FeedBackType);

#endif // EVALUATE_SIZES_ALL_TOGETHER //}

    // Release energy and heavy elements to the surrouding ISM.
    ReleaseEnergyHeavyElements(NExplosion,List);

    // End procedures
    StellarFeedbackEndProcedure(NExplosion,List,FeedBackType);

    return ;
}

#ifdef TASK_TEST_STELLARFEEDBACK //{
void StellarFeedbackTestCalcFeedbackRadius(int NActives, const int IndexList[restrict], const int TypeList[restrict], double Radius[restrict], int Nlist[restrict], int CheckSum[restrict]){
    CalcFeedbackRadius(NActives,IndexList,TypeList);
    for(int i=0;i<NActives;i++){
        Radius[i] = ActiveSNParticle[i].Radius;
        Nlist[i] = ActiveSNParticle[i].Nlist;
#ifdef __CHECK_SUM__
        CheckSum[i] = ActiveSNParticle[i].CheckSum;
#endif // __CHECK_SUM__
        gprintlmpi(ActiveSNParticle[i].Density);
    }

    return ;
}

void StellarFeedbackTestReleaseEnergyHeavyElements(const int NExplosion, const int IndexList[restrict], const int TypeList[restrict]){
    ReleaseEnergyHeavyElements(NExplosion,IndexList);
    StellarFeedbackEndProcedure(NExplosion,IndexList,TypeList);
    return ;
}
#endif // TASK_TEST_STELLARFEEDBACK //}


#endif // USE_CELIB //}

