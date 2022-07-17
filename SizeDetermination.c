#include "config.h"
#include "PlantHydroTree.h"
#include "NeighborSearch.h"
#include "KernelFunctions.h"
#include "HIIregion.h"
#include "StellarFeedback.h"
#include "SizeDetermination.h"

/*! \file SizeDetermination.c
 * \brief The kernel sizes of SPH particles, the sizes of HII regions, and
 * feedback radii are determined in this routine. 
 */

#define UPDATE_SIZE_LOCAL
//#define ADD_PERTURBATION 
#define USE_DEBUG_MODE

#define RadiusFactInc   (1.14) // 1.5 ^ (1.3)
#define RadiusFactDec   (0.79) // 0.75 ^ (1.3)
#define RadiusFactInc_First   (3) // 1.5 ^ (1.3)
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER
static double SmoothedMassConversionFactor;
#endif //USE_SMOOTHED_NEIGHBOR_NUMBER

static int MaxIterationTimes = 20;
static bool OverMaxIterationTimes = false;
static double LocalKernelMax = 0.e0;


/// For HII region
#define __EXPORT_TIMESTEP__
#define __PHOTON_COUNT_BASE__
/// For HII region

/// For SN region
#define __CHECK_SUM__
#define __CHECK_WEIGHT__
/// For SN region


////////////////////////////////////////////////////////////////////////////////////////////////

struct StructCSExport{
    int       Type;
    double    Radius;  // Radius
    double    Pos[3];  // Position.
    int       Leaf;
    bool      ExtraIterationFlag;
#ifdef USE_DEBUG_MODE
     unsigned long int GlobalID;
#endif // USE_DEBUG_MODE

    // HII region
    // double    Radius;  // Radius.
#ifdef __EXPORT_TIMESTEP__
    int       k_star;
#endif // __EXPORT_TIMESTEP__
};

struct StructCSImport{
    int       Type;
    double    SmoothedNumber;    // Mass.
    // double    Rho;
    int       Nlist;   // Nlist.
    int       Leaf;
    bool      ExtraIterationFlag;
#ifdef USE_DEBUG_MODE
    unsigned long int GlobalID;
#endif // USE_DEBUG_MODE
    // HII region
#ifdef __PHOTON_COUNT_BASE__ //{
    double    PhotonCount;    // Total absorbed photon number per sec.
    double    PhotonCountDistanceMin;
#else // __PHOTON_COUNT_BASE__ //} //{
    double    Mass;    // Mass.
    double    MassDistanceMin;
#endif // __PHOTON_COUNT_BASE__ //}
#ifdef SET_SNII_TEMPERATURE //{
    double    GasMass;
#endif //SET_SNII_TEMPERATURE //}
    double    DistanceMin;
    double    Density;
#if VISCOSITY_TYPE == 1 //{
    double    B[3][3];
#endif // VISCOSITY_TYPE == 1 //}
#ifdef USE_SINK_PARTICLE //{
    double    PotentialMin;
    double    MassTotal;
    double    VCOM[3];
#endif // USE_SINK_PARTICLE //}
#ifdef __CHECK_SUM__ //{
    int       CheckSum;
#endif //__CHECK_SUM__ //}
};

static int CS_NContactedDomains;
static int *CS_ContactedDomainID;

static inline void CS_AllocateContactedDomainID(void){
    CS_ContactedDomainID = malloc(sizeof(int)*MPIGetNumProcs());
    return ;
}

static inline bool CS_CheckLocalExternalDomainsContacted(const int MyID, const int ExtID){

    for(int k=0;k<3;k++){
        if((EdgesForHydro[MyID].PosMax[k] < EdgesForHydro[ExtID].PosMin[k])||
           (EdgesForHydro[MyID].PosMin[k] > EdgesForHydro[ExtID].PosMax[k]))  return false;
    }
    return true;

}

/*
 * This function checkes that how many contacted domains around the local
 * domain. Comparing the local domain edge to the external domain, 
 */
static inline void CS_CheckContactedDomain(void){
    int NProcs = MPIGetNumProcs();
    int MyID = MPIGetMyID();
    CS_NContactedDomains = 0;
    for(int i=0;i<NProcs-1;i++){
        int NodeID = CommunicationTable[i].SendRank;
        if(CS_CheckLocalExternalDomainsContacted(MyID,NodeID)){
            CS_ContactedDomainID[CS_NContactedDomains] = i;
            CS_NContactedDomains ++;
        }
    }
    return ;
}


struct StructCSActiveHydroParticle{
    double Rho;
#ifdef USE_SINK_PARTICLE //{
    double PotentialMin;
    double MassTotal;
    double VCOM[3];
#endif // USE_SINK_PARTICLE //}
#if VISCOSITY_TYPE == 1 //{
    double B[3][3];
#endif // VISCOSITY_TYPE == 1 //}
    int CacheIndex;
    int ExtraIteration;
}; 

struct StructCSActiveHIIParticle{
#ifdef __PHOTON_COUNT_BASE__ //{
    double PhotonCount;
#else //__PHOTON_COUNT_BASE__ //} //{
    double Mass;
#endif // __PHOTON_COUNT_BASE__ //}
    double LyAphoton;  // The number of Lyman continum photons [s^{-1}].
    // double Radius;

    double DistanceMin;
#ifdef __PHOTON_COUNT_BASE__ //{
    double PhotonCountDistanceMin;
#else //__PHOTON_COUNT_BASE__ //} //{
    double MassDistanceMin;
#endif // __PHOTON_COUNT_BASE__ //}
    bool HIIRegion;
};

struct StructCSActiveSNParticle{
    double Density;    // The mass weighted nomalization factor.
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
    int Type; // 1=TypeII, 2=TypeIa, 3=AGB
    int Count;
    double InitialMass;
    double Metallicity;
    int IterationCount;
};

static int CSExportFlagsMaxAllocated = 0;
struct StructActiveParticle{
    int Type;
    int Index; // NBCache[Index].
    int Nlist; // Number of neighbors.
    bool LocalUpdateFlags; // Flag for the local update.
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
    double SmoothedNumber; // Smoothed mas..
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
    double Pos[3];
    double Radius; // 2xRadius is the support scale.
    double Rvalue;
    double Lvalue;
#ifdef USE_DEBUG_MODE //{
    unsigned long long int GlobalID;
#endif // USE_DEBUG_MODE //}
    union {
        struct StructCSActiveHydroParticle Hydro; 
        struct StructCSActiveHIIParticle   HII; 
        struct StructCSActiveSNParticle    SN; 
    } Body;
} *ActiveParticle;

static int HydroEntry(struct StructActiveParticle AP[], bool CalcSizeExportFlags[], const int NProcs){
    
    int NActives = 0;
    int RootNodeID = 0; 
    int NumberofLeaves = HydroNode[RootNodeID].NumberofLeaves;
    int header = HydroNode[RootNodeID].Leaves;
    for(int i=0;i<NumberofLeaves;i++){
        int leaf = header + i; 
        if(HydroRoot.Leaves[leaf] < 0) continue;
        if(NBCache[leaf].Active){
            int Index = NBCache[leaf].Leaf;
            CalcSizeExportFlags[NActives*NProcs+NProcs-1] = true;
            AP[NActives].Index = Index;
            AP[NActives].Nlist = Phydro[Index]->Nlist;
            AP[NActives].Type = CS_TypeHydro;
            AP[NActives].Radius = NBCache[leaf].Kernel;
            /*
            if(NBCache[leaf].Kernel == 0){
                fprintf(stderr,"Special alart! Kernel size = 0\n");
                dlprintlmpi(PhydroBody(Index)->GlobalID);
                dprintlmpi(leaf);
                gprintlmpi(AP[NActives].Radius);
                gprintlmpi(Phydro[Index]->Kernel);
                gprintlmpi(Phydro[Index]->KernelPred);
                fflush(NULL);
            }
            */
            AP[NActives].Rvalue = AP[NActives].Lvalue = 0.e0;
            AP[NActives].Pos[0] = NBCache[leaf].Pos[0];
            AP[NActives].Pos[1] = NBCache[leaf].Pos[1];
            AP[NActives].Pos[2] = NBCache[leaf].Pos[2];
            // AP[NActives].Pos[0] = PhydroBody(Index)->PosP[0];
            // AP[NActives].Pos[1] = PhydroBody(Index)->PosP[1];
            // AP[NActives].Pos[2] = PhydroBody(Index)->PosP[2];
            AP[NActives].Body.Hydro.CacheIndex = leaf;
#ifdef USE_DEBUG_MODE //{
            AP[NActives].GlobalID = PhydroBody(Index)->GlobalID;
#endif // USE_DEBUG_MODE //}
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER
            AP[NActives].Body.Hydro.SmoothedNumber = Phydro[NBCache[leaf].Leaf]->SmoothedNumber;
#endif //ifdef USE_SMOOTHED_NEIGHBOR_NUMBER
            NActives ++;
        }
    }
    return NActives;
}

static int HIIregionEntry(struct StructActiveParticle AP[], bool CalcSizeExportFlags[], const int NProcs){

    int NActiveYoungStars = 0;
    for(int i=0;i<Pall.Nstars;i++){
        if(Pstar[i]->HIIflag){
            int Offset = NActiveYoungStars*NProcs;
            CalcSizeExportFlags[Offset+NProcs-1] = true;

            AP[NActiveYoungStars].Index = i;
            AP[NActiveYoungStars].Nlist = 0;
            AP[NActiveYoungStars].Type = CS_TypeHII;
            AP[NActiveYoungStars].Rvalue = AP[NActiveYoungStars].Lvalue = 0.e0;
            AP[NActiveYoungStars].Pos[0] = PstarBody(i)->PosP[0];
            AP[NActiveYoungStars].Pos[1] = PstarBody(i)->PosP[1];
            AP[NActiveYoungStars].Pos[2] = PstarBody(i)->PosP[2];
            // Set the initial Radius.
            if(Pstar[i]->StromgrenRadius>0.e0){ // Reuse old one.
                AP[NActiveYoungStars].Radius = Pstar[i]->StromgrenRadius;
            } else{ // Set an initial guess.
                AP[NActiveYoungStars].Radius = 0.25*PstarBody(i)->Eps;
            }

            // fprintf(stderr,"-- %d %d %g %d:%s\n",NActiveYoungStars,i,AP[NActiveYoungStars].Radius,__LINE__,__FUNCTION__);
            // fflush(NULL);

#ifdef __PHOTON_COUNT_BASE__
            AP[NActiveYoungStars].Body.HII.PhotonCount = 0.e0;
#else //__PHOTON_COUNT_BASE__
            AP[NActiveYoungStars].Body.HII.Mass = 0.e0;
#endif // __PHOTON_COUNT_BASE__

#ifdef PRESERVE_SNII_EVENTRATE
#if 0
            double MassInMsun = Pstar[i]->InitialMass*Pall.UnitMass/MSUN_CGS;
            double prob = SNIINumber*MassInMsun;
            if(prob >= 1.0){
                AP[NActiveYoungStars].Body.HII.LyAphoton = (PstarBody(i)->Mass*Pall.UnitMass/MSUN_CGS)*
                    ReturnNumberofLyPhoton(i,0);
            } else {
                if(Pstar[i]->TypeIIProb){
                    AP[NActiveYoungStars].Body.HII.LyAphoton = (1.0/SNIINumber)*ReturnNumberofLyPhoton(i,0);
                } 
            }
#else
            AP[NActiveYoungStars].Body.HII.LyAphoton = (PstarBody(i)->Mass*Pall.UnitMass/MSUN_CGS)*
                ReturnNumberofLyPhoton(i,0);
#endif 
#else
            AP[NActiveYoungStars].Body.HII.LyAphoton = (PstarBody(i)->Mass*Pall.UnitMass/MSUN_CGS)*
                ReturnNumberofLyPhoton(i,0);
#endif
            AP[NActiveYoungStars].Body.HII.HIIRegion = true;

#ifdef USE_DEBUG_MODE //{
            AP[NActiveYoungStars].GlobalID = PstarBody(i)->GlobalID;
#endif // USE_DEBUG_MODE //}

            NActiveYoungStars ++;
        }
    }
    return NActiveYoungStars;
}

static inline int __attribute__((always_inline)) CheckEventTime(const int Index){

#ifdef USE_CELIB //{
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
#endif // USE_CELIB_NSM //}

#ifdef USE_CELIB_ECSN //{
    if(Pall.TCurrent > Pstar[Index]->EventTimeECSN && Pstar[Index]->ECSNCount == 0){
        return CELibFeedbackType_ECSN;
    }
#endif // USE_CELIB_ECSN //}
#ifdef USE_CELIB_HN //{
    if(Pall.TCurrent > Pstar[Index]->EventTimeHN && Pstar[Index]->HNCount == 0){
        return CELibFeedbackType_HN;
    }
#endif // USE_CELIB_HN //}




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
#endif // USE_CELIB //}
}

static int StellarFeedbackEntry(struct StructActiveParticle AP[], bool CalcSizeExportFlags[], const int NProcs){

    int NExplosoin = 0;
    for(int i=0;i<Pall.Nstars;i++){
        if(PstarBody(i)->Active){
            int CurrentFeedbackType = CheckEventTime(i);
            if(CurrentFeedbackType != NONE){
                int Offset = NExplosoin*NProcs;
                CalcSizeExportFlags[Offset+NProcs-1] = ON;
                AP[NExplosoin].Index = i;
                AP[NExplosoin].Nlist = 0;
                AP[NExplosoin].Type = CS_TypeSN;
                AP[NExplosoin].Rvalue = AP[NExplosoin].Lvalue = 0.e0;
                AP[NExplosoin].Pos[0] = PstarPos(i)[0];
                AP[NExplosoin].Pos[1] = PstarPos(i)[1]; 
                AP[NExplosoin].Pos[2] = PstarPos(i)[2]; 
#ifdef CELib //}
                AP[NExplosoin].Radius = StellarFeedbackRadiusInitialGuess(i);
#else // CELib //}//{
                AP[NExplosoin].Radius = 2.0*PstarBody(i)->Eps; 
#endif // CELib //}
                AP[NExplosoin].Body.SN.Type = CurrentFeedbackType;
#ifdef USE_CELIB //{
                if(CurrentFeedbackType == CELibFeedbackType_AGB){
                    AP[NExplosoin].Body.SN.Count = Pstar[i]->AGBCount;
#ifdef USE_CELIB_NSM //{
                } else if(CurrentFeedbackType == CELibFeedbackType_NSM){
                    AP[NExplosoin].Body.SN.Count = Pstar[i]->NSMCount;
#endif // USE_CELIB_NSM //}
#ifdef USE_CELIB_ECSN //{
                } else if(CurrentFeedbackType == CELibFeedbackType_ECSN){
                    AP[NExplosoin].Body.SN.Count = Pstar[i]->ECSNCount;
#endif // USE_CELIB_ECSN //}
#ifdef USE_CELIB_HN //{
                } else if(CurrentFeedbackType == CELibFeedbackType_HN){
                    AP[NExplosoin].Body.SN.Count = Pstar[i]->HNCount;
#endif // USE_CELIB_HN //}
                } else {                        
                    AP[NExplosoin].Body.SN.Count = Pstar[i]->SNIaCount;
                }
#endif // USE_CELIB //}
                AP[NExplosoin].Body.SN.InitialMass = Pstar[i]->InitialMass;
                AP[NExplosoin].Body.SN.Metallicity = Pstar[i]->Z;
                AP[NExplosoin].Body.SN.IterationCount = 0;
                AP[NExplosoin].Body.SN.Type = CurrentFeedbackType;
#ifdef USE_DEBUG_MODE //{
                AP[NExplosoin].GlobalID = PstarBody(i)->GlobalID;
#endif // USE_DEBUG_MODE //}

                NExplosoin ++;
            }
        }
    }

    return NExplosoin;
}


static inline bool __attribute__((always_inline)) CS_OverlapDomain(double Pos[restrict], const double h, const int NodeID){ 

    double Dist2 = 0.e0;
    for(int k=0;k<3;k++){
        if(Pos[k] < EdgesForHydro[NodeID].PosMin[k]) 
            Dist2 += SQ(EdgesForHydro[NodeID].PosMin[k]-Pos[k]);
        if(Pos[k] > EdgesForHydro[NodeID].PosMax[k])
            Dist2 += SQ(EdgesForHydro[NodeID].PosMax[k]-Pos[k]);
    }
    return (Dist2 < SQ(h));
}

/*
 * This function return number of particles which should export to the node ID of [Index].
 */
static inline int __attribute__((always_inline)) CS_CheckExportFlags(const int NodeIndex, const int NProcs, bool CSExportFlags[restrict], const int NActives, struct StructActiveParticle AP[restrict]){

    if((Pall.Nhydro+Pall.Nstars) == 0)
        return 0;

    int ExportNodeID = CommunicationTable[NodeIndex].SendRank;

    // Node By Node Comparison
    //double BoxCenter[] = {GravityNode[0].Pos[0],GravityNode[0].Pos[1],GravityNode[0].Pos[2]};
    double BoxCenter[] = {HydroNode[0].Pos[0],HydroNode[0].Pos[1],HydroNode[0].Pos[2]}; // need check
    if(!CS_OverlapDomain(BoxCenter,LocalKernelMax,ExportNodeID)){
        return 0;
    }

    int NExport = 0;
    for(int i=0;i<NActives;i++){
        int Offset = i*NProcs;
        if(CSExportFlags[Offset+NProcs-1]){
            if(CS_OverlapDomain(AP[i].Pos,2.0*AP[i].Radius,ExportNodeID)){
                CSExportFlags[Offset+NodeIndex] = true;
                NExport ++;
            }
        }
    }

    return NExport;
}

struct StructCSHydroLocalInfo {
    int Nlist;
#if VISCOSITY_TYPE == 1 //{
    double B[3][3];
#endif // VISCOSITY_TYPE == 1 //}
#ifdef USE_SINK_PARTICLE //{
    double PotentialMin;
    double MassTotal;
    double VCOM[3];
#endif // USE_SINK_PARTICLE //}
};

static struct StructCSHydroLocalInfo CS_GetNumberofNeighbors(double Pos[restrict], const double Kernel, int Neighbors[restrict]){

    struct StructCSHydroLocalInfo TempHydroLocalInfo = {
        .Nlist = 0,
#if VISCOSITY_TYPE == 1 //{
        .B[0][0] = 0.e0,
        .B[0][1] = 0.e0,
        .B[0][2] = 0.e0,
        .B[1][0] = 0.e0,
        .B[1][1] = 0.e0,
        .B[1][2] = 0.e0,
        .B[2][0] = 0.e0,
        .B[2][1] = 0.e0,
        .B[2][2] = 0.e0,
#endif // VISCOSITY_TYPE //}
#ifdef USE_SINK_PARTICLE //{
        .PotentialMin = 0.e0,
        .MassTotal = 0.e0,
        .VCOM[0] = 0.e0,
        .VCOM[1] = 0.e0,
        .VCOM[2] = 0.e0,
#endif // USE_SINK_PARTICLE //}
    };

    int Iteration = 0;
    int RootNodeID = 0;
    int CurrentNodeID = HydroNode[RootNodeID].Children;
    do {
        CurrentNodeID = GetNeighborsIterativeApproach(CurrentNodeID,Pos,2.e0*Kernel,&(TempHydroLocalInfo.Nlist),Neighbors);

#if VISCOSITY_TYPE == 1 //{
        double InvKerneli = 1.e0/Kernel;
        for( int i=0;i<TempHydroLocalInfo.Nlist;i++){
            int index = Neighbors[i];

            double xij[3],vij[3];
#ifdef PERIODIC_RUN //{
            xij[0] = PeriodicDistance(Pos[0],PhydroPosP(index)[0],0);
            xij[1] = PeriodicDistance(Pos[1],PhydroPosP(index)[1],1);
            xij[2] = PeriodicDistance(Pos[2],PhydroPosP(index)[2],2);
#else // PERIODIC_RUN //}//{
            xij[0] = Pos[0]-PhydroPosP(index)[0];
            xij[1] = Pos[1]-PhydroPosP(index)[1];
            xij[2] = Pos[2]-PhydroPosP(index)[2];
#endif // PERIODIC_RUN //}

            double r2 = NORM2(xij);
            double r = sqrt(r2); 
            double w = SPHKernel(r,InvKerneli);

#ifdef USE_DISPH //{
            double Volumej = Phydro[index]->Mass*Phydro[index]->UPred/Phydro[index]->EnergyDensityPred;
#elif defined(USE_SPSPH) // USE_DISPH//USE_SPSPH //}//{
            double Volumej = Phydro[index]->ZwPred/Phydro[index]->PseudoDensityPred;
#else // USE_SPSPH //}//{
            double Volumej = Phydro[index]->Mass/Phydro[index]->RhoPred;
#endif // USE_SPSPH //}
            double Factor = w*Volumej;
            TempHydroLocalInfo.B[0][0] += Factor*xij[0]*xij[0];
            TempHydroLocalInfo.B[0][1] += Factor*xij[0]*xij[1];
            TempHydroLocalInfo.B[0][2] += Factor*xij[0]*xij[2];
            TempHydroLocalInfo.B[1][0] += Factor*xij[1]*xij[0];
            TempHydroLocalInfo.B[1][1] += Factor*xij[1]*xij[1];
            TempHydroLocalInfo.B[1][2] += Factor*xij[1]*xij[2];
            TempHydroLocalInfo.B[2][0] += Factor*xij[2]*xij[0];
            TempHydroLocalInfo.B[2][1] += Factor*xij[2]*xij[1];
            TempHydroLocalInfo.B[2][2] += Factor*xij[2]*xij[2];
        }
#endif // VISCOSITY_TYPE == 1 //}

#ifdef USE_SINK_PARTICLE //{
        if((TempHydroLocalInfo.Nlist>0)&&(Iteration==0)){
            TempHydroLocalInfo.PotentialMin = PhydroBody(Neighbors[0])->Pot;
        }

        for(int i=0;i<TempHydroLocalInfo.Nlist;i++){
            int index = Neighbors[i];

            TempHydroLocalInfo.MassTotal += PhydroBody(index)->Mass;
            TempHydroLocalInfo.VCOM[0] += PhydroBody(index)->Mass*PhydroBody(index)->Vel[0];
            TempHydroLocalInfo.VCOM[1] += PhydroBody(index)->Mass*PhydroBody(index)->Vel[1];
            TempHydroLocalInfo.VCOM[2] += PhydroBody(index)->Mass*PhydroBody(index)->Vel[2];

            // Evaluate minimum distance.
            double CurrentPotentialMin = PhydroBody(index)->Pot;
            if(TempHydroLocalInfo.PotentialMin >= CurrentPotentialMin){
                TempHydroLocalInfo.PotentialMin = CurrentPotentialMin;
            }
        }
#endif // USE_SINK_PARTICLE //}
        Iteration ++;
    }while(CurrentNodeID != RootNodeID);

    return TempHydroLocalInfo;
}

#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
static int CS_GetSmoothedNumberofNeighbors(double Pos[restrict], const double Kernel, int Neighbors[restrict], double *SmoothedNumber){

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

struct StructCSHIILocalInfo {
    int Nlist;
#ifdef __PHOTON_COUNT_BASE__ //{
    double PhotonCount;
    double PhotonCountDistanceMin;
#else //__PHOTON_COUNT_BASE__ //} //{
    double Mass;
    double MassDistanceMin;
#endif // __PHOTON_COUNT_BASE__ //}
    double DistanceMin;
};

/*
 * This function counts (1) the total number of photons which can be absorbed in
 * the neighbor particles, or (2) the total mass of the neighbor particles.
 */
static struct StructCSHIILocalInfo CS_ReturnHIILocalInfo(double Pos[restrict], const double Radius, int Neighbors[restrict]){

    struct StructCSHIILocalInfo TempHIILocalInfo = {
        .Nlist = 0,
        .DistanceMin = 0.0,
        //.DistMin = 10000,
#ifdef __PHOTON_COUNT_BASE__ //{
        .PhotonCount = 0.e0,
        .PhotonCountDistanceMin = 0.0,
#else //__PHOTON_COUNT_BASE__ //} //{
        .Mass = 0.e0,
        .MassDistanceMin = 0.0,
#endif // __PHOTON_COUNT_BASE__ //}
    };

    int Iteration = 0;
    int RootNodeID = 0;
    int CurrentNodeID = HydroNode[RootNodeID].Children;
    do {
        CurrentNodeID = GetNeighborsIterativeApproach(CurrentNodeID,Pos,
                2.e0*Radius,&(TempHIILocalInfo.Nlist),Neighbors);


        if((TempHIILocalInfo.Nlist>0)&&(Iteration==0)){
            TempHIILocalInfo.DistanceMin = DISTANCE2(PhydroBody(Neighbors[0])->PosP,Pos);
        }

        for(int i=0;i<TempHIILocalInfo.Nlist;i++){
            int index = Neighbors[i];
            if(Phydro[index]->UPred*Pall.ConvertUtoT > 1.e+4) continue;

#ifdef __PHOTON_COUNT_BASE__ //{ // fiducial model.
            const double InvProtonMass = 1.0/PROTON_MASS_CGS;
            double X = (1.0-HeliumAbandance-Phydro[index]->Z);
            double n_H = (X*Phydro[index]->RhoPred*Pall.UnitMass*InvProtonMass)/CUBE(Pall.UnitLength);
            double n_e = n_H;
            double r3 = 3*PhydroBody(index)->Mass/(4*M_PI*Phydro[index]->RhoPred);
            double CurrentPhotonCount = (4.0*M_PI/3.0)*n_H*n_e*ReturnRecombinationConstantAlpha()*r3*CUBE(Pall.UnitLength);

            TempHIILocalInfo.PhotonCount += CurrentPhotonCount;
#else //__PHOTON_COUNT_BASE__ //} //{
            double CurrentMass = (1.0-HeliumAbandance-Phydro[index]->Z)*PhydroBody(index)->Mass;
            TempHIILocalInfo.Mass += CurrentMass;
#endif // __PHOTON_COUNT_BASE__ //}

            // Evaluate minimum distance.
            double TmpDistanceMin = DISTANCE2(PhydroBody(index)->PosP,Pos);
            if(TempHIILocalInfo.DistanceMin >= TmpDistanceMin){
                TempHIILocalInfo.DistanceMin = TmpDistanceMin;
#ifdef __PHOTON_COUNT_BASE__ //{
                TempHIILocalInfo.PhotonCountDistanceMin = CurrentPhotonCount;
                // gprintlmpi(CurrentPhotonCount);
#else //__PHOTON_COUNT_BASE__ //} //{
                TempHIILocalInfo.MassDistanceMin = CurrentMass;
#endif // __PHOTON_COUNT_BASE__ //}
            }
        }
        Iteration ++;
    }while(CurrentNodeID != RootNodeID);

    return TempHIILocalInfo;
}

struct StructCSStellarFeedbackLocalInfo{
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


/*
 * This function returns a structure which involves various values (the neighbor
 * number, the smootehd neighbor number, the density, the local gas mass, the
 * minimum distance, and the globalID of the closest particle).
 */
struct StructCSStellarFeedbackLocalInfo CS_RetrunStellarFeedbackLocalInfo(double Pos[restrict], const double Radius, int Neighbors[restrict]){

    //static int Neighbors[MaxNeighborSize];
    struct StructCSStellarFeedbackLocalInfo TempStellarFeedbackLocalInfo = {
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
        //double w = KernelStellarFeedback(r,InvRadiusi);
        double w = SPHKernel(r,InvRadiusi);

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


static inline void __attribute__((always_inline)) CS_OverwriteNeighborInfo(
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER  //{
        const double Nlist, 
#else 
        const int Nlist, 
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER  //}
        const int leaf){

    int Index = NBCache[leaf].Leaf;
    Phydro[Index]->Nlist = Nlist;
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER
    Phydro[Index]->SmoothedNumber = Nlist;
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER
    Phydro[Index]->Kernel = 
    Phydro[Index]->KernelPred = NBCache[leaf].Kernel;

    return ;
}

static inline bool __attribute__((always_inline)) CS_CheckNeighborNumberAndUpdateHydroRadius_i(const int NProcs, bool CSExportFlags_i[restrict], struct StructActiveParticle *AP_i){ 
    // Index can be removed.

    int NBmin = Pall.Ns-Pall.Npm;
    int NBmax = Pall.Ns+Pall.Npm;
#ifdef USE_MAXIMUM_KERNEL_SIZE
#ifdef MAXIMUM_KERNEL_SIZE
    double MaximumKernelSize = Pall.AdaptiveSofteningFactor*MAXIMUM_KERNEL_SIZE*KPC_CGS/Pall.UnitLength;
#else
#error Set MAXIMUM_KERNEL_SIZE
#endif
#endif

#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
    double Nlist = AP_i->SmoothedNumber*
        SmoothedMassConversionFactor*CUBE(AP_i->Radius);
#else // USE_SMOOTHED_NEIGHBOR_NUMBER
    int Nlist = AP_i->Nlist;
#endif  // USE_SMOOTHED_NEIGHBOR_NUMBER //}

    if(((NBmin)<=Nlist)&&(Nlist<=(NBmax))){
        CSExportFlags_i[NProcs-1] = false;
        // OverwriteNeighborInfo(Nlist,leaf);
        return true;
    }else if((Nlist>=(0.8*NBmin))&&(Nlist<=(2*NBmax))
            &&(AP_i->Rvalue>0.e0)&&(AP_i->Lvalue>0.e0)){
        if(AP_i->Rvalue-AP_i->Lvalue < 1.e-6*AP_i->Lvalue){
            CSExportFlags_i[NProcs-1] = false;
            // OverwriteNeighborInfo(Nlist,leaf);
            return true;
        }
    }
#ifdef USE_MINIMUM_KERNEL_SIZE
    else if (((NBmax)<Nlist)&&(AP_i->Radius<=0.5*PhydroBody(AP_i->Index)->Eps*Pall.AdaptiveSofteningFactor)){
        AP_i->Radius = 0.5*PhydroBody(AP_i->Index)->Eps*Pall.AdaptiveSofteningFactor;
        CSExportFlags_i[NProcs-1] = false;
        // OverwriteNeighborInfo(Nlist,leaf);
        return true;
    }
#endif
#ifdef USE_MAXIMUM_KERNEL_SIZE
    else if (((NBmin)>Nlist)&&(AP_i->Radius>MaximumKernelSize)){
        AP_i->Radius = MaximumKernelSize;
        CSExportFlags_i[NProcs-1] = false;
        //OverwriteNeighborInfo(Nlist,leaf);
        return true;
    }
#endif
    if(CSExportFlags_i[NProcs-1]){
        if(Nlist<NBmin){
            AP_i->Lvalue = fmax(AP_i->Lvalue,AP_i->Radius);
        } else if(Nlist>NBmax){
            if(AP_i->Rvalue > 0.e0){
                AP_i->Rvalue = fmin(AP_i->Rvalue,AP_i->Radius);
            }else{
                AP_i->Rvalue = AP_i->Radius;
            }
        }

        if((AP_i->Lvalue>0.e0)&&(AP_i->Rvalue>0.e0)){
            AP_i->Radius = cbrt(0.5*(CUBE(AP_i->Lvalue)+CUBE(AP_i->Rvalue)));
        }else{
            if((AP_i->Rvalue == 0.e0)&&(AP_i->Lvalue > 0.e0)){
                AP_i->Radius *= RadiusFactInc;
            }else if((AP_i->Rvalue > 0.e0)&&(AP_i->Lvalue == 0.e0)){
                AP_i->Radius *= RadiusFactDec;
            }
        }
    }
    
    return false;
}

static inline bool __attribute__((always_inline)) CS_CheckLocalMassAndUpdateHIIRadius_i(const int NProcs, bool CSExportFlags_i[NProcs], struct StructActiveParticle *AP_i){ 

#define ConvergenceFactor  (1.e-2)

#ifdef __PHOTON_COUNT_BASE__
    double Qest = AP_i->Body.HII.PhotonCount;
#else // __PHOTON_COUNT_BASE__
    double fact = (3.0/(4.0*M_PI))*(ReturnAlpha()/(SQ(PROTON_MASS_CGS)));
    double Qest = fact*(SQ((AP_i->Body.HII.Mass)*Pall.UnitMass)/(CUBE((2.0*(AP_i->Radius))*Pall.UnitLength)));
#endif // __PHOTON_COUNT_BASE__
    double Qmin = (1.0-ConvergenceFactor)*AP_i->Body.HII.LyAphoton;
    double Qmax = (1.0+ConvergenceFactor)*AP_i->Body.HII.LyAphoton;


    if(((Qmin)<=Qest)&&(Qest<=(Qmax))){
        CSExportFlags_i[NProcs-1] = false;
#ifdef __PHOTON_COUNT_BASE__
        AP_i->Body.HII.PhotonCount = 1.0;
#else // __PHOTON_COUNT_BASE__
        AP_i->Body.HII.Mass = 1.0;
#endif // __PHOTON_COUNT_BASE__
        return true;
    }else if((Qest<=(Qmax))&&(AP_i->Rvalue>0.e0)&&(AP_i->Lvalue>0.e0)){
        if(AP_i->Rvalue-AP_i->Lvalue < 1.e-3*AP_i->Lvalue){
            CSExportFlags_i[NProcs-1] = false;
#ifdef __PHOTON_COUNT_BASE__ //{
            AP_i->Body.HII.PhotonCount = 2.0;
#else // __PHOTON_COUNT_BASE__ // { //}
            AP_i->Body.HII.Mass = 2.0;
#endif // __PHOTON_COUNT_BASE__ //}
            return true;
        }
    } 
#ifdef __PHOTON_COUNT_BASE__ //{
        else if((AP_i->Body.HII.PhotonCountDistanceMin > AP_i->Body.HII.LyAphoton))
#else // __PHOTON_COUNT_BASE__ // { //}
#error USE PHOTON COUNT BASE MODE
#endif // __PHOTON_COUNT_BASE__ //}
        {
            CSExportFlags_i[NProcs-1] = false;
#ifdef __PHOTON_COUNT_BASE__ //{
            AP_i->Body.HII.PhotonCount = 4.0;
#else // __PHOTON_COUNT_BASE__ // { //}
            AP_i->Body.HII.Mass = 4.0;
#endif // __PHOTON_COUNT_BASE__ //}
            // dprintlmpi(1);
            // gprintlmpi(AP_i->Body.HII.PhotonCountDistanceMin);
            // gprintlmpi(AP_i->Body.HII.LyAphoton);
            AP_i->Body.HII.HIIRegion = false;
            return true;
    }


    if(CSExportFlags_i[NProcs-1]){
        if(Qest<Qmin){
            AP_i->Lvalue = fmax(AP_i->Lvalue,AP_i->Radius);
        } else if(Qest>Qmax){
            if((AP_i->Nlist) == 1){
                AP_i->Radius *= cbrt(Qmax/Qest);
                CSExportFlags_i[NProcs-1] = false;
#ifdef __PHOTON_COUNT_BASE__
                AP_i->Body.HII.PhotonCount = 3.0;
#else // __PHOTON_COUNT_BASE__
                AP_i->Body.HII.Mass = 3.0;
#endif // __PHOTON_COUNT_BASE__
                return true;
            }else{
                if(AP_i->Rvalue > 0.e0){
                    AP_i->Rvalue = fmin(AP_i->Rvalue,AP_i->Radius);
                }else{
                    AP_i->Rvalue = AP_i->Radius;
                }
            }
        }

        if((AP_i->Lvalue>0.e0)&&(AP_i->Rvalue>0.e0)){
            AP_i->Radius = cbrt(0.5*(CUBE(AP_i->Lvalue)+CUBE(AP_i->Rvalue)));
        }else{
            if((AP_i->Rvalue == 0.e0)&&(AP_i->Lvalue > 0.e0)){
                AP_i->Radius *= RadiusFactInc_First;
            }else if((AP_i->Rvalue > 0.e0)&&(AP_i->Lvalue == 0.e0)){
                AP_i->Radius *= RadiusFactDec;
            }
        }
    }
    return false;
}

//#define StellarFeedbackRadiusFactInc   (1.14) // 1.5 ^ (1.3)
#define StellarFeedbackRadiusFactInc   (1.2596) //  2^(0.333)
#define StellarFeedbackRadiusFactDec   (0.79) // 0.75 ^ (1.3)
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER
static double SmoothedMassConversionFactor;
#endif //USE_SMOOTHED_NEIGHBOR_NUMBER
/*
 * This function checks the convergence of the feedback radius, using bisetion
 * method. If the radius is convergered, this function returns a "true" flag.
 * If not, it returns a "false" flag.
 */
static inline bool __attribute__((always_inline)) CS_CheckNeighborNumberAndUpdateFeedbackRadius_i(const int NProcs, bool CSExportFlags_i[restrict], struct StructActiveParticle *AP_i){ 

#ifdef USE_SN_INPUT_PARTICLE_NUMBER //{
    int NBmin = SN_INPUT_PARTICLE_NUMBER-SN_INPUT_PARTICLE_NUMBER_MARGIN;
    int NBmax = SN_INPUT_PARTICLE_NUMBER+SN_INPUT_PARTICLE_NUMBER_MARGIN;
#else // USE_SN_INPUT_PARTICLE_NUMBER
    int NBmin = Pall.Ns-Pall.Npm;
    int NBmax = Pall.Ns+Pall.Npm;
#endif // USE_SN_INPUT_PARTICLE_NUMBER //}

#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
    double Nlist = AP_i->SmoothedNumber*
        SmoothedMassConversionFactor*CUBE(AP_i->Radius);
#else // USE_SMOOTHED_NEIGHBOR_NUMBER
    int Nlist = AP_i->Nlist;
#endif  // USE_SMOOTHED_NEIGHBOR_NUMBER //}

    // Here, the convergence condition is satisfied. 
    if(((NBmin)<=Nlist)&&(Nlist<=(NBmax))){
        CSExportFlags_i[NProcs-1] = false;
        return true;
    }else if(((0.8*NBmin)<=Nlist)&&(Nlist<=(2.0*NBmax))&&(AP_i->Rvalue>0.e0)&&(AP_i->Lvalue>0.e0)){
        if(AP_i->Rvalue-AP_i->Lvalue < 1.e-6*AP_i->Lvalue){
            CSExportFlags_i[NProcs-1] = false;
            return true;
        }
    }
    if(CSExportFlags_i[NProcs-1]){
        if(Nlist<NBmin){
            AP_i->Lvalue = fmax(AP_i->Lvalue,AP_i->Radius);
        } else if(Nlist>NBmax){
            if(AP_i->Rvalue > 0.e0){
                AP_i->Rvalue = fmin(AP_i->Rvalue,AP_i->Radius);
            }else{
                AP_i->Rvalue = AP_i->Radius;
            }
        }

        if((AP_i->Lvalue>0.e0)&&(AP_i->Rvalue>0.e0)){
            AP_i->Radius = cbrt(0.5*(CUBE(AP_i->Lvalue)+CUBE(AP_i->Rvalue)));
        }else{
            if((AP_i->Rvalue == 0.e0)&&(AP_i->Lvalue > 0.e0)){
                AP_i->Radius *= StellarFeedbackRadiusFactInc;
            }else if((AP_i->Rvalue > 0.e0)&&(AP_i->Lvalue == 0.e0)){
                AP_i->Radius *= StellarFeedbackRadiusFactDec;
            }
        }
    }
    return false;
}

static inline bool __attribute__((always_inline)) CS_CheckConvergence_i(const int NProcs, bool CSExportFlags_i[restrict], struct StructActiveParticle *AP_i, const int Type){ 
    if(AP_i->Type == CS_TypeHydro){
        return CS_CheckNeighborNumberAndUpdateHydroRadius_i(NProcs,CSExportFlags_i,AP_i); 
    } else if(AP_i->Type == CS_TypeHII){
        return CS_CheckLocalMassAndUpdateHIIRadius_i(NProcs,CSExportFlags_i,AP_i); 
    } else if(AP_i->Type == CS_TypeSN){
        return CS_CheckNeighborNumberAndUpdateFeedbackRadius_i(NProcs,CSExportFlags_i,AP_i); 
    }
    return false;
}


static void CS_UpdateRadiusAndOthers_i(struct StructActiveParticle *AP_i, int Neighbors[restrict]){

    if(AP_i->Type == CS_TypeHydro){
#if 0
        AP_i->Nlist = 
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
            CS_GetSmoothedNumberofNeighbors(AP_i->Pos,AP_i->Radius,
                    Neighbors,&(AP_i->SmoothedNumber));
#else // USE_SMOOTHED_NEIGHBOR_NUMBER 
            CS_GetNumberofNeighbors(AP_i->Pos,AP_i->Radius,Neighbors);
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
        // if(AP_i->GlobalID == 6){
            // fprintf(stderr," -- %d %g -- \n",AP_i->Nlist,AP_i->Radius);
        // }
#endif 
        struct StructCSHydroLocalInfo TemporalData = 
            CS_GetNumberofNeighbors(AP_i->Pos,AP_i->Radius,Neighbors);
        AP_i->Nlist = TemporalData.Nlist;
#if VISCOSITY_TYPE == 1 //{
        AP_i->Body.Hydro.B[0][0] = TemporalData.B[0][0];
        AP_i->Body.Hydro.B[0][1] = TemporalData.B[0][1];
        AP_i->Body.Hydro.B[0][2] = TemporalData.B[0][2];
        AP_i->Body.Hydro.B[1][0] = TemporalData.B[1][0];
        AP_i->Body.Hydro.B[1][1] = TemporalData.B[1][1];
        AP_i->Body.Hydro.B[1][2] = TemporalData.B[1][2];
        AP_i->Body.Hydro.B[2][0] = TemporalData.B[2][0];
        AP_i->Body.Hydro.B[2][1] = TemporalData.B[2][1];
        AP_i->Body.Hydro.B[2][2] = TemporalData.B[2][2];
#endif // VISCOSITY_TYPE //}

#ifdef USE_SINK_PARTICLE //{
        AP_i->Body.Hydro.PotentialMin = TemporalData.PotentialMin;
        AP_i->Body.Hydro.MassTotal = TemporalData.MassTotal;
        AP_i->Body.Hydro.VCOM[0] = TemporalData.VCOM[0];
        AP_i->Body.Hydro.VCOM[1] = TemporalData.VCOM[1];
        AP_i->Body.Hydro.VCOM[2] = TemporalData.VCOM[2];
#endif // USE_SINK_PARTICLE //}

    } else if(AP_i->Type == CS_TypeHII){
        struct StructCSHIILocalInfo TemporalData =
            CS_ReturnHIILocalInfo(AP_i->Pos,AP_i->Radius,Neighbors);

        AP_i->Nlist = TemporalData.Nlist;
#ifdef __PHOTON_COUNT_BASE__
        AP_i->Body.HII.PhotonCount = TemporalData.PhotonCount;
        AP_i->Body.HII.PhotonCountDistanceMin = TemporalData.PhotonCountDistanceMin;
#else //__PHOTON_COUNT_BASE__
        AP_i->Body.HII.Mass = TemporalData.Mass;
        AP_i->Body.HII.MassDistanceMin = TemporalData.MassDistanceMin;
#endif // __PHOTON_COUNT_BASE__
        AP_i->Body.HII.DistanceMin = TemporalData.DistanceMin;
    } else if(AP_i->Type == CS_TypeSN){
        struct StructCSStellarFeedbackLocalInfo TemporalData = 
            CS_RetrunStellarFeedbackLocalInfo(AP_i->Pos,AP_i->Radius,Neighbors);

        AP_i->Nlist = TemporalData.Nlist;
        AP_i->Body.SN.Density = TemporalData.Density;
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
        AP_i->Body.SN.SmoothedNumber = TemporalData.SmoothedNumber;
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
#ifdef SET_SNII_TEMPERATURE //{
        AP_i->Body.SN.GasMass = TemporalData.GasMass;
#endif //SET_SNII_TEMPERATURE //}
#ifdef MAXIMUM_ENERGY_INPUT //{
        if(TemporalData.Nlist > 0){
            AP_i->Body.SN.DistanceMin = TemporalData.DistanceMin;
            AP_i->Body.SN.DistanceMinGlobalID = TemporalData.DistanceMinGlobalID;
        }
#endif // MAXIMUM_ENERGY_INPUT //}
#ifdef __CHECK_SUM__ //{
        AP_i->Body.SN.CheckSum = TemporalData.CheckSum;
#endif //__CHECK_SUM__ //}
    }
    return ;
}


static void CS_UpdateRadiusAndOthersImported_i(struct StructCSImport *CSImportSend_i, struct StructCSExport *CSExportRecv_i, int Neighbors[restrict]){

    if(CSExportRecv_i->Type == CS_TypeHydro){

#if 0
        CSImportSend_i->SmoothedNumber = 0.e0;
        CSImportSend_i->Nlist = 
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
            CS_GetSmoothedNumberofNeighbors(CSExportRecv_i->Pos,
                    CSExportRecv_i->Radius,Neighbors,&(CSExportRecv_i->SmoothedNumber));
#else // USE_SMOOTHED_NEIGHBOR_NUMBER 
            CS_GetNumberofNeighbors(CSExportRecv_i->Pos,CSExportRecv_i->Radius,Neighbors);
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}

        // if(CSExportRecv_i->GlobalID == 6){
            // fprintf(stderr," -i %d %g i- \n",CSImportSend_i->Nlist,CSExportRecv_i->Radius);
        // }
#endif
        struct StructCSHydroLocalInfo TemporalData = 
            CS_GetNumberofNeighbors(CSExportRecv_i->Pos,CSExportRecv_i->Radius,Neighbors);
        CSImportSend_i->Nlist = TemporalData.Nlist;

#if VISCOSITY_TYPE == 1 //{
        CSImportSend_i->B[0][0] = TemporalData.B[0][0];
        CSImportSend_i->B[0][1] = TemporalData.B[0][1];
        CSImportSend_i->B[0][2] = TemporalData.B[0][2];
        CSImportSend_i->B[1][0] = TemporalData.B[1][0];
        CSImportSend_i->B[1][1] = TemporalData.B[1][1];
        CSImportSend_i->B[1][2] = TemporalData.B[1][2];
        CSImportSend_i->B[2][0] = TemporalData.B[2][0];
        CSImportSend_i->B[2][1] = TemporalData.B[2][1];
        CSImportSend_i->B[2][2] = TemporalData.B[2][2];
#endif // VISCOSITY_TYPE //}

#ifdef USE_SINK_PARTICLE //{
        CSImportSend_i->PotentialMin = TemporalData.PotentialMin;
        CSImportSend_i->MassTotal = TemporalData.MassTotal;
        CSImportSend_i->VCOM[0] = TemporalData.VCOM[0];
        CSImportSend_i->VCOM[1] = TemporalData.VCOM[1];
        CSImportSend_i->VCOM[2] = TemporalData.VCOM[2];
#endif // USE_SINK_PARTICLE //}

    } else if(CSExportRecv_i->Type == CS_TypeHII){
        struct StructCSHIILocalInfo TemporalData =
            CS_ReturnHIILocalInfo(CSExportRecv_i->Pos,CSExportRecv_i->Radius,Neighbors);

        CSImportSend_i->Nlist = TemporalData.Nlist;
#ifdef __PHOTON_COUNT_BASE__
        CSImportSend_i->PhotonCount = TemporalData.PhotonCount;
        CSImportSend_i->PhotonCountDistanceMin = TemporalData.PhotonCountDistanceMin;
#else //__PHOTON_COUNT_BASE__
        CSImportSend_i->Mass = TemporalData.Mass;
        CSImportSend_i->MassDistMin = TemporalData.MassDistanceMin;
#endif // __PHOTON_COUNT_BASE__
        CSImportSend_i->DistanceMin = TemporalData.DistanceMin;

    } else if(CSExportRecv_i->Type == CS_TypeSN){
        struct StructCSStellarFeedbackLocalInfo TemporalData = 
            CS_RetrunStellarFeedbackLocalInfo(CSExportRecv_i->Pos,CSExportRecv_i->Radius,Neighbors);

        CSImportSend_i->Nlist = TemporalData.Nlist;
        CSImportSend_i->Density = TemporalData.Density;
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
        CSImportSend_i->SmoothedNumber = TemporalData.SmoothedNumber;
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
#ifdef SET_SNII_TEMPERATURE //{
        CSImportSend_i->GasMass = TemporalData.GasMass;
#endif //SET_SNII_TEMPERATURE //}
#ifdef MAXIMUM_ENERGY_INPUT //{
        if(TemporalData.Nlist > 0){
            CSImportSend_i->DistanceMin = TemporalData.DistanceMin;
            CSImportSend_i->DistanceMinGlobalID = TemporalData.DistanceMinGlobalID;
        }
#endif // MAXIMUM_ENERGY_INPUT //}
#ifdef __CHECK_SUM__ //{
        CSImportSend_i->CheckSum = TemporalData.CheckSum;
#endif //__CHECK_SUM__ //}
    }

    CSImportSend_i->Leaf = CSExportRecv_i->Leaf;
    CSImportSend_i->Type = CSExportRecv_i->Type;
#ifdef USE_DEBUG_MODE //{
    CSImportSend_i->GlobalID = CSExportRecv_i->GlobalID;
#endif // USE_DEBUG_MODE //}
    return ;
}


static inline bool  __attribute__((always_inline)) CS_CheckInLocalDomain(double Pos[], double Extent /* Radiusx2 */, const int MyID){
    for(int k=0;k<3;k++){
        if(Pos[k]+Extent > EdgesForHydro[MyID].PosMax[k]) return false;
        if(Pos[k]-Extent < EdgesForHydro[MyID].PosMin[k]) return false;
    }
    return true;
}

static inline bool __attribute__((always_inline)) CS_OverlapDomainKernel(double Pos[restrict], const double h /* 2*Radius */, const int NodeID){ 

    double Dist2 = 0.e0;
    for(int k=0;k<3;k++){
        if(Pos[k] < EdgesForHydro[NodeID].PosMin[k]) 
            Dist2 += SQ(EdgesForHydro[NodeID].PosMin[k]-Pos[k]);
        if(Pos[k] > EdgesForHydro[NodeID].PosMax[k])
            Dist2 += SQ(EdgesForHydro[NodeID].PosMax[k]-Pos[k]);
    }
    return (Dist2 < SQ(h));
}

static inline void __attribute__((always_inline)) CS_UpdateRadiusLocal(struct StructActiveParticle *AP_i, const int MyID, const int NProcs, bool CSExportFlags_i[restrict]){

    if(CS_CheckConvergence_i(NProcs,CSExportFlags_i,AP_i,AP_i->Type) == true)
        return;

    int Neighbors[MaxNeighborSize];
    int counter = 0;
    do{
        if(!CS_CheckInLocalDomain(AP_i->Pos,2.0*AP_i->Radius,MyID)) return;
        for(int i=0;i<CS_NContactedDomains;i++){
            int NodeID = CS_ContactedDomainID[i];
            if(CS_OverlapDomainKernel(AP_i->Pos,2.0*AP_i->Radius,NodeID)) return;
        }

        CS_UpdateRadiusAndOthers_i(AP_i,Neighbors);

        counter ++;
        if(counter > 10) return ;
    }while(CS_CheckConvergence_i(NProcs,CSExportFlags_i,AP_i,AP_i->Type) == false);

    return;
}


static inline int __attribute__((always_inline)) CS_CheckConvergence(const int NActives, const int NProcs, bool CSExportFlags[restrict], struct StructActiveParticle AP[restrict]){ 


    int NLocalActiveLeaves = 0;
    for(int i=0;i<NActives;i++){
        int Offset = i*NProcs;
        if(CSExportFlags[Offset+NProcs-1]){ 
#ifdef UPDATE_SIZE_LOCAL //{
        if(AP[i].LocalUpdateFlags == false){ 
#endif // UPDATE_SIZE_LOCAL //} 
            if(CS_CheckConvergence_i(NProcs,CSExportFlags+Offset,AP+i,AP[i].Type) == false){
                NLocalActiveLeaves ++;
            }
#ifdef UPDATE_SIZE_LOCAL //{
        } else {
            NLocalActiveLeaves ++;
        }
#endif // UPDATE_SIZE_LOCAL //}
        }
    }

    return NLocalActiveLeaves;
}

#ifdef ADD_PERTURBATION //{
#define _ResetTiming_ 1
static void ResetLRvalues(const int NActives, const int Niteration, const int NProcs, bool CSExportFlags[restrict], struct StructActiveParticle AP[restrict]){ 

    if((Niteration+1)%(_ResetTiming_*MaxIterationTimes) == 0){
        for(int i=0;i<NActives;i++){
            int Offset = i*NProcs;
            if(CSExportFlags[Offset+NProcs-1]){ 
                AP[i].Rvalue = AP[i].Lvalue = 0.e0;
            }
        }
    }
    return ;
}
#endif // ADD_PERTURBATION //}

static void CalcMomentInv(double B[3][3]){

#if (DIMENSION == 2)
    B[2][2] = 1.0;
#endif // DIMENSION == 2

#if (DIMENSION == 1)
    B[1][1] = B[2][2] = 1.0;
#endif // DIMENSION == 2

    double C[3][3];
    double det = B[0][0]*B[1][1]*B[2][2] + B[1][0]*B[2][1]*B[0][2] + B[2][0]*B[0][1]*B[1][2]
                -B[0][0]*B[2][1]*B[1][2] - B[2][0]*B[1][1]*B[0][2] - B[1][0]*B[0][1]*B[2][2];
    det += 1.e-16;

    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            C[i][j] = B[i][j];
        }
    }

    B[0][0] = (C[1][1]*C[2][2]-C[1][2]*C[2][1])/det;
    B[0][1] = (C[0][2]*C[2][1]-C[0][1]*C[2][2])/det;
    B[0][2] = (C[0][1]*C[1][2]-C[0][2]*C[1][1])/det;
    B[1][0] = (C[1][2]*C[2][0]-C[1][0]*C[2][2])/det;
    B[1][1] = (C[0][0]*C[2][2]-C[0][2]*C[2][0])/det;
    B[1][2] = (C[0][2]*C[1][0]-C[0][0]*C[1][2])/det;
    B[2][0] = (C[1][0]*C[2][1]-C[1][1]*C[2][0])/det;
    B[2][1] = (C[0][1]*C[2][0]-C[0][0]*C[2][1])/det;
    B[2][2] = (C[0][0]*C[1][1]-C[0][1]*C[1][0])/det;

    return ;
}

extern struct StructActiveSNParticle *ActiveSNParticle;

static void RestoreAllData(const int NActives, struct StructActiveParticle AP[restrict]){

    int counter_SN = 0;
    for(int i=0;i<NActives;i++){
        if(AP[i].Type == CS_TypeHydro){
            int Index = AP[i].Index;
            Phydro[Index]->Kernel = Phydro[Index]->KernelPred = AP[i].Radius;
            int CacheIndex = AP[i].Body.Hydro.CacheIndex;
            NBCache[CacheIndex].Kernel = AP[i].Radius;
#if VISCOSITY_TYPE == 1 //{
            // fprintf(stderr,"AAAA %d  %g %g %g\n",AP[i].Nlist,AP[i].Body.Hydro.B[0][0],AP[i].Body.Hydro.B[0][1],AP[i].Body.Hydro.B[0][2]);
            CalcMomentInv(AP[i].Body.Hydro.B);
            Phydro[Index]->Bxx = AP[i].Body.Hydro.B[0][0];
            Phydro[Index]->Bxy = AP[i].Body.Hydro.B[0][1];
            Phydro[Index]->Bxz = AP[i].Body.Hydro.B[0][2];
            Phydro[Index]->Byx = AP[i].Body.Hydro.B[1][0];
            Phydro[Index]->Byy = AP[i].Body.Hydro.B[1][1];
            Phydro[Index]->Byz = AP[i].Body.Hydro.B[1][2];
            Phydro[Index]->Bzx = AP[i].Body.Hydro.B[2][0];
            Phydro[Index]->Bzy = AP[i].Body.Hydro.B[2][1];
            Phydro[Index]->Bzz = AP[i].Body.Hydro.B[2][2];
            // fprintf(stderr,"BBBB %g %g %g\n",Phydro[Index]->Bxx,Phydro[Index]->Bxy,Phydro[Index]->Bxz);
#endif // VISCOSITY_TYPE == 1 //}
        } else if(AP[i].Type == CS_TypeHII){
            int Index = AP[i].Index;
#ifdef __PHOTON_COUNT_BASE__
            Pstar[Index]->Density = AP[i].Body.HII.PhotonCount;
#else // __PHOTON_COUNT_BASE__
            Pstar[Index]->Density = AP[i].Body.HII.Mass;
#endif // __PHOTON_COUNT_BASE__
            if(AP[i].Body.HII.HIIRegion){
                Pstar[Index]->StromgrenRadius = 2.0*AP[i].Radius;
            } else {
                Pstar[Index]->StromgrenRadius = NONE;
            }
        } else if(AP[i].Type == CS_TypeSN){
            // Use the following function.
            //void CalcSizeSetSNInfo(struct StructActiveSNParticle ActiveSNParticle[]);
#if 0
            ActiveSNParticle[counter_SN].Pos[0] = AP[i].Pos[0];
            ActiveSNParticle[counter_SN].Pos[1] = AP[i].Pos[1];
            ActiveSNParticle[counter_SN].Pos[2] = AP[i].Pos[2];
            ActiveSNParticle[counter_SN].Radius = AP[i].Radius;

            ActiveSNParticle[counter_SN].Type = AP[i].Body.SN.Type;
            ActiveSNParticle[counter_SN].Count = AP[i].Body.SN.Count;
            ActiveSNParticle[counter_SN].InitialMass = AP[i].Body.SN.InitialMass;
            ActiveSNParticle[counter_SN].Metallicity = AP[i].Body.SN.Metallicity;
            ActiveSNParticle[counter_SN].IterationCount = 0;

            ActiveSNParticle[counter_SN].Nlist = AP[i].Nlist;
            ActiveSNParticle[counter_SN].Density = AP[i].Body.SN.Density;
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
            ActiveSNParticle[counter_SN].SmoothedNumber = AP[i].SmoothedNumber;
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
#ifdef SET_SNII_TEMPERATURE //{
            ActiveSNParticle[counter_SN].GasMass = AP[i].Body.SN.GasMass;
#endif //SET_SNII_TEMPERATURE //}
#ifdef MAXIMUM_ENERGY_INPUT //{
            ActiveSNParticle[counter_SN].DistanceMin = AP[i].Body.SN.DistanceMin;
            ActiveSNParticle[counter_SN].DistanceMinGlobalID = AP[i].Body.SN.DistanceMinGlobalID;
#endif // MAXIMUM_ENERGY_INPUT //}
#ifdef __CHECK_SUM__ //{
            ActiveSNParticle[counter_SN].CheckSum = AP[i].Body.SN.CheckSum;
#endif //__CHECK_SUM__ //}

#endif
            counter_SN ++;
        }
    }

    return ;
}

int CSEntryNumbers[3];
int CSEntryOffset[3];
int GlobalCSEntryNumbers[3];
int ReturnCalcSizeElementNumber(const int Type, const bool Global){
    if(Global == true){
        return GlobalCSEntryNumbers[Type];
    } else {
        return CSEntryNumbers[Type];
    }
}

void CalcSizeGetHydroInfo_i(const int Index, double *PotentialMin, double VCOM[]){

#ifdef USE_SINK_PARTICLE //{
    *PotentialMin = ActiveParticle[Index].Body.Hydro.PotentialMin;
    
    VCOM[0] = ActiveParticle[Index].Body.Hydro.VCOM[0];
    VCOM[1] = ActiveParticle[Index].Body.Hydro.VCOM[1];
    VCOM[2] = ActiveParticle[Index].Body.Hydro.VCOM[2];
#endif // USE_SINK_PARTICLE //}

    return ;
}

void CalcSizeSetSNInfo(struct StructActiveSNParticle ActiveSNParticle[]){

    int Offset = CSEntryOffset[CS_TypeSN];
    for(int i=0;i<CSEntryNumbers[CS_TypeSN];i++){
        int Index = Offset+i;
        ActiveSNParticle[i].Index = ActiveParticle[Index].Index;

        ActiveSNParticle[i].Pos[0] = ActiveParticle[Index].Pos[0];
        ActiveSNParticle[i].Pos[1] = ActiveParticle[Index].Pos[1];
        ActiveSNParticle[i].Pos[2] = ActiveParticle[Index].Pos[2];
        ActiveSNParticle[i].Radius = ActiveParticle[Index].Radius;

        ActiveSNParticle[i].Type = ActiveParticle[Index].Body.SN.Type;
        ActiveSNParticle[i].Count = ActiveParticle[Index].Body.SN.Count;
        ActiveSNParticle[i].InitialMass = ActiveParticle[Index].Body.SN.InitialMass;
        ActiveSNParticle[i].Metallicity = ActiveParticle[Index].Body.SN.Metallicity;
        ActiveSNParticle[i].IterationCount = 0;

        ActiveSNParticle[i].Nlist = ActiveParticle[Index].Nlist;
        ActiveSNParticle[i].Density = ActiveParticle[Index].Body.SN.Density;
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
        ActiveSNParticle[i].SmoothedNumber = ActiveParticle[Index].SmoothedNumber;
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
#ifdef SET_SNII_TEMPERATURE //{
        ActiveSNParticle[i].GasMass = ActiveParticle[Index].Body.SN.GasMass;
#endif //SET_SNII_TEMPERATURE //}
#ifdef MAXIMUM_ENERGY_INPUT //{
        ActiveSNParticle[i].DistanceMin = ActiveParticle[Index].Body.SN.DistanceMin;
        ActiveSNParticle[i].DistanceMinGlobalID = ActiveParticle[Index].Body.SN.DistanceMinGlobalID;
#endif // MAXIMUM_ENERGY_INPUT //}
#ifdef __CHECK_SUM__ //{
        ActiveSNParticle[i].CheckSum = ActiveParticle[Index].Body.SN.CheckSum;
#endif //__CHECK_SUM__ //}
    }

    return ;
}

static void FinalProcedure(const int NActives, struct StructActiveParticle AP[restrict]){

    for(int i=0;i<NActives;i++){
        if(AP[i].Type == CS_TypeHydro){
#ifdef USE_SINK_PARTICLE //{
            double InvMass = 1.0/AP[i].Body.Hydro.MassTotal;
            AP[i].Body.Hydro.VCOM[0] *= InvMass;
            AP[i].Body.Hydro.VCOM[1] *= InvMass;
            AP[i].Body.Hydro.VCOM[2] *= InvMass;
#endif // USE_SINK_PARTICLE //}
        } else if(AP[i].Type == CS_TypeHII){
        } else if(AP[i].Type == CS_TypeSN){
        }
    }

    return ;
}


static bool first = true;
void CalcSize(void){

    double TimingResultThisRoutine = GetElapsedTime();

#ifdef EVALUATE_KERNEL_BY_ITERATION
    if(first){
        CS_AllocateContactedDomainID();
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER
        SmoothedMassConversionFactor = (4.0*M_PI/3.0)*8.0;
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER
        first = false;
    }

    OverMaxIterationTimes = false;

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();
    int Neighbors[MaxNeighborSize];

    static bool *CSExportFlags;

    int NLocalActiveYoungStars = CheckHIIflags(0);
#ifdef USE_CELIB //{
    int NLocalActiveFeedback = CountFeedbackNumber();
#else // USE_CELIB //}//{
    int NLocalActiveFeedback = 0;
#endif //USE_CELIB //}

    int MaxEntry = MAX(Pall.Nhydro+(NLocalActiveYoungStars+NLocalActiveFeedback),NAdditionUnit);
    if(CSExportFlagsMaxAllocated < MaxEntry){
        if(CSExportFlagsMaxAllocated > 0){
            free(CSExportFlags);
            free(ActiveParticle);
        }
        CSExportFlagsMaxAllocated = MaxEntry;
        CSExportFlags = malloc(sizeof(bool)*CSExportFlagsMaxAllocated*NProcs);
        ActiveParticle = malloc(sizeof(struct StructActiveParticle)*CSExportFlagsMaxAllocated);
    }

    
    CSEntryNumbers[CS_TypeHydro] = CSEntryNumbers[CS_TypeHII] = CSEntryNumbers[CS_TypeSN] = 0;

    int NActives = 0;
    CSEntryNumbers[CS_TypeHydro] = HydroEntry(ActiveParticle,CSExportFlags,NProcs);
    CSEntryOffset[CS_TypeHydro] = NActives;
    NActives += CSEntryNumbers[CS_TypeHydro];
    CSEntryNumbers[CS_TypeSN] = StellarFeedbackEntry(ActiveParticle+NActives,CSExportFlags+NActives*NProcs,NProcs);
    CSEntryOffset[CS_TypeSN] = NActives;
    NActives += CSEntryNumbers[CS_TypeSN];
    CSEntryNumbers[CS_TypeHII] += HIIregionEntry(ActiveParticle+NActives,CSExportFlags+NActives*NProcs,NProcs);
    CSEntryOffset[CS_TypeHII] = NActives;
    NActives += CSEntryNumbers[CS_TypeHII];
    MPI_Allreduce(CSEntryNumbers,GlobalCSEntryNumbers,3,MPI_INT,MPI_SUM,MPI_COMM_WORLD);


    int BitMask = 0x01; 
    int NExportThisTime[NProcs-1];
    int NImportThisTime[NProcs-1];
    int NExportThisTimeNew[NProcs];
    int NImportThisTimeNew[NProcs];

    struct StructCSExport *CSExportSend[NProcs-1];
    struct StructCSExport *CSExportRecv = NULL;
    struct StructCSImport *CSImportSend = NULL;
    struct StructCSImport *CSImportRecv[NProcs-1];
    MPI_Status mpi_status_Export_Send[NProcs-1];
    MPI_Request mpi_request_Export_Send[NProcs-1];
    MPI_Status mpi_status_Export_Recv[NProcs-1];
    MPI_Request mpi_request_Export_Recv[NProcs-1];

#ifdef PRINT_LOG_KERNEL_ITERATION
    if(MPIGetMyID()==MPI_ROOT_RANK)
        fprintf(stderr,"Kernel Iteration");
#endif // PRINT_LOG_KERNEL_ITERATION

    int NActiveLeaves;
    MPI_Allreduce(&NActives,&NActiveLeaves,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

    CS_CheckContactedDomain();
    double BoxCenter[] = {HydroNode[0].Pos[0],HydroNode[0].Pos[1],HydroNode[0].Pos[2]}; // need check
    //double BoxCenter[] = {GravityNode[0].Pos[0],GravityNode[0].Pos[1],GravityNode[0].Pos[2]};

    int Niteration = 0;
    do{
#ifdef PRINT_LOG_KERNEL_ITERATION
        if(MPIGetMyID()==MPI_ROOT_RANK)
            fprintf(stderr,":[%d] = %d ",Niteration,NActiveLeaves);
#endif // PRINT_LOG_KERNEL_ITERATION

        LocalKernelMax = 0.e0;  
        for(int i=0;i<NActives;i++){  // Clear Entries
            int Offset = i*NProcs;
            if(CSExportFlags[Offset+NProcs-1]){
                ActiveParticle[i].Nlist = 0;
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
                ActiveParticle[i].SmoothedNumber = 0.e0;
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
                for(int k=0;k<NProcs-1;k++)
                    CSExportFlags[Offset+k] = false;
                LocalKernelMax = fmax(LocalKernelMax,DISTANCE(BoxCenter,ActiveParticle[i].Pos)+2.0*ActiveParticle[i].Radius);
                if(ActiveParticle[i].Type == CS_TypeHII){
#if VISCOSITY_TYPE == 1 //{
                    ActiveParticle[i].Body.Hydro.B[0][0] = 
                    ActiveParticle[i].Body.Hydro.B[0][1] = 
                    ActiveParticle[i].Body.Hydro.B[0][2] = 
                    ActiveParticle[i].Body.Hydro.B[1][0] = 
                    ActiveParticle[i].Body.Hydro.B[1][1] = 
                    ActiveParticle[i].Body.Hydro.B[1][2] = 
                    ActiveParticle[i].Body.Hydro.B[2][0] = 
                    ActiveParticle[i].Body.Hydro.B[2][1] = 
                    ActiveParticle[i].Body.Hydro.B[2][2] = 0.e0;
#endif // VISCOSITY_TYPE == 1 //}
#ifdef USE_SINK_PARTICLE //{
                    ActiveParticle[i].Body.Hydro.PotentialMin = 0.e0;
                    ActiveParticle[i].Body.Hydro.MassTotal = 0.e0;
                    ActiveParticle[i].Body.Hydro.VCOM[0] = 
                    ActiveParticle[i].Body.Hydro.VCOM[1] = 
                    ActiveParticle[i].Body.Hydro.VCOM[2] = 0.e0;
#endif // USE_SINK_PARTICLE //}
                } else 
                if(ActiveParticle[i].Type == CS_TypeHII){
                    ActiveParticle[i].Body.HII.DistanceMin = 0.e0;
#ifdef __PHOTON_COUNT_BASE__
                    ActiveParticle[i].Body.HII.PhotonCount = ActiveParticle[i].Body.HII.PhotonCountDistanceMin = 0.e0;
#else //__PHOTON_COUNT_BASE__
                    ActiveParticle[i].Body.HII.Mass = ActiveParticle[i].Body.HII.MassDistanceMin = 0.e0;
#endif // __PHOTON_COUNT_BASE__
                } else if(ActiveParticle[i].Type == CS_TypeSN){
                    ActiveParticle[i].Body.SN.Density = 0.0;
#ifdef SET_SNII_TEMPERATURE //{
                    ActiveParticle[i].Body.SN.GasMass = 0.e0;
#endif //SET_SNII_TEMPERATURE //}
#ifdef MAXIMUM_ENERGY_INPUT //{
                    ActiveParticle[i].Body.SN.DistanceMin = 2.0*ActiveParticle[i].Radius;
#endif // MAXIMUM_ENERGY_INPUT //}
                }
            }
        }

        int NExportMaxThisTime = 0;
        for(int i=0;i<NProcs-1;i++){
            NExportThisTime[i] = CS_CheckExportFlags(i,NProcs,CSExportFlags,NActives,ActiveParticle);

            CheckSizeofBufferExportSendIndex(NExportThisTime[i],sizeof(struct StructCSExport),i);
            CheckSizeofBufferImportRecvIndex(NExportThisTime[i],sizeof(struct StructCSImport),i);
            CSExportSend[i] = BufferExportSend[i];
            CSImportRecv[i] = BufferImportRecv[i];

            int NExport = 0;
            if(NExportThisTime[i] > 0){
                for(int k=0;k<NActives;k++){
                    int Offset = k*NProcs;
                    if(CSExportFlags[Offset+NProcs-1]){ 
                        if(CSExportFlags[Offset+i]&BitMask){ 
                            CSExportSend[i][NExport].Type = ActiveParticle[k].Type;
                            CSExportSend[i][NExport].Pos[0] = ActiveParticle[k].Pos[0];
                            CSExportSend[i][NExport].Pos[1] = ActiveParticle[k].Pos[1];
                            CSExportSend[i][NExport].Pos[2] = ActiveParticle[k].Pos[2];
                            CSExportSend[i][NExport].Radius = ActiveParticle[k].Radius;
                            CSExportSend[i][NExport].Leaf = k;
#ifdef USE_DEBUG_MODE //{
                            CSExportSend[i][NExport].GlobalID = ActiveParticle[k].GlobalID;
#endif // USE_DEBUG_MODE //}
                            NExport ++;
                        }
                    }
                }
            }
            NExportThisTime[i] = NExport;
            NExportMaxThisTime = MAX(NExportMaxThisTime,NExport);
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

        CheckSizeofBufferExportRecv(NImportAll,sizeof(struct StructCSExport));
        CheckSizeofBufferImportSend(NImportAll,sizeof(struct StructCSImport));
        CSExportRecv = BufferExportRecv;
        CSImportSend = BufferImportSend; 

        NImport = 0;

        int counter_send = 0;
        int counter_recv = 0;

        int SendFlag,RecvFlag;
        for(int i=0;i<NProcs-1;i++){
            if(NExportThisTime[i]>0){
                MPI_Isend(CSExportSend[i],
                    NExportThisTime[i]*sizeof(struct StructCSExport),
                        MPI_BYTE,CommunicationTable[i].SendRank,TAG_SPH_KERNEL_EXPORT+i,
                            MPI_COMM_WORLD,mpi_request_Export_Send+counter_send);
                MPI_Test(mpi_request_Export_Send+counter_send,&SendFlag,MPI_STATUS_IGNORE);
                counter_send ++;
            }
            if(NImportThisTime[i]>0){
                MPI_Irecv(CSExportRecv+NImport,
                    NImportThisTime[i]*sizeof(struct StructCSExport),
                        MPI_BYTE,CommunicationTable[i].RecvRank,TAG_SPH_KERNEL_EXPORT+i,
                            MPI_COMM_WORLD,mpi_request_Export_Recv+counter_recv);
                MPI_Test(mpi_request_Export_Recv+counter_recv,&RecvFlag,MPI_STATUS_IGNORE);
                counter_recv ++;
                NImport += NImportThisTime[i];
            }
        }


        double TimeNBS = GetElapsedTime();
        for(int i=0;i<NActives;i++){  // Check local
            int Offset = i*NProcs;
            if(CSExportFlags[Offset+NProcs-1]){
                CS_UpdateRadiusAndOthers_i(ActiveParticle+i,Neighbors);
                
#ifdef UPDATE_SIZE_LOCAL //{
                /// Insert Local Update Routine here.
                ActiveParticle[i].LocalUpdateFlags = false;
                int IsLocal = 0;
                for(int k=0;k<NProcs-1;k++){
                    if(CSExportFlags[Offset+k]&BitMask){
                        IsLocal ++;
                    }
                }

                if(IsLocal == 0){
                    ActiveParticle[i].LocalUpdateFlags = true;
                    CS_UpdateRadiusLocal(ActiveParticle+i,MyID,NProcs,CSExportFlags+Offset);
                }
#endif // UPDATE_SIZE_LOCAL // }
            }
        }
        TimingResults.HydroKernelNeighborSearchThisStep += GetElapsedTime()-TimeNBS;

        double TimeComm = GetElapsedTime();
        MPI_Waitall(counter_send,mpi_request_Export_Send,mpi_status_Export_Send);
        MPI_Waitall(counter_recv,mpi_request_Export_Recv,mpi_status_Export_Recv);
        TimingResults.HydroKernelCommThisStep += GetElapsedTime()-TimeComm;


        TimeNBS = GetElapsedTime();
        for(int i=0;i<NImportAll;i++){ // Imported
            CS_UpdateRadiusAndOthersImported_i(CSImportSend+i,CSExportRecv+i,Neighbors);
        }
        TimingResults.HydroKernelNeighborSearchThisStep += GetElapsedTime()-TimeNBS;

        NImportAll = 0;
        int NImportAllNew = 0;
        for(int i=0;i<NProcs-1;i++){
            NImportThisTimeNew[i] = 0;
            for(int k=0;k<NImportThisTime[i];k++){
                if(CSImportSend[NImportAll].Nlist > 0){
                    CSImportSend[NImportAllNew] = CSImportSend[NImportAll];
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
                MPI_Isend(CSImportSend+NImport,
                    NImportThisTimeNew[i]*sizeof(struct StructCSImport),
                        MPI_BYTE,CommunicationTable[i].SendRank,TAG_SPH_KERNEL_IMPORT+i,
                            MPI_COMM_WORLD,mpi_request_Export_Send+counter_send);
                MPI_Test(mpi_request_Export_Send+counter_send,&SendFlag,MPI_STATUS_IGNORE);
                counter_send ++;
            }
            if(NExportThisTimeNew[i]>0){
                MPI_Irecv(CSImportRecv[i],
                    NExportThisTimeNew[i]*sizeof(struct StructCSImport),
                        MPI_BYTE,CommunicationTable[i].RecvRank,TAG_SPH_KERNEL_IMPORT+i,
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
                int leaf = CSImportRecv[i][k].Leaf;
                if(CSImportRecv[i][k].Type == CS_TypeHydro){
                    ActiveParticle[leaf].Nlist += CSImportRecv[i][k].Nlist;
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
                    ActiveParticle[leaf].SmoothedNumber += CSImportRecv[i][k].SmoothedNumber;
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}

#if VISCOSITY_TYPE == 1 //{
                    ActiveParticle[leaf].Body.Hydro.B[0][0] += CSImportRecv[i][k].B[0][0];
                    ActiveParticle[leaf].Body.Hydro.B[0][1] += CSImportRecv[i][k].B[0][1];
                    ActiveParticle[leaf].Body.Hydro.B[0][2] += CSImportRecv[i][k].B[0][2];
                    ActiveParticle[leaf].Body.Hydro.B[1][0] += CSImportRecv[i][k].B[1][0];
                    ActiveParticle[leaf].Body.Hydro.B[1][1] += CSImportRecv[i][k].B[1][1];
                    ActiveParticle[leaf].Body.Hydro.B[1][2] += CSImportRecv[i][k].B[1][2];
                    ActiveParticle[leaf].Body.Hydro.B[2][0] += CSImportRecv[i][k].B[2][0];
                    ActiveParticle[leaf].Body.Hydro.B[2][1] += CSImportRecv[i][k].B[2][1];
                    ActiveParticle[leaf].Body.Hydro.B[2][2] += CSImportRecv[i][k].B[2][2];
#endif // VISCOSITY_TYPE == 1 //}

#ifdef USE_SINK_PARTICLE //{
                    ActiveParticle[leaf].Body.Hydro.PotentialMin = fmin(ActiveParticle[leaf].Body.Hydro.PotentialMin,CSImportRecv[i][k].PotentialMin);
                    ActiveParticle[leaf].Body.Hydro.MassTotal += CSImportRecv[i][k].MassTotal;
                    ActiveParticle[leaf].Body.Hydro.VCOM[0] += CSImportRecv[i][k].VCOM[0];
                    ActiveParticle[leaf].Body.Hydro.VCOM[1] += CSImportRecv[i][k].VCOM[1];
                    ActiveParticle[leaf].Body.Hydro.VCOM[2] += CSImportRecv[i][k].VCOM[2];
#endif // USE_SINK_PARTICLE //}
                } else if(CSImportRecv[i][k].Type == CS_TypeHII){
                    ActiveParticle[leaf].Nlist += CSImportRecv[i][k].Nlist;
#ifdef __PHOTON_COUNT_BASE__ //{
                    ActiveParticle[leaf].Body.HII.PhotonCount += CSImportRecv[i][k].PhotonCount;
#else // __PHOTON_COUNT_BASE__ //}//{
                    ActiveParticle[leaf].Body.HII.Mass += CSImportRecv[i][k].Mass;
#endif //__PHOTON_COUNT_BASE__ //}
                    if((ActiveParticle[leaf].Body.HII.DistanceMin > CSImportRecv[i][k].DistanceMin)&&(CSImportRecv[i][k].Nlist > 0)){
                        ActiveParticle[leaf].Body.HII.DistanceMin = CSImportRecv[i][k].DistanceMin;
#ifdef __PHOTON_COUNT_BASE__ //{
                        ActiveParticle[leaf].Body.HII.PhotonCountDistanceMin = CSImportRecv[i][k].PhotonCountDistanceMin;
#else // __PHOTON_COUNT_BASE__ //}//{
                        ActiveParticle[leaf].Body.HII.MassDistanceMin = CSImportRecv[i][k].MassDistanceMin;
#endif //__PHOTON_COUNT_BASE__ //}
                    }
                } else if(CSImportRecv[i][k].Type == CS_TypeSN){
                    ActiveParticle[leaf].Nlist += CSImportRecv[i][k].Nlist;
                    ActiveParticle[leaf].Body.SN.Density += CSImportRecv[i][k].Density;
                
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
                    ActiveParticle[leaf].Body.SN.SmoothedNumber += CSImportRecv[i][k].SmoothedNumber;
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
#ifdef SET_SNII_TEMPERATURE //{
                    ActiveParticle[leaf].Body.SN.GasMass += CSImportRecv[i][k].GasMass;
#endif //SET_SNII_TEMPERATURE //}
#ifdef MAXIMUM_ENERGY_INPUT //{
                    if(ActiveParticle[leaf].DistanceMin > CSImportRecv[i][k].DistanceMin){
                        ActiveParticle[leaf].DistanceMin = CSImportRecv[i][k].DistanceMin;
                        ActiveParticle[leaf].Body.SN.DistanceMinGlobalID = CSImportRecv[i][k].DistanceMinGlobalID;
                    }
#endif // MAXIMUM_ENERGY_INPUT //}
#ifdef __CHECK_SUM__ //{
                    ActiveParticle[leaf].Body.SN.CheckSum += CSImportRecv[i][k].CheckSum;
#endif //__CHECK_SUM__ //}
                }
            }
        }

#ifdef ADD_PERTURBATION //{
        // assert(Niteration < 1000);
        if(Niteration > 10*MaxIterationTimes) break;
        if(Niteration > MaxIterationTimes)
            OverMaxIterationTimes = true;
        ResetLRvalues(NActives,Niteration,NProcs,CSExportFlags,ActiveParticle);
#endif // ADD_PERTURBATION //}
        int NLocalActiveLeaves = CS_CheckConvergence(NActives,NProcs,CSExportFlags,ActiveParticle); 
        MPI_Allreduce(&NLocalActiveLeaves,&NActiveLeaves,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

        Niteration ++;
#ifndef ADD_PERTURBATION //{
        if(Niteration > 10*MaxIterationTimes)
            break;
#endif // ADD_PERTURBATION //}
    } while (0<NActiveLeaves);

#else // EVALUATE_KERNEL_BY_ITERATION
    for(int i=0;i<Pall.Nhydro;i++){
        if(Phydro[i]->Active){
            Phydro[i]->Kernel = KERNEL_FACTOR*
                pow(Phydro[i]->Mass/Phydro[i]->Rho,1.0/((double)DIMENSION));
            Phydro[i]->Rho = 0.e0;
        }
    }
#endif // EVALUATE_KERNEL_BY_ITERATION

    // Restore all data for Phydro/Pstar
    RestoreAllData(NActives,ActiveParticle);

    // Final procedure for VCOM
    FinalProcedure(NActives,ActiveParticle);


    PlantHydroTreeKernelMaxUpdate();
#ifdef EVALUATE_KERNEL_BY_ITERATION
#ifdef PRINT_LOG_KERNEL_ITERATION
    if(MPIGetMyID()==MPI_ROOT_RANK)
        fprintf(stderr,"\n");
#else // PRINT_LOG_KERNEL_ITERATION
    if(MPIGetMyID()==MPI_ROOT_RANK)
        fprintf(stderr,"%d iterations for kernel determination.\n",Niteration);
#endif // PRINT_LOG_KERNEL_ITERATION
#endif // EVALUATE_KERNEL_BY_ITERATION

    TimingResults.HydroKernelThisStep += GetElapsedTime()-TimingResultThisRoutine;

    return;
}
