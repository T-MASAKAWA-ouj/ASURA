#include "config.h"
#include "HydroMisc.h"
#include "PlantHydroTreeImported.h"
#include "NeighborSearch.h"
#include "KernelFunctions.h"
#include "ThermalConductivity.h"

#define SuppressTensileInstability (ON)

static struct StructHydroAccExport (*HydroAccExportRecvLog);
static int NumberofHydroAccExportRecvLog;
static int NumberofHydroAccExportRecvLogAllocated = 0;


#if defined(BAROTROPIC_EOS_RUN)
/*
 * Bate et al. (2003) MNRAS 339, 577-599
 * P = K \rho^\eta, where K = cs(T=10K)^2
 * if(\rho < \rho_c) \eta = 1
 * else              \eta = 7/5
 */
static inline double __attribute__((always_inline)) ReturnPoverRho2ForBarotropicRun(const double rho){
    /*
     * This function returns the value of P/\rho^2.
     * P/\rho^2 = K \rho^(\eta-2)
     */
    double K = SQ(SOUND_VELOCITY_FOR_BAROTROPIC_RUN*Pall.UnitTime/Pall.UnitLength);
    if(rho*Pall.ConvertDensityToCGS < CRITICAL_DENSITY_FOR_BAROTROPIC_RUN){
        return (K/rho);
    } else {
        return (K*pow(rho,-0.6));
    }
}

static inline double __attribute__((always_inline)) ReturnSoundSpeedForBarotoropicRun(const double rho){
    /*
     * P = K \rho^\eta, where K = cs(T=10K)^2
     * Cs(\rho)^2 = K \eta \rho^(\eta-1)
     */
    double K = SQ(SOUND_VELOCITY_FOR_BAROTROPIC_RUN*Pall.UnitTime/Pall.UnitLength);
    if(rho*Pall.ConvertDensityToCGS < CRITICAL_DENSITY_FOR_BAROTROPIC_RUN){
        return (K);
    } else {
        return (1.4*K*pow(rho,0.4));
    }
}
#endif

static inline double __attribute__((always_inline)) ReturnCS(const double Rhoi, const double Ui){

#if defined(ISOTHERMAL_EOS_RUN)
    return Pall.CS;
#elif defined(BAROTROPIC_EOS_RUN)
    return ReturnSoundSpeedForBarotoropicRun(Rhoi);
#else
    return sqrt(Pall.GGm1*Ui);
#endif
}

static inline double __attribute__((always_inline)) ReturnCSij(const double Rhoj, const double Uj, double csi){
#if defined(ISOTHERMAL_EOS_RUN)
    return Pall.CS;
#elif defined(BAROTROPIC_EOS_RUN)
    return 0.5*(csi+ReturnSoundSpeedForBarotoropicRun(Rhoj));
#else
    return 0.5*(csi+sqrt(Pall.GGm1*Uj));
#endif
}

static inline double __attribute__((always_inline)) ReturnPoverRho2(const double Rho, const double U){

#if defined(ISOTHERMAL_EOS_RUN)
    return (SQ(Pall.CS)/Rho);
#elif defined(BAROTROPIC_EOS_RUN)
    return ReturnPoverRho2ForBarotropicRun(Rho);
#else
    return Pall.Gm1*(U/Rho);
#endif
}


static inline double __attribute__((always_inline)) DomainDistanceSQR(double Pos[restrict], const int NodeID){

    double Dist2 = 0.e0;
    for(int k=0;k<3;k++){
        if(Pos[k] < EdgesForActiveHydro[NodeID].PosMin[k]){
            Dist2 += SQ(EdgesForActiveHydro[NodeID].PosMin[k]-Pos[k]);
        } if(Pos[k] > EdgesForActiveHydro[NodeID].PosMax[k]){
            Dist2 += SQ(EdgesForActiveHydro[NodeID].PosMax[k]-Pos[k]);
        }
    }
    return (Dist2);
}

static inline bool __attribute__((always_inline)) OverlapDomainGather(double Pos[restrict], const double h, const int NodeID){

    double Dist2 = 0.e0;
    for(int k=0;k<3;k++){
        if(Pos[k] < EdgesForActiveHydro[NodeID].PosMin[k]) 
            Dist2 += SQ(EdgesForActiveHydro[NodeID].PosMin[k]-Pos[k]);
        if(Pos[k] > EdgesForActiveHydro[NodeID].PosMax[k])
            Dist2 += SQ(EdgesForActiveHydro[NodeID].PosMax[k]-Pos[k]);
    }
    return (Dist2 < SQ(h));
}

static inline int __attribute__((always_inline)) CheckHydroAccExportFlags(const int Index, const int NProcs, 
        bool HydroAccExportFlags[][NProcs-1], bool HydroInteractionFlags[][NProcs-1]){

    if(Pall.Nhydro == 0) return 0;

    int NodeID = CommunicationTable[Index].SendRank;
    int NExport = 0;

    // Particles which are overlaped other domains are listed on the export particle list.
    int RootNodeID = 0;
    int CurrentNodeID = HydroNode[RootNodeID].Children;
	while(CurrentNodeID != RootNodeID){
        if( SQ(HydroNode[CurrentNodeID].KernelMax) < DomainDistanceSQR(HydroNode[CurrentNodeID].Pos,NodeID) ){
			CurrentNodeID = HydroNode[CurrentNodeID].Next;
		} else if(HydroNode[CurrentNodeID].Children != NONE){
			CurrentNodeID = HydroNode[CurrentNodeID].Children;
		} else {
			int NumberofLeaves = HydroNode[CurrentNodeID].NumberofLeaves;
            int header = HydroNode[CurrentNodeID].Leaves;
            for(int k=0;k<NumberofLeaves;k++){
                int leaf = header+k;
                // if(HydroRoot.Leaves[leaf] < 0) continue;
                if((HydroInteractionFlags[NBCache[leaf].Leaf][Index] == true)||
                    OverlapDomainGather(NBCache[leaf].Pos,2.0*NBCache[leaf].Kernel,NodeID)){
                    HydroAccExportFlags[NBCache[leaf].Leaf][Index] = true;
                    NExport ++;
                }
			}
			CurrentNodeID = HydroNode[CurrentNodeID].Next;
		}
	}

    // Particles of which export flags are on are listed on the export particle list.
    int NExportExtra = 0;
    for(int k=0;k<Pall.Nhydro;k++){
        if(HydroAccExportFlags[k][Index] == false){
            if(HydroInteractionFlags[k][Index] == true){
                HydroAccExportFlags[k][Index] = true;
                NExportExtra ++;
            }
        }
    }

#if 0
    NExport = 0;
    NExportExtra = 0;
    for(int k=0;k<Pall.Nhydro;k++){
        HydroAccExportFlags[k][Index] = true;
        NExportExtra ++;
    }
#endif

	return (NExport+NExportExtra);
}

static inline int __attribute__((always_inline)) CheckHydroAccExportFlagsModified(const int Index, const int NProcs, 
        bool ExportFlags[restrict], bool InteractionFlags[restrict]){

    if(Pall.Nhydro == 0) return 0;

    int NodeID = CommunicationTable[Index].SendRank;
    int NExport = 0;

    // Particles which are overlaped other domains are listed on the export particle list.
    int RootNodeID = 0;
    int CurrentNodeID = HydroNode[RootNodeID].Children;
	while(CurrentNodeID != RootNodeID){
        if( SQ(HydroNode[CurrentNodeID].KernelMax) < DomainDistanceSQR(HydroNode[CurrentNodeID].Pos,NodeID) ){
			CurrentNodeID = HydroNode[CurrentNodeID].Next;
		} else if(HydroNode[CurrentNodeID].Children != NONE){
		//if(HydroNode[CurrentNodeID].Children != NONE){
			CurrentNodeID = HydroNode[CurrentNodeID].Children;
		} else {
			int NumberofLeaves = HydroNode[CurrentNodeID].NumberofLeaves;
            int header = HydroNode[CurrentNodeID].Leaves;
            for(int k=0;k<NumberofLeaves;k++){
                int leaf = header+k;
                // if(HydroRoot.Leaves[leaf] < 0) continue;
                if((InteractionFlags[NBCache[leaf].Leaf] == true)||
                    OverlapDomainGather(NBCache[leaf].Pos,2.0*NBCache[leaf].Kernel,NodeID)){
                    ExportFlags[NBCache[leaf].Leaf] = true;
                    NExport ++;
                }
			}
			CurrentNodeID = HydroNode[CurrentNodeID].Next;
		}
	}

    // Particles of which export flags are on are listed on the export particle list.
    int NExportExtra = 0;
    for(int k=0;k<Pall.Nhydro;k++){
        //if((ExportFlags[k] == false)&&(InteractionFlags[k] == true)){
        if((!ExportFlags[k])&&(InteractionFlags[k])){
            ExportFlags[k] = true;
            NExportExtra ++;
            //ExportFlags[k] = InteractionFlags[k];
            //NExportExtra += (int)InteractionFlags[k];
        }
    }

	return (NExport+NExportExtra);
}

/*
 * Artificial Viscosity term.
 * Lattanzio et al. (1985)
 */
static double __attribute__((always_inline))  ReturnArtificialViscocityGM(const double vx, 
    const double r2, const double csi, const double csj,
    const double Rhoi, const double Rhoj, const double Kerneli, const double Kernelj, 
    const double Fi, const double Fj
#ifdef USE_VARIABLE_ALPHA //{
    , const double Alphai, const double Alphaj
#endif // USE_VARIABLE_ALPHA //}
    ){

    double visc = 0.e0;
    if(vx < 0.e0){ 
        double KernelMean = 0.5*(Kerneli+Kernelj);
        double mu = KernelMean*vx/(r2+Pall.HydroEta2*SQ(KernelMean));
        double cs = 0.5*(csi+csj);

#ifdef USE_VARIABLE_ALPHA //{
        double Alpha = 0.5*(Alphai+Alphaj);
        double Beta = 2.0*Alpha;
        visc = (-Alpha*cs*mu + Beta*SQ(mu))/(0.5*(Rhoi+Rhoj));
#else // USE_VARIABLE_ALPHA //}//{
        visc = (-Pall.HydroAlpha*cs*mu + Pall.HydroBeta*SQ(mu))/(0.5*(Rhoi+Rhoj));
#endif // USE_VARIABLE_ALPHA //}
#if DIMENSION > 1 //{
        visc *= 0.5*(Fi+Fj);
#endif // DIMENSION //}
    } else {
        visc = 0.e0;
    }

    return visc;
}

/*
 * Artificial Viscosity term.
 * Monaghan (1997)
 */
static double __attribute__((always_inline)) ReturnArtificialViscocityReimann(const double vx, 
    const double r, const double csi, const double csj, const double Rhoi, const double Rhoj, 
    const double Fi, const double Fj, double *Vsig
#ifdef USE_VARIABLE_ALPHA //{
    , const double Alphai, const double Alphaj
#endif // USE_VARIABLE_ALPHA //}
    ){

    double visc = 0.e0;
    *Vsig = fmax(csi+csj,*Vsig);
    if(vx < 0.e0){ 
        double wij = vx/r;
        double vsig = csi+csj-SIGNAL_VELOCITY_BETA*wij; 
#ifdef USE_VARIABLE_ALPHA //{
        double Alpha = 0.5*(Alphai+Alphaj);
        visc = -0.5*Alpha*vsig*wij/(0.5*(Rhoi+Rhoj));
#else // USE_VARIABLE_ALPHA //}/{
        visc = -0.5*Pall.HydroAlpha*vsig*wij/(0.5*(Rhoi+Rhoj));
#endif // USE_VARIABLE_ALPHA //}
#if DIMENSION > 1 //{
        visc *= 0.5*(Fi+Fj);
#endif // DIMENSION //}
        *Vsig = fmax(vsig,*Vsig);
    } else {
        visc = 0.e0;
    }

    return visc;
}

#if VISCOSITY_TYPE == 1 //{
#ifdef USE_VARIABLE_ALPHA 
#define SetParameterVariableAlphaHosono(_Alphai) \
    , _Alphai
#else // USE_VARIABLE_ALPHA
#define SetParameterVariableAlphaHosono(_Alphai) 
#endif // USE_VARIABLE_ALPHA
static double __attribute__((always_inline)) ReturnArtificialViscocityHosono(const double DivVi, 
    const double csi, const double Rhoi, const double Kerneli, const double Fi
#ifdef USE_VARIABLE_ALPHA //{
    , const double Alphai
#endif // USE_VARIABLE_ALPHA //}
    ){

#ifdef USE_VARIABLE_ALPHA //{
    double ViscAlpha = Alphai;
#else // USE_VARIABLE_ALPHA //}/{
    double ViscAlpha = Pall.HydroAlpha;
#endif // USE_VARIABLE_ALPHA //}
    
    double qAVi = (DivVi < 0.0) ? (-ViscAlpha*Rhoi*csi*Kerneli*DivVi
            + 2.0*ViscAlpha*Rhoi*SQ(Kerneli*DivVi)) : 0.0;
#if (DIMENSION > 1)
    qAVi *= Fi;
#endif
    return qAVi;
}
#endif // VISCOSITY_TYPE == 1 //}


#ifdef USE_VARIABLE_ALPHA 
#define SetParameterVariableAlpha(_Alphai, _Alphaj) \
    , _Alphai, _Alphaj
#else // USE_VARIABLE_ALPHA
#define SetParameterVariableAlpha(_Alphai, _Alphaj) 
#endif // USE_VARIABLE_ALPHA

static struct StructHydroAccImport InitializeHydroAccImport = {
    .Du = 0.e0,
    .HydroAcc = {0.e0,0.e0,0.e0},
    .Vsig = 0.e0,
    .Nlist = 0,
    .Leaf = 0,
#ifdef USE_SPSPH //{ 
    .DZw = 0.e0,
#endif // USE_SPSPH //}
#ifdef USE_CELIB // {
#ifdef USE_METAL_DIFFUSION //{
    .dZ[0] = 0.e0,
    .dZ[1] = 0.e0,
    .dZ[2] = 0.e0,
    .dZ[3] = 0.e0,
    .dZ[4] = 0.e0,
    .dZ[5] = 0.e0,
    .dZ[6] = 0.e0,
    .dZ[7] = 0.e0,
    .dZ[8] = 0.e0,
    .dZ[9] = 0.e0,
    .dZ[10] = 0.e0,
    .dZ[11] = 0.e0,
#endif // USE_METAL_DIFFUSION //}
#endif // USE_CELIB // }

};

/*
 * This function evaluates the force (m_i dv_i/dt) and the time derivative of the
 * internal energy (m_i du_i/dt).
 */
static struct StructHydroAccImport ReturnStructureDuDtHydroAccEngine(const int index, const int mode,
    struct StructHydroAccExport HydroAccExportRecv[restrict]){

    double Pos[3] = {PhydroPosP(index)[0],PhydroPosP(index)[1],PhydroPosP(index)[2]};
    double Vel[3] = {Phydro[index]->VelP[0],Phydro[index]->VelP[1],Phydro[index]->VelP[2]};
    double Kerneli = Phydro[index]->KernelPred;
    double Rhoi = Phydro[index]->RhoPred;
    double Fi = Phydro[index]->F;
    double Massi = Phydro[index]->Mass;
#ifdef USE_VARIABLE_ALPHA //{
    double Alphai = Phydro[index]->Alpha;
#else // USE_VARIABLE_ALPHA //}//{
    double Alphai = Pall.HydroAlpha;
#endif // USE_VARIABLE_ALPHA //}

#ifdef USE_DISPH //{
    double p_over_fund2_i = Pall.Gm1/Phydro[index]->EnergyDensityPred;
    double Weighti = Massi*Phydro[index]->UPred;
    double Volumei = Weighti/Phydro[index]->EnergyDensityPred;
#elif defined(USE_SPSPH) //}//{
    // P = (g-1)*\rho u
    // P = (g-1)*(my/z) *u
    // P/y^2 = (g-1)*m*u/(z*y)
    double p_over_fund2_i = (Pall.Gm1*Massi*Phydro[index]->UPred/(Phydro[index]->PseudoDensityPred*Phydro[index]->ZwPred));
    double Weighti = Phydro[index]->ZwPred;
    double PseudoDensityi = Phydro[index]->PseudoDensityPred;
    double Ddifi = Phydro[index]->Ddif;
    double Volumei = Weighti/PseudoDensityi;
#else // USE_DISPH //}//{
#if defined(ISOTHERMAL_EOS_RUN)
    double p_over_fund2_i = (SQ(Pall.CS)/Phydro[index]->RhoPred);
#elif defined(BAROTROPIC_EOS_RUN)
    double p_over_fund2_i = ReturnPoverRho2ForBarotropicRun(Phydro[index]->RhoPred);
#else
    double p_over_fund2_i = Pall.Gm1*(Phydro[index]->UPred/Phydro[index]->RhoPred);
#endif
    double Weighti = Massi;
    double Volumei = Weighti/Phydro[index]->RhoPred;
#endif // USE_DISPH //}
    double InvWeighti = 1.0/Weighti;

#ifdef USE_GRAD_H //{
    double Gradhi = Phydro[index]->Gradh;
#ifdef USE_GRAD_N //{
    double fi = Phydro[index]->fij;
#endif // USE_GRAD_N //}
#endif // USE_GRAD_H //}

    double InvKerneli = 1.e0/Kerneli;
#if defined(ISOTHERMAL_EOS_RUN)
    double csi = Pall.CS;
#elif defined(BAROTROPIC_EOS_RUN)
    double csi = ReturnSoundSpeedForBarotoropicRun(Phydro[index]->RhoPred);
#else
    double csi = sqrt(Pall.GGm1*Phydro[index]->UPred);
#endif

#ifdef USE_CELIB //{
#ifdef USE_METAL_DIFFUSION //{
#if DIFFUSION_TYPE==0 //{
    double DZdiff = Phydro[index]->DZdiff;
#elif DIFFUSION_TYPE==1 //}/{
    double Sxy = Phydro[index]->Sxy;
    double Sxz = Phydro[index]->Sxz;
    double Syx = Phydro[index]->Syx;
    double Syz = Phydro[index]->Syz;
    double Szx = Phydro[index]->Szx;
    double Szy = Phydro[index]->Szy;
#endif  // DIFFUSION_TYPE //}
    double Z[CELibYield_Number-1];
    for(int l=0;l<CELibYield_Number-1;l++){
        Z[l] = Phydro[index]->Elements[l+1]/Massi;
    }
#endif // USE_METAL_DIFFUSION //}
#endif // USE_CELIB //{

#if VISCOSITY_TYPE == 1 //{ 
    double Divi = Phydro[index]->Div;
#endif // VISCOSITY_TYPE == 1 //}

#ifdef USE_THERMAL_CONDUCTIVITY //{
    double Temperaturei = Pall.ConvertUtoT*Phydro[index]->UPred;
    double Kappai = CalcThermalConductivity(index);
    double InvGamma = 1.0/Pall.Gamma;
#endif // USE_THERMAL_CONDUCTIVITY //}

#if (defined(USE_SPSPH) || defined(USE_THERMAL_CONDUCTIVITY) || defined(USE_METAL_DIFFUSION)) //{
    double dti = Phydro[index]->dt_hydro;
#endif //(defined(USE_SPSPH) || defined(USE_THERMAL_CONDUCTIVITY) || defined(USE_METAL_DIFFUSION)) //}

    int Neighbors[MaxNeighborSize];
    struct StructHydroAccImport TempHydroAccImport = InitializeHydroAccImport;

    double TimeNeighborSearch = GetElapsedTime();
    int Nlist;
    if(mode == 0){
        Nlist = GetNeighborsPairsLimited(Pos,2.e0*Kerneli,Neighbors);
    } else {
        Nlist = GetNeighborsPairsLimitedImported(Pos,2.e0*Kerneli,Neighbors);
    }
    TimingResults.HydroAccNeighborSearchThisStep += GetElapsedTime()-TimeNeighborSearch;


    struct StructHydroAccPreFetch{
        double Mass;
        double Kernel;
        double Rho;
        double F;
#ifdef USE_VARIABLE_ALPHA //{
        double Alpha;
#endif // USE_VARIABLE_ALPHA //}
        double p_over_fund2_j;
        double Weightj;
        double Volumej;
        double csj;
#ifdef USE_SPSPH //{
        double Ddifj;
        double PseudoDensityj;
#endif // USE_SPSPH //}
#ifdef USE_GRAD_H //{
        double Gradh;
#ifdef USE_GRAD_N //{
        double fj;
#endif // USE_GRAD_N //}
#endif // USE_GRAD_H //}
#ifdef USE_THERMAL_CONDUCTIVITY //{
        double Kappa;
        double Temperature;
#endif // USE_THERMAL_CONDUCTIVITY //}
        double ActiveFactor; // Use for diffusion terms.
        double Pos[3];
        double Vel[3];
#ifdef USE_CELIB //{
#ifdef USE_METAL_DIFFUSION //{
#if DIFFUSION_TYPE==0 //{
        double DZdiff;
#elif DIFFUSION_TYPE==1 //}/{
        double Sxy;
        double Sxz;
        double Syx;
        double Syz;
        double Szx;
        double Szy;
#endif  // DIFFUSION_TYPE //}
        double Z[CELibYield_Number];
        double dt;
#endif // USE_METAL_DIFFUSION //}
#endif // USE_CELIB //}
#if VISCOSITY_TYPE == 1 //{
        double Div;
#endif // VISCOSITY_TYPE == 1 //}
    } HydroAccPreFetch[Nlist];

    for(int k=0;k<Nlist;k++){
        int leaf = Neighbors[k];
        if(mode == 0){
            HydroAccPreFetch[k].Mass = Phydro[leaf]->Mass;
            HydroAccPreFetch[k].Pos[0] = PhydroPosP(leaf)[0];
            HydroAccPreFetch[k].Pos[1] = PhydroPosP(leaf)[1];
            HydroAccPreFetch[k].Pos[2] = PhydroPosP(leaf)[2];
            HydroAccPreFetch[k].Vel[0] = Phydro[leaf]->VelP[0];
            HydroAccPreFetch[k].Vel[1] = Phydro[leaf]->VelP[1];
            HydroAccPreFetch[k].Vel[2] = Phydro[leaf]->VelP[2];
            HydroAccPreFetch[k].Rho = Phydro[leaf]->RhoPred;
            HydroAccPreFetch[k].Kernel = Phydro[leaf]->KernelPred;

#ifdef USE_DISPH //{
            HydroAccPreFetch[k].p_over_fund2_j = Pall.Gm1/Phydro[leaf]->EnergyDensityPred;
            HydroAccPreFetch[k].Weightj = Phydro[leaf]->Mass*Phydro[leaf]->UPred;
            HydroAccPreFetch[k].Volumej = HydroAccPreFetch[k].Weightj/Phydro[leaf]->EnergyDensityPred;
#elif defined(USE_SPSPH) //}//{
            HydroAccPreFetch[k].p_over_fund2_j = Pall.Gm1*Phydro[leaf]->Mass*Phydro[leaf]->UPred/(Phydro[leaf]->PseudoDensityPred*Phydro[leaf]->ZwPred);
            HydroAccPreFetch[k].Weightj = Phydro[leaf]->ZwPred;
            HydroAccPreFetch[k].PseudoDensityj = Phydro[leaf]->PseudoDensityPred;
            HydroAccPreFetch[k].Ddifj = Phydro[leaf]->Ddif;
            HydroAccPreFetch[k].Volumej = HydroAccPreFetch[k].Weightj/HydroAccPreFetch[k].PseudoDensityj;
#else // USE_DISPH //}//{
#if defined(ISOTHERMAL_EOS_RUN) //{
            HydroAccPreFetch[k].p_over_fund2_j = (SQ(Pall.CS)/Phydro[leaf]->RhoPred);
#elif defined(BAROTROPIC_EOS_RUN) //{//}
            HydroAccPreFetch[k].p_over_fund2_j = ReturnPoverRho2ForBarotropicRun(Phydro[leaf]->RhoPred);
#else //{//}
            HydroAccPreFetch[k].p_over_fund2_j = Pall.Gm1*(Phydro[leaf]->UPred/Phydro[leaf]->RhoPred);
#endif //}
            HydroAccPreFetch[k].Weightj = Phydro[leaf]->Mass;
            HydroAccPreFetch[k].Volumej = HydroAccPreFetch[k].Weightj/Phydro[leaf]->RhoPred;
#endif // USE_DISPH //}

#if defined(ISOTHERMAL_EOS_RUN) //{
            HydroAccPreFetch[k].csj = Pall.CS;
#elif defined(BAROTROPIC_EOS_RUN) //}//{
            HydroAccPreFetch[k].csj = ReturnSoundSpeedForBarotoropicRun(Phydro[leaf]->RhoPred);
#else //}//{
            HydroAccPreFetch[k].csj = sqrt(Pall.GGm1*Phydro[leaf]->UPred);
#endif //}

#ifdef USE_GRAD_H //{
            HydroAccPreFetch[k].Gradh = Phydro[leaf]->Gradh;
#ifdef USE_GRAD_N //{
            HydroAccPreFetch[k].fj = Phydro[leaf]->fij;
#endif // USE_GRAD_N //}
#endif // USE_GRAD_H //}

            HydroAccPreFetch[k].F = Phydro[leaf]->F;
#ifdef USE_VARIABLE_ALPHA //{
            HydroAccPreFetch[k].Alpha = Phydro[leaf]->Alpha;
#endif // USE_VARIABLE_ALPHA //}
#ifdef USE_THERMAL_CONDUCTIVITY //{
            HydroAccPreFetch[k].Temperature = Pall.ConvertUtoT*Phydro[leaf]->UPred;
            HydroAccPreFetch[k].Kappa = CalcThermalConductivity(leaf);
#endif // USE_THERMAL_CONDUCTIVITY //}
#ifdef USE_CELIB //{
#ifdef USE_METAL_DIFFUSION //{
#if DIFFUSION_TYPE==0 //{
            HydroAccPreFetch[k].DZdiff = Phydro[leaf]->DZdiff;
#elif DIFFUSION_TYPE==1 //}/{
            HydroAccPreFetch[k].Sxy = Phydro[leaf]->Sxy;
            HydroAccPreFetch[k].Sxz = Phydro[leaf]->Sxz;
            HydroAccPreFetch[k].Syx = Phydro[leaf]->Syx;
            HydroAccPreFetch[k].Syz = Phydro[leaf]->Syz;
            HydroAccPreFetch[k].Szx = Phydro[leaf]->Szx;
            HydroAccPreFetch[k].Szy = Phydro[leaf]->Szy;
#endif  // DIFFUSION_TYPE //}
            for(int l=0;l<CELibYield_Number-1;l++){
                HydroAccPreFetch[k].Z[l] = Phydro[leaf]->Elements[l+1]/HydroAccPreFetch[k].Mass;
            }
            HydroAccPreFetch[k].dt = Phydro[leaf]->dt_hydro;
#endif // USE_METAL_DIFFUSION //}
#endif // USE_CELIB //{
#if VISCOSITY_TYPE == 1 //{
            HydroAccPreFetch[k].Div = Phydro[leaf]->Div;
#endif // VISCOSITY_TYPE == 1 //}
#if (defined(USE_SPSPH) || defined(USE_THERMAL_CONDUCTIVITY) || defined(USE_METAL_DIFFUSION)) //{
            if((Phydro[leaf]->Active)&&(Phydro[index]->dt_hydro>0.e0)){
                HydroAccPreFetch[k].ActiveFactor = fmax(Phydro[leaf]->dt_hydro,dti)/dti;
            } else {
                HydroAccPreFetch[k].ActiveFactor = 0.e0;
            }
#endif // (defined(USE_SPSPH) || defined(USE_THERMAL_CONDUCTIVITY) || defined(USE_METAL_DIFFUSION)) //}
        } else {
            HydroAccPreFetch[k].Mass = HydroAccExportRecv[leaf].Mass;
            HydroAccPreFetch[k].Pos[0] = HydroAccExportRecv[leaf].Pos[0];
            HydroAccPreFetch[k].Pos[1] = HydroAccExportRecv[leaf].Pos[1];
            HydroAccPreFetch[k].Pos[2] = HydroAccExportRecv[leaf].Pos[2];
            HydroAccPreFetch[k].Vel[0] = HydroAccExportRecv[leaf].Vel[0];
            HydroAccPreFetch[k].Vel[1] = HydroAccExportRecv[leaf].Vel[1];
            HydroAccPreFetch[k].Vel[2] = HydroAccExportRecv[leaf].Vel[2];
            HydroAccPreFetch[k].Rho = HydroAccExportRecv[leaf].Rho;
            HydroAccPreFetch[k].Kernel = HydroAccExportRecv[leaf].Kernel;

#ifdef USE_DISPH //{
            HydroAccPreFetch[k].p_over_fund2_j = Pall.Gm1/HydroAccExportRecv[leaf].EnergyDensity;
            HydroAccPreFetch[k].Weightj = HydroAccExportRecv[leaf].Mass*HydroAccExportRecv[leaf].U;
            HydroAccPreFetch[k].Volumej = HydroAccPreFetch[k].Weightj/HydroAccExportRecv[leaf].EnergyDensity;
#elif defined(USE_SPSPH) //}//{
            HydroAccPreFetch[k].p_over_fund2_j = Pall.Gm1*HydroAccExportRecv[leaf].Mass*HydroAccExportRecv[leaf].U
                                                /(HydroAccExportRecv[leaf].PseudoDensity*HydroAccExportRecv[leaf].Zw);
            HydroAccPreFetch[k].Weightj = HydroAccExportRecv[leaf].Zw;
            HydroAccPreFetch[k].PseudoDensityj = HydroAccExportRecv[leaf].PseudoDensity;
            HydroAccPreFetch[k].Ddifj = HydroAccExportRecv[leaf].Ddif;
            HydroAccPreFetch[k].Volumej = HydroAccPreFetch[k].Weightj/HydroAccPreFetch[k].PseudoDensityj;
#else // USE_DISPH //}//{
#if defined(ISOTHERMAL_EOS_RUN) //{
            HydroAccPreFetch[k].p_over_fund2_j = (SQ(Pall.CS)/HydroAccExportRecv[leaf].Rho);
#elif defined(BAROTROPIC_EOS_RUN) //{//}
            oydroAccPreFetch[k].p_over_fund2_j = ReturnPoverRho2ForBarotropicRun(HydroAccExportRecv[leaf].Rho);
#else //{//}
            HydroAccPreFetch[k].p_over_fund2_j = Pall.Gm1*(HydroAccExportRecv[leaf].U/HydroAccExportRecv[leaf].Rho);
#endif //}
            HydroAccPreFetch[k].Weightj = HydroAccExportRecv[leaf].Mass;
            HydroAccPreFetch[k].Volumej = HydroAccPreFetch[k].Weightj/HydroAccExportRecv[leaf].Rho;
#endif // USE_DISPH //}

#if defined(ISOTHERMAL_EOS_RUN) //{
            HydroAccPreFetch[k].csj = Pall.CS;
#elif defined(BAROTROPIC_EOS_RUN) //}//{
            HydroAccPreFetch[k].csj = ReturnSoundSpeedForBarotoropicRun(HydroAccExportRecv[leaf].Rho);
#else //}//{
            HydroAccPreFetch[k].csj = sqrt(Pall.GGm1*HydroAccExportRecv[leaf].U);
#endif //}

#ifdef USE_GRAD_H //{
            HydroAccPreFetch[k].Gradh = HydroAccExportRecv[leaf].Gradh;
#ifdef USE_GRAD_N //{
            HydroAccPreFetch[k].fj = HydroAccExportRecv[leaf].fij;
#endif // USE_GRAD_N //}
#endif // USE_GRAD_H //}

            HydroAccPreFetch[k].F = HydroAccExportRecv[leaf].F;
#ifdef USE_VARIABLE_ALPHA //{
            HydroAccPreFetch[k].Alpha = HydroAccExportRecv[leaf].Alpha;
#endif // USE_VARIABLE_ALPHA //}
#ifdef USE_THERMAL_CONDUCTIVITY //{
            HydroAccPreFetch[k].Temperature = Pall.ConvertUtoT*HydroAccExportRecv[leaf].U;
            HydroAccPreFetch[k].Kappa = CalcThermalConductivityByTemperature(HydroAccPreFetch[k].Temperature);
#endif // USE_THERMAL_CONDUCTIVITY //}
#ifdef USE_CELIB //{
#ifdef USE_METAL_DIFFUSION //{
#if DIFFUSION_TYPE==0 //{
            HydroAccPreFetch[k].DZdiff = HydroAccExportRecv[leaf].DZdiff;
#elif DIFFUSION_TYPE==1 //}/{
            HydroAccPreFetch[k].Sxy = HydroAccExportRecv[leaf].Sxy;
            HydroAccPreFetch[k].Sxz = HydroAccExportRecv[leaf].Sxz;
            HydroAccPreFetch[k].Syx = HydroAccExportRecv[leaf].Syx;
            HydroAccPreFetch[k].Syz = HydroAccExportRecv[leaf].Syz;
            HydroAccPreFetch[k].Szx = HydroAccExportRecv[leaf].Szx;
            HydroAccPreFetch[k].Szy = HydroAccExportRecv[leaf].Szy;
#endif  // DIFFUSION_TYPE //}
            for(int l=0;l<CELibYield_Number-1;l++){
                HydroAccPreFetch[k].Z[l] = HydroAccExportRecv[leaf].Elements[l]/HydroAccPreFetch[k].Mass;
            }
            HydroAccPreFetch[k].dt = HydroAccExportRecv[leaf].dt;
#endif // USE_METAL_DIFFUSION //}
#endif // USE_CELIB //{
#if VISCOSITY_TYPE == 1 //{
            HydroAccPreFetch[k].Div = HydroAccExportRecv[leaf].Div;
#endif // VISCOSITY_TYPE == 1 //}
#if (defined(USE_SPSPH) || defined(USE_THERMAL_CONDUCTIVITY) || defined(USE_METAL_DIFFUSION)) //{
            if((HydroAccExportRecv[leaf].Active)&&(Phydro[index]->dt_hydro>0.e0)){
                HydroAccPreFetch[k].ActiveFactor = fmax(HydroAccExportRecv[leaf].dt,dti)/dti;
            } else {
                HydroAccPreFetch[k].ActiveFactor = 0.e0;
            }
#endif // (defined(USE_SPSPH) || defined(USE_THERMAL_CONDUCTIVITY) || defined(USE_METAL_DIFFUSION)) //}
        }
    }

    int NlistActive = 0;
	for(int k=0;k<Nlist;k++){
        double xij[3];
#ifdef PERIODIC_RUN //{
        xij[0] = PeriodicDistance(Pos[0],HydroAccPreFetch[k].Pos[0],0);
        xij[1] = PeriodicDistance(Pos[1],HydroAccPreFetch[k].Pos[1],1);
        xij[2] = PeriodicDistance(Pos[2],HydroAccPreFetch[k].Pos[2],2);
#else // PERIODIC_RUN //}//{
        xij[0] = Pos[0]-HydroAccPreFetch[k].Pos[0];
        xij[1] = Pos[1]-HydroAccPreFetch[k].Pos[1];
        xij[2] = Pos[2]-HydroAccPreFetch[k].Pos[2];
#endif // PERIODIC_RUN //}
        double r2 = NORM2(xij);

        if(r2 == 0.e0) continue;

        double r = sqrt(r2);

        double Kernelj = HydroAccPreFetch[k].Kernel;
        double InvKernelj = 1.e0/Kernelj;

        double dwi = dSPHKernel(r,InvKerneli);
        double dwj = dSPHKernel(r,InvKernelj);
        double dw = 0.5*(dwi+dwj);

        double vij[3];
        vij[0] = Vel[0]-HydroAccPreFetch[k].Vel[0];
        vij[1] = Vel[1]-HydroAccPreFetch[k].Vel[1];
        vij[2] = Vel[2]-HydroAccPreFetch[k].Vel[2];

        double vx = DOT_PRODUCT(xij,vij);

#if VISCOSITY_TYPE == 0 //{
#if (HydroViscNormal) //{
        double visc = ReturnArtificialViscocityGM(vx,r2,csi,HydroAccPreFetch[k].csj,
            Rhoi,HydroAccPreFetch[k].Rho,Kerneli,HydroAccPreFetch[k].Kernel,
            Fi,HydroAccPreFetch[k].F
            SetParameterVariableAlpha(Alphai,HydroAccPreFetch[k].Alpha))*dw;
#elif (HydroViscVelSig) // Visc //}//{
        double visc = ReturnArtificialViscocityReimann(vx,r,csi,HydroAccPreFetch[k].csj,
            Rhoi,HydroAccPreFetch[k].Rho,Fi,HydroAccPreFetch[k].F,&(TempHydroAccImport.Vsig)
            SetParameterVariableAlpha(Alphai,HydroAccPreFetch[k].Alpha))*dw;
#endif // Visc //}

        double tmpj_visc = Massi*HydroAccPreFetch[k].Mass*visc;

#ifdef USE_GRAD_H //{
#ifdef USE_GRAD_N //{
        double fij = 1.0-fi/HydroAccPreFetch[k].Weightj;
        double fji = 1.0-HydroAccPreFetch[k].fj*InvWeighti;
#else //USE_GRAD_N //}//{
        double fij = Gradhi;
        double fji = HydroAccPreFetch[k].Gradh;
#endif //USE_GRAD_N //}

        double tmpj = Weighti*HydroAccPreFetch[k].Weightj
                      *(fij*p_over_fund2_i*dwi+fji*HydroAccPreFetch[k].p_over_fund2_j*dwj)
                     +tmpj_visc;
        double tmpj_energy = Weighti*HydroAccPreFetch[k].Weightj*fij*p_over_fund2_i*vx*dwi
                                +0.5*tmpj_visc*vx;
#else // USE_GRAD_H //}//{
        double tmpj = Weighti*HydroAccPreFetch[k].Weightj
                      *(p_over_fund2_i+HydroAccPreFetch[k].p_over_fund2_j)*dw
                     +tmpj_visc;
        double tmpj_energy = (Weighti*HydroAccPreFetch[k].Weightj*p_over_fund2_i)*vx*dw
                                +0.5*tmpj_visc*vx;
#endif //USE_GRAD_H //}

        TempHydroAccImport.Du += tmpj_energy;
        TempHydroAccImport.HydroAcc[0] -= tmpj*xij[0];
        TempHydroAccImport.HydroAcc[1] -= tmpj*xij[1];
        TempHydroAccImport.HydroAcc[2] -= tmpj*xij[2];

#elif VISCOSITY_TYPE == 1 //} //{

        double visc_i = ReturnArtificialViscocityHosono(Divi,csi,Rhoi,Kerneli,Fi
                            SetParameterVariableAlphaHosono(Alphai));
        double visc_j = ReturnArtificialViscocityHosono(HydroAccPreFetch[k].Div,
                HydroAccPreFetch[k].csj,HydroAccPreFetch[k].Rho,HydroAccPreFetch[k].Kernel,HydroAccPreFetch[k].F
                SetParameterVariableAlphaHosono(HydroAccPreFetch[k].Alpha));

        double vsig = csi+HydroAccPreFetch[k].csj;
        if(vx < 0.e0){ 
            double wij = vx/r;
            vsig += -SIGNAL_VELOCITY_BETA*wij; 
        }
        TempHydroAccImport.Vsig = fmax(vsig,TempHydroAccImport.Vsig);

#ifdef USE_GRAD_H //{
#ifdef USE_GRAD_N //{
        double fij = 1.0-fi/HydroAccPreFetch[k].Weightj;
        double fji = 1.0-HydroAccPreFetch[k].fj*InvWeighti;
#else //USE_GRAD_N //}//{
        double fij = Gradhi;
        double fji = HydroAccPreFetch[k].Gradh;
#endif //USE_GRAD_N //}

        double tmpj = Weighti*HydroAccPreFetch[k].Weightj
                      *(fij*p_over_fund2_i*dwi+fji*HydroAccPreFetch[k].p_over_fund2_j*dwj);
        double tmpj_energy = Weighti*HydroAccPreFetch[k].Weightj*fij*p_over_fund2_i*vx*dwi;
#else // USE_GRAD_H //}//{
        double tmpj = Weighti*HydroAccPreFetch[k].Weightj
                      *(p_over_fund2_i+HydroAccPreFetch[k].p_over_fund2_j)*dw;
        double tmpj_energy = (Weighti*HydroAccPreFetch[k].Weightj*p_over_fund2_i)*vx*dw;
#endif //USE_GRAD_H //}

        double tmpj_visc = Weighti*HydroAccPreFetch[k].Weightj*
            (visc_i*fij/SQ(Weighti)*dwi+visc_j*fji/SQ(HydroAccPreFetch[k].Weightj)*dwj);
        double tmpj_energy_visc = +visc_i*(Volumei*HydroAccPreFetch[k].Weightj)*fij*vx*dwi;

        
        TempHydroAccImport.Du += tmpj_energy+tmpj_energy_visc;
        TempHydroAccImport.HydroAcc[0] -= (tmpj+tmpj_visc)*xij[0];
        TempHydroAccImport.HydroAcc[1] -= (tmpj+tmpj_visc)*xij[1];
        TempHydroAccImport.HydroAcc[2] -= (tmpj+tmpj_visc)*xij[2];
#endif // VISCOSITY_TYPE == 1 //{ 

#ifdef USE_SPSPH //{
        TempHydroAccImport.DZw += HydroAccPreFetch[k].ActiveFactor*
            4.0*(Ddifi*HydroAccPreFetch[k].Ddifj)/(Ddifi+HydroAccPreFetch[k].Ddifj)
            *((Weighti*HydroAccPreFetch[k].Weightj)/(PseudoDensityi*HydroAccPreFetch[k].PseudoDensityj))
            *(PseudoDensityi-HydroAccPreFetch[k].PseudoDensityj)*dw;
#endif // USE_SPSPH //}
#ifdef USE_THERMAL_CONDUCTIVITY //{
        //double tmpj_tc = InvGamma*(4*Kappai*HydroAccPreFetch[k].Kappa/(Kappai+HydroAccPreFetch[k].Kappa))
        double tmpj_tc = (4*Kappai*HydroAccPreFetch[k].Kappa/(Kappai+HydroAccPreFetch[k].Kappa))
                    *(Volumei*HydroAccPreFetch[k].Volumej)*(Temperaturei-HydroAccPreFetch[k].Temperature)*dw;
        TempHydroAccImport.Du += HydroAccPreFetch[k].ActiveFactor*tmpj_tc;
#endif // USE_THERMAL_CONDUCTIVITY //}

#ifdef USE_CELIB //{
#ifdef USE_METAL_DIFFUSION //{

#if DIFFUSION_TYPE == 0 //{
        double DZidff_i = DZdiff;
        double DZidff_j = HydroAccPreFetch[k].DZdiff;
#elif DIFFUSION_TYPE == 1 //}//{
        double Sij = sqrt(2.0*(
                    SQ(Sxy+HydroAccPreFetch[k].Sxy)+SQ(Sxz+HydroAccPreFetch[k].Sxz)+
                    SQ(Sxy+HydroAccPreFetch[k].Syx)+SQ(Sxz+HydroAccPreFetch[k].Syz)+
                    SQ(Sxy+HydroAccPreFetch[k].Szx)+SQ(Sxz+HydroAccPreFetch[k].Szy)));
        double DZdiff_i = TURBULENT_DIFFUSION_COEF*Sij*SQ(Kerneli);
        if(dti*DZdiff_i>0.1*SQ(Kerneli)){
            DZdiff_i = 0.1*SQ(Kerneli)/dti;
        }
        double DZdiff_j = TURBULENT_DIFFUSION_COEF*Sij*SQ(HydroAccPreFetch[k].Kernel);
        if(HydroAccPreFetch[k].dt*DZdiff_j>0.1*SQ(HydroAccPreFetch[k].Kernel)){
            DZdiff_j = 0.1*SQ(HydroAccPreFetch[k].Kernel)/HydroAccPreFetch[k].dt;
        }
#endif // DIFFUSION_TYPE //}

#if 0
        if(DZdiff_i*DZdiff_j>0.0){
            //double MassFactor = 2.0*HydroAccPreFetch[k].Mass/(Massi+HydroAccPreFetch[k].Mass);
            double MassFactor = 1.0;
            double Fact = HydroAccPreFetch[k].ActiveFactor*MassFactor*
                4.0*(DZdiff_i*DZdiff_j)/(DZdiff_i+DZdiff_j)*(Volumei*HydroAccPreFetch[k].Volumej)*dw;
            for(int l=0;l<11;l++){
                TempHydroAccImport.dZ[l] += Fact*(Z[l]-HydroAccPreFetch[k].Z[l]);

                if(Z[l]-HydroAccPreFetch[k].Z[l]>0.0){
                    assert(Fact*(Z[l]-HydroAccPreFetch[k].Z[l])<0.e0);
                }
                if(Z[l]-HydroAccPreFetch[k].Z[l]<0.0){
                    assert(Fact*(Z[l]-HydroAccPreFetch[k].Z[l])>0.e0);
                }

                if(fabs(Z[l]-HydroAccPreFetch[k].Z[l])/Z[l]>1.e-2){
                    fprintf(stderr,"Something wrong, think different. %g, l %d\n",(Z[l]-HydroAccPreFetch[k].Z[l])/Z[l],l);
                    fprintf(stderr,"ZZ %g %g / II = %d %d\n",Z[l],HydroAccPreFetch[k].Z[l],index,Neighbors[k]);
                    fprintf(stderr,"Fact %g\n",Fact);
                    fflush(NULL);
                    assert(1==0);
                }
                    //HydroAccPreFetch[k].ActiveFactor*MassFactor*
                    //4.0*(DZdiff*HydroAccPreFetch[k].DZdiff)/(DZdiff+HydroAccPreFetch[k].DZdiff)
                        //*(Volumei*HydroAccPreFetch[k].Volumej)*(Z[l]-HydroAccPreFetch[k].Z[l])*dw;
            }
        }
#else
        if((DZdiff_i*DZdiff_j>0.0)&&(HydroAccPreFetch[k].ActiveFactor>0)){
            double Fact = HydroAccPreFetch[k].ActiveFactor*4.0*(DZdiff_i*DZdiff_j)/(DZdiff_i+DZdiff_j)*
                HydroAccPreFetch[k].Mass*(2.0/(Rhoi+HydroAccPreFetch[k].Rho))*dw;
            for(int l=0;l<CELibYield_Number-1;l++){
                double LocalFact = Fact;
                if(fabs(Z[l]-HydroAccPreFetch[k].Z[l])/(Z[l]+HydroAccPreFetch[k].Z[l])<1.e-10){
                    LocalFact = 0.e0;
                } else {
#if 0
                    gprintlmpi(Fact*(Z[l]-HydroAccPreFetch[k].Z[l]));
                    gprintlmpi(Fact);
                    gprintlmpi((Z[l]-HydroAccPreFetch[k].Z[l]));
                    gprintlmpi(HydroAccPreFetch[k].ActiveFactor);
                    gprintlmpi(DZdiff_i*DZdiff_j);
                    gprintlmpi(HydroAccPreFetch[k].Mass);
                    gprintlmpi(dw);
                    gprintlmpi(Rhoi+HydroAccPreFetch[k].Rho);
                    gprintlmpi(fabs(Z[l]-HydroAccPreFetch[k].Z[l])/(Z[l]+HydroAccPreFetch[k].Z[l]));
                    gprintlmpi(HydroAccPreFetch[k].ActiveFactor*4.0*(DZdiff_i*DZdiff_j)/(DZdiff_i+DZdiff_j)*
                            HydroAccPreFetch[k].Mass*(2.0/(Rhoi+HydroAccPreFetch[k].Rho))*dw);
#endif
                }
                TempHydroAccImport.dZ[l] += LocalFact*(Z[l]-HydroAccPreFetch[k].Z[l]);
#if 0
                if(Fact*(Z[l]-HydroAccPreFetch[k].Z[l])>0.e0){
                    fprintf(stderr,"%ld | %1.15g, %g %1.15g %1.15g\n",PhydroBody(index)->GlobalID,
                            Fact*(Z[l]-HydroAccPreFetch[k].Z[l]),
                            Fact,Z[l],HydroAccPreFetch[k].Z[l]);
                    fflush(NULL);
                }
#endif
            }
#endif // USE_METAL_DIFFUSION //}
        }
#endif

#endif // USE_METAL_DIFFUSION //}

        NlistActive ++; 
    }
    TempHydroAccImport.Nlist = NlistActive;

    return (TempHydroAccImport);
}



static bool NExportImportLogFirst = true;
static int *NExportLog;
static int *NImportLog;
static bool *ExportFlagsLog; 

void CalcDuDtAcc(void){

    double TimingResultThisRoutine = GetElapsedTime();

    int NProcs = MPIGetNumProcs();
    MPI_Status  mpi_status;

    if(NExportImportLogFirst == true){
        NExportLog = malloc(sizeof(int)*NProcs);
        NImportLog = malloc(sizeof(int)*NProcs);
        NExportImportLogFirst = false;
    }

    static int HydroAccExportFlagsMaxAllocated = 0;
    static bool (*HydroAccExportFlags)[NProcs-1];
    static int *ActiveIndexList;
    if(HydroAccExportFlagsMaxAllocated < MAX(Pall.Nhydro,NAdditionUnit)){
        if(HydroAccExportFlagsMaxAllocated > 0){
            free(ActiveIndexList);
            free(HydroAccExportFlags);
            free(ExportFlagsLog);
        }
        HydroAccExportFlagsMaxAllocated = (int)(MAX(ForAngelsShare*Pall.Nhydro,NAdditionUnit));
        HydroAccExportFlags = malloc(sizeof(bool)*HydroAccExportFlagsMaxAllocated*(NProcs-1)+1);
        ExportFlagsLog = malloc(sizeof(bool)*HydroAccExportFlagsMaxAllocated*(NProcs-1)+1);
        ActiveIndexList = malloc(sizeof(int)*HydroAccExportFlagsMaxAllocated);
    }

    int NActives = 0;
    int RootNodeID = 0;
    int NumberofLeaves = HydroNode[RootNodeID].NumberofLeaves;
    int header = HydroNode[RootNodeID].Leaves;
    for(int i=0;i<NumberofLeaves;i++){
        int leaf = header + i;
        //if(HydroRoot.Leaves[leaf] < 0) continue;
        if(NBCache[leaf].Active){
            ActiveIndexList[NActives] = NBCache[leaf].Leaf;
            NActives ++;
        }
    }

    memset(HydroAccExportFlags,0,Pall.Nhydro*(NProcs-1)*sizeof(bool));

    int NExportThisTime[NProcs];
    int NImportThisTime[NProcs];

    struct StructHydroAccExport *HydroAccExportSend[NProcs];
    //bool (*HydroInteractionFlags)[NProcs-1];
    bool (*HydroInteractionFlags)[Pall.Nhydro];
    HydroInteractionFlags = (void *)BufferHydroInteractionFlags;

    bool ExportFlags[Pall.Nhydro];


    int NExportAll = 0;
    for(int i=0;i<NProcs-1;i++){
 
        memset(ExportFlags,0,Pall.Nhydro*sizeof(bool));
        NExportThisTime[i] = CheckHydroAccExportFlagsModified(i,NProcs,ExportFlags,HydroInteractionFlags[i]);
        NExportLog[i] = NExportThisTime[i];

        CheckSizeofBufferExportSendIndex(NExportThisTime[i],sizeof(struct StructHydroAccExport),i);
        HydroAccExportSend[i] = BufferExportSend[i];

        int Offset = i*Pall.Nhydro;
        int NExport = 0;
        if(NExportThisTime[i] > 0){
        for(int k=0;k<Pall.Nhydro;k++){
            ExportFlagsLog[Offset+k] = ExportFlags[k];
            if(ExportFlags[k]){
                HydroAccExportSend[i][NExport].Mass = Phydro[k]->Mass;
                HydroAccExportSend[i][NExport].Kernel = Phydro[k]->KernelPred;
                HydroAccExportSend[i][NExport].Pos[0] = PhydroPosP(k)[0];
                HydroAccExportSend[i][NExport].Pos[1] = PhydroPosP(k)[1];
                HydroAccExportSend[i][NExport].Pos[2] = PhydroPosP(k)[2];
                HydroAccExportSend[i][NExport].Vel[0] = Phydro[k]->VelP[0];
                HydroAccExportSend[i][NExport].Vel[1] = Phydro[k]->VelP[1];
                HydroAccExportSend[i][NExport].Vel[2] = Phydro[k]->VelP[2];
                HydroAccExportSend[i][NExport].Rho = Phydro[k]->RhoPred;
                HydroAccExportSend[i][NExport].U = Phydro[k]->UPred;
                HydroAccExportSend[i][NExport].F = Phydro[k]->F;
#ifdef USE_VARIABLE_ALPHA //{
                HydroAccExportSend[i][NExport].Alpha = Phydro[k]->Alpha;
#endif // USE_VARIABLE_ALPHA //}
#ifdef USE_DISPH //{
                HydroAccExportSend[i][NExport].EnergyDensity = Phydro[k]->EnergyDensityPred;
#endif // USE_DISPH //}
#ifdef USE_SPSPH //{
                HydroAccExportSend[i][NExport].PseudoDensity = Phydro[k]->PseudoDensityPred;
                HydroAccExportSend[i][NExport].Zw = Phydro[k]->ZwPred;
                HydroAccExportSend[i][NExport].Ddif = Phydro[k]->Ddif;
#endif // USE_SPSPH //}
#ifdef USE_GRAD_H //{
                HydroAccExportSend[i][NExport].Gradh = Phydro[k]->Gradh;
#ifdef USE_GRAD_N //{
                HydroAccExportSend[i][NExport].fij = Phydro[k]->fij;
#endif // USE_GRAD_N //}
#endif // USE_GRAD_H //}
#if (defined(USE_SPSPH) || defined(USE_THERMAL_CONDUCTIVITY) || defined(USE_METAL_DIFFUSION)) //{
                HydroAccExportSend[i][NExport].Active = Phydro[k]->Active;
                HydroAccExportSend[i][NExport].dt = Phydro[k]->dt_hydro;
#endif //(defined(USE_SPSPH) || defined(USE_THERMAL_CONDUCTIVITY) || defined(USE_METAL_DIFFUSION)) //}
                HydroAccExportSend[i][NExport].Leaf = k;

#ifdef USE_CELIB //{
#if USE_METAL_DIFFUSION //{
#if DIFFUSION_TYPE==0 //{
                HydroAccExportSend[i][NExport].DZdiff = Phydro[k]->DZdiff;
#elif DIFFUSION_TYPE==1 //}/{
                HydroAccExportSend[i][NExport].Sxy = Phydro[k]->Sxy;
                HydroAccExportSend[i][NExport].Sxz = Phydro[k]->Sxz;
                HydroAccExportSend[i][NExport].Syx = Phydro[k]->Syx;
                HydroAccExportSend[i][NExport].Syz = Phydro[k]->Syz;
                HydroAccExportSend[i][NExport].Szx = Phydro[k]->Szx;
                HydroAccExportSend[i][NExport].Szy = Phydro[k]->Szy;
#endif  // DIFFUSION_TYPE //}
                for(int l=0;l<CELibYield_Number-1;l++)
                    HydroAccExportSend[i][NExport].Elements[l] = Phydro[k]->Elements[l+1];
#endif // USE_METAL_DIFFUSION //}
#endif // USE_CELIB //}
#if VISCOSITY_TYPE == 1 //{
                HydroAccExportSend[i][NExport].Div = Phydro[k]->Div;
#endif // VISCOSITY_TYPE == 1 //}
                NExport ++;
            }
        }
        }
        NExportAll += NExport;
    }

    double TimeImbarance = GetElapsedTime();
    int NImportAll = 0;
    int NImportThisTime2[NProcs];
    int NExportThisTime2[NProcs];
    NImportThisTime2[MPIGetMyID()] = 0;
    for(int i=0;i<NProcs-1;i++){
        NExportThisTime2[CommunicationTable[i].SendRank] = NExportThisTime[i];
    }
    MPI_Alltoall(NExportThisTime2,1,MPI_INT,NImportThisTime2,1,MPI_INT,MPI_COMM_WORLD);
    for(int i=0;i<NProcs-1;i++){
        NImportThisTime[i] = NImportThisTime2[CommunicationTable[i].RecvRank];
        NImportAll += NImportThisTime[i];
        NImportLog[i] = NImportThisTime[i];
    }
    TimingResults.HydroImbaranceThisStep += GetElapsedTime() - TimeImbarance;

    struct StructHydroAccExport (*HydroAccExportRecv);
    CheckSizeofBufferExportRecv(NImportAll,sizeof(struct StructHydroAccExport));
    HydroAccExportRecv = BufferExportRecv;

    MPI_Status mpi_status_Export_Send[NProcs-1];
    MPI_Request mpi_request_Export_Send[NProcs-1];
    MPI_Status mpi_status_Export_Recv[NProcs-1];
    MPI_Request mpi_request_Export_Recv[NProcs-1];

    int NImport = 0;
    int counter_send = 0;
    int counter_recv = 0;
    int SendFlag,RecvFlag;
    for(int i=0;i<NProcs-1;i++){
        if(NExportThisTime[i]>0){
            MPI_Isend(HydroAccExportSend[i],
                NExportThisTime[i]*sizeof(struct StructHydroAccExport),
                    MPI_BYTE,CommunicationTable[i].SendRank,TAG_SPH_ACC_EXPORT+i,
                        MPI_COMM_WORLD,mpi_request_Export_Send+counter_send);
            MPI_Test(mpi_request_Export_Send+counter_send,&SendFlag,MPI_STATUS_IGNORE);
            counter_send ++;
        }
        if(NImportThisTime[i]>0){
            MPI_Irecv(HydroAccExportRecv+NImport,
                NImportThisTime[i]*sizeof(struct StructHydroAccExport),
                    MPI_BYTE,CommunicationTable[i].RecvRank,TAG_SPH_ACC_EXPORT+i,
                        MPI_COMM_WORLD,mpi_request_Export_Recv+counter_recv);
            MPI_Test(mpi_request_Export_Recv+counter_recv,&RecvFlag,MPI_STATUS_IGNORE);
            counter_recv ++;
        }
        NImport += NImportThisTime[i];
    }

    double TimeComm = GetElapsedTime();
    MPI_Waitall(counter_send,mpi_request_Export_Send,mpi_status_Export_Send);
    MPI_Waitall(counter_recv,mpi_request_Export_Recv,mpi_status_Export_Recv);
    TimingResults.HydroAccCommThisStep += GetElapsedTime()-TimeComm;

    if(NImportAll>0)
        PlantHydroTreeImported(NImportAll,HydroAccExportRecv);

    for(int k=0;k<NActives;k++){ // calculation for local
        struct StructHydroAccImport TempHydroAccImport;
        int leaf = ActiveIndexList[k];

        TempHydroAccImport = ReturnStructureDuDtHydroAccEngine(leaf,0,HydroAccExportRecv);

        Phydro[leaf]->Du = TempHydroAccImport.Du;
        Phydro[leaf]->HydroAcc[0] = TempHydroAccImport.HydroAcc[0];
        Phydro[leaf]->HydroAcc[1] = TempHydroAccImport.HydroAcc[1];
        Phydro[leaf]->HydroAcc[2] = TempHydroAccImport.HydroAcc[2];
        Phydro[leaf]->Vsig = TempHydroAccImport.Vsig;
        Phydro[leaf]->Nlist = TempHydroAccImport.Nlist;
#ifdef USE_SPSPH //{
        Phydro[leaf]->DZw = TempHydroAccImport.DZw;
#endif // USE_SPSPH //}
#ifdef USE_CELIB //{
#ifdef USE_METAL_DIFFUSION //{
        double Z[CELibYield_Number-1];
        for(int l=0;l<CELibYield_Number-1;l++){
            Z[l] = Phydro[leaf]->Elements[l+1]/Phydro[leaf]->Mass + TempHydroAccImport.dZ[l]*Phydro[leaf]->dt_hydro;
#if 0
            if(TempHydroAccImport.dZ[l]*Phydro[leaf]->dt_hydro/Z[l] > 1.e-4){
                fprintf(stderr,"Something wrong, think different. %g %d\n",TempHydroAccImport.dZ[l],l);
                fprintf(stderr,"IDs %d %d\n",leaf,k);
                fflush(NULL);
                assert(1==0);
            }
#endif
        }
#endif // USE_METAL_DIFFUSION //}
#endif // USE_CELIB //}

        if(NImportAll>0){
            //TempHydroAccImport = ReturnStructureDuDtHydroAccEngineImported(leaf,HydroAccExportRecv);
            TempHydroAccImport = ReturnStructureDuDtHydroAccEngine(leaf,1,HydroAccExportRecv);

            Phydro[leaf]->Du += TempHydroAccImport.Du;
            Phydro[leaf]->HydroAcc[0] += TempHydroAccImport.HydroAcc[0];
            Phydro[leaf]->HydroAcc[1] += TempHydroAccImport.HydroAcc[1];
            Phydro[leaf]->HydroAcc[2] += TempHydroAccImport.HydroAcc[2];
            Phydro[leaf]->Vsig = fmax(Phydro[leaf]->Vsig,TempHydroAccImport.Vsig);
            Phydro[leaf]->Nlist += TempHydroAccImport.Nlist;
#ifdef USE_SPSPH //{
            Phydro[leaf]->DZw += TempHydroAccImport.DZw;
#endif // USE_SPSPH //}

#ifdef USE_CELIB //{
#ifdef USE_METAL_DIFFUSION //{
            for(int l=0;l<CELibYield_Number-1;l++){
                Z[l] += TempHydroAccImport.dZ[l]*Phydro[leaf]->dt_hydro;
            }
#endif // USE_METAL_DIFFUSION //}
#endif // USE_CELIB //}

        }

        double InvMass = 1.0/Phydro[leaf]->Mass;
        Phydro[leaf]->Du *= InvMass;
        Phydro[leaf]->HydroAcc[0] *= InvMass;
        Phydro[leaf]->HydroAcc[1] *= InvMass;
        Phydro[leaf]->HydroAcc[2] *= InvMass;
#ifdef USE_CELIB //{
#ifdef USE_METAL_DIFFUSION //{
        double MassZ = 0.e0;
        for(int l=0;l<CELibYield_Number-1;l++){
            Phydro[leaf]->Elements[l+1] = Phydro[leaf]->Mass*Z[l];
            MassZ += Phydro[leaf]->Elements[l+1];
        }
        Phydro[leaf]->Elements[0] = Phydro[leaf]->Mass-MassZ;
        double LightLements = Phydro[leaf]->Elements[CELibYield_H]+Phydro[leaf]->Elements[CELibYield_He];
        Phydro[leaf]->Z = fmax((Phydro[leaf]->Mass-LightLements)/Phydro[leaf]->Mass,0.e0);
#if 0
        assert(Phydro[leaf]->Elements[0]>=0.0);
        int flag = 0;    
        for(int l=0;l<12;l++){
            if(Phydro[leaf]->Elements[l]<=0.0){
                fprintf(stderr,"Problem[%02d] %d %d %g\n",MPIGetMyID(),leaf,l,Phydro[leaf]->Elements[l]);
                flag ++;
            }
            //assert(Phydro[leaf]->Elements[l]>=0.0);
        }
        if(flag > 0){
            for(int l=0;l<12;l++){
                fprintf(stderr,"Trouble[%02d] %d %d %g\n",MPIGetMyID(),leaf,l,Phydro[leaf]->Elements[l]);
            }
        }
        assert(flag == 0);
#endif
#endif // USE_METAL_DIFFUSION //}
#endif // USE_CELIB //}
    }


    // Copy StructHydroAcc to temporary structures.
    if(NumberofHydroAccExportRecvLogAllocated < NImportAll){
        if(NumberofHydroAccExportRecvLogAllocated > 0){
            free(HydroAccExportRecvLog);
        }
        NumberofHydroAccExportRecvLogAllocated = NImportAll;
        HydroAccExportRecvLog = malloc(sizeof(struct StructHydroAccExport)*NImportAll);
    }
    NumberofHydroAccExportRecvLog = NImportAll;
    for(int i=0;i<NImportAll;i++){
        HydroAccExportRecvLog[i] = HydroAccExportRecv[i];
    }

    TimingResults.HydroAccThisStep += GetElapsedTime()-TimingResultThisRoutine;

    return;
}

static struct StructHydroAccImport ReturnStructureDuDtHydroAccCorrectionEngine(const int index, const int mode){

    double Pos[3] = {PhydroPos(index)[0],PhydroPos(index)[1],PhydroPos(index)[2]};
    double Vel[3] = {Phydro[index]->VelP[0],Phydro[index]->VelP[1],Phydro[index]->VelP[2]};
    double Kerneli = Phydro[index]->KernelPred;
    double Rhoi = Phydro[index]->RhoPred;
    double Fi = Phydro[index]->F;
    double Massi = Phydro[index]->Mass;
#ifdef USE_VARIABLE_ALPHA
    double Alphai = Phydro[index]->Alpha;
#else // USE_VARIABLE_ALPHA
    double Alphai = Pall.HydroAlpha;
#endif // USE_VARIABLE_ALPHA

#ifdef USE_DISPH //{
    double p_over_fund2_i = Pall.Gm1/Phydro[index]->EnergyDensityPred;
    double Weighti = Massi*Phydro[index]->UPred;
    double Volumei = Weighti/Phydro[index]->EnergyDensityPred;
#elif defined(USE_SPSPH) //}//{
    // P = (g-1)*\rho u = (g-1)*my/z *u
    // P/y^2 = (g-1)*m*u/(z*y)
    double p_over_fund2_i = (Pall.Gm1*Massi*Phydro[index]->UPred/(Phydro[index]->PseudoDensityPred*Phydro[index]->ZwPred));
    double Weighti = Phydro[index]->ZwPred;
    double PseudoDensityi = Phydro[index]->PseudoDensityPred;
    double Ddifi = Phydro[index]->Ddif;
    double Volumei = Weighti/PseudoDensityi;
#else // USE_DISPH //}//{
    double p_over_fund2_i = Pall.Gm1*(Phydro[index]->UPred/Phydro[index]->RhoPred);
    double Weighti = Massi;
    double Volumei = Weighti/Phydro[index]->RhoPred;
#endif // USE_DISPH //}
    double InvWeighti = 1.0/Weighti;

#ifdef USE_GRAD_H //{
    double Gradhi = Phydro[index]->Gradh;
#ifdef USE_GRAD_N //{
    double fi = Phydro[index]->fij;
#endif // USE_GRAD_N //}
#endif // USE_GRAD_H //}


#ifdef USE_THERMAL_CONDUCTIVITY //{
    double Temperaturei = Pall.ConvertUtoT*Phydro[index]->UPred;
    double Kappai = CalcThermalConductivity(index);
#endif // USE_THERMAL_CONDUCTIVITY //}

#if (defined(USE_SPSPH) || defined(USE_THERMAL_CONDUCTIVITY) || defined(USE_METAL_DIFFUSION)) //{
    double dti = Phydro[index]->dt_hydro;
#endif //(defined(USE_SPSPH) || defined(USE_THERMAL_CONDUCTIVITY) || defined(USE_METAL_DIFFUSION)) //}

    double InvKerneli = 1.e0/Kerneli;
    double csi = sqrt(Pall.GGm1*Phydro[index]->UPred);
    int Neighbors[MaxNeighborSize];
    struct StructHydroAccImport TempHydroAccImport = InitializeHydroAccImport;

    double TimeNeighborSearch = GetElapsedTime();
    int Nlist;
    if(mode == 0){
        Nlist = GetNeighborsPairsLimited(Pos,2.e0*Kerneli,Neighbors);
    } else {
        Nlist = GetNeighborsPairsLimitedImported(Pos,2.e0*Kerneli,Neighbors);
    }
    TimingResults.HydroAccNeighborSearchThisStep += GetElapsedTime()-TimeNeighborSearch;

    struct StructHydroAccPreFetch{
        double Mass;
        double Kernel;
        double Rho;
        //double U;
        double F;
#ifdef USE_VARIABLE_ALPHA
        double Alpha;
#endif // USE_VARIABLE_ALPHA
        double p_over_fund2_j;
        double Weightj;
        double Volumej;
        double csj;
#ifdef USE_SPSPH //{
        double PseudoDensityj;
        double Ddifj;
#endif // USE_SPSPH //}
#ifdef USE_GRAD_H //{
        double Gradh;
#ifdef USE_GRAD_N //{
        double fj;
#endif // USE_GRAD_N //}
#endif // USE_GRAD_H //}
#ifdef USE_THERMAL_CONDUCTIVITY //{
        double Kappa;
        double Temperature;
#endif // USE_THERMAL_CONDUCTIVITY //}
        double ActiveFactor; // Use for diffusion terms.
        double Pos[3];
        double Vel[3];
    } HydroAccPreFetch[Nlist];
    for(int k=0;k<Nlist;k++){
        int leaf = Neighbors[k];
        if(mode == 0){
            HydroAccPreFetch[k].Mass = PhydroMass(leaf);
            HydroAccPreFetch[k].Pos[0] = PhydroPosP(leaf)[0];
            HydroAccPreFetch[k].Pos[1] = PhydroPosP(leaf)[1];
            HydroAccPreFetch[k].Pos[2] = PhydroPosP(leaf)[2];
            HydroAccPreFetch[k].Vel[0] = Phydro[leaf]->VelP[0];
            HydroAccPreFetch[k].Vel[1] = Phydro[leaf]->VelP[1];
            HydroAccPreFetch[k].Vel[2] = Phydro[leaf]->VelP[2];
            HydroAccPreFetch[k].Rho = Phydro[leaf]->RhoPred;
            HydroAccPreFetch[k].Kernel = Phydro[leaf]->KernelPred;
#ifdef USE_DISPH //{
            HydroAccPreFetch[k].p_over_fund2_j = Pall.Gm1/Phydro[leaf]->EnergyDensityPred;
            HydroAccPreFetch[k].Weightj = Phydro[leaf]->Mass*Phydro[leaf]->UPred;
            HydroAccPreFetch[k].Volumej = HydroAccPreFetch[k].Weightj/Phydro[leaf]->EnergyDensityPred;
#elif defined(USE_SPSPH) //}//{
            HydroAccPreFetch[k].p_over_fund2_j = Pall.Gm1*Phydro[leaf]->Mass*Phydro[leaf]->UPred/(Phydro[leaf]->PseudoDensityPred*Phydro[leaf]->ZwPred);
            HydroAccPreFetch[k].Weightj = Phydro[leaf]->ZwPred;
            HydroAccPreFetch[k].PseudoDensityj = Phydro[leaf]->PseudoDensityPred;
            HydroAccPreFetch[k].Ddifj = Phydro[leaf]->Ddif;
            HydroAccPreFetch[k].Volumej = HydroAccPreFetch[k].Weightj/HydroAccPreFetch[k].PseudoDensityj;
#else // USE_DISPH //}//{
            HydroAccPreFetch[k].p_over_fund2_j = Pall.Gm1*Phydro[leaf]->UPred/Phydro[leaf]->RhoPred;
            HydroAccPreFetch[k].Weightj = Phydro[leaf]->Mass;
            HydroAccPreFetch[k].Volumej = HydroAccPreFetch[k].Weightj/Phydro[leaf]->RhoPred;
#endif // USE_DISPH //}
            HydroAccPreFetch[k].csj = sqrt(Pall.GGm1*Phydro[leaf]->UPred);

#ifdef USE_GRAD_H //{
            HydroAccPreFetch[k].Gradh = Phydro[leaf]->Gradh;
#ifdef USE_GRAD_N //{
            HydroAccPreFetch[k].fj = Phydro[leaf]->fij;
#endif // USE_GRAD_N //}
#endif // USE_GRAD_H //}
            HydroAccPreFetch[k].F = Phydro[leaf]->F;
#ifdef USE_VARIABLE_ALPHA
            HydroAccPreFetch[k].Alpha = Phydro[leaf]->Alpha;
#endif // USE_VARIABLE_ALPHA
#ifdef USE_THERMAL_CONDUCTIVITY //{
            HydroAccPreFetch[k].Temperature = Pall.ConvertUtoT*Phydro[index]->UPred;
            HydroAccPreFetch[k].Kappa = CalcThermalConductivity(index);
#endif // USE_THERMAL_CONDUCTIVITY //}
#if (defined(USE_SPSPH) || defined(USE_THERMAL_CONDUCTIVITY) || defined(USE_METAL_DIFFUSION)) //{
            if((Phydro[leaf]->Active)&&(Phydro[index]->dt_hydro>0.e0)){
                HydroAccPreFetch[k].ActiveFactor = fmax(Phydro[leaf]->dt_hydro,dti)/dti;
            } else {
                HydroAccPreFetch[k].ActiveFactor = 0.e0;
            }
#endif // (defined(USE_SPSPH) || defined(USE_THERMAL_CONDUCTIVITY) || defined(USE_METAL_DIFFUSION)) //}
        } else {
            HydroAccPreFetch[k].Mass = HydroAccExportRecvLog[leaf].Mass;
            HydroAccPreFetch[k].Pos[0] = HydroAccExportRecvLog[leaf].Pos[0];
            HydroAccPreFetch[k].Pos[1] = HydroAccExportRecvLog[leaf].Pos[1];
            HydroAccPreFetch[k].Pos[2] = HydroAccExportRecvLog[leaf].Pos[2];
            HydroAccPreFetch[k].Vel[0] = HydroAccExportRecvLog[leaf].Vel[0];
            HydroAccPreFetch[k].Vel[1] = HydroAccExportRecvLog[leaf].Vel[1];
            HydroAccPreFetch[k].Vel[2] = HydroAccExportRecvLog[leaf].Vel[2];
            HydroAccPreFetch[k].Rho = HydroAccExportRecvLog[leaf].Rho;
            HydroAccPreFetch[k].Kernel = HydroAccExportRecvLog[leaf].Kernel;
#ifdef USE_DISPH //{
            HydroAccPreFetch[k].p_over_fund2_j = Pall.Gm1/HydroAccExportRecvLog[leaf].EnergyDensity;
            HydroAccPreFetch[k].Weightj = HydroAccExportRecvLog[leaf].Mass*HydroAccExportRecvLog[leaf].U;
            HydroAccPreFetch[k].Volumej = HydroAccPreFetch[k].Weightj/HydroAccExportRecvLog[leaf].EnergyDensity;
#elif defined(USE_SPSPH) //}//{
            HydroAccPreFetch[k].p_over_fund2_j = Pall.Gm1*HydroAccExportRecvLog[leaf].Mass*HydroAccExportRecvLog[leaf].U
                                                /(HydroAccExportRecvLog[leaf].PseudoDensity*HydroAccExportRecvLog[leaf].Zw);
            HydroAccPreFetch[k].Weightj = HydroAccExportRecvLog[leaf].Zw;
            HydroAccPreFetch[k].PseudoDensityj = HydroAccExportRecvLog[leaf].PseudoDensity;
            HydroAccPreFetch[k].Ddifj = HydroAccExportRecvLog[leaf].Ddif;
            HydroAccPreFetch[k].Volumej = HydroAccPreFetch[k].Weightj/HydroAccPreFetch[k].PseudoDensityj;
#else // USE_DISPH //}//{
            HydroAccPreFetch[k].p_over_fund2_j = Pall.Gm1*HydroAccExportRecvLog[leaf].U/HydroAccExportRecvLog[leaf].Rho;
            HydroAccPreFetch[k].Weightj = HydroAccExportRecvLog[leaf].Mass;
            HydroAccPreFetch[k].Volumej = HydroAccPreFetch[k].Weightj/HydroAccExportRecvLog[leaf].Rho;
#endif // USE_DISPH //{
            HydroAccPreFetch[k].csj = sqrt(Pall.GGm1*HydroAccExportRecvLog[leaf].U);

#ifdef USE_GRAD_H //{
            HydroAccPreFetch[k].Gradh = HydroAccExportRecvLog[leaf].Gradh;
#ifdef USE_GRAD_N //{
            HydroAccPreFetch[k].fj = HydroAccExportRecvLog[leaf].fij;
#endif // USE_GRAD_N //}
#endif // USE_GRAD_H //}
            HydroAccPreFetch[k].F = HydroAccExportRecvLog[leaf].F;
#ifdef USE_VARIABLE_ALPHA
            HydroAccPreFetch[k].Alpha = HydroAccExportRecvLog[leaf].Alpha;
#endif // USE_VARIABLE_ALPHA
#ifdef USE_THERMAL_CONDUCTIVITY //{
            HydroAccPreFetch[k].Temperature = Pall.ConvertUtoT*HydroAccExportRecvLog[leaf].U;
            HydroAccPreFetch[k].Kappa = CalcThermalConductivityByTemperature(HydroAccPreFetch[k].Temperature);
#endif // USE_THERMAL_CONDUCTIVITY //}
#if (defined(USE_SPSPH) || defined(USE_THERMAL_CONDUCTIVITY) || defined(USE_METAL_DIFFUSION)) //{
            if((HydroAccExportRecvLog[leaf].Active)&&(Phydro[index]->dt_hydro>0.e0)){
                HydroAccPreFetch[k].ActiveFactor = fmax(HydroAccExportRecvLog[leaf].dt,dti)/dti;
            } else {
                HydroAccPreFetch[k].ActiveFactor = 0.e0;
            }
#endif // (defined(USE_SPSPH) || defined(USE_THERMAL_CONDUCTIVITY) || defined(USE_METAL_DIFFUSION)) //}
        }
    }

    int NlistActive = 0;
	for(int k=0;k<Nlist;k++){
        double xij[3];
#ifdef PERIODIC_RUN //{
        xij[0] = PeriodicDistance(Pos[0],HydroAccPreFetch[k].Pos[0],0);
        xij[1] = PeriodicDistance(Pos[1],HydroAccPreFetch[k].Pos[1],1);
        xij[2] = PeriodicDistance(Pos[2],HydroAccPreFetch[k].Pos[2],2);
#else // PERIODIC_RUN //{//}
        xij[0] = Pos[0]-HydroAccPreFetch[k].Pos[0];
        xij[1] = Pos[1]-HydroAccPreFetch[k].Pos[1];
        xij[2] = Pos[2]-HydroAccPreFetch[k].Pos[2];
#endif // PERIODIC_RUN //}
        double r2 = NORM2(xij);

        if(r2 > 0.e0){
            double r = sqrt(r2);

            double Kernelj = HydroAccPreFetch[k].Kernel;
            double InvKernelj = 1.e0/Kernelj;

            double dwi = dSPHKernel(r,InvKerneli);
            double dwj = dSPHKernel(r,InvKernelj);
            double dw = 0.5*(dwi+dwj);

            double vij[3];
            vij[0] = Vel[0]-HydroAccPreFetch[k].Vel[0];
            vij[1] = Vel[1]-HydroAccPreFetch[k].Vel[1];
            vij[2] = Vel[2]-HydroAccPreFetch[k].Vel[2];

            double vx = DOT_PRODUCT(xij,vij);

#if (HydroViscNormal) //{
            double visc = ReturnArtificialViscocityGM(vx,r2,csi,HydroAccPreFetch[k].csj,
                Rhoi,HydroAccPreFetch[k].Rho,Kerneli,HydroAccPreFetch[k].Kernel,
                Fi,HydroAccPreFetch[k].F
                SetParameterVariableAlpha(Alphai,HydroAccPreFetch[k].Alpha))*dw;
#elif (HydroViscVelSig) //}//{
            double visc = ReturnArtificialViscocityReimann(vx,r,csi,HydroAccPreFetch[k].csj,
                Rhoi,HydroAccPreFetch[k].Rho,Fi,HydroAccPreFetch[k].F,&(TempHydroAccImport.Vsig)
                SetParameterVariableAlpha(Alphai,HydroAccPreFetch[k].Alpha))*dw;
#endif // Visc //}

            double tmpj_visc = Massi*HydroAccPreFetch[k].Mass*visc;

#ifdef USE_GRAD_H //{
#ifdef USE_GRAD_N //{
            double fij = 1.0-fi/HydroAccPreFetch[k].Weightj;
            double fji = 1.0-HydroAccPreFetch[k].fj*InvWeighti;
#else //USE_GRAD_N //}//{
            double fij = Gradhi;
            double fji = HydroAccPreFetch[k].Gradh;
#endif //USE_GRAD_N //}

            double tmpj = Weighti*HydroAccPreFetch[k].Weightj
                          *(fij*p_over_fund2_i*dwi+fji*HydroAccPreFetch[k].p_over_fund2_j*dwj)
                         +tmpj_visc;
            double tmpj_energy = (Weighti*HydroAccPreFetch[k].Weightj*p_over_fund2_i)*fij*vx*dwi
                                    +0.5*tmpj_visc*vx;
#else // USE_GRAD_H //}//{
            double tmpj = Weighti*HydroAccPreFetch[k].Weightj
                          *(p_over_fund2_i+HydroAccPreFetch[k].p_over_fund2_j)*dw
                         +tmpj_visc;
            double tmpj_energy = (Weighti*HydroAccPreFetch[k].Weightj*p_over_fund2_i)*vx*dw
                                    +0.5*tmpj_visc*vx;
#endif //USE_GRAD_H //}

            TempHydroAccImport.Du += tmpj_energy;
            TempHydroAccImport.HydroAcc[0] -= tmpj*xij[0];
            TempHydroAccImport.HydroAcc[1] -= tmpj*xij[1];
            TempHydroAccImport.HydroAcc[2] -= tmpj*xij[2];

#ifdef USE_SPSPH //{
            TempHydroAccImport.DZw += HydroAccPreFetch[k].ActiveFactor*
                4.0*(Ddifi*HydroAccPreFetch[k].Ddifj)/(Ddifi+HydroAccPreFetch[k].Ddifj)
                *((Weighti*HydroAccPreFetch[k].Weightj)/(PseudoDensityi*HydroAccPreFetch[k].PseudoDensityj))
                *(PseudoDensityi-HydroAccPreFetch[k].PseudoDensityj)*dw;
#endif // USE_SPSPH //}
#ifdef USE_THERMAL_CONDUCTIVITY //{
            double tmpj_tc = (4*Kappai*HydroAccPreFetch[k].Kappa/(Kappai+HydroAccPreFetch[k].Kappa))
                        *(Volumei*HydroAccPreFetch[k].Volumej)*(HydroAccPreFetch[k].Temperature-Temperaturei)*dw;
#endif // USE_THERMAL_CONDUCTIVITY //}

            NlistActive ++; 
        }
    }
    TempHydroAccImport.Nlist = NlistActive;

    return (TempHydroAccImport);
}



struct StructHydroAccCorrectEnegyDensityGradh{ 
#ifdef USE_DISPH //{
    double EnergyDensity; 
#endif // USE_DISPH //}
#ifdef USE_SPSPH //{
    double PseudoDensity; 
#endif // USE_SPSPH //}
    double Vsig;
#ifdef USE_GRAD_H //{
    double Gradh; 
#ifdef USE_GRAD_N //{
    double fij; 
#endif // USE_GRAD_N //}
#endif // USE_GRAD_H //}
};

struct StructHydroAccCorrectEnegyDensityGradh ReturnStructureHydroAccCorrectionEnergyDensityGradhEngine(const int index, const int mode){

    struct StructHydroAccCorrectEnegyDensityGradh Temp = {
#ifdef USE_DISPH //{
        .EnergyDensity = 0.e0,
#endif // USE_DISPH //}
#ifdef USE_SPSPH //{
        .PseudoDensity = 0.e0,
#endif // USE_SPSPH //}
        .Vsig = 0.e0,
#ifdef USE_GRAD_H //{
        .Gradh = 0.e0,
#ifdef USE_GRAD_N //{
        .fij = 0.e0, 
#endif // USE_GRAD_N //}
#endif // USE_GRAD_H //}
    };

#if defined(USE_DISPH)||defined(USE_SPSPH) //{

    double Pos[3] = {PhydroPosP(index)[0],PhydroPosP(index)[1],PhydroPosP(index)[2]};
    double Vel[3] = {Phydro[index]->VelP[0],Phydro[index]->VelP[1],Phydro[index]->VelP[2]};
    double Kerneli = Phydro[index]->KernelPred;

    double InvKerneli = 1.e0/Kerneli;
    int Neighbors[MaxNeighborSize];

    double TimeComm = GetElapsedTime();
    int Nlist;
    if(mode == 0){
        Nlist = GetNeighborsLimited(Pos,2.e0*Kerneli,Neighbors);
    } else {
        Nlist = GetNeighborsLimitedImported(Pos,2.e0*Kerneli,Neighbors);
    }
    TimingResults.HydroDensityNeighborSearchThisStep += GetElapsedTime()-TimeComm;

    struct StructHydroDensityPreFetch{
        double Mass;
        double Pos[3];
        double Vel[3];
        double csj;
        double Weight;
    } HydroDensityPreFetch[Nlist];

	for(int k=0;k<Nlist;k++){
		int leaf = Neighbors[k];
        if(mode == 0){
            HydroDensityPreFetch[k].Mass = PhydroMass(leaf);
            // HydroDensityPreFetch[k].Kernel = Phydro[j]->KernelPred;
            HydroDensityPreFetch[k].Pos[0] = PhydroPosP(leaf)[0];
            HydroDensityPreFetch[k].Pos[1] = PhydroPosP(leaf)[1];
            HydroDensityPreFetch[k].Pos[2] = PhydroPosP(leaf)[2];
            HydroDensityPreFetch[k].Vel[0] = Phydro[leaf]->VelP[0];
            HydroDensityPreFetch[k].Vel[1] = Phydro[leaf]->VelP[1];
            HydroDensityPreFetch[k].Vel[2] = Phydro[leaf]->VelP[2];
            //HydroDensityPreFetch[k].cs = sqrt(Pall.GGm1*Phydro[j]->UPred);
#ifdef USE_DISPH //{
            HydroDensityPreFetch[k].Weight = Phydro[leaf]->Mass*Phydro[leaf]->UPred;
#endif // USE_DISPH //}
#ifdef USE_SPSPH //{
            HydroDensityPreFetch[k].Weight = Phydro[leaf]->ZwPred;
#endif // USE_SPSPH //}
        } else {
            HydroDensityPreFetch[k].Mass = HydroAccExportRecvLog[leaf].Mass;
            HydroDensityPreFetch[k].Pos[0] = HydroAccExportRecvLog[leaf].Pos[0];
            HydroDensityPreFetch[k].Pos[1] = HydroAccExportRecvLog[leaf].Pos[1];
            HydroDensityPreFetch[k].Pos[2] = HydroAccExportRecvLog[leaf].Pos[2];
            HydroDensityPreFetch[k].Vel[0] = HydroAccExportRecvLog[leaf].Vel[0];
            HydroDensityPreFetch[k].Vel[1] = HydroAccExportRecvLog[leaf].Vel[1];
            HydroDensityPreFetch[k].Vel[2] = HydroAccExportRecvLog[leaf].Vel[2];
#ifdef USE_DISPH //{
            HydroDensityPreFetch[k].Weight 
                = HydroAccExportRecvLog[leaf].Mass*HydroAccExportRecvLog[leaf].U;
#endif // USE_DISPH //}
#ifdef USE_SPSPH //{
            HydroDensityPreFetch[k].Weight = HydroAccExportRecvLog[leaf].Zw;
#endif // USE_SPSPH //}
        }
    }

	for(int k=0;k<Nlist;k++){
        double xij[3],vij[3];
#ifdef PERIODIC_RUN
        xij[0] = PeriodicDistance(Pos[0],HydroDensityPreFetch[k].Pos[0],0);
        xij[1] = PeriodicDistance(Pos[1],HydroDensityPreFetch[k].Pos[1],1);
        xij[2] = PeriodicDistance(Pos[2],HydroDensityPreFetch[k].Pos[2],2);
#else // PERIODIC_RUN
        xij[0] = Pos[0]-HydroDensityPreFetch[k].Pos[0];
        xij[1] = Pos[1]-HydroDensityPreFetch[k].Pos[1];
        xij[2] = Pos[2]-HydroDensityPreFetch[k].Pos[2];
#endif // PERIODIC_RUN

        double r2 = SQ(xij[0])+SQ(xij[1])+SQ(xij[2]);
        double r = sqrt(r2); 
		double w = SPHKernel(r,InvKerneli);
        double dw = dSPHKernel(r,InvKerneli);
#ifdef USE_DISPH //{
        Temp.EnergyDensity += HydroDensityPreFetch[k].Weight*w;
#endif // USE_DISPH //}
#ifdef USE_SPSPH //{
        Temp.PseudoDensity += HydroDensityPreFetch[k].Weight*w;
#endif // USE_SPSPH //}
#ifdef USE_GRAD_H //{
        Temp.Gradh += HydroDensityPreFetch[k].Weight*dw*r2;
#endif //USE_GRAD_H //}
    }

#endif // defined(USE_DISPH)||defined(USE_SPSPH) //}
    return Temp;
}

// Update U
static void CorrectU(const int NProcs, bool HydroAccExportFlags[][NProcs-1]){

    MPI_Status  mpi_status;

    int NExportThisTime[NProcs];
    int NImportThisTime[NProcs];

    float *UExport[NProcs];
    bool (*HydroInteractionFlags)[NProcs-1];
    HydroInteractionFlags = (void *)BufferHydroInteractionFlags;

    for(int i=0;i<NProcs-1;i++){
        NExportThisTime[i] = CheckHydroAccExportFlags(i,NProcs,HydroAccExportFlags,HydroInteractionFlags);
        UExport[i] = malloc(sizeof(float)*NExportThisTime[i]);

        if(NExportThisTime[i] > 0){
            int NExport = 0;
            for(int k=0;k<Pall.Nhydro;k++){
                if(HydroAccExportFlags[k][i]){
                    if(Phydro[k]->Active){
                        UExport[i][NExport] = Phydro[k]->U;
                    } else {
                        UExport[i][NExport] = Phydro[k]->UPred;
                    }
                    NExport ++;
                }
            }
            assert(NExport == NExportThisTime[i]);
        }
    }

    double TimeImbarance = GetElapsedTime();
    int NImportAll = 0;
    for(int i=0;i<NProcs-1;i++){
        MPI_Sendrecv(NExportThisTime+i,1,MPI_INT,CommunicationTable[i].SendRank,TAG_SPH_ACC_PRECOMM,
            NImportThisTime+i,1,MPI_INT,CommunicationTable[i].RecvRank,TAG_SPH_ACC_PRECOMM,
                MPI_COMM_WORLD,&mpi_status);
        NImportAll += NImportThisTime[i];
    }
    TimingResults.HydroImbaranceThisStep += GetElapsedTime() - TimeImbarance;

    float *UImport;
    UImport = malloc(sizeof(float)*NImportAll);

    MPI_Status mpi_status_Export_Send[NProcs-1];
    MPI_Request mpi_request_Export_Send[NProcs-1];
    MPI_Status mpi_status_Export_Recv[NProcs-1];
    MPI_Request mpi_request_Export_Recv[NProcs-1];

    int NImport = 0;
    for(int i=0;i<NProcs-1;i++){
        MPI_Isend(UExport[i],NExportThisTime[i],
            MPI_FLOAT,CommunicationTable[i].SendRank,TAG_SPH_ACC_EXPORT+i,
                MPI_COMM_WORLD,mpi_request_Export_Send+i);
        MPI_Irecv(UImport+NImport,NImportThisTime[i],
                MPI_FLOAT,CommunicationTable[i].RecvRank,TAG_SPH_ACC_EXPORT+i,
                    MPI_COMM_WORLD,mpi_request_Export_Recv+i);
        NImport += NImportThisTime[i];
    }
    assert(NImport == NImportAll);

    double TimeComm = GetElapsedTime();
    MPI_Waitall(NProcs-1,mpi_request_Export_Recv,mpi_status_Export_Recv);
    MPI_Waitall(NProcs-1,mpi_request_Export_Send,mpi_status_Export_Send);
    TimingResults.HydroAccCommThisStep += GetElapsedTime()-TimeComm;

    // Update
    for(int i=0;i<NImportAll;i++){
        HydroAccExportRecvLog[i].U = UImport[i];
    }
    free(UImport);

    return ;
}

// Update EnergyDensity and Gradh

void CalcDuDtAccEnergyDensityForCorrection(void){

#ifdef COOLING_RUN //{
    double TimingResultThisRoutine = GetElapsedTime();

    int NProcs = MPIGetNumProcs();
    MPI_Status  mpi_status;

    static int HydroAccExportFlagsMaxAllocated = 0;
    static int *ActiveIndexList;
    if(HydroAccExportFlagsMaxAllocated < MAX(Pall.Nhydro,NAdditionUnit)){
        if(HydroAccExportFlagsMaxAllocated > 0){
            free(ActiveIndexList);
        }
        HydroAccExportFlagsMaxAllocated = (int)(MAX(ForAngelsShare*Pall.Nhydro,NAdditionUnit));
        ActiveIndexList = malloc(sizeof(int)*HydroAccExportFlagsMaxAllocated);
    }

    int NActives = 0;
    int RootNodeID = 0;
    int NumberofLeaves = HydroNode[RootNodeID].NumberofLeaves;
    int header = HydroNode[RootNodeID].Leaves;
    for(int i=0;i<NumberofLeaves;i++){
        int leaf = header + i;
        if(NBCache[leaf].Active){
            ActiveIndexList[NActives] = NBCache[leaf].Leaf;
            NActives ++;
        }
    }

    struct StructHydroAccCorrectUk {
#ifndef USE_SPSPH //{
        double  U;
#endif // USE_DISPH //}
#ifdef USE_SPSPH //{
        double  Zw;
#endif // USE_SPSPH //}
    } *UkExport[NProcs];

    // double *UExport[NProcs];
    struct StructHydroAccCorrectEnegyDensityGradh *CorrectEnergyDensityGradh[NProcs];
    bool (*HydroInteractionFlags)[Pall.Nhydro];
    HydroInteractionFlags = (void *)BufferHydroInteractionFlags;


    // Update Particle predictors.
    for(int i=0;i<Pall.Nhydro;i++){
        if(Phydro[i]->Active){
            Phydro[i]->UPred = Phydro[i]->U;
        }
    }


    bool (*Log)[Pall.Nhydro];
    Log = (void *)ExportFlagsLog;

    for(int i=0;i<NProcs-1;i++){

        // UExport[i] = malloc(sizeof(double)*NExportLog[i]);
        UkExport[i] = malloc(sizeof(struct StructHydroAccCorrectUk)*NExportLog[i]);
        CorrectEnergyDensityGradh[i] = malloc(sizeof(struct StructHydroAccCorrectEnegyDensityGradh)*NExportLog[i]);

        int Offset = i*Pall.Nhydro;
        if(NExportLog[i] > 0){
            int NExport = 0;
            for(int k=0;k<Pall.Nhydro;k++){
                if(ExportFlagsLog[Offset+k]){
                    if(Phydro[k]->Active){
#ifndef USE_SPSPH //{
                        UkExport[i][NExport].U = Phydro[k]->UPred;
#endif //USE_DISPH //}
#ifdef USE_SPSPH //{
                        UkExport[i][NExport].Zw = Phydro[k]->ZwPred;
#endif //USE_SPSPH //}
                    } else {
#ifndef USE_SPSPH //{
                        UkExport[i][NExport].U = Phydro[k]->UPred;
#endif //USE_DISPH //}
#ifdef USE_SPSPH //{
                        UkExport[i][NExport].Zw = Phydro[k]->ZwPred;
#endif //USE_SPSPH //}
                    }
                    NExport ++;
                }
            }
            assert(NExport == NExportLog[i]);
        }
    }

    int NImportAll = 0;
    for(int i=0;i<NProcs-1;i++){
        NImportAll += NImportLog[i];
    }

    struct StructHydroAccCorrectUk *UkImport;
    UkImport = malloc(sizeof(struct StructHydroAccCorrectUk)*NImportAll);

    MPI_Status mpi_status_Export_Send[NProcs-1];
    MPI_Request mpi_request_Export_Send[NProcs-1];
    MPI_Status mpi_status_Export_Recv[NProcs-1];
    MPI_Request mpi_request_Export_Recv[NProcs-1];


    size_t StructureSizeUk = sizeof(struct StructHydroAccCorrectUk);
    int NImport = 0;
    int counter_send = 0;
    int counter_recv = 0;
    int SendFlag,RecvFlag;
    for(int i=0;i<NProcs-1;i++){
        if(NExportLog[i] > 0){
            MPI_Isend(UkExport[i],StructureSizeUk*NExportLog[i],
                MPI_BYTE,CommunicationTable[i].SendRank,TAG_SPH_ACC_EXPORT+i,
                    MPI_COMM_WORLD,mpi_request_Export_Send+counter_send);
            MPI_Test(mpi_request_Export_Send+counter_send,&SendFlag,MPI_STATUS_IGNORE);
            counter_send ++;
        } 
        if(NImportLog[i] > 0){
            MPI_Irecv(UkImport+NImport,StructureSizeUk*NImportLog[i],
                    MPI_BYTE,CommunicationTable[i].RecvRank,TAG_SPH_ACC_EXPORT+i,
                        MPI_COMM_WORLD,mpi_request_Export_Recv+counter_recv);
            MPI_Test(mpi_request_Export_Recv+counter_recv,&RecvFlag,MPI_STATUS_IGNORE);
            counter_recv ++;
        }
        NImport += NImportLog[i];
    }
    assert(NImport == NImportAll);

    double TimeComm = GetElapsedTime();
    MPI_Waitall(counter_send,mpi_request_Export_Send,mpi_status_Export_Send);
    MPI_Waitall(counter_recv,mpi_request_Export_Recv,mpi_status_Export_Recv);
    TimingResults.HydroAccCommThisStep += GetElapsedTime()-TimeComm;

    for(int i=0;i<NImportAll;i++){
#ifndef USE_SPSPH //{
        HydroAccExportRecvLog[i].U = UkImport[i].U;
#endif //USE_DISPH //}
#ifdef USE_SPSPH //{
        HydroAccExportRecvLog[i].Zw = UkImport[i].Zw;
#endif //USE_SPSPH //}
    }



#if defined(USE_DISPH)||defined(USE_SPSPH) //{
    // Update EnergyDensity and Gradh
    for(int k=0;k<NActives;k++){
        struct StructHydroAccCorrectEnegyDensityGradh TempHydroAccCorrectEnergyDensityGradh;
        int leaf = ActiveIndexList[k];
        TempHydroAccCorrectEnergyDensityGradh = ReturnStructureHydroAccCorrectionEnergyDensityGradhEngine(leaf,0);
#ifdef USE_DISPH //{
        Phydro[leaf]->EnergyDensity = TempHydroAccCorrectEnergyDensityGradh.EnergyDensity;
#endif //USE_DISPH //}
#ifdef USE_SPSPH //{
        Phydro[leaf]->PseudoDensity = TempHydroAccCorrectEnergyDensityGradh.PseudoDensity;
#endif //USE_SPSPH //}
#ifdef USE_GRAD_H //{
        Phydro[leaf]->Gradh = TempHydroAccCorrectEnergyDensityGradh.Gradh;
#endif // USE_GRAD_H //}

        if(NImportAll>0){
            TempHydroAccCorrectEnergyDensityGradh = ReturnStructureHydroAccCorrectionEnergyDensityGradhEngine(leaf,1);
#ifdef USE_DISPH //{
            Phydro[leaf]->EnergyDensity += TempHydroAccCorrectEnergyDensityGradh.EnergyDensity;
#endif //USE_DISPH //}
#ifdef USE_SPSPH //{
            Phydro[leaf]->PseudoDensity += TempHydroAccCorrectEnergyDensityGradh.PseudoDensity;
#endif //USE_SPSPH //}
#ifdef USE_GRAD_H //{
            Phydro[leaf]->Gradh += TempHydroAccCorrectEnergyDensityGradh.Gradh;
#endif // USE_GRAD_H //}
        }

#ifdef USE_DISPH //{
        Phydro[leaf]->EnergyDensityPred = Phydro[leaf]->EnergyDensity;
        double Fund = Phydro[leaf]->EnergyDensity;
#endif //USE_DISPH //}
#ifdef USE_SPSPH //{
        Phydro[leaf]->PseudoDensityPred = Phydro[leaf]->PseudoDensity;
        double Fund = Phydro[leaf]->PseudoDensity;
#endif //USE_SPSPH //}
        
#ifdef USE_GRAD_H //{
#ifdef USE_GRAD_N //{
        Phydro[leaf]->fij = -(Fund+Phydro[leaf]->Gradh/DIMENSION)*Phydro[leaf]->GradN/Phydro[leaf]->NumberDensity;
#endif // USE_GRAD_N //}
        Phydro[leaf]->Gradh = -DIMENSION*Fund/Phydro[leaf]->Gradh;
#endif // USE_GRAD_H //}
    }



    // Exhange EnergyDensity and Gradh
    for(int i=0;i<NProcs-1;i++){

        int Offset = i*Pall.Nhydro;
        //if(NExportThisTime[i] > 0){
        if(NExportLog[i] > 0){
            int NExport = 0;
            for(int k=0;k<Pall.Nhydro;k++){
                if(ExportFlagsLog[Offset+k]){ 
#ifdef USE_DISPH //{
                    CorrectEnergyDensityGradh[i][NExport].EnergyDensity = Phydro[k]->EnergyDensityPred;
#endif // USE_DISPH //}
#ifdef USE_SPSPH //{
                    CorrectEnergyDensityGradh[i][NExport].PseudoDensity = Phydro[k]->PseudoDensityPred;
#endif // USE_SPSPH //}
#ifdef USE_GRAD_H //{
                    CorrectEnergyDensityGradh[i][NExport].Gradh = Phydro[k]->Gradh;
#ifdef USE_GRAD_N //{
                    CorrectEnergyDensityGradh[i][NExport].fij = Phydro[k]->fij;
#endif // USE_GRAD_N //}
#endif // USE_GRAD_H //}
                    NExport ++;
                }
            }
            assert(NExport == NExportLog[i]);
        }
    }

    struct StructHydroAccCorrectEnegyDensityGradh *CorrectEnergyDensityGradhImport;
    CorrectEnergyDensityGradhImport = malloc(sizeof(struct StructHydroAccCorrectEnegyDensityGradh)*NImportAll);


    NImport = 0;
    size_t StructureSize = sizeof(struct StructHydroAccCorrectEnegyDensityGradh);
    counter_send = counter_recv = 0;
    for(int i=0;i<NProcs-1;i++){
        if(NExportLog[i]>0){
            MPI_Isend(CorrectEnergyDensityGradh[i],StructureSize*NExportLog[i],
                MPI_BYTE,CommunicationTable[i].SendRank,TAG_SPH_ACC_EXPORT+i,
                    MPI_COMM_WORLD,mpi_request_Export_Send+counter_send);
            MPI_Test(mpi_request_Export_Send+counter_send,&SendFlag,MPI_STATUS_IGNORE);
            counter_send ++;
        }
        if(NImportLog[i]>0){
            MPI_Irecv(CorrectEnergyDensityGradhImport+NImport,StructureSize*NImportLog[i],
                    MPI_BYTE,CommunicationTable[i].RecvRank,TAG_SPH_ACC_EXPORT+i,
                        MPI_COMM_WORLD,mpi_request_Export_Recv+counter_recv);
            MPI_Test(mpi_request_Export_Recv+counter_recv,&RecvFlag,MPI_STATUS_IGNORE);
            counter_recv ++;
        }
        NImport += NImportLog[i];
    }
    assert(NImport == NImportAll);

    TimeComm = GetElapsedTime();
    MPI_Waitall(counter_send,mpi_request_Export_Send,mpi_status_Export_Send);
    MPI_Waitall(counter_recv,mpi_request_Export_Recv,mpi_status_Export_Recv);
    TimingResults.HydroAccCommThisStep += GetElapsedTime()-TimeComm;

    // Update imported data 
    for(int i=0;i<NImportAll;i++){
#ifdef USE_DISPH //{
        HydroAccExportRecvLog[i].EnergyDensity = CorrectEnergyDensityGradhImport[i].EnergyDensity;
#endif // USE_DISPH //}
#ifdef USE_SPSPH //{
        HydroAccExportRecvLog[i].PseudoDensity = CorrectEnergyDensityGradhImport[i].PseudoDensity;
#endif // USE_SPSPH //}
#ifdef USE_GRAD_H //{
        HydroAccExportRecvLog[i].Gradh = CorrectEnergyDensityGradhImport[i].Gradh;
#ifdef USE_GRAD_N //{
        HydroAccExportRecvLog[i].fij = CorrectEnergyDensityGradhImport[i].fij;
#endif // USE_GRAD_N //}
#endif // USE_GRAD_H //}
    }
#endif //defined(USE_DISPH)||defined(USE_SPSPH) //}

    // calculation for local
    for(int k=0;k<NActives;k++){
        struct StructHydroAccImport TempHydroAccImport;
        int leaf = ActiveIndexList[k];
        TempHydroAccImport = ReturnStructureDuDtHydroAccCorrectionEngine(leaf,0);

        Phydro[leaf]->Du = TempHydroAccImport.Du;
        Phydro[leaf]->HydroAcc[0] = TempHydroAccImport.HydroAcc[0];
        Phydro[leaf]->HydroAcc[1] = TempHydroAccImport.HydroAcc[1];
        Phydro[leaf]->HydroAcc[2] = TempHydroAccImport.HydroAcc[2];
        Phydro[leaf]->Vsig = TempHydroAccImport.Vsig;
        Phydro[leaf]->Nlist = TempHydroAccImport.Nlist;
#ifdef USE_SPSPH //{
        Phydro[leaf]->DZw = TempHydroAccImport.DZw;
#endif // USE_SPSPH //}

        if(NImportAll>0){
            TempHydroAccImport 
                = ReturnStructureDuDtHydroAccCorrectionEngine(leaf,1);

            Phydro[leaf]->Du += TempHydroAccImport.Du;
            Phydro[leaf]->HydroAcc[0] += TempHydroAccImport.HydroAcc[0];
            Phydro[leaf]->HydroAcc[1] += TempHydroAccImport.HydroAcc[1];
            Phydro[leaf]->HydroAcc[2] += TempHydroAccImport.HydroAcc[2];
            Phydro[leaf]->Vsig = fmax(Phydro[leaf]->Vsig,TempHydroAccImport.Vsig);
            Phydro[leaf]->Nlist += TempHydroAccImport.Nlist;
#ifdef USE_SPSPH //{
            Phydro[leaf]->DZw = TempHydroAccImport.DZw;
#endif // USE_SPSPH //}

        }
        double InvMass = 1.0/Phydro[leaf]->Mass;
        Phydro[leaf]->Du *= InvMass;
        Phydro[leaf]->HydroAcc[0] *= InvMass;
        Phydro[leaf]->HydroAcc[1] *= InvMass;
        Phydro[leaf]->HydroAcc[2] *= InvMass;
    }

    free(UkImport);
#if defined(USE_DISPH)||defined(USE_SPSPH) //{
    free(CorrectEnergyDensityGradhImport);
#endif // USE_DISPH //}

    for(int i=0;i<NProcs-1;i++){
        free(UkExport[i]);
#if defined(USE_DISPH)||defined(USE_SPSPH) //{
        free(CorrectEnergyDensityGradh[i]);
#endif // USE_DISPH //}
    }

    TimingResults.HydroAccThisStep += GetElapsedTime()-TimingResultThisRoutine;
#endif //}

    return;
}


/* 
 * For simulations with sink particles.
 */
 
void CalcSinkPressureForcesToHydro(double Pos[], double Vel[], double Mass, double Kernel, double Rho){

    // double Rhoi = Phydro[index]->RhoPred;
    double Ui = 1.0;
    double Fi = 1.0;
    double Alpha = Pall.HydroAlpha;
    double Kerneli = Kernel;
    double Rhoi = Rho;

    // or Rho = (M/Kernel)^(1/3)

#if defined(ISOTHERMAL_EOS_RUN)
    double p_over_rho2_i = (SQ(Pall.CS)/Rho);
#elif defined(BAROTROPIC_EOS_RUN)
    double p_over_rho2_i = ReturnPoverRho2ForBarotropicRun(Rho);
#else
    double p_over_rho2_i = Pall.Gm1*(Ui/Rho);
#endif

    double InvKerneli = 1.e0/Kerneli;
#if defined(ISOTHERMAL_EOS_RUN)
    double csi = Pall.CS;
#elif defined(BAROTROPIC_EOS_RUN)
    double csi = ReturnSoundSpeedForBarotoropicRun(Rhoi);
#else
    double csi = sqrt(Pall.GGm1*Ui);
#endif
    int Neighbors[MaxNeighborSize];
    //struct StructHydroAccImport TempHydroAccImport = InitializeHydroAccImport;

    double TimeNeighborSearch = GetElapsedTime();
    //int Nlist = GetNeighborsPairsLimitedImported(Pos,2.e0*Kerneli,Neighbors);
    int Nlist = GetNeighborsLimited(Pos,Kerneli,Neighbors);
    TimingResults.HydroAccNeighborSearchThisStep += GetElapsedTime()-TimeNeighborSearch;

    if(Nlist > MaxNeighborSize){
        fprintf(stderr,"The neighbor list is over flow(Nlist = %d). function:%s,line:%d,file:%s\n",
                Nlist,__FUNCTION__,__LINE__,__FILE__);
        MPI_Finalize();
        exit(NeighborListOverFlow);
    }

    struct StructHydroAccPreFetch{
        double Mass;
        double Kernel;
        double Rho;
        double U;
        double F;
        double Pos[3];
        double Vel[3];
    } HydroAccPreFetch[Nlist];

    int counter = 0;
    for(int k=0;k<Nlist;k++){
        int j = Neighbors[k];
        if(Phydro[j]->Active){
            HydroAccPreFetch[counter].Mass = PhydroMass(j);
            HydroAccPreFetch[counter].Rho = Phydro[j]->RhoPred;
            HydroAccPreFetch[counter].Kernel = Phydro[j]->KernelPred;
            HydroAccPreFetch[counter].U = Phydro[j]->UPred;
            HydroAccPreFetch[counter].F = Phydro[j]->F;

            HydroAccPreFetch[counter].Pos[0] = PhydroPosP(j)[0];
            HydroAccPreFetch[counter].Pos[1] = PhydroPosP(j)[1];
            HydroAccPreFetch[counter].Pos[2] = PhydroPosP(j)[2];
            HydroAccPreFetch[counter].Vel[0] = Phydro[j]->VelP[0];
            HydroAccPreFetch[counter].Vel[1] = Phydro[j]->VelP[1];
            HydroAccPreFetch[counter].Vel[2] = Phydro[j]->VelP[2];
            counter ++;
        }
    }
    Nlist = counter;

    /// Calibrate.
    //double A = Rho/(dKernelHydroAccR(0.5*Kernel,2.0/Kernel)*0.5*Kernel);
    
    int NlistActive = 0;
	for(int k=0;k<Nlist;k++){
        double xij[3];
#ifdef PERIODIC_RUN
        xij[0] = PeriodicDistance(Pos[0],HydroAccPreFetch[k].Pos[0],0);
        xij[1] = PeriodicDistance(Pos[1],HydroAccPreFetch[k].Pos[1],1);
        xij[2] = PeriodicDistance(Pos[2],HydroAccPreFetch[k].Pos[2],2);
#else // PERIODIC_RUN
        xij[0] = Pos[0]-HydroAccPreFetch[k].Pos[0];
        xij[1] = Pos[1]-HydroAccPreFetch[k].Pos[1];
        xij[2] = Pos[2]-HydroAccPreFetch[k].Pos[2];
#endif // PERIODIC_RUN
        double r2 = NORM2(xij);

        //if(r2 > 0.e0){
        if(r2 > SQ(TINY*(Kerneli+ HydroAccPreFetch[k].Kernel))){
            double r = sqrt(r2);

            double Kernelj = HydroAccPreFetch[k].Kernel;
            double InvKernelj = 1.e0/Kernelj;
            double InvRhoj = 1.e0/HydroAccPreFetch[k].Rho;

            double dwi = dSPHKernel(r,InvKerneli);
            double dwj = dSPHKernel(r,InvKernelj);
            double dw = 0.5*(dwi+dwj);

            double vij[3];
            vij[0] = Vel[0]-HydroAccPreFetch[k].Vel[0];
            vij[1] = Vel[1]-HydroAccPreFetch[k].Vel[1];
            vij[2] = Vel[2]-HydroAccPreFetch[k].Vel[2];

            double vx = DOT_PRODUCT(xij,vij);

            double visc;
#if (HydroViscNormal)
            if(vx < 0.e0){ 
                double KernelMean = 0.5*(Kerneli+Kernelj);
                double mu = KernelMean*vx/(r2+Pall.HydroEta2*SQ(KernelMean));
#if defined(ISOTHERMAL_EOS_RUN)
                double cs = Pall.CS;
#elif defined(BAROTROPIC_EOS_RUN)
                double cs = 0.5*(csi+ReturnSoundSpeedForBarotoropicRun(HydroAccPreFetch[k].RhoPred));
#else
                double cs = 0.5*(csi+sqrt(Pall.GGm1*HydroAccPreFetch[k].U));
#endif
                visc = (-Pall.HydroAlpha*cs*mu + Pall.HydroBeta*SQ(mu))/(0.5*(Rhoi+HydroAccPreFetch[k].Rho));
                visc *= 0.5*(Fi+HydroAccPreFetch[k].F);
            } else {
                visc = 0.e0;
            }
#elif (HydroViscVelSig)
#if defined(ISOTHERMAL_EOS_RUN)
            double csj = Pall.CS;
#elif defined(BAROTROPIC_EOS_RUN)
            double csj = ReturnSoundSpeedForBarotoropicRun(HydroAccPreFetch[k].Rho);
#else
            double csj = sqrt(Pall.GGm1*HydroAccPreFetch[k].U);
#endif
            Phydro[Neighbors[k]]->Vsig = fmax(csi+csj,Phydro[Neighbors[k]]->Vsig);
            if(vx < 0.e0){ 
                double wij = vx/r;
                double vsig = csi+csj-SIGNAL_VELOCITY_BETA*wij; 
                visc = -0.5*Pall.HydroAlpha*vsig*wij/(0.5*(Rhoi+HydroAccPreFetch[k].Rho));
                visc *= 0.5*(Fi+HydroAccPreFetch[k].F);
                Phydro[Neighbors[k]]->Vsig = fmax(vsig,Phydro[Neighbors[k]]->Vsig);
            } else {
                visc = 0.e0;
            }
#endif
            visc *= dw;

#if defined(ISOTHERMAL_EOS_RUN)
            double p_over_rho2_j = SQ(Pall.CS)*InvRhoj;
#elif defined(BAROTROPIC_EOS_RUN)
            double p_over_rho2_j = ReturnPoverRho2ForBarotropicRun(HydroAccPreFetch[k].Rho);
#else
            double p_over_rho2_j = Pall.Gm1*HydroAccPreFetch[k].U*InvRhoj;
#endif
            //double pgrd = p_over_rho2_i*dwi+p_over_rho2_j*dwj;
            double pgrd = (p_over_rho2_i+p_over_rho2_j)*dw;
            //double pgrd = p_over_rho2_j*dw;

            double tmpj = Mass*(pgrd+visc);
            //double tmpj = A*dw+HydroAccPreFetch[k].Mass*(pgrd+visc);
            Phydro[Neighbors[k]]->HydroAcc[0] += tmpj*xij[0];
            Phydro[Neighbors[k]]->HydroAcc[1] += tmpj*xij[1];
            Phydro[Neighbors[k]]->HydroAcc[2] += tmpj*xij[2];
#if !(defined(ISOTHERMAL_EOS_RUN)||defined(BAROTROPIC_EOS_RUN))
            Phydro[Neighbors[k]]->Du += HydroAccPreFetch[k].Mass*(p_over_rho2_i*dw+0.5e0*visc)*vx;
#endif
            NlistActive ++;
        }
    }
    return ;
}

