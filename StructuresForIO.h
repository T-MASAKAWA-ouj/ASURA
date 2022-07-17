
// This file is generated at 2022-03-15 16:13:14 +0900

struct StructPbodyIOCompact{

// omit StructPbodyptr Next; // <TMP> <LEAN>
// omit void    (*Baryon);   // <TMP> <LEAN>

    unsigned long int GlobalID;
// omit unsigned long long int OrderingKey; // <TMP> <LEAN> Key for ordering.
// omit unsigned long int    NextLeaf;   // <TMP> <LEAN> Next Leaf for Tree (Grav). This indicates NextLeaf in Pbody[NextLeaf].

// omit bool    Use;    // <TMP> <LEAN> If this structure is not use, the flag = false.
    bool Active;

// omit int InteractionList; // <TMP> <LEAN>

    float Pos[3];
// omit double    PosP[3];      // <TMP> <LEAN> The predictor of the particle position.
    float Vel[3];
    float Velh[3];
    float Acc[3];
    float AccOld[3];
    float Pot;
    float Mass;
    float Eps;

    short Type;
// omit short     k;            // <TMP> <LEAN>
    float dt;
// omit double    EraLocal;  // <TMP> <LEAN>

};

struct StructPbodyIOCompact CopyPbodyToTemporalStructureCompact(const int index);

StructPbody CopyTemporalStructureCompactToPbodyCompact(struct StructPbodyIOCompact PbodyIOCompact);

struct StructPbodyIOCompactDouble{

// omit StructPbodyptr Next; // <TMP> <LEAN>
// omit void    (*Baryon);   // <TMP> <LEAN>

    unsigned long int GlobalID;
// omit unsigned long long int OrderingKey; // <TMP> <LEAN> Key for ordering.
// omit unsigned long int    NextLeaf;   // <TMP> <LEAN> Next Leaf for Tree (Grav). This indicates NextLeaf in Pbody[NextLeaf].

// omit bool    Use;    // <TMP> <LEAN> If this structure is not use, the flag = false.
    bool Active;

// omit int InteractionList; // <TMP> <LEAN>

    double Pos[3];
// omit double    PosP[3];      // <TMP> <LEAN> The predictor of the particle position.
    double Vel[3];
    double Velh[3];
    double Acc[3];
    double AccOld[3];
    double Pot;
    double Mass;
    double Eps;

    short Type;
// omit short     k;            // <TMP> <LEAN>
    double dt;
// omit double    EraLocal;  // <TMP> <LEAN>

};

struct StructPbodyIOCompactDouble CopyPbodyToTemporalStructureCompactDouble(const int index);

StructPbody CopyTemporalStructureCompactDoubleToPbodyCompactDouble(struct StructPbodyIOCompactDouble PbodyIOCompactDouble);

struct StructPbodyIOLean{

// omit StructPbodyptr Next; // <TMP> <LEAN>
// omit void    (*Baryon);   // <TMP> <LEAN>

    unsigned long int GlobalID;
// omit unsigned long long int OrderingKey; // <TMP> <LEAN> Key for ordering.
// omit unsigned long int    NextLeaf;   // <TMP> <LEAN> Next Leaf for Tree (Grav). This indicates NextLeaf in Pbody[NextLeaf].

// omit bool    Use;    // <TMP> <LEAN> If this structure is not use, the flag = false.
// omit bool    Active; // <LEAN> This is used for the individual time step.

// omit int InteractionList; // <TMP> <LEAN>

    float Pos[3];
// omit double    PosP[3];      // <TMP> <LEAN> The predictor of the particle position.
    float Vel[3];
// omit double    Velh[3];      // <LEAN> The \"half-step slided\" velocity of the particle. For leapfrog integration.
    float Acc[3];
// omit double    AccOld[3];    // <LEAN> The acceleration of the particle.
    float Pot;
    float Mass;
    float Eps;

    short Type;
// omit short     k;            // <TMP> <LEAN>
// omit double    dt;           // <LEAN>
// omit double    EraLocal;  // <TMP> <LEAN>

};

struct StructPbodyIOLean CopyPbodyToTemporalStructureLean(const int index);

StructPbody CopyTemporalStructureLeanToPbodyLean(struct StructPbodyIOLean PbodyIOLean);

struct StructPhydroIOCompact{

// omit StructPhydroptr Next;   // <TMP> <LEAN>
// omit StructPbodyptr  Body;   // <TMP> <LEAN>

// omit bool Use;                // <TMP> <LEAN>
// omit unsigned int ExportFlag; // <TMP> <LEAN> if you use NProcs>32, increase this bit field.
    unsigned int Nlist;
// omit unsigned long int NextLeaf; // <TMP> <LEAN> Next Leaf for Tree (Neighbor Search).
// omit bool    Active;          // <TMP> <LEAN> if \"OFF\", the particle skips the cooling routine.
// omit bool    CoolingFlag;     // <TMP> <LEAN> if \"OFF\", the particle skips the cooling routine.

    int Leaf;

    bool GravityKickFlag;
    int k_hydro;
    float dt_hydro;
    float EraLocal_hydro;
// omit bool      HydroTimeStepLimiterFlag; // <TMP> <LEAN>
    int k_hydro_localmin;
    int k_hydro_localmin_old;
    float NextUpdateEra;
    float dt_hydro_localmin;

// omit double    Mass;         // <TMP> <LEAN>
// omit double    PosP[3];      // <TMP> <LEAN>
// omit double    VelP[3];      // <TMP> <LEAN>

    float Rho;
// omit double    RhoPred;    // <TMP> <LEAN> Predictor of density.
    float EnergyDensity;
// omit double    EnergyDensityPred; // <TMP> <LEAN> Predictor of the energy density used in DISPH.
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER
// omit double    SmoothedNumber;    // <TMP> <LEAN> Smoothed Number.
#endif
    float Kernel;
// omit double    KernelPred; // <TMP> <LEAN> Predictor of kernel. 0<- rij/Kernel ->2
    float Gradh;
    float NumberDensity;
    float GradN;
    float fij;
    float Div;
    float Rot[3];
    float F;
    float HydroAcc[3];
    float U;
// omit double    UPred;      // <TMP> <LEAN>
    float Du;
    float DuPrev;
    float DuCooling;

    float PseudoDensity;
// omit double    PseudoDensityPred;      // <TMP> The predictor of y
    float Zw;
// omit double    ZwPred;     // <TMP> The predictor of the weight used in SPSPH
    float DZw;
    float Ddif;


    float DQheat;
    float dMass;
    float Z;
    float ZII;
    float ZIa;
    float dZII;
    float dZIa;
    float Vsig;
    float Alpha;
    float DAlpha;

    float DZdiff;
    float Sxy;
    float Sxz;
    float Syx;
    float Syz;
    float Szx;
    float Szy;

#if VISCOSITY_TYPE == 1 //{
    float Bxx;
    float Bxy;
    float Bxz;
    float Byx;
    float Byy;
    float Byz;
    float Bzx;
    float Bzy;
    float Bzz;
#endif // VISCOSITY_TYPE == 1 //}


#ifdef USE_SPAANS2008_COOLING_FUNCTIONS
    float G0;
    float fH2;
#endif // USE_SPAANS2008_COOLING_FUNCTIONS

#if (UseSFModelSpawn)
    short SpawnTimes;
    float SpawnMass;
#endif
#ifdef USE_CELIB
    float Elements[16];
#endif // USE_CELIB
#ifdef USE_PARTICLE_TAG
    int Tag;
#endif

};

struct StructPhydroIOCompact CopyPhydroToTemporalStructureCompact(const int index);
struct StructPhydroIOCompact CopyPhydroToTemporalStructureCompactElement(StructPhydroptr const Ph);
StructPhydro CopyTemporalStructureCompactToPhydroCompact(struct StructPhydroIOCompact PhydroIOCompact);

struct StructPhydroIOCompactDouble{

// omit StructPhydroptr Next;   // <TMP> <LEAN>
// omit StructPbodyptr  Body;   // <TMP> <LEAN>

// omit bool Use;                // <TMP> <LEAN>
// omit unsigned int ExportFlag; // <TMP> <LEAN> if you use NProcs>32, increase this bit field.
    unsigned int Nlist;
// omit unsigned long int NextLeaf; // <TMP> <LEAN> Next Leaf for Tree (Neighbor Search).
// omit bool    Active;          // <TMP> <LEAN> if \"OFF\", the particle skips the cooling routine.
// omit bool    CoolingFlag;     // <TMP> <LEAN> if \"OFF\", the particle skips the cooling routine.

    int Leaf;

    bool GravityKickFlag;
    int k_hydro;
    double dt_hydro;
    double EraLocal_hydro;
// omit bool      HydroTimeStepLimiterFlag; // <TMP> <LEAN>
    int k_hydro_localmin;
    int k_hydro_localmin_old;
    double NextUpdateEra;
    double dt_hydro_localmin;

// omit double    Mass;         // <TMP> <LEAN>
// omit double    PosP[3];      // <TMP> <LEAN>
// omit double    VelP[3];      // <TMP> <LEAN>

    double Rho;
// omit double    RhoPred;    // <TMP> <LEAN> Predictor of density.
    double EnergyDensity;
// omit double    EnergyDensityPred; // <TMP> <LEAN> Predictor of the energy density used in DISPH.
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER
// omit double    SmoothedNumber;    // <TMP> <LEAN> Smoothed Number.
#endif
    double Kernel;
// omit double    KernelPred; // <TMP> <LEAN> Predictor of kernel. 0<- rij/Kernel ->2
    double Gradh;
    double NumberDensity;
    double GradN;
    double fij;
    double Div;
    double Rot[3];
    double F;
    double HydroAcc[3];
    double U;
// omit double    UPred;      // <TMP> <LEAN>
    double Du;
    double DuPrev;
    double DuCooling;

    double PseudoDensity;
// omit double    PseudoDensityPred;      // <TMP> The predictor of y
    double Zw;
// omit double    ZwPred;     // <TMP> The predictor of the weight used in SPSPH
    double DZw;
    double Ddif;


    double DQheat;
    double dMass;
    double Z;
    double ZII;
    double ZIa;
    double dZII;
    double dZIa;
    double Vsig;
    double Alpha;
    double DAlpha;

    double DZdiff;
    double Sxy;
    double Sxz;
    double Syx;
    double Syz;
    double Szx;
    double Szy;

#if VISCOSITY_TYPE == 1 //{
    double Bxx;
    double Bxy;
    double Bxz;
    double Byx;
    double Byy;
    double Byz;
    double Bzx;
    double Bzy;
    double Bzz;
#endif // VISCOSITY_TYPE == 1 //}


#ifdef USE_SPAANS2008_COOLING_FUNCTIONS
    double G0;
    double fH2;
#endif // USE_SPAANS2008_COOLING_FUNCTIONS

#if (UseSFModelSpawn)
    short SpawnTimes;
    double SpawnMass;
#endif
#ifdef USE_CELIB
    double Elements[16];
#endif // USE_CELIB
#ifdef USE_PARTICLE_TAG
    int Tag;
#endif

};

struct StructPhydroIOCompactDouble CopyPhydroToTemporalStructureCompactDouble(const int index);
struct StructPhydroIOCompactDouble CopyPhydroToTemporalStructureCompactDoubleElement(StructPhydroptr const Ph);
StructPhydro CopyTemporalStructureCompactDoubleToPhydroCompactDouble(struct StructPhydroIOCompactDouble PhydroIOCompactDouble);

struct StructPhydroIOLean{

// omit StructPhydroptr Next;   // <TMP> <LEAN>
// omit StructPbodyptr  Body;   // <TMP> <LEAN>

// omit bool Use;                // <TMP> <LEAN>
// omit unsigned int ExportFlag; // <TMP> <LEAN> if you use NProcs>32, increase this bit field.
// omit unsigned int Nlist;      // <LEAN>
// omit unsigned long int NextLeaf; // <TMP> <LEAN> Next Leaf for Tree (Neighbor Search).
// omit bool    Active;          // <TMP> <LEAN> if \"OFF\", the particle skips the cooling routine.
// omit bool    CoolingFlag;     // <TMP> <LEAN> if \"OFF\", the particle skips the cooling routine.

// omit int Leaf; // <LEAN> The ID to the HydroTree leaves.

// omit bool      GravityKickFlag; // <LEAN> The kick flag of gravitataional acceleration.
// omit int       k_hydro;         // <LEAN>
// omit double    dt_hydro;        // <LEAN>
// omit double    EraLocal_hydro;  // <LEAN> The local time in a current block.
// omit bool      HydroTimeStepLimiterFlag; // <TMP> <LEAN>
// omit int       k_hydro_localmin;  // <LEAN>
// omit int       k_hydro_localmin_old;  // <LEAN>
// omit double    NextUpdateEra;     // <LEAN> The next update timing in Era. This variable is used when the HydroTimeStepLimiterFlag is true.
// omit double    dt_hydro_localmin; // <LEAN>

// omit double    Mass;         // <TMP> <LEAN>
// omit double    PosP[3];      // <TMP> <LEAN>
// omit double    VelP[3];      // <TMP> <LEAN>

    float Rho;
// omit double    RhoPred;    // <TMP> <LEAN> Predictor of density.
    float EnergyDensity;
// omit double    EnergyDensityPred; // <TMP> <LEAN> Predictor of the energy density used in DISPH.
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER
// omit double    SmoothedNumber;    // <TMP> <LEAN> Smoothed Number.
#endif
    float Kernel;
// omit double    KernelPred; // <TMP> <LEAN> Predictor of kernel. 0<- rij/Kernel ->2
// omit double    Gradh; // <LEAN> Grad-h term
// omit double    NumberDensity; // <LEAN> Number density
// omit double    GradN;         // <LEAN> Grad-N term
// omit double    fij;           // <LEAN> Part of fij. fij = 1-fij/(Uj)
// omit double    Div;        // <LEAN>
// omit double    Rot[3];     // <LEAN>
// omit double    F;          // <LEAN>
    float HydroAcc[3];
    float U;
// omit double    UPred;      // <TMP> <LEAN>
    float Du;
// omit double    DuPrev;     // <LEAN> Du/Dt.
// omit double    DuCooling;  // <LEAN> The amount of the energy loss by the radiative cooling.

    float PseudoDensity;
    float PseudoDensityPred;
    float Zw;
    float ZwPred;
    float DZw;
    float Ddif;


// omit double    DQheat;     // <LEAN> The amount of the heating energy by SNe, which injects during dt.
// omit double    dMass;      // <LEAN> The additional mass by SNe.
    float Z;
    float ZII;
    float ZIa;
    float dZII;
    float dZIa;
// omit double    Vsig;       // <LEAN> Signal Velocity  or Maximum Velocity difference.
    float Alpha;
// omit double    DAlpha;     // <LEAN>

    float DZdiff;
    float Sxy;
    float Sxz;
    float Syx;
    float Syz;
    float Szx;
    float Szy;

#if VISCOSITY_TYPE == 1 //{
    float Bxx;
    float Bxy;
    float Bxz;
    float Byx;
    float Byy;
    float Byz;
    float Bzx;
    float Bzy;
    float Bzz;
#endif // VISCOSITY_TYPE == 1 //}


#ifdef USE_SPAANS2008_COOLING_FUNCTIONS
    float G0;
    float fH2;
#endif // USE_SPAANS2008_COOLING_FUNCTIONS

#if (UseSFModelSpawn)
// omit short   SpawnTimes; // <LEAN>
// omit double  SpawnMass;  // <LEAN>
#endif
#ifdef USE_CELIB
    float Elements[16];
#endif // USE_CELIB
#ifdef USE_PARTICLE_TAG
    int Tag;
#endif

};

struct StructPhydroIOLean CopyPhydroToTemporalStructureLean(const int index);
struct StructPhydroIOLean CopyPhydroToTemporalStructureLeanElement(StructPhydroptr const Ph);
StructPhydro CopyTemporalStructureLeanToPhydroLean(struct StructPhydroIOLean PhydroIOLean);

struct StructPstarIOCompact{

// omit StructPstarptr  Next;   // <TMP> <LEAN>
// omit StructPbodyptr  Body;   // <TMP> <LEAN>

// omit bool    Use;            // <TMP> <LEAN>
    bool TypeII;
    bool TypeIa;
// omit bool    HIIflag;        // <TMP> <LEAN>
#if defined(PRESERVE_SNII_EVENTRATE)
    bool TypeIIProb;
#endif

    short IMFTYPE;
    short NthChildren;
    unsigned long int ParentGlobalID;

    float InitialMass;
    float Mass;
    float FormationTime;
    float Z;
    float ZII;
    float ZIa;
    float TempForm;
    float RhoForm;
    float StromgrenRadius;
// omit double  Density;        // <TMP> <LEAN>
#ifdef USE_CELIB //{
    int SNIaCount;
    float EventTime;
    int AGBCount;
    float EventTimeAGB;
#ifdef USE_CELIB_NSM //{
    int NSMCount;
    float EventTimeNSM;
#endif // USE_CELIB_NSM //}
#ifdef USE_CELIB_ECSN //{
    int ECSNCount;
    float EventTimeECSN;
#endif // USE_CELIB_ECSN //}
#ifdef USE_CELIB_HN //{
    int HNCount;
    float EventTimeHN;
#endif // USE_CELIB_HN //}


    float Elements[16];
#endif // USE_CELIB //}

};

struct StructPstarIOCompact CopyPstarToTemporalStructureCompact(const int index);
struct StructPstarIOCompact CopyPstarToTemporalStructureCompactElement(StructPstarptr const Ps);

StructPstar CopyTemporalStructureCompactToPstarCompact(struct StructPstarIOCompact PstarIOCompact);

struct StructPstarIOCompactDouble{

// omit StructPstarptr  Next;   // <TMP> <LEAN>
// omit StructPbodyptr  Body;   // <TMP> <LEAN>

// omit bool    Use;            // <TMP> <LEAN>
    bool TypeII;
    bool TypeIa;
// omit bool    HIIflag;        // <TMP> <LEAN>
#if defined(PRESERVE_SNII_EVENTRATE)
    bool TypeIIProb;
#endif

    short IMFTYPE;
    short NthChildren;
    unsigned long int ParentGlobalID;

    double InitialMass;
    double Mass;
    double FormationTime;
    double Z;
    double ZII;
    double ZIa;
    double TempForm;
    double RhoForm;
    double StromgrenRadius;
// omit double  Density;        // <TMP> <LEAN>
#ifdef USE_CELIB //{
    int SNIaCount;
    double EventTime;
    int AGBCount;
    double EventTimeAGB;
#ifdef USE_CELIB_NSM //{
    int NSMCount;
    double EventTimeNSM;
#endif // USE_CELIB_NSM //}
#ifdef USE_CELIB_ECSN //{
    int ECSNCount;
    double EventTimeECSN;
#endif // USE_CELIB_ECSN //}
#ifdef USE_CELIB_HN //{
    int HNCount;
    double EventTimeHN;
#endif // USE_CELIB_HN //}


    double Elements[16];
#endif // USE_CELIB //}

};

struct StructPstarIOCompactDouble CopyPstarToTemporalStructureCompactDouble(const int index);
struct StructPstarIOCompactDouble CopyPstarToTemporalStructureCompactDoubleElement(StructPstarptr const Ps);

StructPstar CopyTemporalStructureCompactDoubleToPstarCompactDouble(struct StructPstarIOCompactDouble PstarIOCompactDouble);

struct StructPstarIOLean{

// omit StructPstarptr  Next;   // <TMP> <LEAN>
// omit StructPbodyptr  Body;   // <TMP> <LEAN>

// omit bool    Use;            // <TMP> <LEAN>
    bool TypeII;
    bool TypeIa;
// omit bool    HIIflag;        // <TMP> <LEAN>
#if defined(PRESERVE_SNII_EVENTRATE)
    bool TypeIIProb;
#endif

    short IMFTYPE;
    short NthChildren;
    unsigned long int ParentGlobalID;

// omit double  InitialMass;    // <LEAN> Initial Mass
    float Mass;
    float FormationTime;
    float Z;
    float ZII;
    float ZIa;
// omit double  TempForm;       // <LEAN> Temperature(t=FormationTime)
// omit double  RhoForm;        // <LEAN> Density(t=FormationTime) // UnitMass/UnitLength^3
// omit double  StromgrenRadius;// <LEAN>
// omit double  Density;        // <TMP> <LEAN>
#ifdef USE_CELIB //{
    int SNIaCount;
    float EventTime;
    int AGBCount;
    float EventTimeAGB;
#ifdef USE_CELIB_NSM //{
    int NSMCount;
    float EventTimeNSM;
#endif // USE_CELIB_NSM //}
#ifdef USE_CELIB_ECSN //{
    int ECSNCount;
    float EventTimeECSN;
#endif // USE_CELIB_ECSN //}
#ifdef USE_CELIB_HN //{
    int HNCount;
    float EventTimeHN;
#endif // USE_CELIB_HN //}


    float Elements[16];
#endif // USE_CELIB //}

};

struct StructPstarIOLean CopyPstarToTemporalStructureLean(const int index);
struct StructPstarIOLean CopyPstarToTemporalStructureLeanElement(StructPstarptr const Ps);

StructPstar CopyTemporalStructureLeanToPstarLean(struct StructPstarIOLean PstarIOLean);

struct StructPsinkIOCompact{

// omit StructPsinkptr  Next;   // <TMP> <LEAN>
// omit StructPbodyptr  Body;   // <TMP> <LEAN>

// omit bool    Use;            // <TMP> <LEAN>
    unsigned long int ParentGlobalID;

// omit double  PosP[3];        // <TMP> <LEAN>
// omit double  VelP[3];        // <TMP> <LEAN>
    float dt_localmin;

    int NumberofAbsorbedParticles;
    float AccretionRadius;
    float MergingDistance;
    float AM[3];
    float FormationTime;
    float Z;
    float ZII;
    float ZIa;
    float AccretionMass;
    float AccretionMassGas;
    float AccretionMassStar;
    float AccretionMassToBH;

#ifdef USE_CELIB
    float Elements[16];
#endif // USE_CELIB

};

struct StructPsinkIOCompact CopyPsinkToTemporalStructureCompact(const int index);
struct StructPsinkIOCompact CopyPsinkToTemporalStructureCompactElement(StructPsinkptr const Psk);

StructPsink CopyTemporalStructureCompactToPsinkCompact(struct StructPsinkIOCompact PsinkIOCompact);

struct StructPsinkIOCompactDouble{

// omit StructPsinkptr  Next;   // <TMP> <LEAN>
// omit StructPbodyptr  Body;   // <TMP> <LEAN>

// omit bool    Use;            // <TMP> <LEAN>
    unsigned long int ParentGlobalID;

// omit double  PosP[3];        // <TMP> <LEAN>
// omit double  VelP[3];        // <TMP> <LEAN>
    double dt_localmin;

    int NumberofAbsorbedParticles;
    double AccretionRadius;
    double MergingDistance;
    double AM[3];
    double FormationTime;
    double Z;
    double ZII;
    double ZIa;
    double AccretionMass;
    double AccretionMassGas;
    double AccretionMassStar;
    double AccretionMassToBH;

#ifdef USE_CELIB
    double Elements[16];
#endif // USE_CELIB

};

struct StructPsinkIOCompactDouble CopyPsinkToTemporalStructureCompactDouble(const int index);
struct StructPsinkIOCompactDouble CopyPsinkToTemporalStructureCompactDoubleElement(StructPsinkptr const Psk);

StructPsink CopyTemporalStructureCompactDoubleToPsinkCompactDouble(struct StructPsinkIOCompactDouble PsinkIOCompactDouble);

struct StructPsinkIOLean{

// omit StructPsinkptr  Next;   // <TMP> <LEAN>
// omit StructPbodyptr  Body;   // <TMP> <LEAN>

// omit bool    Use;            // <TMP> <LEAN>
    unsigned long int ParentGlobalID;

// omit double  PosP[3];        // <TMP> <LEAN>
// omit double  VelP[3];        // <TMP> <LEAN>
// omit double  dt_localmin;    // <LEAN> The local minimum time-steps.

    int NumberofAbsorbedParticles;
    float AccretionRadius;
    float MergingDistance;
    float AM[3];
    float FormationTime;
    float Z;
    float ZII;
    float ZIa;
    float AccretionMass;
    float AccretionMassGas;
    float AccretionMassStar;
    float AccretionMassToBH;

#ifdef USE_CELIB
    float Elements[16];
#endif // USE_CELIB

};

struct StructPsinkIOLean CopyPsinkToTemporalStructureLean(const int index);
struct StructPsinkIOLean CopyPsinkToTemporalStructureLeanElement(StructPsinkptr const Psk);

StructPsink CopyTemporalStructureLeanToPsinkLean(struct StructPsinkIOLean PsinkIOLean);


