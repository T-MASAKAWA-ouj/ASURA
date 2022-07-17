/* ./StellarFeedback.c */ 
#ifndef __STELLARFEEDBACK_H_INCLUDED__ //{
#define __STELLARFEEDBACK_H_INCLUDED__

#define __CHECK_SUM__
#define __CHECK_WEIGHT__

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
    int Type; // 1=TypeII, 2=TypeIa, 3=AGB
    int Count;
    double InitialMass;
    double Metallicity;
    double Rvalue;
    double Lvalue;
    int IterationCount;
};

void InitializeStellarFeedback(void);
int StellarFeedbackGetIMFType(void);
int StellarFeedbackGetIaType(void);
int CountFeedbackNumber(void);
double StellarFeedbackRadiusInitialGuess(const int Index);
double StellarFeedbackGetNextExplosionTime(const int mode, const double Metallicity, const double InitialMass, const int Count);
struct StructStellarFeedbackLocalInfo RetrunStellarFeedbackLocalInfo(double Pos[restrict], const double Radius);
void StellarFeedback(void);

#ifdef TASK_TEST_STELLARFEEDBACK //{
void StellarFeedbackTestCalcFeedbackRadius(int NActives, const int IndexList[restrict], const int TypeList[restrict], double Radius[restrict], int Nlist[restrict], int CheckSum[restrict]);
void StellarFeedbackTestReleaseEnergyHeavyElements(const int NExplosion, const int IndexList[restrict], const int TypeList[restrict]);
#endif // TASK_TEST_STELLARFEEDBACK //}

#endif // __STELLARFEEDBACK_H_INCLUDED__ //}
