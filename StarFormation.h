#ifndef	__STARFORMATION_H_INCLUDED__
#define	__STARFORMATION_H_INCLUDED__
void InitializeStarFormation(void);
void LogStarFormationRate(void);
void StarFormation(void);
double CheckTotalBaryonMass(void);
void ArtificialStarFormation(const double SFR, const double Duration);
#endif // __STARFORMATION_H_INCLUDED__
