/* Cosmological.c */
#ifndef __COSMOLOGICAL_H_INCLUDED__
#define __COSMOLOGICAL_H_INCLUDED__
void UpdateCosmologicalParameters(void);
void UpdateAdaptiveSofteningFactor(void);
double CalcCriticalDensityZ0(void);
double CalcHubbleParameterZ0(void);
double CalcCriticalDensityZ(void);
double CalcHubbleParameterZ(void);
double CalcOverDensityZ(void);
double CalcVirializedDensityZ(void);
void CalcUniformSphereForce(void);

#endif //__COSMOLOGICAL_H_INCLUDED__
