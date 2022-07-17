/* CosmologicalTime.c */
#ifndef __COSMOLOGICALTIME_H_INCLUDED__
#define __COSMOLOGICALTIME_H_INCLUDED__

double CalcCurrentTtoZ(void);
double CalcCurrentZtoT(void);
double CalcTtoZ(const double TCurrent);
double CalcZtoT(const double Redshift);
double CalcHubbleZ(void);
double CalcOmegaMZ(void);
double CalcOmegaLZ(void);

#endif //__COSMOLOGICALTIME_H_INCLUDED__
