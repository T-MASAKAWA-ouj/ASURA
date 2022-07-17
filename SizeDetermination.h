/* ./SizeDetermination.c */
#ifndef __SIZEDETERMINATION_H_INCLUDED__
#define __SIZEDETERMINATION_H_INCLUDED__

enum {
    CS_TypeHydro,
    CS_TypeHII,
    CS_TypeSN,
};


#include "StellarFeedback.h"

void CalcSize(void);
int ReturnCalcSizeElementNumber(const int Type, const bool Global);
void CalcSizeGetHydroInfo_i(const int Index, double *PotentialMin, double VCOM[]);
void CalcSizeSetSNInfo(struct StructActiveSNParticle ActiveSNParticle[]);

#endif // __SIZEDETERMINATION_H_INCLUDED__
