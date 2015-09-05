/* This is REALLY boring */
#include "../common.h"
#include "modules.h"

#ifdef DPP_ZERO
inline void Reacceleration(const struct Gas_Data SphP_i, const double *p, .
		const int nGrid, double *Dpp)
{
	for (int i = 0; i < nGrid; i++)
		Dpp[i] = 0;

    return ;
}
#endif


#ifdef Q_ZERO
inline void Injection(const struct Gas_Data SphP_i, const double *p, 
		const int nGrid, double *Q)
{
		for (int i = 0; i < nGrid; i++)
		Q[i] = 0;

    return 0;
}
#endif

#ifdef T_ZERO
void Escape(const struct Gas_Data SphP_i, const double *p, const int nGrid, 
		double *T)
{
	for (int i = 0; i < nGrid; i++)
		T[i] = DBL_MAX-1;

    return ;
}
#endif

#ifdef N0_ZERO
inline double Initial_Spectrum(size_t ipart, double p, double time,double Dpp)
{
    return -DBL_MAX;
}
#endif


