#include "../common.h"
#include "modules.h"

#ifdef DPP_EXPANSION
inline double Reacceleration(size_t ipart,double p,double time)
{
    return DBL_MIN;
}
#endif

#ifdef HP_EXPANSION
static float a1, b1,  drhodt;
static size_t kpart = 1;

double Cooling(size_t ipart,double p,double time,double Dpp)
{
    float a_time = t2a(time); // time in cosmological expansion factor

    if (ipart != kpart || a_time == Snap.LastTime) {
        kpart = ipart;

        float dt = (Snap.LastTime - Snap.Time); // interpolation in a not t
        
        float last_rho = LastP[ipart].Rho;
        float rho = P[ipart].Rho;

        a1 = (last_rho - rho) /dt; 
        b1 = rho - a1*Snap.Time;

        drhodt = (last_rho - rho) / (a2t(Snap.Time) - a2t(Snap.LastTime));
    }

    float rho = (a1*a_time + b1);

    return - p/3.0*drhodt/rho; 
}

#endif

#ifdef N0_POWERLAW
inline double Initial_Spectrum(size_t ipart, double p,double time,double Dpp)
{
    double pp = p/(m_e*c);

    double result = pow(pp, -2)*(1 - pp / 1e3) *(1 - 1e2 /pp);

    if (result < 0 || !isfinite(result)) // Think of a better function 
        result = 0;

    return log10( result ) ;
}
#endif
