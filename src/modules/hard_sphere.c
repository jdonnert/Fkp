/* Test Cases from Park & Petrosian 1995/6 */
#include "../common.h"
#include "modules.h"

#ifdef DPP_HARD_SPHERE
double Reacceleration(size_t ipart, double p, double time)
{
#if HP_HARD_SPHERE == 5

    return p*p*p;

#else

    return p*p;

#endif
}
#endif

#ifdef HP_HARD_SPHERE
double Cooling(size_t ipart, double p,double time,double Dpp)
{
#if HP_HARD_SPHERE == 0         // B=-1, r=2 

    return - p - p*p;    

#elif HP_HARD_SPHERE == 1       // B=1, r=0 

    return - p + 1;               

#elif HP_HARD_SPHERE == 2       // B=1, r=2 

    return - p + p*p;    

#elif HP_HARD_SPHERE == 3       // B=-1, r=0 

    return - p - 1;              

#elif HP_HARD_SPHERE == 4       // B=0, r=0
    
    return -p;

#elif HP_HARD_SPHERE == 5       // time dependent solution
    
    return -p*p;

#endif


}
#endif // HP_HARD_SPHERE

#ifdef Q_HARD_SPHERE
double Injection(size_t ipart, double p,double time,double Dpp)
{
    const double sigma = 0.01, p0 = 0.1; // Gaussian 

#if HP_HARD_SPHERE == 5

    return 0;

#else 

    return 1/sigma/sqrt(2*pi) * exp(-0.5 * pow((p-p0)/sigma,2));

#endif
}
#endif

#ifdef T_HARD_SPHERE
double Escape(size_t ipart,double p,double time,double Dpp)
{

#if HP_HARD_SPHERE == 4

    return p;

#else

    return 1;

#endif
}
#endif

#ifdef N0_HARD_SPHERE
double Initial_Spectrum(size_t ipart,double p,double time,double Dpp)
{
#if HP_HARD_SPHERE == 5
    const double sigma = 0.01, p0 = 0.1; // Gaussian 

    return log10(1/sigma/sqrt(2*pi) * exp(-0.5 * pow((p-p0)/sigma,2)));

#else 

    return -DBL_MAX;
    
#endif

}
#endif

