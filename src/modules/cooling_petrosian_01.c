/* Petrosian 2001, Sarazin 1988 p 79, Schlickeiser 2002 */

#include "../common.h"
#include "modules.h"

#ifdef HP_PETROSIAN_01

static float a1, b1, a2, b2, a3, b3, dNthDt;
static size_t kpart = 1;

static inline float cross_section_bremsstrahlung(float);
#define r0 (e*e/m_e/c/c)
#define Sigma 5.670373e-8 // Stefan-Boltzmann-constant

static const double fac_IC = 4.0/9.0 * r0*r0 *p2(Bcmb0);
static const double fac_Syn = 4.0/9.0 * r0*r0;
static const double fac_Coul = 4*pi*r0*r0 * m_e*c*c;
static const double fac_Brems = 8 * finestruct *r0*r0*m_e*c*c;

inline double Cooling(size_t ipart, double p,double time,double Dpp)
{
    float a_time = t2a(time); // time in cosmological expansion factor

    if (ipart != kpart || a_time == Snap.LastTime) {
        kpart = ipart;

        float dt = (Snap.LastTime - Snap.Time); // interpolation in a not t
        
        float sqrt_last_nth = sqrt(number_density_cgs(LastP[ipart].Rho));
        float sqrt_nth = sqrt(number_density_cgs(P[ipart].Rho));

        dNthDt = (number_density_cgs(P[ipart].Rho) 
                - number_density_cgs(LastP[ipart].Rho)) / dt;

        a1 = (sqrt_last_nth - sqrt_nth) /dt; 
        b1 = sqrt_nth - a1*Snap.Time;

        a2 = (length3(LastSphP[ipart].Bfld)
            -length3(SphP[ipart].Bfld))/dt;
        b2 = length3(SphP[ipart].Bfld) - a2*Snap.Time;

        a3 = (temperature_cgs(LastSphP[ipart].U) 
                -  temperature_cgs(SphP[ipart].U) ) /dt;
        b3 = temperature_cgs(SphP[ipart].U) - a3*Snap.Time;
    }

    float sqrt_nth = (a1*a_time + b1);  
    float nth = sqrt_nth * sqrt_nth; 
    float bfld = (a2*a_time + b2);
    float temp = (a3*a_time + b3);
  
    float pp = p/(m_e*c);

    float gamma = pp - 1;
    float gamma2 = gamma * gamma;
    float beta2 = 1 - 1/gamma;
    float beta = sqrt(beta2);
    
    float coulomb_log = 37.8 + log(temp*1e-8 / (sqrt_nth * 1e-3 ) ); // T>4e5

    double loss_IC = fac_IC * beta2*gamma2;
    double loss_Synchro = fac_Syn * beta2*gamma2 * bfld*bfld;

    double loss_Coul = fac_Coul / beta * nth * coulomb_log ;
    /* valid only for complete ionisation */
    double loss_Brems = fac_Brems * nth *gamma* (H_frac + 3*He_frac) 
        * ( log(2*gamma) - 1.0/3.0 );

    double losses = loss_IC + loss_Synchro + loss_Coul + loss_Brems; 

    return -2.0*p*Dpp - p/3.0*dNthDt/nth + losses; 
}

#endif // HP_PETROSIAN_01

