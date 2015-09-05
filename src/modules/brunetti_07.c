/* Brunetti & Lazarian 2007, Cassano & Brunetti 2006 */

#include "../common.h"
#include "modules.h"

#if defined(DPP_BRUNETTI_07) 

static inline float calc_Dpp_Brunetti_07(size_t, float, float); 

static double a, b;
static size_t jpart = 1;

double Reacceleration(size_t ipart, double p, double time)
{   
    float a_time = t2a(time); // cosmological expansion factor

    if ( a_time == Snap.LastTime || ipart != jpart ) {  

        a = (SphP[ipart].Dpp - LastSphP[ipart].Dpp) 
            / (Snap.Time - Snap.LastTime);

        b = LastSphP[ipart].Dpp - a*Snap.LastTime;
 
        jpart = ipart;
    }
    
   return (a*a_time + b) *p*p;
}


void Calculate_Dpp()
{
    const double eta_t = 0.20; // Fraction of turb energy in magnetos waves
    const double turb_limiter = 0.8; // limit vturb to this fraction of csound
    
    rprintf("Calculating Dpp/p^2, (Brunetti & Lazarian 2007);\n"
        "   eta_t = <%g>  \n   turb limiter = <%g>\n\n"
        , eta_t, turb_limiter);

	#pragma omp parallel
    for (size_t ipart=0; ipart<ThisTask.Npart[0]; ipart++) {

        SphP[ipart].Dpp = calc_Dpp_Brunetti_07(ipart, eta_t, turb_limiter);
        
        Assert(isfinite(SphP[ipart].Dpp),"Found non-finite Dpp");
    }

    return ;
}

static inline float calc_Dpp_Brunetti_07(size_t ipart, float eta_t, 
        float turb_limiter) 
{
    double scale = 2*P[ipart].Hsml*Unit.Length;
    
    double vturb2 = SphP[ipart].VRms*Unit.Vel * SphP[ipart].VRms*Unit.Vel;

    double csound2 = adiabatic_index * (adiabatic_index-1) 
            * SphP[ipart].U * Unit.Vel*Unit.Vel;
    
    double vlimit2 = turb_limiter * turb_limiter * csound2;
    if (vturb2 > vlimit2)   //  Limit vturb to fraction of sound speed
        vturb2 = vlimit2;

    double result = 9e-8*eta_t*vturb2*vturb2/(scale*csound2);
        
    return result;
}
#endif // DPP_BRUNETTI_07

#ifdef HP_BRUNETTI_07
static double a1, b1, a2, b2;
static size_t kpart = 1;

inline double Cooling(size_t ipart,double p,double time,double Dpp)
{
    float a_time = t2a(time);       // cosmological expansion factor

    if ( ipart != kpart || a_time == Snap.LastTime ){
        kpart = ipart;

        double dt = (Snap.LastTime - Snap.Time); // interp in a not t
                
        a1 = (number_density_cgs(LastP[ipart].Rho) 
            -number_density_cgs(P[ipart].Rho)) /dt; 
        b1 = number_density_cgs(P[ipart].Rho) - a1*Snap.Time;

        a2 = (length3(LastSphP[ipart].Bfld)
            -length3(SphP[ipart].Bfld))/dt;
        b2 = length3(SphP[ipart].Bfld) - a2*Snap.Time;
    }
    
    double nth = (a1*a_time + b1);
    double bfld = (a2*a_time + b2);

    double loss_ic = 3.3e-29 * nth * ( 1 + log( p/(m_e*c) /nth ) / 75.0 );
    
    double loss_rad = 4.8e-4 * p*p * ( pow(bfld/BCMB0,2) + 1);

#ifdef TIMEVARIABLECMB 
	loss_rad +=  4.8e-4 * p*p * ( pow(a_time,-4) - 1 );
#endif
    
	return -2.0/p * Dpp + loss_rad + loss_ic;
}
#endif // HP_BRUNETTI_07

#ifdef Q_CASSANO_05
void Injection(const size_t ipart, const double *p, const int nGrid, 
		const double *Dpp, const double time, double *Q)
{
	for (int i = 0; i < nGrid; i++) {
	    
		double pp = p[i]/(m_e*c);

      	Q[i] = pow(pp, -2) * sqrt(1 - pp / 1e-4/ 1e8)
        	* (1 - 498820 * 1e-4/pp);

	    if (Q[i] < 0 || !isfinite(Q[i]))  
    	    Q[i] = 0; // Think of a better function
   	}

    return ;
}
#endif // Q_CASSANO_05

