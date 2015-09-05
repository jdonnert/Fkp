/* Petrosian 2001, Sarazin 1988 p 79, Schlickeiser 2002 */

#include "../common.h"
#include "modules.h"

#ifdef DPP_DONNERT_13

#define MACHCUT 0.8

static inline float calc_Dpp_Donnert_13(size_t, float, float); 

void Reacceleration(const struct Gas_Data SphP_i, const double *p, 
		const int nGrid, double *Dpp)
{   
    for (int i = 0; i < nGrid; i++) 
		Dpp[i] = SphP_i.Dpp * p[i] * p[i];
	
   return ;
}

void Calculate_Dpp()
{
    const double eta_t = 0.3; // Fraction of turb energy in magnetos waves
    const double turb_limiter = 0.8; // limit vturb to this fraction of csound
    
    rprintf("Calculating Dpp/p^2, (Brunetti & Lazarian 2007);\n"
        "   eta_t = <%g>  \n   turb limiter = <%g>\n\n"
        , eta_t, turb_limiter);

    for (size_t ipart=0; ipart<ThisTask.Npart[0]; ipart++) {

        SphP[ipart].Dpp = calc_Dpp_Donnert_13(ipart, eta_t, turb_limiter);
        
        Assert(isfinite(SphP[ipart].Dpp),"Found non-finite Dpp");
    }

	rprintf("\n Assuming IC cooling at z=%g  \n\n", REDSHIFT);

    return ;
}

static inline float calc_Dpp_Donnert_13(size_t ipart, float eta_t, 
        float turb_limiter) 
{
	if (SphP[ipart].Mach > MACHCUT)
		return 0;

    double scale = 2*SphP[ipart].Hsml*Unit.Length;
    
    double vturb2 = SphP[ipart].VRms*Unit.Vel * SphP[ipart].VRms*Unit.Vel;

    double csound2 = adiabatic_index * (adiabatic_index-1) 
            * SphP[ipart].U * Unit.Vel*Unit.Vel;
    
    double vlimit2 = turb_limiter * turb_limiter * csound2;
    if (vturb2 > vlimit2)   //  Limit vturb to fraction of sound speed
        vturb2 = vlimit2;

    double result = 1.64e-8*eta_t*vturb2*vturb2/(scale*csound2) * 
		(-0.25 - log(sqrt(csound2)/c));
        
    return result;
}

#endif // DPP_DONNERT_13

#ifdef HP_DONNERT_13

#define r0 (e*e/m_e/c/c)

const double fac_IC = 4.0/9.0 * r0*r0 *p2(BCMB0) * p2(p2(1+REDSHIFT));
const double fac_Syn = 4.0/9.0 * r0*r0;
const double fac_Coul = 4*pi*r0*r0 * m_e*c*c;
const double fac_Brems = 8 * finestruct *r0*r0*m_e*c*c;

void Cooling(const struct Gas_Data SphP_i, const double *q, const int nGrid, 
		const double *Dpp, const int ipart, double *Hp_out)
{
	double nth = number_density_cgs(SphP_i.Rho);
	double sqrt_nth = sqrt(nth);

	double dt = a2t(Snap.Time) - a2t(Snap.LastTime);

    double dNthDt = number_density_cgs(SphP_i.DRho) / dt;

	double bfld = length3(SphP_i.Bfld);
    double temp = temperature_cgs(SphP_i.U);

    for (int i = 0; i < nGrid; i++) {
		
	    double pp = q[i]/(m_e*c);
    
    	double gamma = pp - 1;
	    double gamma2 = gamma * gamma;
    	double beta2 = 1 - 1/gamma2;
	    double beta = sqrt(beta2);
    	
	    double coulomb_log = 37.8 + log(temp*1e-8 / (sqrt_nth * 1e-3 ) );
			// T>4e5

    	double loss_IC = fac_IC * beta2*gamma2;
    	double loss_Synchro = fac_Syn * beta2*gamma2 * bfld*bfld;

    	double loss_Coul = fac_Coul / beta * nth * coulomb_log ;

    	double loss_Brems = fac_Brems * nth *gamma* (H_frac + 3*He_frac) 
        	* ( log(2*gamma) - 1.0/3.0 ); // valid only for complete ionisation

		double losses = loss_IC + loss_Synchro + loss_Coul + loss_Brems;

		Hp_out[i] = -2.0*p[i]*Dpp[i]  + losses  - p[i]/3.0*dNthDt/nth;
	}
    
	return ; 
}

#endif // HP_DONNERT_13

