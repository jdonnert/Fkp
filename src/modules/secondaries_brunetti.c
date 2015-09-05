/* Brunetti & 2005 Secondaries */
#include "../common.h"

#ifdef Q_BRUNETTI_05

#define MIN 0.0			/*Integrals parameter for P123 */
#define MAX 1.0
#define NSTEP 500000.0

#define mMu (m_mu*c*c*1e-6/eV2cgs)	/* [MeV] */
#define mPi (m_pi*c*c*1e-6/eV2cgs)	/* [MeV] */

/* Parameters of injection function */
#define X_CRP Param.X_crp /* 1% per unit time */

static double P1, P2, P3;

static double a, A_s, Kp_prefac; 

/* for interpolation */
static double a1, a2, b1, b2;
static size_t jpart;

/* Q(p) from Brunetti & Blasi 2005 */
double Injection_Brunetti05(size_t ipart, double p,double time,double Dpp)
{
    double Ee, Ee_GeV, Kp, result;
    double dt,rho, u;

    if (p >= p[Param.NGrid-3]) /* minimize boundary error */
        return(0);

    if ( ipart != jpart || time == Snap.LastTime ){ /* interpolation */
        jpart = ipart;

        dt = (Snap.LastTime - Snap.Time);
                
        a1 = (LastP[ipart].Rho-P[ipart].Rho) /dt; 
        b1 = P[ipart].Rho - a1*Snap.Time;

        a2 = (LastSphP[ipart].U -SphP[ipart].U) /dt; 
        b2 = SphP[ipart].U - a2*Snap.Time;
        
    }

    rho = (a1*time + b1);
    u = (a2*time + b2);

    Ee = p*c;
    Ee_GeV = Ee/GeV2cgs;
    Kp = Kp_prefac  * thermal_energy_density_cgs(u, rho);    

    result = A_s * proton_number_density_cgs(rho)*Kp*pow(Ee,-Param.a_crp)
        * (P1 + P2*log(a * Ee_GeV/6.4) + P3*sqrt(a*Ee_GeV) );
    

    if (result < 0 || !isfinite(result))
        result = 0;

     /* e+ & e- => factor 2 */
     return(2*result);
}

void Init_Brunetti05()
{
    const double E0_proton = Param.p0_crp*m_p*c*c; // minimum of CR protons

	compute_P1P2P3();

	a = 2 * mPi * mPi / (mPi * mPi + mMu * mMu);
	A_s = 2 * c * sigma_pp * pow(a, 1.0 - s);
    
    /* norm rel to thermal energy density  */
    Kp_prefac  = X_CRP * (s-2)/(pow(E0_proton,2-s)); 
    
	return;
}

/* The spectrum is relatively sensitive to these 3 
 * parameters. Brunetti 05  does not define them 
 * accurate enough.
 */
void compute_P1P2P3()
{
    const double s = Param.a_crp;
    const double c1 = 1.22;
    const double c2 = 0.92;

	double m_fac1 = 2 * mPi * mPi / (mMu * mMu + mPi * mPi);
	double m_fac2 = (m_fac1 * m_fac1 - 1) * (1 + (mMu / mPi) * (mMu / mPi))
	    / (1 - (mMu / mPi) * (mMu / mPi));

	double aa = 5.0 / 12 * (1 + 0.2 * m_fac2);
	double bb = -3.0 / 4 * (1 + m_fac2);
	double cc = 1.0 / 3 * (1 + 2.0 * m_fac2);

	double P_fac1 = aa /s/s + bb / (s + 2) / (s + 2) + cc / (s + 3) / (s + 3);
	double P_fac2 = aa / s + bb / (s + 2) + cc / (s + 3);
	double P_fac3 = 1.5 * (aa / (s + 0.5) + bb / (s + 2.5) + cc / (s + 3.5));

	double h = fabs(MIN - MAX / NSTEP);

	double I0 = 0, I1 = 0, I2 = 0;

	for (int i = 0; i < NSTEP; i++) {
		double x = MIN + i * h + 0.5 * h;
		double f1 = c1 * pow(1.0 - x, 3.5) + c2 * exp(-18.0 * x);
		double f2 = pow(x, s - 2.0);

		I0 += h * f1 * f2;
		I1 += h * f1 * f2 * log(x);
		I2 += h * f1 * pow(x, s - 1.5);
	}

	P1 = P_fac1 * I0 - P_fac2 * I1;
	P2 = P_fac2 * I0;
	P3 = P_fac3 * I2;

	rprintf("Using Brunetti's Secondaries:\n"
		 "s    = %1.4f\n"
		 "P1   = %1.6f\n"
		 "P2   = %1.6f\n" 
		 "P3   = %1.6f\n\n"
		 , s, P1, P2, P3);

	return;
}

#undef MIN
#undef MAX
#undef NSTEP
#undef mMu
#undef mPi

#endif // Q_BRUNETTI_05
