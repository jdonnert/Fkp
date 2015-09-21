/* Injection from primary shock acceleration, Kang+ 2007, Donnert+ 2015 
 * We join the spectrum with a Maxwellian at downstream density and Temp. */

#include "../common.h"

#ifdef Q_SHOCK_PRIMARIES

//#define DEBUG

#include "modules.h"

#define BISECTION_ACCURACY 1e-4
#define NSTEP 2048 // err < 5e-3 in norm integral

#define HIGHCUT 1e3     // Upper Energy Cut of PowerLaw 
#define LOWCUT 50       // Lower Energy Cut of PowerLaw 
#define INDEX -2.0		// Spectral index of power-law
#define X_CRE 0.01

static double cre_spectrum_hadronic_background(const double, const double);
static double maxwellian(const double T, const double p, const double nth);
static double injection_efficiency_Kang07(const double mach);
static double cre_energy_integral(const double s, const double p_inj_cgs);
static double cre_spectrum_normalisation(const struct Gas_Data, const double );
static double cre_spectral_index(double mach);
static double cre_spectrum(const double norm, const double p, const double s);
static void downstream_quantities(const struct Gas_Data, double *Temp_dw, 
									 double *Nth_dw, double *epsCRe);

static double Integral = 0;

void Injection(const struct Gas_Data SphP_i, const double *p_cgs, 
		const int nGrid, const double dt, double * Q)
{
#ifdef Q_SHOCK_WITH_HADRONIC_BACKGROUND
    double eps_therm = thermal_energy_density_cgs(SphP_i.U, SphP_i.Rho);
	
	double norm = p2(eps_therm) * X_CRE / Integral; 
		
	for (int i = 0; i < nGrid; i++)
		Q[i] += cre_spectrum_hadronic_background(norm, p_cgs[i]);
#endif

	if (SphP_i.Mach < 1.7)
		return ; 

	double s = cre_spectral_index(SphP_i.Mach);

	double CRe_norm = cre_spectrum_normalisation(SphP_i, s);

	for (int i = 0; i < nGrid; i++) 
		Q[i] = cre_spectrum(CRe_norm, p_cgs[i], s) / dt; 

	return ;
}

static double cre_spectrum(const double norm, const double p, const double s)
{
	return norm * pow(p, -s);
}

static double cre_spectrum_normalisation(const struct Gas_Data SphP_i, 
		const double s)
{
	double temp_dw = 0, nth_dw = 0, epsCRe = 0;

	downstream_quantities(SphP_i, &temp_dw, &nth_dw, &epsCRe);
	
	if (epsCRe == 0)
		return 0;

	const double cre_norm_fac = epsCRe /c /pow(m_e*c, 2-s);

	double cre_norm = 0;

	double left = sqrt(2*m_e*k_B*temp_dw); // peak of maxwellian
	double right = left * 30;	// maxwellian is very steep

	double delta = 1, 
		   p_inj = 0;

	int it = 0;

	while (fabs(delta) > BISECTION_ACCURACY) { // root finding

		p_inj = left + 0.5 * (right-left);

		double n_maxw = maxwellian(temp_dw, p_inj, nth_dw);

		cre_norm = cre_norm_fac / cre_energy_integral(s, p_inj/(m_e*c));

		double n_cre = cre_spectrum(cre_norm, p_inj, s);

		delta = (n_maxw - n_cre) / (n_maxw + n_cre); // weighted difference

#ifdef DEBUG
	if (SphP_i.ID == 1607043) 
		printf("%g %g %g %g %g %g %g %g %d\n", left, p_inj, right, n_maxw, 
				n_cre, delta, cre_norm, cre_energy_integral(s, p_inj), it++);
#endif

		if (delta > 0)
			left = p_inj;
		else
			right = p_inj;
	}

#ifdef DEBUG
	if (SphP_i.ID == 1607043) 
		printf("Found M=%g p_inj=%g CRe_norm=%g \n ",SphP_i.Mach, p_inj/m_e/c, cre_norm);
#endif

	return cre_norm;
}

/* 
 * injected CRe energy density from shock conditions
 */

static void downstream_quantities(const struct Gas_Data SphP_i, double *Temp_dw,
		double *Nth_dw, double *EpsCRe)
{
	const double gamma = adiabatic_index;

	const double mach = SphP_i.Mach;

	const double rho_up = SphP_i.Shock_Density;
	const double pres_up = SphP_i.Shock_Pressure;
	
	double mach2 = mach*mach; // Mo, v.d.Bosch & White, eq. 8.49
											
	double sigma = 1/(1/mach2 + (gamma-1)/(gamma+1) * (1-1/mach2));
			
	double u_up = pres_up /rho_up /(gamma - 1);

	double epsTherm_up = thermal_energy_density_cgs(u_up, rho_up); // cgs !
	
	double rho_dw  = rho_up * sigma;
	
	*Nth_dw = number_density_cgs(rho_dw);

	double pres_dw = pres_up * ( 2*gamma/(gamma + 1) * p2(mach) 
		  - (gamma - 1)/(gamma + 1) ); 	// Mo, v.d.Bosch & White, eq. 8.50

	double u_dw = pres_dw /rho_dw /(gamma - 1);

	*Temp_dw = temperature_cgs(u_dw);

	double epsTherm_dw = thermal_energy_density_cgs(u_dw, rho_dw);  // cgs

	double eta_cre = Param.Xi_ep * injection_efficiency_Kang07(mach); 

	*EpsCRe = eta_cre * (epsTherm_dw - epsTherm_up * pow(sigma, gamma)); 
		// cgs, Ensslin+ 2007

if (SphP_i.ID == 1607043) 
	printf("EpsCRe = %g, epsdw=%g epsup=%g sigma=%g mach=%g "
			"rho_up=%g pres_up=%g ID=%d \n", 
			*EpsCRe, epsTherm_dw,epsTherm_up, sigma, mach, rho_up,
			pres_up, SphP_i.ID);

	if (*EpsCRe < 0 || !isfinite(*EpsCRe)) {

		double vshock = SphP_i.Shock_Velocity;

		printf("EpsCRe = %g, epsdw=%g epsup=%g sigma=%g mach=%g vshock=%g "
			"rho_up=%g pres_up=%g ID=%d \n", 
			*EpsCRe, epsTherm_dw,epsTherm_up, sigma, mach, vshock, rho_up,
			pres_up, SphP_i.ID);

		*EpsCRe = 0;
	}

	return ;
}

/* Maxwell Boltzmann distribution in CGS & dp */
static double maxwellian(const double temp, const double p, const double nth)
{
	const double arg = 2*m_e*k_B*temp;
	
	return nth * 4*pi*p*p * pow(pi*arg, -1.5) * exp(-(p*p)/arg);
}

/* Kang, Ryu, Cen, Ostriker 2007 */
static double injection_efficiency_Kang07(const double mach)
{
	const double coeff[5] = { 5.46, -9.78, 4.17, -0.334, 0.570 };
	
	if (mach <= 2)
		return 1.96e-3 * (mach*mach - 1) ; // eq A2
	
	const double mach4 = p2(mach)*p2(mach);

	double eta = 0;
	
	for (int i = 0; i < 5; i++ )
		eta += coeff[i] * pow(mach-1, i) / mach4; // eq A5
	
	return eta;
}

/* 
 * integration & pmin in units of m_e*c 
 */

static double cre_energy_integral(const double s, const double pmin)
{
	const double pmax = 1e50; // ~= inf

	const double step = log10(pmax/pmin) / (NSTEP-1);

	double p_last = pmin;
	double f_last = pow(pmin, -s) * (sqrt(pmin*pmin + 1) - 1);

	double integral = 0;

	for (int i = 1; i < NSTEP; i++) {
	
		double p = pmin * pow(10, i*step);

		double dp = p - p_last;
		
		double f = pow(p, -s) * (sqrt(p*p + 1) - 1);

		integral += 0.5 * dp * (f + f_last); // trapezoidal rule

		p_last = p;
		f_last = f;
	}

	return integral;
}

/* Standard DSA, e.g. Kang+ 2012 */
static double cre_spectral_index(double mach)
{
	if (mach <= 1.01) // formula diverges, return maximum 
		return DBL_MAX;
	
	const double m2 = p2(mach);

	double s = 2 * (m2 + 1) / (m2 - 1);

	return s;
}

#ifdef Q_SHOCK_WITH_HADRONIC_BACKGROUND

static double cre_spectrum_hadronic_background(const double norm, 
		const double p)
{
	const double gamma = p/(m_e*c);

    if (gamma < LOWCUT || gamma > HIGHCUT)
	    return 0;

	return norm *  pow(gamma , INDEX) * sqrt(1 - (gamma) /HIGHCUT) 
	    * (1 - LOWCUT/(gamma));
}

void Init_Hadronic_Background()
{
	const double H = HIGHCUT, L = LOWCUT;

    double x = H;

    Integral = -1.0/(8*p2(x)*p2(H)) * (p2(x)*(L-4*H)*log(1-sqrt(1-x/H)) 
		+ p2(x) * (4*H-L)*log(sqrt(1-x/H)+1)
		- 2*H*sqrt(1-x/H) * (-2*H*L + 4*H*x + L*x));

    x = L;

    Integral +=  -1.0/(8*p2(x)*p2(H)) * (p2(x)*(L-4*H)*log(1-sqrt(1-x/H)) 
		+ p2(x) * (4*H-L)*log(sqrt(1-x/H)+1)
		- 2*H*sqrt(1-x/H) * (-2*H*L + 4*H*x + L*x));
	
	printf("Shock Primaries & Hadronic Background Q(p):"
			" Low=%g, High=%g, Norm=%g \n", L, H, 1/Integral);

    return;
	
}
#endif


#undef NSTEP
#undef BISECTION_ACCURACY 

#endif // Q_SHOCK_PRIMARIES
