/* Simple Secondary Injection */
#include "../common.h"

#ifdef Q_SIMPLE_SECONDARIES

#define HIGHCUT 1e3     // Upper Energy Cut of PowerLaw 
#define LOWCUT 50       // Lower Energy Cut of PowerLaw 
#define INDEX -2.0		// Spectral index of power-law
#define X_CRE 1         // CRe fraction in thermal energy density per Gyr 

static double norm;
#pragma omp threadprivate(norm)

void Injection(const struct Gas_Data SphP_i, const double *p_cgs, 
		const int nGrid, double *Q)
{
	for (int i = 0; i < nGrid; i++) {
	
    	double gamma = p[i]/(m_e*c);
		
    	if (gamma < LOWCUT || gamma > HIGHCUT) {

			Q[i] = 0;
		
			continue;
		}

		double result = norm * pow(gamma , INDEX) * sqrt(1 - (gamma) /HIGHCUT) 
	   	 	* (1 - LOWCUT/(gamma));
		
   		double eps_therm = thermal_energy_density_cgs(SphP_i.U, SphP_i.Rho);

		Q[i] = X_CRE * pow(eps_therm,2) * result;
	}

    return ;
}

void Init_Simple_Secondaries()
{
    const double H = HIGHCUT, L = LOWCUT;

    double x = H;
    double integral = -1.0/(8*p2(x)*p2(H)) * (p2(x)*(L-4*H)*log(1-sqrt(1-x/H)) 
		+ p2(x) * (4*H-L)*log(sqrt(1-x/H)+1)
		- 2*H*sqrt(1-x/H) * (-2*H*L + 4*H*x + L*x));

    x = L;
    integral +=  -1.0/(8*p2(x)*p2(H)) * (p2(x)*(L-4*H)*log(1-sqrt(1-x/H)) 
		+ p2(x) * (4*H-L)*log(sqrt(1-x/H)+1)
		- 2*H*sqrt(1-x/H) * (-2*H*L + 4*H*x + L*x));
	
	#pragma omp parallel
    norm = 1.0/integral;

	rprintf("Q(p): Low=%g, High=%g, Norm=%g \n", L, H, norm);

    return;
}

#undef INDEX
#undef HIGHCUT
#undef LOWCUT
#undef X_CRE

#endif
