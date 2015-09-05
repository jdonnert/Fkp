#include "common.h"

void set_units()
{
    Unit.Time = Unit.Length/Unit.Vel;

#if defined(I_GADGET2) || defined(I_NONE)
    rprintf("Setting System of Units <GADGET> \n"
            "   Unit Length = %g cm \n"
            "   Unit Time   = %g sec\n"
            "   Unit Mass   = %g g  \n"
            ,Unit.Length,Unit.Time, Unit.Mass);
#endif

    return;
}

/*Default Units:  cm/kpc, g/(1e10 M_sol), cm/sec/(km/sec), */
double density_cgs(float rho) 
{   
#if defined(I_GADGET2) || defined(I_NONE)
	return ((double)(rho) * Unit.Mass /
	    (Unit.Length * Unit.Length * Unit.Length)) * p2(Cosmo.h);
#endif
#if defined(I_TXT)
    return(rho);
#endif
}

double number_density_cgs(float rho)
{   
#if defined(I_GADGET2) || defined(I_NONE) 
    return(density_cgs(rho) * n2ne / (u_mol * m_p));
#endif
#if defined(I_TXT) 
    return(0);
#endif
}

double proton_number_density_cgs(float rho)
{
#if defined(I_GADGET2) || defined(I_NONE)
    return(density_cgs(rho) * H_frac / m_p);
#endif
#if defined(I_TXT) 
    return(0);
#endif
}

double temperature_cgs(float u)
{		
#if defined(I_GADGET2) || defined(I_NONE)
    return (2.0 / 3.0 * u * Unit.Vel * Unit.Vel * m_p *
		mean_mol_weight / k_B);
#endif
#if defined(I_TXT) 
    return(0);
#endif
}

double thermal_energy_density_cgs(float u, float rho)
{   
#if defined(I_GADGET2) || defined(I_NONE)
    return density_cgs(rho) / (m_p * u_mol) * k_B * temperature_cgs(u);
#endif
#if defined(I_TXT) 
    return 0 ;
#endif
}


