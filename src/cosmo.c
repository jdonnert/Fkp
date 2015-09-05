#include "common.h"

struct cosmology Cosmo;

void set_cosmology(int cosmology)
{
    switch (cosmology){ // Parameters 
    case 0: 
        Cosmo = null;
        a2t = &a2t_null;
        t2a = &t2a_null;
        break;
    case 1:    
        Cosmo = WMAP3;
        a2t = &a2t_LCDM;
        t2a = &t2a_LCDM;
        break;
    case 2:    
        Cosmo = deSitter;
        a2t = &a2t_deSitter;
        t2a = &t2a_deSitter;
        break;
    default:
        Assert(0, "Cosmology not found");
        break;
    }

    if ( ! ThisTask.Rank ) 
		printf("\nSetting Cosmology : %s \n\n",Cosmo.name);

    return;
}

/* Einstein - deSitter Universe */
struct cosmology deSitter = {"deSitter", 0.7, 0.17, 1, 0,
    0, 0, 1, 0, 4.35e17};

double a2t_deSitter(double a)
{
    return pow(a ,3./2.) * deSitter.t0;
}

double t2a_deSitter(double t)
{
    return pow(t/deSitter.t0 , 2./3.);
}

/* NULL Cosmology = nocomoving */
struct cosmology null = {"NULL   ", 1, 0, 0, 0, 0, 0, 1, 0, 0};

double a2t_null(double a)
{   
    return a * Unit.Time;
}

double t2a_null(double t)
{
    return t / Unit.Time;
}

/* Lambda CDM Universe using WMAP3 cosmology*/
struct cosmology WMAP3 = { "WMAP3  ", 0.732, 0.0444, 0.266,
    0.732,-1.0, 0.079, 0.94E-26, 0.772, 13.73E9*yr2sec};

double a2t_LCDM(double a)
{
    return -1;
}

double t2a_LCDM(double t)
{
    return -1;
}
