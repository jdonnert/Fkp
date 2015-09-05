void set_cosmology(int);

double NoComov(float);

struct cosmology {	
	char name[10];		    /* Identifier */
	float h;		        /* hubble parameter */
	float Omega_b;		    /* baryon overdensity */
	float Omega_m;		    /* matter overdensity */
	float Omega_l;		    /* dark energy overdensity */
	float w;		        /* equation of state parameter */
	float tau;		        /* optical depth of reionization */
	double rho0;		    /* critical density */
	float sigma8;		    /* galaxy fluctuation amplitude */
	float t0;		        /* age of the universe */
} Cosmo, WMAP3,null,deSitter;	

/* Time conversions */
double (*a2t)(double); 
double (*t2a)(double);

double a2t_deSitter(double);
double t2a_deSitter(double);

double a2t_null(double);
double t2a_null(double);

double a2t_LCDM(double);
double t2a_LCDM(double);

