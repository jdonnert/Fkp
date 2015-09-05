/*physical constants cgs*/
#define c			GSL_CONST_CGSM_SPEED_OF_LIGHT
#define e			GSL_CONST_CGSM_ELECTRON_CHARGE*GSL_CONST_CGSM_SPEED_OF_LIGHT
#define h_planck	GSL_CONST_CGSM_PLANCKS_CONSTANT_H
#define hbar 		GSL_CONST_CGSM_PLANCKS_CONSTANT_HBAR
#define k_B 		GSL_CONST_CGSM_BOLTZMANN
#define m_p 		GSL_CONST_CGSM_MASS_PROTON
#define m_pi0		(134.9766/c/c*1e6*eV2cgs)	/* Eidelman 04 [g]*/
#define m_pi		(139.57018/c/c*1e6*eV2cgs)
#define m_mu		(105.65837/c/c*1e6*eV2cgs)
#define m_e			GSL_CONST_CGSM_MASS_ELECTRON
#define sigma_T	    GSL_CONST_CGSM_THOMSON_CROSS_SECTION
#define sigma_pp	(32*1e-3*barn2cgs)	
#define finestruct	(1/137.036)

/* unit conversions */
#define barn2cgs	GSL_CONST_CGSM_BARN
#define eV2cgs		GSL_CONST_CGSM_ELECTRON_VOLT
#define GeV2cgs	    (eV2cgs*1e9)
#define Msol2cgs    ((double)(1.98892e33))
#define kpc2cgs 	((double)(3.08568025e21))
#define pc2cgs 	    ((double)(3.08568025e18))
#define yr2sec		31556926
#define Gyr2sec     ((double)(1e9*yr2sec))
#define rad2deg     57.29577791

/* other */
#define Tcmb		(2.728*(1+Snap.Redshift))	
#define BCMB0       3.24516e-6  /* @ z=0 */
#define Bcmb		(Bcmb0*p2(1+Snap.Redshift))	
#define H_frac		0.76	
#define He_frac	    (1.0-H_frac)	
#define u_mol		(4.0/(5.0*H_frac+3.0))	
#define n2ne 		(H_frac+0.5*He_frac)/(2.0*H_frac+0.75*He_frac)
#define yHelium	    (He_frac / (4.0 *H_frac))
#define mean_mol_weight (1.0+4.0*yHelium)/(1.0+3.0*yHelium+1.0)
#define adiabatic_index (5.0/3.0)
#define hubbletime     4.35e17     

