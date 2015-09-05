void Reacceleration(const struct Gas_Data, const double *, const int, double*);
void Cooling(const struct Gas_Data, const double *,const int, const double *, 
		const int ,double *);
void Injection(const struct Gas_Data, const double *, const int, const double,
		double *);
void Escape(const struct Gas_Data,const double *,const int, double *);
double Initial_Spectrum(size_t, double, double, double);

/* Init functions for modules */
#ifdef COMPUTE_DPP
void Calculate_Dpp();
#endif

#ifdef DPP_BRUNETTI_07
void Init_Brunetti_07();
#endif

#ifdef Q_SIMPLE_SECONDARIES
void Init_Simple_Secondaries();
#endif

#ifdef Q_BRUNETTI_05
void Init_Brunetti05();
#endif


#ifdef Q_SHOCK_WITH_HADRONIC_BACKGROUND
void Init_Hadronic_Background();
#endif
