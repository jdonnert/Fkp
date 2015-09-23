#include "common.h"
#include "modules/modules.h"
#include "compress.h"

#define nGrid N_SPEC_BINS

static void set_interpolation_coefficients(const int ipart,double *a,double *b);
static void interpolate(const double *a, const double *b, const double a_time, 
	 const int ipart, struct Gas_Data *SphP_i);

void solver()
{
    const double aFac_low = 1/(log10(p[LowIdx]/p[LowIdx+1]));
    const double bFac_low = log10(p[LowIdx+1]);
    const double aFac_high = 1/(log10(p[HighIdx]/p[HighIdx-1]));
    const double bFac_high = log10(p[HighIdx-1]);

    double log_p[N_SPEC_BINS] = { 0 };

    for (int i = 0; i < nGrid; i++)
        log_p[i] = log10(p[i]);

    const double tmin = a2t(Snap.LastTime);  // integration in cgs 
    const double tmax = a2t(Snap.Time);
    
    rprintf("Evolving Spectra for <%g> Myr (%s)\n", 
            (tmax-tmin) /yr2sec/1e6,Cosmo.name); fflush(stdout);

	const int counter_const = floor(ThisTask.Npart[0]/Omp.NThreads/10);

	#pragma omp parallel for \
		schedule(dynamic, ThisTask.Npart[0]/Omp.NThreads/64)
    for (int ipart = 0; ipart < ThisTask.Npart[0]; ipart++) {

		if (ipart % counter_const == 0)
			rprintf(".");

		double np[nGrid] = { -DBL_MAX };

		double a[32] = { 0 }, b[32] = { 0 };

		set_interpolation_coefficients(ipart, a, b);
		
#if defined(COMPRESSION_INTERNAL) 
        Uncompress(LastSphP[ipart].Ncre_total,LastSphP[ipart].CReSpectrum,np);
#else
        for (int i = LowIdx; i < HighIdx + 1; i++) 
            np[i] = pow(10, LastSphP[ipart].Ncre[i]);
#endif

#ifndef VARIABLE_TIMESTEPS 
        float dt = timesteps[ipart];    // timesteps in cgs already
#else
        float dt = TimeBase * (1ULL << TimeBins[ipart]);
#endif

		float a_dt = t2a(dt);

        double t = tmin; 
		double a_time = Snap.LastTime; // cosmological expansion factor

        double ntot = 0;

		while ( t < tmax ) { 
			
			struct Gas_Data SphP_i = { 0 };

			interpolate(a, b, a_time, ipart, &SphP_i);
            
            if ( (np[LowIdx] > 0) && (LowIdx > 2) ) { // boundary extrapolation
                
                double a = (log10(np[LowIdx]/np[LowIdx+1]) ) * aFac_low;
                
                double b = log10(np[LowIdx+1]) - a * bFac_low;
                
                for (int i = 1; i < LowIdx; i++)
                    np[i] = pow(10, a*log_p[i]+ b);
            }

            if ( (np[HighIdx] > 0) && (HighIdx < nGrid-2) ) {
                
                double a = (log10(np[HighIdx]/np[HighIdx-1]) ) * aFac_high;
                
                double b = log10(np[HighIdx-1]) - a*bFac_high;
                
                for (int i = HighIdx; i < nGrid; i++)
                    np[i] = pow(10, a*log_p[i]+ b);
            }

			double  Dpp[nGrid] = { 0 }, Hp[nGrid] = { 0 }, Q[nGrid] = { 0 }, 
					T[nGrid] = { 0 };

			Reacceleration(SphP_i, q, nGrid, Dpp); //get coefficients
			
			Cooling(SphP_i, q, nGrid, Dpp, ipart, Hp);
			
			Injection(SphP_i, p, nGrid, dt, Q);
			
			Escape(SphP_i,  p, nGrid, T);
			   
//if (P[ipart].ID == 1607043)
 //  for (int i = 0; i < nGrid; i++)
//		printf("%d p=%g np=%g Dpp=%g Hp=%g Q=%g M=%g dt=%g\n", 
	//			i, p[i], np[i], Dpp[i], Hp[i], Q[i], SphP_i.Mach, dt);

            /* no-flux boundary conditions at edge of comp domain */
            Dpp[0] = Hp[0] = 0;
            Dpp[nGrid-1] = Hp[nGrid-1] = 0;
    
            /* Park & Petrosian 1995, Chang & Cooper 1970 */
			double Wp[nGrid] = { 0 }, Wm[nGrid] = { 0 };

            for (int i = 1; i < nGrid; i++) { // w, Wp, Wm @ i+1/2 
                
                double w = Hp[i]/Dpp[i]*dq[i];

				w = fmin(700, w); // overflow ?, Hp dominant
				
				w = fsign(w) * fmax(1e-8, fabs(w)); // underflow ?, Dpp dom.

                double exp_w  = exp(w);
                
                Wp[i] = 1 / ( 1.0 - 1.0/exp_w ); // forward weight

                Wm[i] = 1 / ( exp_w - 1.0 ); // backward weight
            }

            Wp[nGrid-1] = Wm[nGrid-1] = 0;

            double A[nGrid] = { 0 }, B[nGrid] = { 0 }, C[nGrid] = { 0 }, 
				O[nGrid] = { 0 };

            for (int i = 1; i < nGrid; i++) {

                A[i] = -dt/dp[i] * Hp[i] * Wp[i];

                C[i] = -dt/dp[i] * Hp[i-1] * Wm[i-1];

                B[i] = 1.0 + dt/T[i] + dt/dp[i] * (Hp[i-1] * Wp[i-1] 
												 + Hp[i] * Wm[i]) ;
                O[i] = np[i] + dt * Q[i];
            }
            
			double E[nGrid] = { 0 }, F[nGrid] = { 0 };

            for (int i = 1; i < nGrid-1; i++) { // TriDiagonal MAtrix Solver 

                E[i] = (A[i])/(B[i] - C[i]*E[i-1]);
                
                F[i] = (O[i] - C[i]*F[i-1])/(B[i] - C[i]*E[i-1]);
            }
			
            np[nGrid-1] = O[nGrid-1];

            for (int i = 1; i < nGrid; i++) {  // Update Spectrum 

                int j = nGrid - 1 - i;     

                np[j] = F[j] - E[j]*np[j+1];

                ntot += np[j] * dp[j];

            }

            t += dt;
    		a_time += a_dt; 
        }
        
        if (ntot < 1e-100) { // prevent underflow 

            ntot = 0;
            
			memset(np, 0, nGrid*sizeof(*np));
        }
        
        for (int i = 0; i < LowIdx; i++) // clean boundaries 
            np[i] = 0;

        for (int i = HighIdx; i < N_SPEC_BINS; i++)
            np[i] = 0;
         
#ifdef COMPRESSION_INTERNAL  

        Compress(np, &SphP[ipart].Ncre_total, SphP[ipart].CReSpectrum);

#else

        for (int i = 0; i < nGrid; i++) // save spectrum
           SphP[ipart].Ncre[i] = log10(np[i]);

#endif // COMPRESSION_INTERNAL

    } // for(ipart)
    
	MPI_Barrier(MPI_COMM_WORLD);
    
	return;
}

/* 
 * These two functions time interpolate a particle between two snapshots.
 * Because of finite precision, the interpolation from or two the bracketing
 * values can cause problems. Thus we copy the particle in this case
 */

static void set_interpolation_coefficients(const int ipart, double *a, 
		double *b)
{
	const double a_dt = Snap.Time - Snap.LastTime;

#ifdef Q_SHOCK_PRIMARIES

	if (!isfinite(SphP[ipart].Shock_Density)) // correct shock finder
		SphP[ipart].Shock_Density = 0;

	if (!isfinite(LastSphP[ipart].Shock_Density))
		LastSphP[ipart].Shock_Density = 0;

	if (!isfinite(SphP[ipart].Shock_Velocity))
		SphP[ipart].Shock_Velocity = 0;

	if (!isfinite(LastSphP[ipart].Shock_Velocity))
		LastSphP[ipart].Shock_Velocity = 0;

#endif // Q_SHOCK_PRIMARIES

	if (!isfinite(SphP[ipart].Mach))
		SphP[ipart].Mach = 0;

	if (!isfinite(LastSphP[ipart].Mach))
		LastSphP[ipart].Mach = 0;
	
//	double f_shear = fabs(SphP[ipart].DivVel) 
//		/ (fabs(SphP[ipart].DivVel) + fabs(SphP[ipart].CurlVel));

//	if (f_shear < 0.9) // filter criteria
//		SphP[ipart].Mach = 0;

//	if (SphP[ipart].Mach > 3) {
	
//		double c_s = adiabatic_index * (adiabatic_index-1) * SphP[ipart].U;

//		SphP[ipart].Mach = 2*SphP[ipart].Hsml*fabs(SphP[ipart].DivVel)
//			/c_s;
//	}
	
	int i = 0;

	a[i] = (SphP[ipart].Rho - LastSphP[ipart].Rho)/a_dt;
	b[i] = LastSphP[ipart].Rho - a[i++] * Snap.LastTime;
	
	a[i] = (SphP[ipart].Hsml - LastSphP[ipart].Hsml)/a_dt;
	b[i] = LastSphP[ipart].Hsml - a[i++] * Snap.LastTime;

	a[i] = (SphP[ipart].U - LastSphP[ipart].U)/a_dt;
	b[i] = LastSphP[ipart].U - a[i++] * Snap.LastTime;

	a[i] = (SphP[ipart].Bfld[0] - LastSphP[ipart].Bfld[0])/a_dt;
	b[i] = LastSphP[ipart].Bfld[0] - a[i++] * Snap.LastTime;
	
	a[i] = (SphP[ipart].Bfld[1] - LastSphP[ipart].Bfld[1])/a_dt;
	b[i] = LastSphP[ipart].Bfld[1] - a[i++] * Snap.LastTime;

	a[i] = (SphP[ipart].Bfld[2] - LastSphP[ipart].Bfld[2])/a_dt;
	b[i] = LastSphP[ipart].Bfld[2] - a[i++] * Snap.LastTime;

	a[i] = (SphP[ipart].VRms - LastSphP[ipart].VRms)/a_dt;
	b[i] = LastSphP[ipart].VRms - a[i++] * Snap.LastTime;

	a[i] = (SphP[ipart].Dpp - LastSphP[ipart].Dpp)/a_dt;
	b[i] = LastSphP[ipart].Dpp - a[i++] * Snap.LastTime;

	a[i] = (SphP[ipart].DRho - LastSphP[ipart].DRho)/a_dt;
	b[i] = LastSphP[ipart].DRho - a[i++] * Snap.LastTime;

#ifdef Q_SHOCK_PRIMARIES

	a[i] = (SphP[ipart].Mach - LastSphP[ipart].Mach)/a_dt;
	b[i] = LastSphP[ipart].Mach - a[i++] * Snap.LastTime;
	
	a[i] = (SphP[ipart].Shock_Velocity - LastSphP[ipart].Shock_Velocity)/a_dt;
	b[i] = LastSphP[ipart].Shock_Velocity - a[i++] * Snap.LastTime;

	a[i] = (SphP[ipart].Shock_Compression_Ratio 
			- LastSphP[ipart].Shock_Compression_Ratio)/a_dt;
	b[i] = LastSphP[ipart].Shock_Compression_Ratio - a[i++] * Snap.LastTime;

	a[i] = (SphP[ipart].Shock_Pressure - LastSphP[ipart].Shock_Pressure)/a_dt;
	b[i] = LastSphP[ipart].Shock_Pressure - a[i++] * Snap.LastTime;

	a[i] = (SphP[ipart].Shock_Density - LastSphP[ipart].Shock_Density)/a_dt;
	b[i] = LastSphP[ipart].Shock_Density - a[i++] * Snap.LastTime;
	
	a[i] = (SphP[ipart].DivVel - LastSphP[ipart].DivVel)/a_dt;
	b[i] = LastSphP[ipart].DivVel - a[i++] * Snap.LastTime;

	a[i] = (SphP[ipart].CurlVel - LastSphP[ipart].CurlVel)/a_dt;
	b[i] = LastSphP[ipart].CurlVel - a[i++] * Snap.LastTime;

	a[i] = (SphP[ipart].Vel[0] - LastSphP[ipart].Vel[0])/a_dt;
	b[i] = LastSphP[ipart].Vel[0] - a[i++] * Snap.LastTime;

	a[i] = (SphP[ipart].Vel[1] - LastSphP[ipart].Vel[1])/a_dt;
	b[i] = LastSphP[ipart].Vel[1] - a[i++] * Snap.LastTime;

	a[i] = (SphP[ipart].Vel[2] - LastSphP[ipart].Vel[2])/a_dt;
	b[i] = LastSphP[ipart].Vel[2] - a[i++] * Snap.LastTime;

	a[i] = (SphP[ipart].Shock_Normal[0] - LastSphP[ipart].Shock_Normal[0])/a_dt;
	b[i] = LastSphP[ipart].Shock_Normal[0] - a[i++] * Snap.LastTime;

	a[i] = (SphP[ipart].Shock_Normal[1] - LastSphP[ipart].Shock_Normal[1])/a_dt;
	b[i] = LastSphP[ipart].Shock_Normal[1] - a[i++] * Snap.LastTime;
	
	a[i] = (SphP[ipart].Shock_Normal[2] - LastSphP[ipart].Shock_Normal[2])/a_dt;
	b[i] = LastSphP[ipart].Shock_Normal[2] - a[i++] * Snap.LastTime;

#endif // Q_SHOCK_PRIMARIES

	return ;
}

static void interpolate(const double *a, const double *b, const double a_time, 
		const int ipart, struct Gas_Data *SphP_i)
{
	if (a_time == Snap.LastTime) { // avoid roundoff with interpolation to 0
	
		memcpy(SphP_i, &LastSphP[ipart], sizeof(*LastSphP));

		return;
	}

	if (a_time == Snap.Time) {
	
		memcpy(SphP_i, &SphP[ipart], sizeof(*SphP));

		return;
	}

	int i = 0;

	SphP_i->ID = P[ipart].ID; // SphP carries ID with it ...

	SphP_i->Rho = a[i] * a_time + b[i++];
	SphP_i->Hsml = a[i] * a_time + b[i++];
	SphP_i->U = a[i] * a_time + b[i++];
	SphP_i->Bfld[0] = a[i] * a_time + b[i++];
	SphP_i->Bfld[1] = a[i] * a_time + b[i++];
	SphP_i->Bfld[2] = a[i] * a_time + b[i++];
	SphP_i->VRms = a[i] * a_time + b[i++];
	SphP_i->Dpp = fmax(0, a[i] * a_time + b[i++]); 

	SphP_i->DRho = a[i] * a_time + b[i++];

#ifdef Q_SHOCK_PRIMARIES

	SphP_i->Mach = a[i] * a_time + b[i++];
	SphP_i->Shock_Velocity = a[i] * a_time + b[i++];
	SphP_i->Shock_Compression_Ratio = a[i] * a_time + b[i++];
	SphP_i->Shock_Pressure = a[i] * a_time + b[i++];

	SphP_i->Shock_Density = a[i] * a_time + b[i++];
	SphP_i->DivVel = a[i] * a_time + b[i++];
	SphP_i->CurlVel = a[i] * a_time + b[i++];
	
	SphP_i->Vel[0] = a[i] * a_time + b[i++];
	SphP_i->Vel[1] = a[i] * a_time + b[i++];
	SphP_i->Vel[2] = a[i] * a_time + b[i++];

	SphP_i->Shock_Normal[0] = a[i] * a_time + b[i++];
	SphP_i->Shock_Normal[1] = a[i] * a_time + b[i++];
	SphP_i->Shock_Normal[2] = a[i] * a_time + b[i++];

#endif


/*	printf("rho=%g h=%g u=%g b0=%g b1=%g b2=%g \n"
			"vrms=%g dpp=%g  M=%g SHVel=%g sigma=%g shpr=%g shrho=%g \n",
			SphP_i->Rho, SphP_i->Hsml, SphP_i->U, SphP_i->Bfld[0], 
			SphP_i->Bfld[1], SphP_i->Bfld[2], SphP_i->VRms, SphP_i->Dpp, 
			SphP_i->Mach, SphP_i->Shock_Velocity, 
			SphP_i->Shock_Compression_Ratio, SphP_i->Shock_Pressure, 
			SphP_i->Shock_Density ); */

	return ;
}

#undef nGrid
