#include "common.h"
#include "modules/modules.h"

#ifndef VARIABLE_TIMESTEPS
#ifndef I_NONE
void calc_timesteps()
{
    const int npart = ThisTask.Npart[0];

#pragma omp parallel for
	for (int i = 0; i < npart; i++)
		SphP[i].DRho = SphP[i].Rho - LastSphP[i].Rho;

    double dtmax = a2t(Snap.Time) - a2t(Snap.LastTime);

    if (timesteps == NULL)
        timesteps = Malloc(npart*sizeof(*timesteps));
    
    for (size_t ipart=0; ipart<npart; ipart++)
        timesteps[ipart] = 0.5*yr2sec*1e6;
    
    rprintf("Setting global timestep <%4.2g> Myr\n\n"
                ,timesteps[0]/yr2sec/1e6);
    return;
}    
#endif // !I_NONE 

#ifdef I_NONE
void calc_timesteps()
{
    size_t ipart;
    double dtmax;

    if (timesteps == NULL)
        timesteps = Malloc(ThisTask.Npart[0]*sizeof(*timesteps));

    dtmax = Snap.Time - Snap.LastTime;
    
    for ( ipart=0; ipart<ThisTask.Npart[0]; ipart++)
        timesteps[ipart] = dtmax / 1e6;

    printf("Setting global timestep %g \n\n", dtmax / 5e2);
    
    return;
}       
#endif // I_NONE
#endif // !VARIABLE_TIMESTEPS

#ifdef VARIABLE_TIMESTEPS
void calc_timesteps()
{
    const int npart = ThisTask.Npart[0];

#pragma omp parallel for
	for (int i = 0; i < npart; i++)
		SphP[i].DRho = SphP[i].Rho - LastSphP[i].Rho;

    const double t2Myr = 1.0/yr2sec/1e6;
    double thisTime = a2t(Snap.Time); // [sec]
    double lastTime = a2t(Snap.LastTime); // [sec]

    double dtMax = thisTime - lastTime; // [sec]
    TimeBase =  dtMax / (1UL<<62); // dyn. range >1e9

    if (TimeBins == NULL) 
        TimeBins = Malloc(npart * sizeof(*TimeBins));
    
	int bin_max = 0;
    int bin_min = INT_MAX;
    uint64_t bin_counter[63] = { 0 };

	rprintf("Finding Timesteps ... ");

#pragma omp parallel for reduction(max:bin_max) reduction(min:bin_min)
    for (int ipart = 0; ipart < npart; ipart++) {

        int imax = HighIdx; // find largest occupied bin

        while (imax > LowIdx) {   

            imax--;

            if (LastSphP[ipart].Ncre[imax] > 1) 
                break;
        }

        if (imax == 0)
            imax = N_SPEC_BINS*0.5;
            
        double dpp = fmax(SphP[ipart].Dpp, LastSphP[ipart].Dpp);
        double dt_Dpp = 1/(4.0*dpp);

		double hp_last[N_SPEC_BINS] = { 0 };
		double hp_this[N_SPEC_BINS] = { 0 };
		double dpp_zero[N_SPEC_BINS] = { 0 }; // just a dummy

		Cooling(LastSphP[ipart], p, N_SPEC_BINS, dpp_zero, ipart, hp_last);
		Cooling(SphP[ipart], p, N_SPEC_BINS, dpp_zero, ipart,  hp_this);
        
		double hp_low = fmax(hp_last[LowIdx], hp_this[LowIdx]);
        double hp_high = fmax(hp_last[imax], hp_this[imax]);

        double dt_Hp = fmin(p[LowIdx]/hp_low, p[imax]/hp_high);
        double dt = fmin(dt_Dpp, dt_Hp);

        TimeBins[ipart] = 62 - max(0 , ceil(log2(dtMax/dt))) - 2;
        
        Assert(TimeBins[ipart] > 0, 
				"TimeBin out of dynamic range : ipart=%d dtDpp=%g,"
				" dtHp=%g dt=%g imax=%d %g %g ", ipart, dt_Dpp, dt_Hp, dt, 
				imax, hp_high, hp_low);

		#pragma omp atomic
        bin_counter[TimeBins[ipart]]++; // do some statistics

        if (TimeBins[ipart] > bin_max)
            bin_max = TimeBins[ipart];
        
        if (TimeBins[ipart] < bin_min)
            bin_min = TimeBins[ipart];
    }

    /* print statistics */
    int global_bin_max = 0;
    MPI_Reduce(&bin_max, &global_bin_max, 1, MPI_INT, MPI_MAX, 0, 
            MPI_COMM_WORLD);

    int global_bin_min = 0;
    MPI_Reduce(&bin_min, &global_bin_min, 1, MPI_INT, MPI_MIN, 0, 
            MPI_COMM_WORLD);

    long long global_counter[63] = { 0 };
    MPI_Reduce(bin_counter, global_counter, 63, MPI_LONG_LONG, MPI_SUM, 0, 
            MPI_COMM_WORLD);
    
    rprintf("\nSnapshot %d -> %d, a=%g -> %g \n\n", Snap.LastSnapNum, 
            Snap.SnapNum, Snap.LastTime, Snap.Time);
    
    rprintf("   TimeBin        Npart     dt [Myr]  \n");

    for (int i = global_bin_max; i > global_bin_min-1; i--)
        rprintf("    %3d      %10lli     %8.2g \n", i, 
				global_counter[i], TimeBase * (1ULL<<i)*t2Myr );

    rprintf("\n");
    
    return;
}

#endif // VARIABLE_TIMESTEPS

