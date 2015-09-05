#include "common.h"
#include "io/io.h"
#include "modules/modules.h"

void print_compile_time_settings();
void preamble(int, char**);

int main(int argc, char *argv[])
{
    preamble(argc, argv);

    setup();
    
	Snap.Time = read_input();

    Assert(Snap.Time < Param.EndTime, "First snapshot later than EndTime "
			"%g <> %g", Snap.Time, Param.EndTime);

    sortP_global();    

    Init_spectrum();
    
    start_timing(CPU_MAIN);

    for (;;) { 

        Shift_Particle_Data();
            
        Snap.Time = read_input();
    
        if (Snap.Time > Param.EndTime)
            break; // End program 

        sortP_global();    
        
		calc_timesteps(); 

        solver();
        
		write_spectra();
    }
    
    rprintf("\nEnd of Simulation reached \n\n");

    stop_timing(CPU_MAIN);

    finish_timing();

    MPI_Finalize();

    return EXIT_SUCCESS;
}

void preamble(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask.Rank);
	MPI_Comm_size(MPI_COMM_WORLD, &ThisTask.NTask);

	#pragma omp parallel
    {
    Omp.ThreadID = omp_get_thread_num();
    Omp.NThreads = omp_get_num_threads();
    }

  	rprintf("\n*** This is P-Fkp, a parallel Fokker Planck Solver ***\n");
	
	#pragma omp parallel
	rprintf("\nRunning on <%d> nodes  and <%d> threads per node\n\n", 
				ThisTask.NTask, Omp.NThreads);

    print_compile_time_settings();

    if ( !(argc == 2 || argc == 4) && !ThisTask.Rank) {
        printf("\nUsage   : ./P-Fkp $PARAMETERFILE [Flag] [SnapNum]\n"
            "Flag    : 1 = Resume, 2 = Convert Spectrum File\n"
            "SnapNum : Snapshot number to restart/convert from\n");

        Assert(0, "Program Call incomplete");
    }

    init_timing();

    read_param_file(argv[1]);

    if (argc == 4) {
        
        Param.Flag_Start = atoi(argv[2]);

        Snap.SnapNum = atoi(argv[3]);

        sprintf(Param.Input_File, "snap_%03d", Snap.SnapNum);

        if (Param.Flag_Start == 1)
            rprintf("\nResuming simulation at Snapshot %d \n",Snap.SnapNum);
        else if (Param.Flag_Start == 2)
            rprintf("\nConverting Spectrum File No: %d \n",Snap.SnapNum);
    }

    return ;
}
