/* Timings & Benchmarking */

#include "common.h"

char TimeMarkDescr[CPU_LASTENTRY][20] = {
    "Solver         ",
    "Main           ",
};

void init_timing()
{
	long i;

	for (i = 0; i < CPU_LASTENTRY; i++) {
		Cpu.Total[i] = 0;
		Cpu.Start[i] = -1;
	}
	Cpu.Zero = get_current_time();
	return;
}

void finish_timing()
{
	long mark = 0;
	double myTime = 0, timeTotal = 0;

	if (ThisTask.Rank == 0) {
		fprintf(stdout, "\n Timing Statistics [s] on <%d> CPUs : \n"
			"----------------------------\n", ThisTask.NTask);
	}

	for (mark = 0; mark < CPU_LASTENTRY; mark++) {
		myTime = Cpu.Total[mark];
		timeTotal += Cpu.Total[mark];
		if (ThisTask.Rank == 0)
			fprintf(stdout, "%s	%g		\n",
				TimeMarkDescr[mark], myTime);
	}
	if (ThisTask.Rank == 0) {
		fprintf(stdout, "----------------------------\n"
			"Total		%g		 	        \n",
			timeTotal);
	}

	if (ThisTask.Rank == 0)
		fprintf(stdout, "\n");

	return;
}

void start_timing(enum TimeMarks mark)
{
	if (Cpu.Start[mark] > 0) {
		fprintf(stderr, "Timing for Mark %i already running \n",
			(int)mark);
	} else {
		Cpu.Start[mark] = get_current_time();
	}
	return;
}

void stop_timing(enum TimeMarks mark)
{
	if (Cpu.Start[mark] < 0) {
		fprintf(stderr, "Timing for Mark %i wasn't started \n",
			(int)mark);
	} else {
		Cpu.Total[mark] += get_current_time() - Cpu.Start[mark];
	}
	return;
}

double get_current_time()
{
	return MPI_Wtime();
}
