/* PROFILING
 * We define marks in this enum,
 * which are used to let the start
 * and stop timing functions keep
 * track of what is measured.
 */

enum TimeMarks {
	CPU_SOLVER,
    CPU_MAIN,
	CPU_LASTENTRY		/* Keep this entry at the end */
} mark;

struct Timings {
	double Zero;		/* Absolute Zeropoint */
	double Total[CPU_LASTENTRY];	/* Total Time for mark */
	double Start[CPU_LASTENTRY];	/* Start time for mark */
} Cpu;

void init_timing();
void finish_timing();
void start_timing(enum TimeMarks mark);
void stop_timing(enum TimeMarks mark);
double get_current_time();
