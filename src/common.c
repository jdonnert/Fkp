#include "common.h"

struct Parameters Param;
struct ParallelInfos ThisTask;
struct SnapProperties Snap = {0};
struct Particle_Data *P=NULL, *LastP=NULL;
struct Gas_Data *SphP=NULL, *LastSphP=NULL;

struct ParallelInfos ThisTask = { 0, 0, 0, {0} };
#pragma omp threadprivate(Omp)
struct OpenMP_infos Omp = { 0 };

short *TimeBins = NULL;
float TimeBase = 0; 
float *timesteps = NULL;

size_t LowIdx = 0, HighIdx = 0;

double p[N_SPEC_BINS] = { 0 }, q[N_SPEC_BINS] = { 0 };
double dp[N_SPEC_BINS] = { 0 }, dq[N_SPEC_BINS] = { 0 }; 

#pragma omp threadprivate(Flag)
int Flag = 0, Ipart = 0;


/* 
 * Memory Management 
 */

void *malloc_info(size_t size, const char* file
        , const char* func, const int line)
{
    char errmsg[MAXLINELENGTH];

	void *result = malloc(size);

    sprintf(errmsg, "Allocation failed, <%zu> bytes at %s in %s line %d \n"
        ,size, file, func, line);
    Assert(result != NULL, errmsg);

	return (result);
}

void *realloc_info(void *ptr, size_t size, const char* file
        , const char* func, const int line)
{
    char errmsg[MAXLINELENGTH];

	void *result = realloc(ptr, size);

    sprintf(errmsg, "Reallocation failed, <%zu> bytes at %s in %s line %d \n"
        ,size, file, func, line);
    Assert(result != NULL && size != 0, errmsg);

	return (result);
}

void free_info(void *ptr) 
{
    Assert(ptr != NULL, "Can't free a NULL pointer");
        
    free(ptr);
    
    return;
}

/* 
 * Error Handling, we use variable arguments to be able
 * to print more informative error messages 
 */

void Assert_Info(const char *func, const char *file, int line,
		int64_t expr, const char *errmsg, ...)
{
    if (expr)
        return;

	va_list varArgList;

	va_start(varArgList, errmsg);

	/* we fucked up, tell them */
    fprintf(stderr, 
			"\nERROR Task %d: In file %s, function %s(), line %d :\n\n	", 
			ThisTask.Rank, file, func, line);

	vfprintf(stderr, errmsg, varArgList); 
	
	fprintf(stderr, "\n\n"); 
	
	fflush(stderr);

	va_end(varArgList);

    MPI_Abort(MPI_COMM_WORLD, -1); // finish him ...

    exit(EXIT_FAILURE); // ... fatality

    return;
}

/* 
 * Particle Metrics 
 */

extern void Reallocate_P(long long partTotal, 
        long long nPart[N_PART_TYPES], int sign)
{
	if (partTotal == 0)
		return;

	ThisTask.PartTotal += sign * partTotal;

	for (int type = 0; type < N_PART_TYPES; type++) {

		ThisTask.Npart[type] += sign * nPart[type];
		
        Assert(ThisTask.Npart[type] >= 0, 
                "Can't allocate negative particles");
	}

	P = Realloc((void *)P, sizeof(*P)*ThisTask.PartTotal);

	SphP = Realloc((void *)SphP, sizeof(*SphP)*ThisTask.Npart[0]);

	return;
}

extern void Shift_Particle_Data()
{
	rprintf("Shifting Particle data ...");

    Snap.LastSnapNum = Snap.SnapNum;
    Snap.LastTime = Snap.Time;
    Snap.LastRedshift = Snap.Redshift;
    
    struct Particle_Data *tmp_P = LastP;    // save the pointer, can be NULL
    LastP = P;
    P = tmp_P;

    struct Gas_Data *tmp_SphP = LastSphP;
    LastSphP = SphP;
    SphP = tmp_SphP;

    if (P != NULL)    // clear memory
        memset(P, 0, ThisTask.Npart[0] * sizeof(*P));
    
    if (SphP != NULL)
        memset(SphP, 0, ThisTask.Npart[0] * sizeof(*SphP));
   
    ThisTask.Npart[0] = 0;  // these are set at readin 
    ThisTask.PartTotal = 0;
    Snap.Npart[0] = 0;
    Snap.PartTotal = 0;

	rprintf("done\n");

    return;
}

