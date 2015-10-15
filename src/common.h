/* Main Variables Header. Every thing global dump in here.
 * */

#include <stdlib.h>             // system       
#include <stdio.h>
#include <stdint.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <limits.h>

#include <mpi.h>
#include <omp.h>

/* gsl */
#include <gsl/gsl_math.h>
#include <gsl/gsl_const_cgsm.h>
#include <gsl/gsl_const_num.h>

/* Where all the defines are */
#include "config.h"

/* We need those everywhere */
#include "unit.h"
#include "cosmo.h"
#include "constants.h"
#include "timing.h"
#include "sort.h"

/* code parameters */
#define MAXLINELENGTH 256L

#define N_SPEC_BINS 128L
#define SPECSIZE_BYTES 60L

#define N_PART_TYPES 6L 

/* mathematical constants */
#define pi 			M_PI
#define sqrt2		M_SQRT2
#define sqrt3       1.73205081

/* MACROS */
#define rprintf if(!ThisTask.Rank && !Omp.ThreadID) printf

#define Assert(...) Assert_Info(__func__, __FILE__, __LINE__, __VA_ARGS__)
#define Malloc(x) malloc_info(x, __FUNCTION__, __FILE__, __LINE__)
#define Realloc(x,y) realloc_info(x, y, __FUNCTION__, __FILE__, __LINE__)
#define Free(x) free_info(x)

#define min(a,b) ((a)>(b)?(b):(a))
#define max(a,b) ((a)<(b)?(b):(a))

#define p2(a) ((a)*(a)) 
#define p3(a) ((a)*(a)*(a))

#define length3(a) sqrt(p2(a[0]) + p2(a[1]) + p2(a[2]))
#define fsignf(a) ((0.0f < a) - (a < 0.0f))
#define fsign(a) ((0.0 < a) - (a < 0.0))

/* 
 * MAIN VARIABLES 
 * */
extern struct ParallelInfos {
	int Rank;		/* Rank of local Processor */
	int NTask;		/* Number of Processors */
	long long PartTotal;	/* Total Npart on this Processor */
	long long Npart[N_PART_TYPES];	/* Npart of species stored locally */
} ThisTask;

extern struct OpenMP_infos{
    int NThreads;       // Number of openMP threads 
    int ThreadID;       // Thread ID of this thread 
}Omp;
#pragma omp threadprivate(Omp)

extern struct Parameters {	/* parameters from par file */
    int N_IOTasks;		    /* number of read tasks */
    int Cosmology;          /* cosmology to use */
    char Input_File[MAXLINELENGTH]; 
    char Output_File[MAXLINELENGTH];
    double EndTime;            /* redshift to finish at */
    double Pmin;            /* minimum of momentum grid */
    double Pmax;            /* maximum of momentum grid */
    double Plow;            /* maximum of lower boundary region */
    double Phigh;           /* minimum of upper boundary region */
    double Dpp_spec_idx;    /* spectral index for power law Dpp */
    double a_crp;           /* spectral index in sec injection Q */
    double p0_crp;          /* spectral cutoff in sec injection Q */
    double X_crp;           /* fraction rel to thermain in sec Q */
    int Flag_Start;         // run mode
#ifdef Q_SHOCK_PRIMARIES
	double Xi_ep;			// electron to proton fraction
#endif
} Param;

extern struct SnapProperties {		
	int SnapNum;	
    int LastSnapNum;
	double Boxsize;						
	double Redshift;
    double LastRedshift;
	double Time;
    double LastTime;
	long long Npart[N_PART_TYPES];	    /* over all processes */
	long long PartTotal;	            /* over all processes */
	double Masstab[N_PART_TYPES];
    int Have_Arepo;
} Snap;

extern struct Particle_Data { /*struct holding data */
    //float Pos[3];             /* Position vector */
    //float Vel[3];             /* Velocity vector */
	unsigned long ID;       /* Unique ID */
	int Type;               /* Type */
	float Mass;             /* Mass */
} *P, *LastP;

extern struct Gas_Data {
	int ID;
    float Vel[3];             /* Velocity vector */
    float Hsml;             /* Smoothing Length */
	float Rho;				// Density 
	float DRho;				// change in density
	float U;				/* Internal Energy, add below ! */
	float Bfld[3];	        /* Magnetic Field */
	float VRms;				/* RMS around mean  */
    float DivVel;           /* Velocity Divergence */
    float CurlVel;          /* Velocity Curl */
#ifdef READ_MORE_INFO
	float VTurb;			/* Turbulent Velocity */
	float VBulk[3];	        /* Mean Velocity */
#endif
	float TNgb;				// True Number of Neighbours 
    double Dpp;              // Dpp/p^(Dpp_spec_idx) 
#ifdef COMPRESSION_INTERNAL
    float Ncre_total;        // CR electron total number density
    char CReSpectrum[SPECSIZE_BYTES];// CR electron spectrum 
#else
    float Ncre[N_SPEC_BINS]; // CR electron spectrum 
#endif
	float Mach;
#ifdef Q_SHOCK_PRIMARIES
	float Shock_Velocity;	// upwind
	float Shock_Compression_Ratio; 
	float Shock_Pressure; 	// upwind
	float Shock_Density; 	// upwind
	float Shock_Normal[3];
#endif
} *SphP, *LastSphP;

extern double p[N_SPEC_BINS];            // momentum grid 
extern double q[N_SPEC_BINS];            // momentum grid at midpoint
extern double dp[N_SPEC_BINS];           // cellsize @ gridpoint 
extern double dq[N_SPEC_BINS];           // cellsize @ midpoint 

#ifdef VARIABLE_TIMESTEPS
extern short *TimeBins;       /* dt = TimeBase * (1ULL<<TimeBins[ipart]) */
extern float TimeBase;          /* = ( Time1 - Time0 ) / (1ULL<<62) */
#endif

extern float *timesteps;

extern size_t LowIdx,HighIdx;/* boundary indizes */

extern int Flag, Ipart;

/* functions */
void setup();
void Init_spectrum();

void solver();

void calc_timesteps();

void sortP_global();

#ifdef COMPRESSION
void Compress(const int, const float*, float *, char *);
void Uncompress(const float, const char *, double*);
void Setup_compression();
#endif

/* These replace my_malloc etc via macro */
void *malloc_info(size_t size,
        const char* file, const char* func, const int line);
void *realloc_info(void *ptr, size_t size, 
        const char* file, const char* func, const int line);
void Assert_Info(const char *, const char *, int, int64_t, const char *, ...);
void free_info(void *ptr);

extern void Reallocate_P(long long, long long *, int);
extern void Shift_Particle_Data();

