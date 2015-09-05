/* Common input functions
 * */
#include "../common.h"
#include "io.h"

#define REAL 1
#define STRING 2
#define INT 3

int id[MAXTAGS] = {0};
void *addr[MAXTAGS] = {NULL};
char tag[MAXTAGS][50],comment[MAXTAGS][50];

/* Reads a number of tags from an ascii file
 * the comment sign is  %
 * */
void read_param_file(char *filename)
{

	FILE *fd=NULL;
	char buf[MAXLINELENGTH], buf1[MAXLINELENGTH];
	char buf2[MAXLINELENGTH],buf3[2 * MAXLINELENGTH];
	int tagDone[MAXTAGS]= {0};
	int i, j, nt=0;
	
    /* read parameter file on process 0 */
    if (ThisTask.Rank == 0) {

        strcpy(tag[nt], "NIOTasks");
		strcpy(comment[nt], "Number of files read in parallel");
		addr[nt] = &Param.N_IOTasks;
		id[nt++] = INT;
        
        strcpy(tag[nt], "Cosmo");
		strcpy(comment[nt], "Cosmology Flag");
		addr[nt] = &Param.Cosmology;
		id[nt++] = INT;
        
        strcpy(tag[nt], "Input_File");
		strcpy(comment[nt], "Input File Name");
		addr[nt] = &Param.Input_File;
		id[nt++] = STRING;
        
        strcpy(tag[nt], "Output_File");
		strcpy(comment[nt], "Output File Name");
		addr[nt] = &Param.Output_File;
		id[nt++] = STRING;

        strcpy(tag[nt], "EndTime");
		strcpy(comment[nt], "Minimum Redshift");
		addr[nt] = &Param.EndTime;
		id[nt++] = REAL;

        strcpy(tag[nt], "Pmin");
		strcpy(comment[nt], "minimum of momentum grid");
		addr[nt] = &Param.Pmin;
		id[nt++] = REAL;

        strcpy(tag[nt], "Pmax");
		strcpy(comment[nt], "maximum of momentum grid");
		addr[nt] = &Param.Pmax;
		id[nt++] = REAL;
        
        strcpy(tag[nt], "Plow");
		strcpy(comment[nt], "minimum of lower boundary region");
		addr[nt] = &Param.Plow;
		id[nt++] = REAL;

		strcpy(tag[nt], "Phigh");
		strcpy(comment[nt], "maximum of lower boundary region");
		addr[nt] = &Param.Phigh;
		id[nt++] = REAL;

#ifdef Q_SHOCK_PRIMARIES
        strcpy(tag[nt], "Xi_ep");
		strcpy(comment[nt], "Shock injection  e-p ratio");
		addr[nt] = &Param.Xi_ep;
		id[nt++] = REAL;
#endif

        strcpy(tag[nt], "Dpp_spec_idx");
		strcpy(comment[nt], "Spectral index for power law Dpp");
		addr[nt] = &Param.Dpp_spec_idx;
		id[nt++] = INT;

		strcpy(tag[nt], "X_crp");
		strcpy(comment[nt], "Norm. fraction rel. to thermal density");
		addr[nt] = &Param.X_crp;
		id[nt++] = REAL;

        strcpy(tag[nt], "a_crp");
		strcpy(comment[nt], "spectral index in secondary injection");
		addr[nt] = &Param.a_crp;
		id[nt++] = REAL;

        strcpy(tag[nt], "p0_crp");
		strcpy(comment[nt], "spectral cutoff in secondary injection");
		addr[nt] = &Param.p0_crp;
		id[nt++] = REAL;

        strcpy(tag[nt], "UnitLength_in_cm");
		strcpy(comment[nt], "[cm] Unit Length");
		addr[nt] = &Unit.Length;
		id[nt++] = REAL;

		strcpy(tag[nt], "UnitMass_in_g");
		strcpy(comment[nt], "[g] Unit Mass");
		addr[nt] = &Unit.Mass;
		id[nt++] = REAL;

		strcpy(tag[nt], "UnitVelocity_in_cm_per_s");
		strcpy(comment[nt], "[cm/s] Unit Vel");
		addr[nt] = &Unit.Vel;
		id[nt++] = REAL;

        id[nt] = -1111;		// Last entry 

        if ((fd = fopen(filename, "r"))) {
		    printf("\nReading Parameter file : %s \n", filename);

			while (fgets(buf, MAXLINELENGTH, fd)) {
                if (sscanf(buf, "%s%s%s",buf1,buf2,buf3) < 2)
						continue;
				if (buf1[0] == '%')
					continue;

				for (i = 0, j = -1; i < nt; i++)
					if ((strcmp(buf1, tag[i]) == 0) 
					 && (tagDone[i] != 1) ){
						j = i;
						tagDone[i] = 1;
						break;
					}

                if (j >= 0) {
					switch (id[j]) {
					case REAL:
						*((double *)addr[j]) = atof(buf2);
						break;
					case STRING:
						strcpy((char *)addr[j], buf2);
						break;
					case INT:
						*((int *)addr[j]) = atoi(buf2);
						break;
					}
				}
			}
		fclose(fd);
		} else 
			Assert(0, "\nParameter file %s not found.\n",filename);

		for (i = 0; i < nt; i++) 
            Assert(tagDone[i],
					"Value for tag '%s' missing in parameter file '%s'.\n",
					tag[i], filename);
    }
    
	MPI_Bcast(&Param, sizeof(Param), MPI_BYTE, 0, MPI_COMM_WORLD);
    
	MPI_Bcast(&Unit, sizeof(Unit), MPI_BYTE, 0, MPI_COMM_WORLD);

    /* consistency checks on input values */
    Assert(Param.Pmin <= Param.Plow, "Pmin has to smaller than Plow");
    Assert(Param.Phigh <= Param.Pmax, "Phigh has to smaller than Pmax");
    Assert(Param.N_IOTasks <= ThisTask.NTask, "NIOTasks not <= NTask");

    return ;
}

#undef REAL
#undef STRING
#undef INT

#ifdef I_NONE   // Dummy Input for Code tests 

#define OUTPUTFREQUENCY 0.1

static int call_count = 0;

double read_input()
{   
    long long npart[N_PART_TYPES]= { 0 };

    npart[0] = 1;
    Reallocate_P(1,npart,1);

    Snap.SnapNum = call_count;
    Snap.Time = OUTPUTFREQUENCY * call_count;
    Snap.Redshift = 1/Snap.Time ;
    Snap.Boxsize = 0;
    Snap.Npart[0] = 1;
    Snap.PartTotal = 1;
    Snap.Masstab[0] = 0;

    call_count++;

    rprintf("Input at Time <%g>, Max <%g> \n", Snap.Time, Param.EndTime);

	int ipart = 0; 

	P[ipart].Rho = 1e-3 /number_density_cgs(1) ; // nth = 1e-3
	SphP[ipart].U = 4.75e9; // T= 1e8
	SphP[ipart].VRms = 400; // km / s
	P[ipart].Hsml = 150; // kpc
	SphP[ipart].Bfld[0] = 1e-6;
	SphP[ipart].Bfld[1] = 0;
	SphP[ipart].Bfld[2] = 0;

	printf("DUMMY INPUT: \n"
			" rho 	= %g \n"
			" U 	= %g \n"
			" Vrms 	= %g \n"
			" Hsml 	= %g \n"
			" B		= %g \n",
			P[ipart].Rho, SphP[ipart].U, SphP[ipart].VRms, P[ipart].Hsml, 
			length3(SphP[ipart].Bfld));
			
    return Snap.Time;
}
    

#endif

/*Determine number of files to read
 * */
int Find_files(char *fname)
{
	char buf[MAXLINELENGTH];
	int n_files = 0;
	FILE *fd = NULL;

    if (!(fd = fopen(fname, "r"))) {

		for (;;) {
		
            sprintf(buf, "%s.%i", fname, n_files);
			
            if (!(fd = fopen(buf, "r"))) 
				break;
			
			fclose(fd);
			
            n_files++;

            Assert(n_files<1000, "Found more than 1000 files");
		}

		if (n_files == 0) 
            Assert(0,"Can't open file <%s> or <%s> !", fname, buf);
	    else 
			rprintf(" \nFound <%i> file(s) ! \n\n",	n_files);
		
	} else 
		n_files = 1;

	return  n_files;
}
