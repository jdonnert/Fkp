#include "../common.h"

#ifdef I_TXT // Read & write ASCII Data

static double *z = NULL, *a = NULL
    ,*Dpp = NULL,  *nth = NULL, *bfld = NULL;

static size_t nLines = 0;

double read_input()
{
    size_t i;
    long long npart[N_PART_TYPES]= {0};
    char buf[MAXLINELENGTH];
    FILE *fp;
    char *stat;

    Assert(ThisTask.NTask == 1, "ASCII Input only with NTask = 1");

    if ( z == NULL ){ /* Read File */
        
        rprintf("ASCII Input from <%s>\n",Param.Input_File);

        npart[0] = 1;
        Reallocate_P(1,npart,1);
    
        fp = fopen(Param.Input_File,"r");

        while ( fgets(buf,MAXLINELENGTH,fp) )
            nLines++;
    
        rprintf("Read <%zu> lines\n\n",nLines);

        z = my_malloc(nLines * sizeof(*z));
        a = my_malloc(nLines * sizeof(*a));
        Dpp = my_malloc(nLines * sizeof(*Dpp));
        nth = my_malloc(nLines * sizeof(*nth));
        bfld = my_malloc(nLines * sizeof(*bfld));
    
        rewind(fp);

        for ( i=0; i<nLines; i++ ){
            stat = fgets(buf,MAXLINELENGTH,fp);
        
            sscanf(buf,"%12le %12le %12le %12le\n"
                    ,z+i,Dpp+i, nth+i, bfld+i);

            a[i] = 1/(1+z[i]);
        }
        
        fclose(fp);

        Snap.SnapNum = -1;
    }
    
    if ( ++Snap.SnapNum == nLines){
        free(z); free(a); 
        free(Dpp); free(nth); free(bfld);

        printf("\nNo Snapshots left.\n\n");

        return(Param.EndTime - 1);
    }
   
    i = Snap.SnapNum; 
    Snap.Redshift = z[i];
    Snap.Time = a2t(a[i]);
    Snap.Boxsize = 0;
    Snap.Npart[0] = 1;
    Snap.PartTotal = 1;
    Snap.Masstab[0] = 0;

    SphP[0].Dpp = Dpp[i];
    P[0].Rho = nth[i];
    SphP[0].Bfld[0] = bfld[i] /sqrt3;
    SphP[0].Bfld[1] = bfld[i] /sqrt3;
    SphP[0].Bfld[2] = bfld[i] /sqrt3;
    
    if (!ThisTask.Rank) 
        printf("New Snap at z = %1.3f = %1.2e s\n \n"
                ,Snap.Redshift, Snap.Time);

    return(Snap.Redshift);
}
#endif

#ifdef O_TXT
void write_output()
{
    char fname[MAXLINELENGTH];
    size_t ipart = 0;
    FILE *fp = NULL;

	const long ngrid = N_SPEC_BINS;

    if (!ThisTask.Rank) {

        sprintf(fname, "%s_%03i_%ld.txt"
                ,Param.Output_File,Snap.SnapNum,ngrid);

        printf("\nWriting ASCII output to <%s>\n\n",fname);

        fp = fopen(fname, "w");

        Assert(fp != NULL, "Can't open output file");

        ipart = 0; // output only one particle 
        
        fprintf(fp, "# %g %g %g %g\n",SphP[ipart].Dpp, P[ipart].Rho
                    ,length3(SphP[ipart].Bfld), SphP[ipart].VRms);

        double *np = my_malloc(N_SPEC_BINS * sizeof(*np));

#ifdef COMPRESSION_INTERNAL
        Uncompress(SphP[ipart].Ncre_total, SphP[ipart].CReSpectrum,  np);
#else
        for (int i = 0; i < N_SPEC_BINS; i++)
            np[i] = pow(10,SphP[ipart].Ncre[i]);
#endif

        for (int i = 0; i < N_SPEC_BINS; i++)
            fprintf(fp, "%g %g\n",p[i], log10(np[i]));

        fclose(fp); 
    }
    MPI_Barrier(MPI_COMM_WORLD);

    return;
}

#endif
