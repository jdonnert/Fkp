#include "common.h"
#include "modules/modules.h"
#include "io/io.h"

void Init_spectrum();
void init_modules();

/* 
 * Set up grid, module pointers etc
 */

void setup()
{
    int i;
    char name[MAXLINELENGTH];

    set_cosmology(Param.Cosmology);

    set_units();

    const int ngrid = N_SPEC_BINS;

    Param.Pmin *= m_e*c;    // convert to cgs
    Param.Pmax *= m_e*c;

    Param.Plow *= m_e*c;
    Param.Phigh *= m_e*c;

    const double pstep = log10(Param.Pmax/Param.Pmin)/(ngrid-1);

    for ( i=0; i<ngrid; i++ )
        p[i] = Param.Pmin * pow(10,i*pstep); 

    for ( i=0; i<ngrid-1; i++ )
        q[i] = 0.5*(p[i+1] + p[i]);

    for ( i=1; i<ngrid-1; i++ ) {
        dp[i] = 0.5 * (p[i+1] - p[i-1]);
        dq[i] = p[i+1] - p[i];
    }

    double pout = Param.Pmin*pow(10,-pstep);
    dp[0] = 0.5 * (p[1] - pout);
    dq[0] = p[1] - p[0];
    
    i = ngrid-1;
    pout = Param.Pmin*pow(10,(i+1)*pstep); 
    q[i] = 0.5*(pout + p[i]);
    dp[i] = 0.5 * (pout - p[i-1]);
    dq[i] = pout - p[i];
    
    /* find open boundary indices */
    for ( i=0; i<ngrid-1; i++ ) { 

        if ( p[i] <= Param.Plow && p[i+1] > Param.Plow)
           LowIdx = i;

        if ( p[i] < Param.Phigh && p[i+1] >= Param.Phigh)
           HighIdx = i+1;
    }
    
    rprintf("\nGrid spans <%d> gridcells \n\n"
        "Boundary Conditions: \n"
        "  pmin = %1.1e g*cm/s \n"
        "       = %1.1e mec    \n"
        "  pmax = %1.1e g*cm/s \n"
        "       = %1.1e mec    \n"
        "  plow = %1.1e g*cm/s \n"
        "       = %7zu cells   \n"
        "  phigh= %1.1e g*cm/s    \n"
        "       = %7zu cells \n\n",
        ngrid, 
        Param.Pmin, Param.Pmin/m_e/c,
        Param.Pmax, Param.Pmax/m_e/c,
        Param.Plow, LowIdx,
        Param.Phigh, ngrid-HighIdx);

#ifdef COMPRESSION
    Setup_compression();
#endif // COMPRESSION

    init_modules();

    return;
}

void Init_spectrum()
{
    const size_t nGas = ThisTask.Npart[0];

    switch (Param.Flag_Start) {

        case 0:

            rprintf("Initialising Spectra ...");

#pragma omp parallel for
            for (size_t ipart = 0; ipart < nGas; ipart++) {
#ifdef COMPRESSION_INTERNAL
                double np[N_SPEC_BINS] = { 0 };

                for (int i = 0; i < N_SPEC_BINS; i++) 
                    np[i] = pow(10, Initial_Spectrum(ipart,p[i],0,0));

                Compress(np, &SphP[ipart].Ncre_total, SphP[ipart].CReSpectrum);
#else
                for (int i=0; i<N_SPEC_BINS; i++) 
                    SphP[ipart].Ncre[i] = Initial_Spectrum(ipart,p[i],0,0);
#endif
            }
            
        break ;
    
        case 1:

            read_spectra();

        break;

        case 2:

#ifdef COMPRESSION
            rprintf("\nConverting read spectra to COMPRESSED\n");
#else
            rprintf("\nConverting read spectra to UNCOMPRESSED\n");
#endif
            read_spectra();
            
            sprintf(Param.Output_File, "%s_new", Param.Output_File);
		    
            write_spectra();

            rprintf("Conversion Complete. Bye ! \n\n");

            Free(P); 
            Free(SphP);

            MPI_Finalize();

            exit(EXIT_SUCCESS);

        break;

        default:
            Assert(0, "Start Flag %d not handled", Param.Flag_Start);
    }

    rprintf("done\n\n");

    return;
}

void init_modules()
{

#ifdef Q_SIMPLE_SECONDARIES
    Init_Simple_Secondaries();
#endif

#ifdef Q_BRUNETTI_05
    Init_Brunetti05();
#endif

#ifdef Q_SHOCK_WITH_HADRONIC_BACKGROUND
	Init_Hadronic_Background();
#endif

    return;
}
