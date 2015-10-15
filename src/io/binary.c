#include "../common.h"

#ifdef O_BINARY

#include "../compress.h"
#include "io.h"

#define MAXFILESIZE (512*1024*1024) // 512 MB

#ifdef COMPRESSION
const size_t SizeSpectrumBytes = (SPECSIZE_BYTES + sizeof(float));
#else
const size_t SizeSpectrumBytes = (N_SPEC_BINS * sizeof(float));
#endif // COMPRESSION

static struct CR_Spectrum_Header {
    	uint64_t SnapNum;
    	uint64_t Nbins;
    	uint64_t Nall;
    	uint64_t Npart;
    	uint64_t Nfiles;
    	uint64_t StartID;
    	double Pmin;
    	double Pmax;
    	double Plow;
    	double Phigh;
		long FlagCompressed;
		long SpecSizeBytes;
		char Fill[84];
} Head;

static void write_file_parallel(int,int,int);
static void read_file_parallel(int,int,int,int);
static void write_header(FILE *,char*, uint64_t);
static void read_header(FILE *,char*);
static size_t get_num_files();

void write_spectra()
{
    const int nTask = ThisTask.NTask;
	const int nGroups = get_num_files();

    MPI_Barrier(MPI_COMM_WORLD);

    int groupSize = nTask/nGroups; // construct write groups 
    
    if (nTask % nGroups)
        groupSize++;

    const int groupMaster = (ThisTask.Rank / groupSize) * groupSize;
    const int groupLast = fmin(groupMaster + groupSize - 1, nTask-1);
	
	size_t nBytes = Snap.Npart[0] * SizeSpectrumBytes;

#ifndef COMPRESSION
    
	rprintf("\nParallel Write of %zu MB in %ld files smaller than %d MB :\n", 
			nBytes/1024/1024, nGroups, MAXFILESIZE/1024/1024);

#else
    
	rprintf("\nParallel Compressed Write of %zu MB in %d files smaller "
			"than %d MB :\n", nBytes/1024/1024, nGroups, MAXFILESIZE/1024/1024);

#endif // COMPRESSION

    write_file_parallel(groupSize, groupMaster, groupLast);

    MPI_Barrier(MPI_COMM_WORLD);

    rprintf("done \n\n"); fflush(stdout);

    return;
}

static size_t get_num_files()
{
	double bunchSize = SizeSpectrumBytes * Snap.Npart[0];

	int nGroups = ceil( fmax(bunchSize/MAXFILESIZE, Param.N_IOTasks) );

	if (nGroups != Param.N_IOTasks) {

		rprintf("Increasing Number of files to %d \n", nGroups);

		Param.N_IOTasks = nGroups;
	}

	return nGroups;
}

static void write_file_parallel(int groupSize, int groupMaster, int groupLast)
{
    int i;
    MPI_Comm groupComm;

    /* setup group communication */
    int thisGroup = groupMaster/groupSize;
	const int groupRank = ThisTask.Rank-groupMaster;
    
    MPI_Comm_split(MPI_COMM_WORLD, thisGroup, ThisTask.Rank, &groupComm);
    
    /* construct filename */
    char fname[MAXLINELENGTH];

    if (Param.N_IOTasks < 2 || ThisTask.NTask == 1)
        sprintf(fname,"%s_%03i",Param.Output_File, Snap.SnapNum);
    else 
        sprintf(fname,"%s_%03i.%d", Param.Output_File, Snap.SnapNum, 
				thisGroup);

    /* find number of particles in file */
	int groupSizeCorr = groupLast - groupMaster + 1; // correct for end

    long long groupNpart[groupSizeCorr];

    MPI_Allgather(ThisTask.Npart, 1, MPI_LONG_LONG, groupNpart, 1
            ,MPI_LONG_LONG, groupComm);

    long long nFile = 0;

    for (i = 0; i<groupSizeCorr; i++)
        nFile += groupNpart[i];

    /* start writing */
    char *recvbuf = NULL;
    FILE *fp = NULL;

    if (ThisTask.Rank == groupMaster) {

        printf("%zu spectra in %s from %d-%d \n"
                ,nFile, fname, groupMaster, groupLast);

        fp = fopen(fname, "w");

        Assert(fp != NULL, "Can't open file for writing.\n");
        
        write_header(fp, fname, nFile);

    	size_t nBytes = nFile * SizeSpectrumBytes;

		recvbuf = Malloc(nBytes);
    }

    /* treat sendbuf */
    size_t nBytes = ThisTask.Npart[0] * SizeSpectrumBytes;
    
	char *sendbuf = Malloc(nBytes);
    
	//#pragma omp parallel for
	for (int ipart = 0; ipart < ThisTask.Npart[0]; ipart++) { // fill sendbuf

#ifdef COMPRESSION

		char creSpectrum[SPECSIZE_BYTES] = { 0 };
		float nCRe = 0;

		Ipart = P[ipart].ID; // for DEBUG

if (Ipart == 119455+1) {
	for (int i = 0; i < 128; i++)
printf("%d %g %d %g %g %g \n", ipart, nCRe, i, p[i], SphP[ipart].Ncre[i] ,SphP[ipart].Mach );

		Compress(ipart, &SphP[ipart].Ncre[0], &nCRe, creSpectrum);

	//nCRe = 14041981;
}

		size_t idx = ipart * (sizeof(nCRe) + SPECSIZE_BYTES);

		memcpy(&sendbuf[idx], &nCRe, sizeof(nCRe));

		idx += sizeof(nCRe);
		
		memcpy(&sendbuf[idx], creSpectrum, SPECSIZE_BYTES);

#else // !COMPRESSION

		size_t idx = ipart * N_SPEC_BINS*sizeof(float);
		
		memcpy(&sendbuf[idx], SphP[ipart].Ncre, N_SPEC_BINS*sizeof(float));

#endif // !COMPRESSION
	}

    /* gatherv book keeping */
    int recvcnts[groupSize];

    for (i = 0; i<groupSize; i++) 
        recvcnts[i] = groupNpart[i] * SizeSpectrumBytes;

    int offsets[groupSize];
    offsets[0] = 0;

    for (i = 1; i<groupSize; i++)
        offsets[i] = offsets[i-1] + recvcnts[i-1];

    /* communicate spectra */
    MPI_Gatherv(sendbuf, recvcnts[groupRank], MPI_BYTE, 
				recvbuf, recvcnts, offsets, MPI_BYTE, 0, groupComm) ;

    /* Write on groupMaster */
    nBytes = nFile * SizeSpectrumBytes;
    
	if (ThisTask.Rank == groupMaster)  { 
    
		fwrite(recvbuf, nBytes , 1, fp);
        
		fclose(fp);
    }

	if (recvbuf != NULL)
    	Free(recvbuf); 

	Free(sendbuf);

	MPI_Comm_free(&groupComm);

    return;
}

void read_spectra()
{
#ifdef COMPRESSION_INTERNAL
    Assert(0, "Reading Spectra into internalt compression not implemented");
#endif

    const int nTask = ThisTask.NTask;

    MPI_Barrier(MPI_COMM_WORLD);

    int nFiles = 0;

    if (!ThisTask.Rank) { // find nFiles

        char fname[MAXLINELENGTH];
        
        sprintf(fname,"%s_%03i", Param.Output_File, Snap.SnapNum);

        nFiles = Find_files(fname);

        Assert(nTask >=  nFiles, "More files (%d) than tasks", nFiles );
    }

	MPI_Bcast(&nFiles, 1, MPI_INT, 0, MPI_COMM_WORLD);  

    rprintf("Reading CR spectra from %d files : \n\n", nFiles);
    
    const int nGroups = nFiles;

    int groupSize = nTask/nGroups; // construct read groups 
    
    if (nTask % nGroups)
        groupSize++;

    const int groupMaster = (ThisTask.Rank / groupSize) * groupSize;
    const int groupLast = fmin(groupMaster + groupSize - 1, nTask-1);
 
    read_file_parallel(nFiles, groupSize, groupMaster, groupLast);

    return ;
}

static void read_file_parallel(int nFiles, int groupSize, int groupMaster, 
		int groupLast)
{
	const int groupRank = ThisTask.Rank-groupMaster;
    const int nPart = ThisTask.Npart[0];
    
    int thisGroup = groupMaster/groupSize; // setup group communication
    
    MPI_Comm groupComm;

    MPI_Comm_split(MPI_COMM_WORLD, thisGroup, ThisTask.Rank, &groupComm);

    long long nPartGroup = 0;

    MPI_Reduce(&ThisTask.Npart[0], &nPartGroup, 1, MPI_LONG_LONG, MPI_SUM,
           0, groupComm);
    
    void *readBuf = NULL;  
    
    if (ThisTask.Rank == groupMaster) {
        
        char fname[MAXLINELENGTH]; // construct filename 

        if (nFiles == 1)
            sprintf(fname,"%s_%03i",Param.Output_File, Snap.SnapNum);
        else 
            sprintf(fname,"%s_%03i.%d", Param.Output_File, Snap.SnapNum, 
                thisGroup);

        FILE *fp = fopen(fname, "r"); 

        Assert(fp != NULL, "Can't open file for reading. %s\n", fname);
    
        read_header(fp, fname);

        Assert(Head.Npart == nPartGroup, "Inconsistent particle numbers"
        "Ngroup=%lld, Nfile=%lld ", nPartGroup, Head.Npart);

        printf("Reading %lld spectra from %s to tasks %d to %d\n", nPartGroup, 
                fname, groupMaster, groupLast);

        size_t nBytes = Head.SpecSizeBytes * Head.Npart; 

        readBuf = Malloc(nBytes);   
         
        fread(readBuf, nBytes, 1, fp);

        fclose(fp);
    }

    MPI_Bcast(&Head, sizeof(Head), MPI_BYTE, 0, groupComm);  
                                                         
    Assert(Head.StartID <= P[0].ID, // Check ID range
            "Not all needed IDs found in file StartFileID=%d  FirstTaskID=%d",
            Head.StartID, P[0].ID);

    Assert(Head.StartID+Head.Npart-1 >= P[nPart-1].ID, 
           "Not all needed IDs found in file, LastFileID=%d LastTaskID=%d",
                Head.StartID+Head.Npart-1, P[nPart-1].ID);

    const int nBytes = Head.SpecSizeBytes * nPart; 

    int nBytesGroup[groupSize]; // find transfer sizes

    MPI_Allgather(&nBytes, 1, MPI_INT, nBytesGroup, 1, MPI_INT, groupComm);

    int offsets[groupSize]; // find displacements
    
    offsets[0] = 0;

    for (int i = 1; i < groupSize; i++) 
        offsets[i] = offsets[i-1] + nBytesGroup[i-1]; 

    char *commBuf = Malloc(nBytes); // allocate comm buf

    MPI_Scatterv(readBuf, nBytesGroup, offsets, MPI_BYTE, commBuf, nBytes, 
            MPI_BYTE, 0, groupComm);

    if (ThisTask.Rank == groupMaster)
        Free(readBuf);


    if (Head.FlagCompressed == 0) { // uncompressed file

        float *src = (float*) commBuf;

		#pragma omp parallel for private(src)
        for (int ipart = 0; ipart < nPart; ipart++) 
            for (int j = 0; j < N_SPEC_BINS; j++) 
                SphP[ipart].Ncre[j] = src[ipart*N_SPEC_BINS + j];

    } else { // commBuf is compressed

        Assert(SizeSpectrumBytes == Head.SpecSizeBytes, 
                "Compressed spec size differ: code=%zu, file=%zu",
                SizeSpectrumBytes, Head.SpecSizeBytes);

        for (int ipart = 0; ipart < nPart; ipart++) {
                
        	char creSpectrum[SPECSIZE_BYTES] = {0};
	        float nCRe = 0;

			size_t idx = ipart * (SizeSpectrumBytes);

            memcpy(&nCRe, &commBuf[idx], sizeof(float));

            idx += sizeof(float);

            memcpy(creSpectrum, &commBuf[idx], SPECSIZE_BYTES);

    		double np[N_SPEC_BINS] = { -FLT_MAX };

			Uncompress(nCRe, creSpectrum, np);

			for (int i = 0; i < N_SPEC_BINS; i++)
				SphP[ipart].Ncre[i] = log10(np[i]);
        }
        
    }

    Free(commBuf);
    
	MPI_Comm_free(&groupComm);

    MPI_Barrier(MPI_COMM_WORLD);
    
    return ;
}

static void write_header(FILE *fp, char *fname, uint64_t npart_file)
{
    Head.SnapNum = Snap.SnapNum;
    Head.Nbins = N_SPEC_BINS;
    Head.Nall = Snap.Npart[0];
    Head.Npart = npart_file;
    Head.Nfiles = Param.N_IOTasks;
    Head.StartID = P[0].ID;
    Head.Pmin = Param.Pmin/(m_e*c);
    Head.Pmax = Param.Pmax/(m_e*c);
    Head.Plow = Param.Plow/(m_e*c);
    Head.Phigh = Param.Phigh/(m_e*c);
#ifdef COMPRESSION
    Head.FlagCompressed = 1;
    Head.SpecSizeBytes = SPECSIZE_BYTES + sizeof(float);
#else 
    Head.FlagCompressed = 0;
    Head.SpecSizeBytes = N_SPEC_BINS * sizeof(float);
#endif

    size_t nItem = fwrite(&Head, sizeof(Head), 1, fp);

    Assert(nItem == 1, "Problem Writing file header ! \n");

    return;
}

static void read_header(FILE *fp, char *fname)
{
    size_t nItem = fread(&Head, sizeof(Head), 1, fp);

    Assert(nItem == 1, "Problem Reading file header ! \n");
    
    Assert(Head.Nbins == N_SPEC_BINS, 
            "N_SPEC_BINS inconsistent %d : %d", N_SPEC_BINS, Head.Nbins);
    Assert(Snap.SnapNum == Head.SnapNum, 
            "Snapshot numbers inconsistent %d : %d", 
            Snap.SnapNum, Head.SnapNum);
    Assert(Head.Nall == Snap.Npart[0], "Npart[0] inconsistent %g : %g",
            Snap.Npart[0],Head.Nall);
    Assert(Head.Pmin == Param.Pmin/(m_e*c),"Pmin inconsistent %g : %g",
       Param.Pmin/(m_e*c),Head.Pmin);
    Assert(Head.Pmax == Param.Pmax/(m_e*c),"Pmax inconsistent %g : %g",
        Param.Pmax/(m_e*c),Head.Pmax);
    Assert(Head.Plow == Param.Plow/(m_e*c),"Plow inconsistent %g : %g",
        Param.Plow/(m_e*c),Head.Plow);
    Assert(Head.Phigh == Param.Phigh/(m_e*c),"Phigh inconsistent %g : %g",
        Param.Phigh/(m_e*c),Head.Phigh);

#ifndef COMPRESSION
    Assert(Head.FlagCompressed == 0, 
            "Can't read compressed data, recompile");
#endif

    return ;
}
#endif
