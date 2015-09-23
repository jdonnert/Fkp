#include <gsl/gsl_sort.h>

#include "common.h"

void sortP_local();

/* Here we do a global sort for ID. As we know how many IDs (nID) every
 * task carries, we do an All to All communication, send what belongs to 
 * the other task, append what we recieve to the end of the particle
 * structure and move the end of the particle array to close the sent gap. 
 */

void sortP_global()
{
    const int nTask = ThisTask.NTask;
    const int myRank = ThisTask.Rank;
    const int lastTask = ThisTask.NTask-1;

    rprintf("Starting global sort ..."); fflush(stdout);

    /* remove non SPH particles from memory */
    long long nNonSph = 0, nPart[N_PART_TYPES] = {0};

	for (int i = 1; i < N_PART_TYPES; i++) {

        nNonSph += ThisTask.Npart[i];
        
		nPart[i]+= ThisTask.Npart[i];
    }

    if (nNonSph > 0)
        Reallocate_P(nNonSph, nPart, -1);

	Assert(ThisTask.PartTotal == ThisTask.Npart[0],
			"Problem in particle removal");

    /* sort P+SphP by ID */
	sortP_local();

    /* exchange particles */
    const size_t nID = floor(Snap.Npart[0] / nTask); // not true for last
    const size_t bitmax = pow(2, ceil(log2(nTask))); // hypercube

    for (int bitmask = 1; bitmask < bitmax; bitmask++) { 
    
        MPI_Barrier(MPI_COMM_WORLD);

		int sendTask = myRank;
        int recvTask = myRank ^ bitmask; // talk to everyone
        
        int tag = sendTask & recvTask;

        if (recvTask >= nTask || recvTask == myRank)
            continue;

        /* IDs req. in recvTask */
        size_t recvIDrange[2] = { recvTask*nID+1, (recvTask+1)*nID }; 

        if ( recvTask == lastTask ) // treat last CPU if nPart%nTask != 0 
            recvIDrange[1] = Snap.Npart[0];

        /* find IDs to be send  */
        size_t nSend = 0, iSend = 0;

        for (size_t ipart = 0; ipart < ThisTask.Npart[0]; ipart++) {

            if (P[ipart].ID < recvIDrange[0]) 
				continue;

            if (P[ipart].ID > recvIDrange[1])
				break;
            
            if (nSend == 0)
                iSend = ipart;

            nSend++;
        }

        /* communicate no of IDs to recieve */
    	MPI_Status status;

        size_t nRecv = 0;

        MPI_Sendrecv(&nSend, sizeof(nSend), MPI_BYTE, recvTask, tag,
                     &nRecv, sizeof(nRecv), MPI_BYTE, recvTask, tag,
                     MPI_COMM_WORLD, &status);

        /* add space at end of P and SphP */
        memset(nPart, 0, N_PART_TYPES*sizeof(*nPart));
		nPart[0] = nRecv;

        Reallocate_P(nRecv, nPart, +1);

        /* communicate particles - two way*/
        size_t iRecv = ThisTask.Npart[0] - nRecv;

        MPI_Sendrecv(&(P[iSend]), nSend*sizeof(*P), MPI_BYTE, recvTask, tag,
                     &(P[iRecv]), nRecv*sizeof(*P), MPI_BYTE, recvTask, tag,
                     MPI_COMM_WORLD, &status);

        MPI_Sendrecv(&(SphP[iSend]),nSend*sizeof(*SphP),MPI_BYTE,recvTask,tag, 
                     &(SphP[iRecv]),nRecv*sizeof(*SphP),MPI_BYTE,recvTask,tag,
                     MPI_COMM_WORLD, &status);

        /* move particles to close send gap */
        size_t iLast = iSend + nSend;
        size_t nBytes = (ThisTask.Npart[0] - iLast) * sizeof(*P);

        memmove(&(P[iSend]), &(P[iLast]), nBytes);
        
        nBytes = (ThisTask.Npart[0] - iLast)*sizeof(*SphP);

        memmove(&(SphP[iSend]), &(SphP[iLast]),nBytes );
        
		/* remove space at the end of P & SphP */
        memset(nPart, 0, N_PART_TYPES*sizeof(*nPart));

        nPart[0] = nSend;
        
        Reallocate_P(nSend,nPart,-1);
    }
  	
	sortP_local(); 

    if (LastP == NULL)
		goto skip_test;
    
	for (int ipart = 0; ipart < ThisTask.Npart[0]; ipart++ ) {

            ptrdiff_t delta = P[ipart].ID - LastP[ipart].ID;
            
			Assert(delta == 0, "Global Sort not successful %d %d %d", 
					ipart,  P[ipart].ID, LastP[ipart].ID);
    }

	skip_test:;
        
    rprintf("done\n");

    MPI_Barrier(MPI_COMM_WORLD);

	#pragma omp parallel for
    for (int ipart = 0; ipart < ThisTask.Npart[0] -1; ipart++ ) 
		Assert(P[ipart+1].ID - P[ipart].ID == 1, 
				"Local Sort not successful, ipart=%d ID=%d NextID=%d npart=%d ",
				ipart, P[ipart].ID, P[ipart+1].ID, ThisTask.Npart[0] );
	
    return;
}

int compare_ids(const void *a, const void *b)
{
	const int *x = (const int*) a;
	const int *y = (const int*) b;
	
	return (*x > *y) - (*x < *y);
}

/* 
 * memory efficient out of place sorting of both particle structures 
 * Knowing where idx starts in memory we can reconstruct ipart 
 */

void sortP_local() 
{
    const int nPart = ThisTask.Npart[0];

    size_t *idx = Malloc(nPart*sizeof(*idx));
    int *ids = Malloc(nPart*sizeof(*ids));

    for (size_t ipart = 0; ipart < nPart; ipart++) 
        ids[ipart] = P[ipart].ID;

	#pragma omp parallel
	Qsort_Index(1, idx, ids, nPart, sizeof(*ids), compare_ids);
    
	for (int i = 0; i < nPart; i++) {

        if (idx[i] == i)
            continue;

		size_t dest = i;

        struct Particle_Data Ptmp = P[dest];
        struct Gas_Data SphPtmp = SphP[dest];

		size_t src = idx[i];

        for (;;) {

			P[dest] = P[src];
            SphP[dest] = SphP[src];
            idx[dest] = dest;

			dest = src;

			src = idx[dest];

            if (src == i) 
                break;
        }

		P[dest] = Ptmp;
        SphP[dest] = SphPtmp;

        idx[dest] = dest;
    }

    Free(idx); Free(ids);

	//#pragma omp parallel for
    //for (int ipart = 0; ipart < ThisTask.Npart[0]-1; ipart++ ) 
//		Assert(P[ipart+1].ID - P[ipart].ID == 1, 
//				"Local Sort not successful, ipart=%d ID=%d NextID=%d ",
//				ipart, P[ipart].ID, P[ipart+1].ID );

    return;
}

