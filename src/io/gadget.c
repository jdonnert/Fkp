#include "../common.h"
#include "io.h"
#include "gadget.h"
#include "../modules/modules.h"
#include <limits.h>

#define int4bytes int

#define SKIP  {if(Fread(&blksize,sizeof(int),1,fd)){swap_Nbyte((char*)&blksize,1,4);}}

#define UINT sizeof(unsigned int)
#define INT sizeof(long)
#define FLOAT sizeof(float)
#define DBL sizeof(double)

#ifdef I_GADGET2

int4bytes blksize, swap = 0;
enum iofields blocknr;

static long n_files = 0;
static char fbase[MAXLINELENGTH];

static void swap_Nbyte(char *, int, int);
static size_t Fread(void *, size_t, size_t, FILE *);
static int find_block(FILE *, char *);
static void read_file(char *, int, int);
static int read_gadget_block(void *, char *, FILE *, size_t);
static int read_gadget_head(FILE *);
static void set_block_prop(enum iofields);
static void empty_comm_buffer(enum iofields, void *, size_t, size_t);

/* Main Reading :
 * This routine tries first to put one snapshot per task
 * Then nGroups Masters distribute to their groups
 * All files left are then read at the same time.
 * */
double read_input()
{
	long rest_files, file_number, nGroups
        , groupMaster, groupLast, groupSize, group;
	char fname[MAXLINELENGTH], buf[MAXLINELENGTH];

	Snap.Have_Arepo = 0;

    nGroups = Param.N_IOTasks;

    /* Construct Filename */
    if ( !n_files ){    /* First Time */

        Snap.SnapNum = 
            guess_snapnum(Param.Input_File, fbase);
        Snap.LastSnapNum = -1;
        sprintf(fname, "%s_%03i",fbase, Snap.SnapNum);
    } else {            /* next time */

        Snap.LastSnapNum = Snap.SnapNum;
        Snap.SnapNum++;
        sprintf(fname,"%s_%03i",fbase, Snap.SnapNum);
    }

    n_files = rest_files = Find_files(fname);
    
    if (!n_files) /* simulation done */
        return DBL_MAX;

    Assert(n_files, "Input file not found");

    /* Start parallel reading */
	if (nGroups > n_files) {
		nGroups = n_files;
	    rprintf("\nReducing Number of Read Tasks to Number of Files\n");
	}

	groupSize = ThisTask.NTask / nGroups;
	
    if (ThisTask.NTask % nGroups) 
		groupSize++;

	groupMaster = (ThisTask.Rank / groupSize) * groupSize;
	groupLast = groupMaster + groupSize - 1;
	
    if (groupLast > ThisTask.NTask - 1) 
		groupLast = ThisTask.NTask - 1;

    while (rest_files) { 
        /* read 1 to all CPUs*/
		if (n_files == 1) {	
			sprintf(buf, "%s", fname);
			
            read_file(buf, 0, ThisTask.NTask - 1);
			
            MPI_Barrier(MPI_COMM_WORLD);
			
            rest_files--;

         /* Read in nIO blocks */   
		} else if (rest_files >= ThisTask.NTask) { 
			file_number = ThisTask.Rank 
                + (rest_files - ThisTask.NTask);
			
            sprintf(buf, "%s.%li", fname, file_number);
			
            for (group = 0; group < groupSize; group++) {
				if (ThisTask.Rank == (groupMaster + group)) 
					read_file(buf, ThisTask.Rank,
						  ThisTask.Rank);
				
                MPI_Barrier(MPI_COMM_WORLD);
			}

			rest_files -= ThisTask.NTask;

        /* Read and distribute by master in nIO groups */
		} else if (rest_files >= nGroups) { 
			file_number = groupMaster / groupSize
                + (rest_files - nGroups);
			
            sprintf(buf, "%s.%li", fname, file_number);
			
            read_file(buf, groupMaster, groupLast);
			
            MPI_Barrier(MPI_COMM_WORLD);
			
            rest_files -= nGroups;

        /* reduce nIO to rest files */
		} else { 
			groupSize = ThisTask.NTask / rest_files;
			if (ThisTask.NTask % nGroups) 
				groupSize++;
			groupMaster = (ThisTask.Rank / groupSize)
                * groupSize;
			
            file_number = groupMaster / groupSize;
			
            sprintf(buf, "%s.%li", fname,file_number);
			
            read_file(buf, groupMaster, groupLast);
			
            MPI_Barrier(MPI_COMM_WORLD);
			
            rest_files -= groupSize;
		}
	}

	rprintf("\nTime of Snapshot : z=%g, a=%4g, t=%g Myr (%s) \n\n"
            , Snap.Redshift, Snap.Time, a2t(Snap.Time)/yr2sec/1e6, Cosmo.name);
	
    if ( Snap.Have_Arepo )
        rprintf("Read an AREPO Snapshot ! \n\n");

#ifdef COMPUTE_DPP  
    Calculate_Dpp();
#endif

	return Snap.Time;
}

/* Reads and Distributes a file 
 * */
void read_file(char *fname, int ReadTask, int LastTask)
{
	long i = 0, j = 0, task = 0;
	long target, src, nTask, blockExist;
	long long nRead[N_PART_TYPES] = { 0 }, nReadTot = 0;
	long long nSend[N_PART_TYPES] = { 0 }, ntot = 0;
	void *comm_buf = NULL;
	size_t nBytes, byteOffset, partOffset, bufOffset;
	FILE *fd = NULL;
	int tag = 0;
	MPI_Status status;
	
    if (ThisTask.Rank == ReadTask) {
		fd = fopen(fname, "r");
		ntot = read_gadget_head(fd);

		printf("\nReading file <%s> on Task <%i-%i> \n"
            "   Sph   <%9lli>  DM     <%9lli>    \n"
            "   Disk  <%9lli>  Bulge  <%9lli>    \n"
            "   Star  <%9lli>  Bndry  <%9lli>    \n"
            "   Total <%9lli> \n\n",
			fname, ReadTask, LastTask,
			Header.Npart[0], Header.Npart[1], Header.Npart[2],
			Header.Npart[3], Header.Npart[4], Header.Npart[5], 
            ntot); fflush(stdout);

		for (task = ReadTask + 1; task <= LastTask; task++) {

			MPI_Ssend(&Header, sizeof(Header), MPI_BYTE, task, tag,
				  MPI_COMM_WORLD);
		}
	} else 
		MPI_Recv(&Header, sizeof(Header), MPI_BYTE, ReadTask, tag,
			 MPI_COMM_WORLD, &status);
    
    /* Set Snapshot Properties */
    Snap.Boxsize = Header.BoxSize;
    Snap.Redshift = Header.Redshift;
    Snap.Time = Header.Time;

	Snap.PartTotal = 0;
	
	for (i = 0; i < N_PART_TYPES; i++) {

		Snap.Npart[i] = Header.Nall[i];
		Snap.PartTotal += Header.Nall[i];
        Snap.Masstab[i] = Header.Mass[i] * 1 / Cosmo.h;
	}

    MPI_Bcast(&Snap, sizeof(Snap),
			MPI_BYTE, 0, MPI_COMM_WORLD);

    /* Set Comoving units for Gadget */
	Comv2phys.Length = 1;	
	Comv2phys.Mass = 1;
	Comv2phys.Vel = 1;
    
    /* Determine particle distribution over CPUs */
	nTask = LastTask - ReadTask + 1;
	for (i = 0; i < N_PART_TYPES; i++) {
		for (j = ThisTask.Rank - ReadTask; j < Header.Npart[i];
		     j += nTask) {
			nRead[i]++; /* npart[i] ThisTask reads */
			nReadTot++; /* ntot ThisTask reads */
		}
	}

	Reallocate_P(nReadTot, nRead, +1);

	/*Shift collisionless particles if multiple files */

	src = ThisTask.Npart[0] - nRead[0];	
	target = src + nReadTot;
	nBytes = (ThisTask.PartTotal - nReadTot - src) * sizeof(*P);
	memmove(&P[target], &P[src], nBytes);

	/* Read blocks  */
	for (blocknr = IO_ID; blocknr < IO_LASTENTRY; blocknr++) {	
        
		set_block_prop((enum iofields)blocknr); // Set Block Properties + Read
        
        byteOffset = 0;
    
    	if (ThisTask.Rank == ReadTask) {

            nBytes = Block.Ntot * Block.Bytes_per_element;

	    	comm_buf = Malloc(nBytes);
   
   			blockExist = read_gadget_block(comm_buf, Block.Label, fd,
   					  Block.Data_type);

            if (blockExist)
                printf("    %s  %zu MB\n", Block.Name, nBytes/1024/1024);

            for (task = ReadTask + 1; task < LastTask + 1; task++)
                MPI_Ssend(&blockExist, 1,MPI_LONG, task, tag, MPI_COMM_WORLD);

        } else 
            MPI_Recv(&blockExist, 1,MPI_LONG, ReadTask,tag
                    ,MPI_COMM_WORLD, &status);     
        
        
        /* Communicate Blockdata */
        if (ThisTask.Rank == ReadTask && blockExist){
   			for (i = 0; i < N_PART_TYPES; i++) {
                
                if (!Block.Npart[i])    
                    continue;

   				byteOffset += nRead[i]  /* Keep some on ReadTask */
                    * Block.Bytes_per_element; 

   				for (task = ReadTask + 1; task < LastTask + 1;
   				     task++) {
   				
                    MPI_Recv(nSend, N_PART_TYPES,
   						 MPI_LONG_LONG, task, tag,
   						 MPI_COMM_WORLD, &status);
   					
                    nBytes = nSend[i] * Block.Bytes_per_element;
   					
                    MPI_Ssend(comm_buf + byteOffset, nBytes,
   						  MPI_BYTE, task, tag,
   						  MPI_COMM_WORLD);

                    byteOffset += nBytes;
   				}
   			}
   		} else if (blockExist){
   			comm_buf =
   			    Malloc(Block.Bytes_per_element * nReadTot);

   			for (i = 0; i < N_PART_TYPES; i++) {

                if (!Block.Npart[i])    
                    continue;

   				MPI_Ssend(nRead, N_PART_TYPES, MPI_LONG_LONG,
   					  ReadTask, tag, MPI_COMM_WORLD);
   				
                nBytes = nRead[i] * Block.Bytes_per_element;
   				
                MPI_Recv(comm_buf + byteOffset, nBytes,
   					 MPI_BYTE, ReadTask, tag,
   					 MPI_COMM_WORLD, &status);

   				byteOffset += nBytes;
   			}
   		}
    
        /* Readout commbuffer */
   		partOffset = ThisTask.Npart[0] - nRead[0];
   		bufOffset = 0;
   		for (i = 0; i < N_PART_TYPES; i++) {

   			if (Block.Npart[i] && blockExist) {

   				for (j = 0; j < nRead[i]; j++) {
   				
                    empty_comm_buffer((enum iofields)
   							  blocknr, comm_buf,partOffset + j,
   							  (bufOffset + j) * Block.
   							  Val_per_element);
   				}
   
   				if ( blocknr == IO_ID) /* set type */
   					for (j = 0; j < nRead[i]; j++) 
   						P[partOffset+j].Type = i;

   			} else if( blocknr == IO_MASS ){
                for (j = 0; j < nRead[i]; j++)/* masses from Header */
		            P[partOffset+j].Mass = Snap.Masstab[i];
            }  

            if (Block.Npart[i]){
                if (ThisTask.Rank == ReadTask) {
   					bufOffset += Block.Npart[i];
   				} else {
   					bufOffset += nRead[i];
   				}
   				partOffset += nRead[i];
            }
   		}
        if( blockExist || ThisTask.Rank == ReadTask)
	        free(comm_buf);
	}

	if (ThisTask.Rank == ReadTask)
		fclose(fd);

    /* treat AREPO VOL -> HSML */
    const double fourpi3 = 4*pi/3.0;
    if (Snap.Have_Arepo)
        for (j = 0; j < nRead[0]; j++) 
            SphP[j].Hsml = pow(SphP[j].Hsml/fourpi3,1.0/3.0);

	return;
}

/* Basic routine to read data from a file 
 * */
size_t Fread(void *ptr, size_t size, size_t nmemb, FILE * stream)
{
	size_t nRead;
	if ((nRead = fread(ptr, size, nmemb, stream)) != nmemb) {
        Assert(feof(stream), "I/O error");

        nRead = 0;	// EOF reached 
	}
	return nRead;
}

/* Routine to swap ENDIAN
 * */
void swap_Nbyte(char *data, int n, int m)
{
	int i, j;
	char old_data[16];

	if (swap > 0) {
		for (j = 0; j < n; j++) {
			memcpy(&old_data[0], &data[j * m], m);
			for (i = 0; i < m; i++) {
				data[j * m + i] = old_data[m - i - 1];
			}
		}
	}
	return;
}

/* Find block label
 * */
int find_block(FILE * fd, char *label)
{
	int4bytes blocksize = 0;
	char blocklabel[5] = { "    " };

	rewind(fd);

	while (!feof(fd) && blocksize == 0) {
		SKIP;
		if (blksize == 134217728) { /* 256 in other endianess */
			fprintf(stdout, "Enabling ENDIAN swapping !\n");
			swap = 1 - swap;
			swap_Nbyte((char *)&blksize, 1, 4);
		}

        Assert(blksize==8, "incorrect F77 Header found");

		if (Fread(blocklabel, 4 * sizeof(char), 1, fd)) {
		    Fread(&blocksize, sizeof(int4bytes), 1, fd);
			swap_Nbyte((char *)&blocksize, 1, 4);
			SKIP;
			if (strncmp(label, blocklabel,5) != 0) {
				fseek(fd, blocksize, 1);
				blocksize = 0;
			}
		} else {
			blocksize = 8;
			break;
		}
	}
	return (blocksize - 8);
}

/* Read the header information, returns total nr of particles in file
 * */
int read_gadget_head(FILE * fd)
{
	int blocksize, dummysize, i;
	unsigned int Npart[N_PART_TYPES], Nall[N_PART_TYPES],
	    NallHW[N_PART_TYPES];
	long long ntot = 0;

	blocksize = find_block(fd, "HEAD");

    Assert(blocksize>0, "Header not found");
	
    dummysize = blocksize - 2 * N_PART_TYPES * sizeof(int) 
        - 4 * sizeof(long) - 12 * sizeof(double);
	SKIP;

	Fread(Npart, N_PART_TYPES * sizeof(unsigned int), 1, fd);
	swap_Nbyte((char *)Npart, N_PART_TYPES, 4);

	Fread(Header.Mass, N_PART_TYPES * sizeof(double), 1, fd);
	swap_Nbyte((char *)Header.Mass, N_PART_TYPES, 8);

	Fread((void *)&Header.Time, sizeof(double), 1, fd);
	swap_Nbyte((char *)&Header.Time, 1, 8);

	Fread((void *)&Header.Redshift, sizeof(double), 1, fd);
	swap_Nbyte((char *)&Header.Redshift, 1, 8);

	Fread((void *)&Header.FlagSfr, sizeof(int), 1, fd);
	swap_Nbyte((char *)&Header.FlagSfr, 1, 4);

	Fread((void *)&Header.FlagFeedback, sizeof(int), 1, fd);
	swap_Nbyte((char *)&Header.FlagFeedback, 1, 4);

	Fread(Nall, N_PART_TYPES * sizeof(unsigned int), 1, fd);
	swap_Nbyte((char *)Nall, N_PART_TYPES, 4);

	Fread((void *)&Header.FlagCooling, sizeof(int), 1, fd);
	swap_Nbyte((char *)&Header.FlagCooling, 1, 4);

	Fread((void *)&Header.NumFiles, sizeof(int), 1, fd);
	swap_Nbyte((char *)&Header.NumFiles, 1, 4);

	Fread((void *)&Header.BoxSize, sizeof(double), 1, fd);
	swap_Nbyte((char *)&Header.BoxSize, 1, 8);

	Fread((void *)&Header.Omega0, sizeof(double), 1, fd);
	swap_Nbyte((char *)&Header.Omega0, 1, 8);

	Fread((void *)&Header.OmegaLambda, sizeof(double), 1, fd);
	swap_Nbyte((char *)&Header.OmegaLambda, 1, 8);

	Fread((void *)&Header.HubbleParam, sizeof(double), 1, fd);
	swap_Nbyte((char *)&Header.HubbleParam, 1, 8);

	Fread((void *)&Header.FlagAge, sizeof(int), 1, fd);
	swap_Nbyte((char *)&Header.FlagAge, 1, 8);

	Fread((void *)&Header.FlagMetals, sizeof(int), 1, fd);
	swap_Nbyte((char *)&Header.FlagMetals, 1, 8);

	Fread((void *)NallHW, N_PART_TYPES * sizeof(unsigned int), 1,fd);
	swap_Nbyte((char *)NallHW, N_PART_TYPES, 4);

	if (NallHW[0] != 0) 
		printf("Nall HW not well tested !! ");
	
	fseek(fd, dummysize, 1);
	SKIP;

	for (i = 0; i < N_PART_TYPES; i++) { /* HighWord */

		Header.Npart[i] = (long long)(Npart[i] 
                + (((long long)NallHW[i]) << 32));
		Header.Nall[i] = (long long)(Nall[i] 
                + (((long long)NallHW[i]) << 32));
		ntot += Header.Npart[i];
	}

	return (ntot);
}

int read_gadget_block(void *data, char *label, FILE * fd, size_t sizeof_type)
{
	int blocksize = 0;

	blocksize = find_block(fd, label);

	if (blocksize <= 0) {
	
		return (0);
	
	} else {
	
		SKIP;
		
		Fread(data, blocksize, 1, fd);
		
		swap_Nbyte((char *)data, blocksize / sizeof_type, 4);
		
		SKIP;
	}

	return (blocksize);
}

/* Set Block characteristics 
 * "You are all different" - 
 * "We are all different" - 
 * "I'm not !" 
 * 				(life of brian) */
void set_block_prop(enum iofields blocknr)
{
	int i;

	for (i = 0; i < N_PART_TYPES; i++) 
		Block.Npart[i] = 0;

	switch (blocknr) {
/*	case IO_POS:
		Block.Label = "POS ";	// Has to be 4 chars 
		Block.Name = "Coordinates";
		for (i = 0; i < N_PART_TYPES; i++) 
			Block.Npart[i] = Header.Npart[i];	
		Block.Val_per_element = 3;	
		Block.Data_type = FLOAT;	
		Block.Rmv_comoving = Comv2phys.Length; 
		break;*/
	case IO_VEL:
		Block.Label = "VEL ";
		Block.Name = "Velocities";
		for (i = 0; i < N_PART_TYPES; i++) 
			Block.Npart[i] = Header.Npart[i];
		Block.Val_per_element = 3;
		Block.Data_type = FLOAT;
		Block.Rmv_comoving = Comv2phys.Vel;
		break; 
	case IO_ID:
		Block.Label = "ID  ";
		Block.Name = "Particle IDs";
		for (i = 0; i < N_PART_TYPES; i++) 
			Block.Npart[i] = Header.Npart[i];
		Block.Val_per_element = 1;
		Block.Data_type = UINT;
		Block.Rmv_comoving = 1;
		break;
	case IO_U:
		Block.Label = "U   ";
		Block.Name = "Internal Energy";
		Block.Npart[0] = Header.Npart[0];
		Block.Val_per_element = 1;
		Block.Data_type = FLOAT;
		Block.Rmv_comoving = 1;
		break;
	case IO_RHO:
		Block.Label = "RHO ";
		Block.Name = "Density";
		Block.Npart[0] = Header.Npart[0];
		Block.Val_per_element = 1;
		Block.Data_type = FLOAT;
		Block.Rmv_comoving =
		    (Comv2phys.Mass / pow(Comv2phys.Length, 3));
		break;
	case IO_HSML:
		Block.Label = "HSML";
		Block.Name = "Smoothing Length";
		Block.Npart[0] = Header.Npart[0];
		Block.Val_per_element = 1;
		Block.Data_type = FLOAT;
		Block.Rmv_comoving = Comv2phys.Length;
		break;
    case IO_VOL: /* AREPO only */
		Block.Label = "VOL ";
		Block.Name = "AREPO Cell Volume";
		Block.Npart[0] = Header.Npart[0];
		Block.Val_per_element = 1;
		Block.Data_type = FLOAT;
		Block.Rmv_comoving = pow(Comv2phys.Length,3);
		break;      
	case IO_BFLD:
		Block.Label = "BFLD";
		Block.Name = "Magnetic Field";
		Block.Npart[0] = Header.Npart[0];
		Block.Val_per_element = 3;
		Block.Data_type = FLOAT;
		Block.Rmv_comoving = 1;
		break;
	case IO_MASS:
		Block.Label = "MASS";
		Block.Name = "Particle Mass";
		for (i = 0; i < N_PART_TYPES; i++) 
			Block.Npart[i] = Header.Npart[i];
		Block.Val_per_element = 1;
		Block.Data_type = FLOAT;
		Block.Rmv_comoving = Comv2phys.Mass;
		break;
	case IO_MACH:
		Block.Label = "MACH";
		Block.Name = "Mach Number";
		Block.Npart[0] = Header.Npart[0];
		Block.Val_per_element = 1;
		Block.Data_type = FLOAT;
		Block.Rmv_comoving = 1;
		break;
	case IO_VRMS:
#ifndef VGRD
		Block.Label = "VRMS";
#else
		Block.Label = "VGRD";
#endif
		Block.Name = "Local Turbulent Velocity around mean";
		Block.Npart[0] = Header.Npart[0];
		Block.Val_per_element = 1;
		Block.Data_type = FLOAT;
		Block.Rmv_comoving = Comv2phys.Vel;
		break;
#ifdef Q_SHOCK_PRIMARIES
	case IO_SHSP: // shock speed
		Block.Label = "SHSP";
		Block.Name = "Shock Speed";
		Block.Npart[0] = Header.Npart[0];
		Block.Val_per_element = 1;
		Block.Data_type = FLOAT;
		Block.Rmv_comoving = Comv2phys.Vel;
		break;
	case IO_SHCR: // shock compression ratio
		Block.Label = "SHCP";
		Block.Name = "Shock Compression Ratio";
		Block.Npart[0] = Header.Npart[0];
		Block.Val_per_element = 1;
		Block.Data_type = FLOAT;
		Block.Rmv_comoving = 1;
		break;
	case IO_SHRH: // shock rho
		Block.Label = "SHRH";
		Block.Name = "Shock Density";
		Block.Npart[0] = Header.Npart[0];
		Block.Val_per_element = 1;
		Block.Data_type = FLOAT;
		Block.Rmv_comoving = (Comv2phys.Mass / pow(Comv2phys.Length, 3));
		break;
	case IO_SHPR: // shock pressure
		Block.Label = "SHPR";
		Block.Name = "Shock Pressure";
		Block.Npart[0] = Header.Npart[0];
		Block.Val_per_element = 1;
		Block.Data_type = FLOAT;
		Block.Rmv_comoving = 1;	
		break;
	case IO_SHNR: // shock normal
		Block.Label = "SHNR";
		Block.Name = "Shock Pressure";
		Block.Npart[0] = Header.Npart[0];
		Block.Val_per_element = 3;
		Block.Data_type = FLOAT;
		Block.Rmv_comoving = 1;	
		break;
#endif	
#ifdef READ_MORE_INFO
	case IO_VBULK:
		Block.Label = "VBULK";
		Block.Name = "Local Mean Velocity";
		Block.Npart[0] = Header.Npart[0];
		Block.Val_per_element = 3;
		Block.Data_type = FLOAT;
		Block.Rmv_comoving = Comv2phys.Vel;
		break;	
	case IO_VELT:
		Block.Label = "VELT";
		Block.Name = "Local Turbulent Velocity around part";
		Block.Npart[0] = Header.Npart[0];
		Block.Val_per_element = 1;
		Block.Data_type = FLOAT;
		Block.Rmv_comoving = Comv2phys.Vel;
		break;
#endif

    case IO_VDIV:
		Block.Label = "VDIV";
		Block.Name = "Local Velocity Divergence";
		Block.Npart[0] = Header.Npart[0];
		Block.Val_per_element = 1;
		Block.Data_type = FLOAT;
		Block.Rmv_comoving = Comv2phys.Vel/Comv2phys.Length;
		break;
    case IO_VROT:
		Block.Label = "VROT";
		Block.Name = "Local Velocity Curl";
		Block.Npart[0] = Header.Npart[0];
		Block.Val_per_element = 1;
		Block.Data_type = FLOAT;
		Block.Rmv_comoving = Comv2phys.Vel/Comv2phys.Length;
		break;    
	case IO_TNGB:
		Block.Label = "TNGB";
		Block.Name = "True Number Of Neighbours";
		Block.Npart[0] = Header.Npart[0];
		Block.Val_per_element = 1;
		Block.Data_type = FLOAT;
		Block.Rmv_comoving = 1;
		break;
	case IO_DPP:
		Block.Label = "DPP ";
		Block.Name = "MagnetosReaccCoefficient";
		Block.Npart[0] = Header.Npart[0];
		Block.Val_per_element = 1;
		Block.Data_type = FLOAT;
		Block.Rmv_comoving = 1;
		break;
    /*Add above, not below !! */
	case IO_LASTENTRY:
		Block.Label = "LAST";
		Block.Name = "";
		Block.Val_per_element = 0;
		Block.Data_type = 0;
		break;
	}
	Block.Bytes_per_element = Block.Data_type * Block.Val_per_element;
	for (i = Block.Ntot = 0; i < N_PART_TYPES; i++) {
		Block.Ntot += Block.Npart[i];
	}
	return;
}

/*Fill P and SphP with data buffer 'fp'. 
 * */
void
empty_comm_buffer(enum iofields blocknr, void *fp, size_t Pindex,
		  size_t fpIndex)
{
	switch (blocknr) {
	/*case IO_POS:
		P[Pindex].Pos[0] = ((float *)fp)[fpIndex] * Block.Rmv_comoving;
		P[Pindex].Pos[1]= ((float *)fp)[fpIndex+1] * Block.Rmv_comoving;
		P[Pindex].Pos[2] = ((float *)fp)[fpIndex+2] * Block.Rmv_comoving;
		break;*/
	case IO_VEL:
		SphP[Pindex].Vel[0] = ((float *)fp)[fpIndex] * Block.Rmv_comoving;
		SphP[Pindex].Vel[1] = ((float *)fp)[fpIndex+1] * Block.Rmv_comoving;
		SphP[Pindex].Vel[2] = ((float *)fp)[fpIndex+2] * Block.Rmv_comoving;
		break;
	case IO_ID:
		P[Pindex].ID = ((unsigned int *)fp)[fpIndex];
		break;
	case IO_U:
		SphP[Pindex].U = ((float *)fp)[fpIndex] * Block.Rmv_comoving;
		break;
	case IO_RHO:
		SphP[Pindex].Rho = ((float *)fp)[fpIndex] * Block.Rmv_comoving;
		break;
    case IO_VOL: /* Store VOL in HSML for now */
        Snap.Have_Arepo = 1;    
	case IO_HSML:
		SphP[Pindex].Hsml = ((float *)fp)[fpIndex] * Block.Rmv_comoving;
		break;
	case IO_BFLD:
		SphP[Pindex].Bfld[0] = ((float *)fp)[fpIndex] * Block.Rmv_comoving;
		SphP[Pindex].Bfld[1] = ((float *)fp)[fpIndex + 1] * Block.Rmv_comoving;
		SphP[Pindex].Bfld[2] = ((float *)fp)[fpIndex + 2] * Block.Rmv_comoving;
		break;
	case IO_MASS:
		P[Pindex].Mass = ((float *)fp)[fpIndex] * Block.Rmv_comoving;
		break;
	case IO_VRMS:
		SphP[Pindex].VRms = ((float *)fp)[fpIndex] * Block.Rmv_comoving;
		break;
#ifdef READ_MORE_INFO
	case IO_VELT:
		SphP[Pindex].VTurb = ((float *)fp)[fpIndex] * Block.Rmv_comoving;
		break;
	case IO_VBULK:
		SphP[Pindex].VBulk[0] = ((float *)fp)[fpIndex] * Block.Rmv_comoving;
		SphP[Pindex].VBulk[1] = ((float *)fp)[fpIndex + 1] * Block.Rmv_comoving;
		SphP[Pindex].VBulk[2] = ((float *)fp)[fpIndex + 2] * Block.Rmv_comoving;
		break;	
#endif
    case IO_VDIV:
		SphP[Pindex].DivVel = ((float *)fp)[fpIndex] * Block.Rmv_comoving;
		break;
    case IO_VROT:
		SphP[Pindex].CurlVel = ((float *)fp)[fpIndex] * Block.Rmv_comoving;
		break;    
	case IO_TNGB:
		SphP[Pindex].TNgb = ((float *)fp)[fpIndex];
		break;
	case IO_DPP:
		SphP[Pindex].Dpp = ((float *)fp)[fpIndex] * Block.Rmv_comoving;
		break;
#ifdef Q_SHOCK_PRIMARIES
	case IO_MACH:
		SphP[Pindex].Mach = ((float *)fp)[fpIndex];
		break;
	case IO_SHSP: // shock speed
		SphP[Pindex].Shock_Velocity = ((float *)fp)[fpIndex];
		break;
	case IO_SHCR: // shock compression ratio
		SphP[Pindex].Shock_Compression_Ratio = ((float *)fp)[fpIndex];
		break;
	case IO_SHRH: // shock rho
		SphP[Pindex].Shock_Density = ((float *)fp)[fpIndex];
		break;
	case IO_SHPR: // shock pressure
		SphP[Pindex].Shock_Pressure = ((float *)fp)[fpIndex];
		break;
	case IO_SHNR:
		SphP[Pindex].Shock_Normal[0] = ((float *)fp)[fpIndex];
		SphP[Pindex].Shock_Normal[1] = ((float *)fp)[fpIndex+1];
		SphP[Pindex].Shock_Normal[2] = ((float *)fp)[fpIndex+2];
		break;
#endif
		/*Add above, not below !! */
	case IO_LASTENTRY:
		break;
	}
	return;
}


long guess_snapnum(char *fname, char *basename)
{
	long snapnum=0;
	char *token=NULL, *last_token=NULL,file[MAXLINELENGTH];

	strcpy(file, fname); /* protect fname from strtok */

	token = strtok(file,"_"); /* ptr funk */
	while (token != NULL){
        last_token = token;
	    token = strtok(NULL,"_");
        if (token != NULL)
            sprintf(basename,"%s_%s",basename,last_token);
	}
	snapnum = (long)atoi(last_token);

    memmove(basename, basename+1, MAXLINELENGTH-2);

	if (!ThisTask.Rank)
		printf("Guessing snapshot number as : %ld \n\n"
                ,snapnum);

	return(snapnum);
}
#endif

#undef SKIP
#undef int4bytes
#undef UINT
#undef INT
#undef FLOAT
#undef DBL
