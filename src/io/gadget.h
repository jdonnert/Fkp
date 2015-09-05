#ifdef I_GADGET2		// Additional Stuff for GADGET-2 Input 

int swap;

enum iofields {
//	IO_POS,
	IO_VEL,
	IO_ID,
	IO_U,
	IO_RHO,
	IO_HSML,
    IO_VOL, // AREPO only 
	IO_BFLD,
	IO_MASS,
	IO_VRMS,
#ifdef READ_MORE_INFO
	IO_VELT,
	IO_VBULK,
#endif
    IO_VDIV,
    IO_VROT,
	IO_TNGB,
    IO_DPP,
	IO_MACH,
#ifdef Q_SHOCK_PRIMARIES
	IO_SHSP, // shock speed
	IO_SHCR, // shock compression ratio
	IO_SHRH, // shock rho
	IO_SHPR, // shock pressure
	IO_SHNR, // shock normal
#endif
	IO_LASTENTRY		// Keep this entry at the end for termination 
};


long guess_snapnum(char*,char*);

struct Gadget_head { /*Input Header*/
	long long Npart[N_PART_TYPES];
	double Mass[N_PART_TYPES];
	double Time;
	double Redshift;
	long long Nall[N_PART_TYPES];
	long NumFiles;
	double BoxSize;
	double Omega0;
	double OmegaLambda;
	double HubbleParam;
	int FlagSfr;
	int FlagFeedback;
	int FlagCooling;
	int FlagAge;
	int FlagMetals;
} Header;

/* Block specifications */
struct blockdef {
	char *Label;
	char *Name;
	void *DataPtr;
	long long Npart[N_PART_TYPES];
	long long Ntot;
	int Val_per_element;
	size_t Data_type;
	size_t Bytes_per_element;
	double Rmv_comoving;
} Block;
#endif
