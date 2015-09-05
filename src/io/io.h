#define MAXTAGS 300

extern void read_param_file(char *filename);

/* This Function fills the Particle structure
 * P and G, and sets in the Comv structure.
 * It should do so in parallel. It returns 0
 * if there is no more snapshot
 * */
extern double read_input();

/* This functions output one or more files in 
 * F77 binary containing the spectra sorted by
 * ID */

extern void write_spectra();
extern void read_spectra();

/* Parameter File Tags, also used to write fits header */
int id[MAXTAGS];
void *addr[MAXTAGS];
char tag[MAXTAGS][50],comment[MAXTAGS][50];

int Find_files(char *);
