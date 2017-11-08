#include "percolation.h"

int N = 128;	/* Define default lattice size */

int main(int argc, char *argv[])
{
	printf("SHPC4002, PROJECT 2: Emily Hackett, 21489688\n\n");
	
	clock_t start = clock();	/* Begin timing */

	/* CHECK COMMAND LINE ARGUMENTS:
	 * If 3 command line arguments read, assume they are:
	 * 	argv[1]: Occupancy probability, btwn 0 and 1	 
	 * 	argv[2]: Percolation type, 's' or 'b'
	 * 	argv[3]: Spanning type, 0 (col), 1 (row), 2 (both)
	 * If 4 command line arguments read, assume:
	 * 	argv[4]: Lattice size (must be power of 2)
	 */
	if (argc > 3)	{
		float occ_prob=atof(argv[1]);	/* Occupancy probability */
		if (occ_prob<0||occ_prob>1)	{
			fprintf(stderr,"ERROR: Occupancy probability must lie within [0,1]\n");
			exit(EXIT_FAILURE);
		}
		printf("	Occupancy probability = %5.4f\n",occ_prob);
		
		char perc_type=*argv[2];	/* Percolation type, site or bond */
		if (perc_type=='s'||perc_type=='b')	{
			printf("	Percolation type = %c\n",perc_type);
		}
		else	{
			fprintf(stderr,"ERROR: Percolation type must be 's' (site) or 'b' (bond)\n");
			exit(EXIT_FAILURE);
		}
		
		int span_type=atoi(argv[3]);	/* Spanning type */
		if (span_type==0||span_type==1||span_type==0)	{
			printf("	Spanning type = %i\n",span_type);
		}
		else	{
			fprintf(stderr,"ERROR: Spanning type must be 0 (col), 1 (row) or 2 (both)\n");
			exit(EXIT_FAILURE);
		}
		
		if (argc > 4)	{
			N = atoi(argv[4]);	/* Lattice size */
			if (( N & (N - 1)) != 0 )	{	/* Check if power of 2 */
				fprintf(stderr,"ERROR: Lattice size must be a power of 2\n");
				exit(EXIT_FAILURE);
			}
		}
		printf("	Lattice size = %i\n",N);
	}
	else	{
		fprintf(stderr, "Usage:	./main occ_prob perc_type span_type\n");
		exit(EXIT_FAILURE);
	}
	
	int** sites = allocate_lattice(N);	/* Allocate the site lattice */
	int** hbonds = allocate_lattice(N);	/* Allocate horizontal bonds */
	int** vbonds = allocate_lattice(N);	/* Allocate vertical bonds */
	
	/* Initialise the lattice and bonds based on percolation type */
//	if (perc_type == 's') {
//		initialise_lattice(sites,N,occ_prob);	/* Fill lattice according to occ_prob */
//		initialise_bonds(hbonds,vbonds,N,1);	/* Set all bonds equal to 1 */
//	}
//	if (perc_type == 'b') {
//		initialise_lattice(sites,N,occ_prob);	/* Set all sites to occupied */
//		initialise_bonds(hbonds,vbonds,N,occ_prob);	/* Fill bonds according to occ_prob */
//	}
	
	clock_t end=clock();	/* End the timer */
	double time_spent = (double)(end-start)/CLOCKS_PER_SEC;
	
	printf("TOTAL TIME: %10.6f\n",time_spent);
	
	return 0;
}
