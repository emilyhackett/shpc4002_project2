#include "percolation.h"

/* Allocates space for the MxM lattice, checking for memory faults */
int** allocate_lattice(int N)
{
	int i,j;
	int** lattice=malloc(N * sizeof(int*));
	if (lattice == NULL)	{	/* Check memory alloc */
		fprintf(stderr,"ERROR: Lattice allocation unsuccessful\n");
	}
	
	for(i = 0; i < N; i++) {
		lattice[i]=malloc(N * sizeof(int));
		if (lattice[i] == NULL)	{	/* Check memory alloc */
			fprintf(stderr,"ERROR: Lattice allocation unsuccessful\n");
		}
	}
	
	return lattice;
}

/* Fills the lattice site or bonds according to the occupancy probability.
 * 	If the occupancy probability is passed as less than 0, all lattice sites are occupied 
 * 	(such as required for bond percolation) or bonds (for site percolation)
 */
void initialise_lattice(int** lattice, int N, float occ_prob)
{
	int i,j;
	float num;
	
	srand(time(NULL));	/* Initialise random seed */
	
	for(i = 0; i < N; i++) {
		for(j = 0; j < N; j++) {
			
			if (occ_prob > 0) {
			
			num = (double)rand()/(double)RAND_MAX;	/* Random number */
				
				if (num <= occ_prob) {	/* Check against probability */
					lattice[i][j]=1;
				}
				else	{
					lattice[i][j]=0;
				}

			}
			else	{
				/* No random number, simply set to 1 */
				lattice[i][j]=1;
			}
		}
	}

}

/* Print the lattice:
 * 	Each lattice site takes up a 3x3 region
 */
void display_lattice(int** sites, int** hbonds, int** vbonds, int N)
{
	int i,j;
	
	printf("\n%ix%i LATTICE:\n",N,N);
	if (N > 32) {
		printf("Warning! Will only print lattice 32x32 and smaller\n\n");
		return;
	}

	for(i = 0; i < N; i++) {	/* Loop over rows */
		
		for(j = 0; j < N; j++ )	{	/* Loop over columns */

			/* Print site */
			if (sites[i][j] == 1) {
				printf("o");
			}
			else	{
				printf(" ");
			}

			/* Print horizontal bond */
			if (hbonds[i][j] == 1) {
				printf("-");
			}
			else	{
				printf(" ");
			}

		}

		printf("\n");

		for(j = 0; j < N; j++ )	{
			
			/* Print vertical bond */
			if (vbonds[i][j] == 1) {
				printf("| ");
			}
			else {
				printf("  ");
			}
		}

		printf("\n");

	}

	printf("\n");
	
}
