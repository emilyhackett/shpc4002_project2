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

/* Fills the lattice site according to the occupancy probability.
 * 	If the occupancy probability is passed as 1, all lattice sites are occupied 
 * 	(such as required for bond percolation)
 */
void initialise_lattice(int** lattice, int N, float occ_prob)
{
}

/* Fills the horizontal and vertical bonds according to the occupancy probability.
 * 	If the occupancy probability is passed as 1, all bonds exist 
 * 	(such as required for site percolation)
 */
void initialise_bonds(int** hbonds, int** vbonds, int N, float occ_prob)
{
}
