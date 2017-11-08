#include "percolation.h"

/* Returns the wrapped value of an integer */
int wrap(int val, int max)
{
	int upper = max - 1;
	int lower = 0;

	int range = upper - lower + 1;
	val = ((val - lower) % range);

	if (val < 0) {
		return upper + 1 + val;
	}
	else	{
		return lower + val;
	}
}

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
 * 	Will only print the lattice if it is 32x32 size or smaller. 
 * 	Otherwise takes too long and is not particuarly useful
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

/* Depth first search: Recursively checks each lattice site given that the lattice site 
 * is occupied, and that a bond exists that connects it to the previous. 
 */
int depth_first_search(int** sites, int** hbonds, int** vbonds, int N, int row, int col)
{
	int num_nodes;
	int temp_nodes;
	
	sites[row][col] = -1;	/* Mark current site as visited */

	num_nodes = 1;
	
	/* Check that a horizontal bond exists to the next occupied site (left) */
	int left_col = wrap(col + 1, N);	
	if (hbonds[row][col] == 1 && sites[row][left_col] == 1)	{
		temp_nodes = depth_first_search(sites,hbonds,vbonds,N,row,left_col);
		num_nodes = num_nodes + temp_nodes;
	}

	/* Check that a vertical bond exists to the next occupied site (down) */
	int down_row = wrap(row + 1, N);
	if (vbonds[row][col] == 1 && sites[down_row][col] == 1)	{
		temp_nodes = depth_first_search(sites,hbonds,vbonds,N,down_row,col);
		num_nodes = num_nodes + temp_nodes;
	}
	
	/* Check that a horizontal bond exists to the previous occupied site (right) */
	int right_col = wrap(col - 1, N);
	if (hbonds[row][right_col] == 1 && sites[row][right_col] == 1)	{
		temp_nodes = depth_first_search(sites,hbonds,vbonds,N,row,right_col);
		num_nodes = num_nodes + temp_nodes;
	}

	/* Check that a vertical bond exists to the next occupied site (up) */
	int up_row = wrap(row - 1,N);
	if (vbonds[up_row][col] == 1 && sites[up_row][col] == 1)	{
		temp_nodes = depth_first_search(sites,hbonds,vbonds,N,up_row,col);
		num_nodes = num_nodes + temp_nodes;
	}

	return num_nodes;
}

