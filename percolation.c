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
CLUSTER* depth_first_search(int** sites, int** hbonds, int** vbonds, int N, int chunk_size, int row, int col, CLUSTER* tmp)
{	
	sites[row][col] = -1;	/* Mark current site as visited */
	
	/* Check that a horizontal bond exists to the next occupied site (right) */
	int right_col = wrap(col + 1, N);	
	if (hbonds[row][col] == 1 && sites[row][right_col] == 1)	{
		/* Add this node to the cluster */
		tmp->num_nodes = tmp->num_nodes + 1;
		tmp->cols_reached[right_col] = tmp->cols_reached[right_col] + 1;
		tmp->rows_reached[row] = tmp->rows_reached[row] + 1;
		
		/* Continue depth first search */
		tmp = depth_first_search(sites,hbonds,vbonds,N,chunk_size,row,right_col,tmp);
	}

	/* Check that a vertical bond exists to the next occupied site (down) */
	int down_row = wrap(row + 1, N);
	if (vbonds[row][col] == 1 && sites[down_row][col] == 1)	{
		/* Add this node to the cluster */
		tmp->num_nodes = tmp->num_nodes + 1;
		tmp->cols_reached[col] = tmp->cols_reached[col] + 1;
		tmp->rows_reached[down_row] = tmp->rows_reached[down_row] + 1;
		
		/* Continue depth first search */
		tmp = depth_first_search(sites,hbonds,vbonds,N,chunk_size,down_row,col,tmp);
	}
	
	/* Check that a horizontal bond exists to the previous occupied site (left) */
	int left_col = wrap(col - 1, N);
	if (hbonds[row][left_col] == 1 && sites[row][left_col] == 1)	{
		/* Add this node to the cluster */
		tmp->num_nodes = tmp->num_nodes + 1;
		tmp->cols_reached[left_col] = tmp->cols_reached[left_col] + 1;
		tmp->rows_reached[row] = tmp->rows_reached[row] + 1;
		
		/* Continue the depth first search */
		tmp = depth_first_search(sites,hbonds,vbonds,N,chunk_size,row,left_col,tmp);
	}

	/* Check that a vertical bond exists to the next occupied site (up) */
	int up_row = wrap(row - 1,N);
	if (vbonds[up_row][col] == 1 && sites[up_row][col] == 1)	{
		/* Add this node to the cluster */
		tmp->num_nodes = tmp->num_nodes + 1;
		tmp->cols_reached[col] = tmp->cols_reached[col] + 1;
		tmp->rows_reached[up_row] = tmp->rows_reached[up_row] + 1;

		/* Continue the depth first search */
		tmp = depth_first_search(sites,hbonds,vbonds,N,chunk_size,up_row,col,tmp);
	}

	return tmp;
}

/* Allocates and clears a new cluster to be added to a node in the linked list */
CLUSTER* initialise_cluster(int N, int row, int col)
{
	CLUSTER* tmp;

	tmp = malloc( sizeof( CLUSTER ) );	/* Allocate new cluster in memory */
	if (tmp == NULL)	{
		fprintf(stderr,"ERROR: Cluster could not be allocated.\n");
		exit(EXIT_FAILURE);
	}

	tmp->num_nodes = 1;	/* Set the number of nodes to 1 to count start point */

	tmp->cols_reached = malloc(N * sizeof(int));	/* Allocate columns reached array */
	tmp->rows_reached = malloc(N * sizeof(int));	/* Allocate rows reached array */

	int i;
	for (i = 0; i < N; i++ )	{
		tmp->cols_reached[i] = 0;	/* Set values to 0 */
		tmp->rows_reached[i] = 0;
	}

	tmp->cols_reached[col] = 1;	/* Set starting point as visited */
	tmp->rows_reached[row] = 1;
	
	return tmp;
}

/* Checks if the cluster spans all columns or rows (or both) based on spanning type */
int check_spanning(CLUSTER* tmp, int N, int span_type)
{	
	int i;
	int row_sum = 0;
	int col_sum = 0;

	for(i = 0; i < N; i++)	{
		if (tmp->rows_reached[i] != 0)	{
			row_sum = row_sum + 1;
		}

		if (tmp->cols_reached[i] != 0)	{
			col_sum = col_sum + 1;
		}
	}

	if (row_sum == N && col_sum == N && span_type == 2)	{
		return 1;
	}
	else if (row_sum == N && span_type == 0)	{
		return 1;
	}
	else if (col_sum == N && span_type == 1)	{
		return 1;
	}
	else	{
		return 0;
	}
}
