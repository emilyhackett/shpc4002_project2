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
 *
 * MODIFIED -> For parallel sectioning, remove wrapping between top and bottom rows (since 
 * aren't actually connected in the entire cluster.
 */
CLUSTER* depth_first_search(int** sites, int** hbonds, int** vbonds, int N, int* chunk, int row, int col, CLUSTER* tmp)
{	
	sites[row][col] = -1;	/* Mark current site as visited */

	if (row == chunk[0])	{	/* Add column to bottom bounds */
		tmp->bottom_bounds[col]=row;
	}

	if (row == chunk[1])	{	/* Add column to top bounds */
		tmp->top_bounds[col]=row;
	}
	
	/* Check that a horizontal bond exists to the next occupied site (right) */
	int right_col = wrap(col + 1, N);	
	if (hbonds[row][col] == 1 && sites[row][right_col] == 1)	{
		/* Add this node to the cluster */
		tmp->num_nodes = tmp->num_nodes + 1;
		tmp->cols_reached[right_col] = tmp->cols_reached[right_col] + 1;
		tmp->rows_reached[row] = tmp->rows_reached[row] + 1;
		
		/* Continue depth first search */
		tmp = depth_first_search(sites,hbonds,vbonds,N,chunk,row,right_col,tmp);
	}

	/* Check that a vertical bond exists to the next occupied site (down) */
	int down_row = row + 1;	/* Remove wrapping */
	if (down_row >= chunk[0] && down_row <= chunk[1])	{	/* If no wrapping occured */
		if (vbonds[row][col] == 1 && sites[down_row][col] == 1)	{
			/* Add this node to the cluster */
			tmp->num_nodes = tmp->num_nodes + 1;
			tmp->cols_reached[col] = tmp->cols_reached[col] + 1;
			tmp->rows_reached[down_row] = tmp->rows_reached[down_row] + 1;

			/* Continue depth first search */
			tmp = depth_first_search(sites,hbonds,vbonds,N,chunk,down_row,col,tmp);
		}
	}
	
	/* Check that a horizontal bond exists to the previous occupied site (left) */
	int left_col = wrap(col - 1, N);
	if (hbonds[row][left_col] == 1 && sites[row][left_col] == 1)	{
		/* Add this node to the cluster */
		tmp->num_nodes = tmp->num_nodes + 1;
		tmp->cols_reached[left_col] = tmp->cols_reached[left_col] + 1;
		tmp->rows_reached[row] = tmp->rows_reached[row] + 1;
		
		/* Continue the depth first search */
		tmp = depth_first_search(sites,hbonds,vbonds,N,chunk,row,left_col,tmp);
	}

	/* Check that a vertical bond exists to the next occupied site (up) */
	int up_row = row - 1;	/* Remove wrapping */
	if (up_row >= chunk[0] && up_row <= chunk[1])	{	/* If no wrapping occured */
		if (vbonds[up_row][col] == 1 && sites[up_row][col] == 1)	{
			/* Add this node to the cluster */
			tmp->num_nodes = tmp->num_nodes + 1;
			tmp->cols_reached[col] = tmp->cols_reached[col] + 1;
			tmp->rows_reached[up_row] = tmp->rows_reached[up_row] + 1;

			/* Continue the depth first search */
			tmp = depth_first_search(sites,hbonds,vbonds,N,chunk,up_row,col,tmp);
		}
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

	tmp->top_bounds = malloc(N * sizeof(int));	/* Allocate top row boundaries */
	tmp->bottom_bounds = malloc(N * sizeof(int));	/* Allocate bottom row boundaries */

	int i;
	for (i = 0; i < N; i++ )	{
		tmp->cols_reached[i] = 0;	/* Set values to 0 */
		tmp->rows_reached[i] = 0;
		tmp->top_bounds[i]=0;
		tmp->bottom_bounds[i]=0;
	}

	tmp->cols_reached[col] = 1;	/* Set starting point as visited */
	tmp->rows_reached[row] = 1;
	
	return tmp;
}

/* Merges two adjacent cluster lists.
 * 	Returns the pointer to the head node of the new combined list.
 */
NODE* merge_cluster_lists(NODE* bottom_list_head, NODE* top_list_head, int N)
{
	NODE* new_head = NULL;

	NODE* bottom_node = bottom_list_head;

	CLUSTER* bottom_cluster;
	CLUSTER* top_cluster;

	int* control = malloc(N * sizeof(int));
	memset(control, 1, N);
	
	while (bottom_node != NULL)	{
		
		NODE* top_node = top_list_head;	/* Reset to list head */
		
		int changed = 0;

		/* Check if on the boundary */
		if (check_bounds_crossover(bottom_cluster->top_bounds,control,N) == 0)	{
			
		}
		
		else	{	/* Compare against top segment clusters */
			
			while (top_node != NULL)	{
			
				bottom_cluster = bottom_node->data;
				top_cluster = top_node->data;
		
				if (check_bounds_crossover(top_cluster->bottom_bounds,control,N) == 0)	{
					changed = 0;
				}
				else if (check_bounds_crossover(top_cluster->bottom_bounds,bottom_cluster->top_bounds,N) == 1)	{
					changed = 1;
				}

				top_node = top_node->next;	/* Move to next node/cluster */	
				
			}
		}
			bottom_node = bottom_node->next;	/* Move to next node/cluster */
	}	
	
	return new_head;
}

/* Take two bounds array and checks to see if they match up */
int check_bounds_crossover(int* bounds_1, int* bounds_2, int N)
{
	int bounds_sum = 0;
	int i;

	for(i = 0; i < N; i++ )	{
		bounds_sum = bounds_sum + bounds_1[i]*bounds_2[i];
	}

	if (bounds_sum == 0)	{
		return 0;
		}
	else	{
		return 1;
	}
	
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

