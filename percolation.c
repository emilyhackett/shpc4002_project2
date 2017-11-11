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
	
	/* Check that a horizontal bond exists to the next occupied site (right) */
	int right_col = wrap(col + 1, N);	
	if (hbonds[row][col] == 1 && sites[row][right_col] == 1)	{
		/* Add this node to the cluster */
		tmp->num_nodes = tmp->num_nodes + 1;
		tmp->cols_reached[right_col] = tmp->cols_reached[right_col] + 1;
		tmp->rows_reached[row] = tmp->rows_reached[row] + 1;

		/* Add column to top bounds */
		if (row == chunk[0])	{	tmp->top_bounds[col]=row;	}
		/* Add column to bottom bounds */
		if (row == chunk[1])	{	tmp->bottom_bounds[col]=row;	}
		
		/* Continue depth first search */
		tmp = depth_first_search(sites,hbonds,vbonds,N,chunk,row,right_col,tmp);
	}

	/* Check that a vertical bond exists to the next occupied site (down) */
	int down_row = wrap(row + 1, N);	
	if (vbonds[row][col] == 1 && sites[down_row][col] != 0)	{
		/* If we can get to this point, add the current site to the bounds list (bottom) */
		if (row == chunk[1])	{	tmp->bottom_bounds[col] = row + 1;	}

		if (down_row >= chunk[0] && down_row <= chunk[1] && sites[down_row][col] == 1)	{	/* If no wrapping occured */
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
	int up_row = wrap(row - 1, N);	
	if (vbonds[up_row][col] == 1 && sites[up_row][col] != 0)	{
		/* If we can get to this point, add the current site to the bounds list (top) */
		if (row == chunk[0])	{	tmp->top_bounds[col] = row + 1;	}	
	
		if (up_row >= chunk[0] && up_row <= chunk[1] && sites[up_row][col] == 1)	{	/* If no wrapping occured */
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
 * 	
 * 	A -> top segment (lower thread id)
 * 	B -> bottom segment (higher thread id)
 */
NODE* merge_cluster_lists(NODE* head_A, NODE* head_B, int N)
{
	NODE* new_list = NULL;	/* Initialise new list to export */
	NODE* merge_list = NULL;	/* Initialise new list to track merges */

	NODE* current_A = head_A;	/* Current cluster in top segment */
	NODE* current_B	= head_B;	/* Current cluster in bottom segment */

	CLUSTER* cluster_A;	/* Cluster data for A */
	CLUSTER* cluster_B;	/* Cluster data for B */

	NODE* previous_A = NULL;	/* Tracks previous cluster, to remove from linked list */
	NODE* previous_B = NULL;

	int* control = malloc(N * sizeof(int));	
	memset(control, 1, N);	/* Used to see if a cluster lies on the boundary */

	int changed = 1;	/* Tracks changes */
	
	/* Loop over list A (bottom segment) and remove clusters not on the bounds */
	while (current_A != NULL)	{
		cluster_A = current_A->data;

		if (check_bounds_crossover(cluster_A->bottom_bounds,control,N) == 0)	{
			/* No cluster on the boundaries, add to the new list */
			new_list = push(new_list, cluster_A);

			if (previous_A == NULL)	{	/* At the start of the list, tf pop off */
				head_A = pop(head_A, cluster_A);
			}
			else	{
				previous_A->next = current_A->next;	/* Remove current from list */
			}

		}
		else	{
			/* Cluster exists on the bounds, add to merge list */
			merge_list = push(merge_list, cluster_A);
		}

		previous_A = current_A;
		current_A = current_A->next;
		
	}

	/* Loop over list B (top segment) and remove clusters not on the bounds */
	while (current_B != NULL)	{
		cluster_B = current_B->data;

		if (check_bounds_crossover(cluster_B->top_bounds,control,N) == 0)	{
			/* No cluster on the boundaries, add to the new list */
			new_list = push(new_list, cluster_B);

			if (previous_B == NULL)	{	/* At the start of the list, tf pop off */
				head_B = pop(head_B, cluster_B);
			}
			else	{
				previous_B->next = current_B->next;	/* Remove current from list */
			}
		}
		else	{
			/* Cluster exists on the bounds, add to merge list */
			merge_list = push(merge_list, cluster_B);
		}

		previous_B = current_B;
		current_B = current_B->next;
	}

	printf("MERGE LIST:  ");
	display_list(merge_list);

	return new_list;
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

