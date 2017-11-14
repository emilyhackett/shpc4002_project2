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

/* Allocates space for the Nxchunk lattice, checking for memory faults */
int** allocate_lattice(int N, int chunk)
{
	int i,j;

	/* NOTE! Modified allocation to make sure contiguous (for MPI transfer) */
	int* data = malloc(N * chunk * sizeof(int));
	int** lattice=malloc(chunk * sizeof(int*));
	if (lattice == NULL || data == NULL)	{	/* Check memory alloc */
		fprintf(stderr,"ERROR: Lattice allocation unsuccessful\n");
	}
	
	for(i = 0; i < chunk; i++) {
		lattice[i]=&(data[N*i]);
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
void display_lattice(int** sites, int** hbonds, int** vbonds, int N, int chunk_size)
{
	int i,j;
	
	printf("\n%ix%i LATTICE:\n",N,chunk_size);
	if (N > 32) {
		printf("Warning! Will only print lattice 32x32 and smaller\n\n");
		return;
	}

	for(i = 0; i < chunk_size; i++) {	/* Loop over rows */
		
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
CLUSTER* depth_first_search(int** sites, int** hbonds, int** vbonds, int N, int* chunk, int row, int col, CLUSTER* tmp, int rank, int id)
{
	int right_col = wrap(col + 1, N);
	int down_row = row + 1;
	int left_col = wrap(col - 1, N);
        int up_row = row - 1;
	printf("	%i: Thread %i at [%i,%i], will check R[%i,%i] D[%i,%i] L[%i,%i] U[%i,%i]\n",rank,id,row,col,row,right_col,down_row,col,row,left_col,up_row,col);
	
	sites[row][col] = -1;	/* Mark current site as visited */

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	
	/* Check that a horizontal bond exists to the next occupied site (right) */
	if (hbonds[row][col] == 1 && sites[row][right_col] == 1)	{

		/* Add this node to the cluster */
		tmp->num_nodes = tmp->num_nodes + 1;
		tmp->cols_reached[right_col] = tmp->cols_reached[right_col] + 1;
		tmp->rows_reached[row] = tmp->rows_reached[row] + 1;

		/* Continue depth first search */
		tmp = depth_first_search(sites,hbonds,vbonds,N,chunk,row,right_col,tmp,rank,id);
	}

	/* Check that a vertical bond exists to the next occupied site (down) */
	if (down_row <= chunk[1])	{
		if (vbonds[row][col] == 1 && sites[down_row][col] != 0)	{
			/* If we can get to this point, add the current site to the bounds list (bottom) */
			if (row == chunk[1])	{	tmp->bottom_bounds[col] = row + 1;	}

			if (down_row >= chunk[0] && down_row <= chunk[1] && sites[down_row][col] == 1)	{	/* If no wrapping occured */
				/* Add this node to the cluster */
				tmp->num_nodes = tmp->num_nodes + 1;
				tmp->cols_reached[col] = tmp->cols_reached[col] + 1;
				tmp->rows_reached[down_row] = tmp->rows_reached[down_row] + 1;

				/* Continue depth first search */
				tmp = depth_first_search(sites,hbonds,vbonds,N,chunk,down_row,col,tmp,rank,id);
			}
		}
	}
	
	/* Check that a horizontal bond exists to the previous occupied site (left) */
	if (hbonds[row][left_col] == 1 && sites[row][left_col] == 1)	{
		/* Add this node to the cluster */
		tmp->num_nodes = tmp->num_nodes + 1;
		tmp->cols_reached[left_col] = tmp->cols_reached[left_col] + 1;
		tmp->rows_reached[row] = tmp->rows_reached[row] + 1;
		
		/* Continue the depth first search */
		tmp = depth_first_search(sites,hbonds,vbonds,N,chunk,row,left_col,tmp,rank,id);
	}

	/* Check that a vertical bond exists to the next occupied site (up) */
	if (up_row >= chunk[0])	{
		if (vbonds[up_row][col] == 1 && sites[up_row][col] != 0)	{
			/* If we can get to this point, add the current site to the bounds list (top) */
			if (row == chunk[0])	{	tmp->top_bounds[col] = row + 1;	}	
	
			if (up_row >= chunk[0] && up_row <= chunk[1] && sites[up_row][col] == 1)	{	/* If no wrapping occured */
				/* Add this node to the cluster */
				tmp->num_nodes = tmp->num_nodes + 1;
				tmp->cols_reached[col] = tmp->cols_reached[col] + 1;
				tmp->rows_reached[up_row] = tmp->rows_reached[up_row] + 1;
				
				/* Continue the depth first search */
				tmp = depth_first_search(sites,hbonds,vbonds,N,chunk,up_row,col,tmp,rank,id);
			}
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

	tmp->merged = 0;

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
NODE* merge_cluster_lists(NODE* head_A, NODE* head_B, int N, int end_idx_A)
{
	int i;

	NODE* new_list = NULL;	/* Initialise new list to export */
	NODE* merge_list = NULL;	/* Initialise new list of merging clusters */

	NODE* current;
	
	int* control = malloc(N * sizeof(int));	
	for (i = 0; i < N; i++ )	{
		control[i] = 1;
	}
	
	/* Loop over list A (bottom segment) and remove clusters not on the bounds */
	current = head_A;
	while (current != NULL)	{
		if (check_bounds_crossover(current->data->bottom_bounds,control,N) == 0)	{
			/* No cluster on the boundaries, add to the new list */
			new_list = push(new_list, current->data);

		}
		else	{	/* Add to the list to be merged */
			merge_list = push(merge_list, current->data);
		}

		current = current->next;	
	}

	/* Loop over list B (top segment) and remove clusters not on the bounds */
	current = head_B;
	while (current != NULL)	{
		if (check_bounds_crossover(current->data->top_bounds,control,N) == 0)	{
			/* No cluster on the boundaries, add to the new list */
			new_list = push(new_list, current->data);

		}
		else	{	/* Add to the merge list */
			merge_list = push(merge_list, current->data);
		}
		
		current = current->next;
	}

//	printf("MERGE LIST: ");
//	display_list(merge_list);

	current = merge_list;
	while (current != NULL)	{	/* Loop over the merge list */
		
		if (current->data->merged != 1)	{
			//printf("Checking %i node cluster with bounds top: %i %i %i %i, bottom: %i %i %i %i ...\n",current->data->num_nodes,current->data->top_bounds[0],current->data->top_bounds[1],current->data->top_bounds[2],current->data->top_bounds[3],current->data->bottom_bounds[0],current->data->bottom_bounds[1],current->data->bottom_bounds[2],current->data->bottom_bounds[3]);;
		
			CLUSTER* new_cluster = initialise_cluster(N,1,1);	
			new_cluster->cols_reached[1] = 0;	/* Fix allocation */
			new_cluster->rows_reached[1] = 0;
			new_cluster->num_nodes = 0;

			new_cluster = merge(current->data, merge_list, new_cluster, N);

			new_list = push(new_list, new_cluster);
		}

		current=current->next;
	}
	
	return new_list;
}

/* Recursively checks each cluster in a list against the main cluster to determine
 * whether or not they should be merged or not.
 */
CLUSTER* merge(CLUSTER* current, NODE* head, CLUSTER* new_cluster, int N)
{
	NODE* loop = head;

	current->merged = 1;

	/* Add this cluster information to the new cluster */
	new_cluster->num_nodes = new_cluster->num_nodes + current->num_nodes;
	
	int i;
	for (i = 0; i < N; i++)	{
		new_cluster->cols_reached[i] = new_cluster->cols_reached[i] + current->cols_reached[i];
		new_cluster->rows_reached[i] = new_cluster->rows_reached[i] + current->rows_reached[i];
	}		

	while(loop != NULL)	{	/* Loop over list */

		if (loop->data->merged != 1)	{	/* If cluster not merged, check */

			/* Check bounds (if on top) */
			if (loop->data->bottom_row_idx + 1 == current->top_row_idx)	{
				//printf(" ... against %i node cluster with bottom bounds: %i %i %i %i\n",loop->data->num_nodes,loop->data->bottom_bounds[0],loop->data->bottom_bounds[1],loop->data->bottom_bounds[2],loop->data->bottom_bounds[3]);
				if ( check_bounds_crossover(loop->data->bottom_bounds, current->top_bounds, N) == 1)	{

					/* If match up on bounds, merge */
					new_cluster = merge(loop->data, head, new_cluster, N);
				}
			}
			
			/* Check bounds (if on bottom) */			
			if (loop->data->top_row_idx + 1 == current->bottom_row_idx)	{

				//printf(" ... against %i node cluster with top bounds: %i %i %i %i\n",loop->data->num_nodes,loop->data->top_bounds[0],loop->data->top_bounds[1],loop->data->top_bounds[2],loop->data->top_bounds[3]);
				if ( check_bounds_crossover(loop->data->top_bounds, current->bottom_bounds, N) == 1)	{
					
					/* If match up on bounds, merge */
					new_cluster = merge(loop->data, head, new_cluster, N);
				}
			}
		}

		loop = loop->next;
	}

	return new_cluster;

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

