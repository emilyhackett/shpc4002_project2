#include "percolation.h"

/* Adds an element to the linked list, returns the new head node */
struct NODE* push(struct NODE* head, struct CLUSTER* data)
{
	struct NODE* tmp = (struct NODE*)malloc( sizeof(struct NODE) );	/* Allocate new node */

	if (tmp == NULL)	{
		fprintf(stderr,"ERROR: New node in linked list unallocated.\n");
		exit(EXIT_FAILURE);
	}

	tmp->data = data;	/* Place data in new structure */
	tmp->next = head;	/* Set previous list head as the next node */

	head=tmp;		/* Set new node as the new head node */

	return head;
}

/* Pop an element from the top of the linked list, return the new head node */
struct NODE* pop(struct NODE* head, struct CLUSTER* data)
{
	data = head->data;	/* Save the data into the memory location passed */

	head = head->next;	/* Set the head to the next node */

	return head;		/* Return the new head node */
}

/* Display the entire linked list of clusters.
 * 	If bond percolation, will ignore the 'clusters' that have 1 node only 
 */
void display_list(struct NODE* head)
{
	NODE* current;
	CLUSTER* data;
	current = head;
	
	int i=0;
	int max_print_clusters = 30;
	
	if (current == NULL)	{
		printf("No clusters in the list.\n");
		return;
	}
	
	printf("CLUSTERS -->\n");
	while ( current != NULL )	{
		data = current->data;
		printf("  Cluster %i with %i nodes\n",i,data->num_nodes);
		i++;
		current=current->next;

		if (i > max_print_clusters)	{
			printf("Warning! Will only display up to %i clusters.\n",max_print_clusters);
			return;
		}
	}
	printf("\n");
}

/* Traverse the linked list and check if spanning cluster exists, and also
 * find the maximum number of nodes.
 * 	Returns the cluster number of the one with the max number of nodes
 * 	(NOTE! This doesn't necessarily mean that this is the spanning cluster.
 */
int traverse_list(struct NODE* head, int N, int span_type, int* spanning, int* max_num_nodes)
{
	NODE* current;
	CLUSTER* data;
	current = head;

	int span_result;	/* Variable to hold the spanning result */
	int max_result = 0;	/* Initialise maximum number of nodes as 0 */
	int max_cluster_idx;

	int i = 0;	/* Counts how far into the list we are */

	if (current == NULL)	{
		printf("No clusters in the list.\n");
		return -1;
	}
	
	while (current != NULL)	{	/* Loop until the end of the list reached */
		data = current->data;
		i++;
		
		if (span_result != 1)	{	/* If no spanning cluster found, check for it */
			span_result = check_spanning(data,N,span_type);
		}

		if (data->num_nodes > max_result)	{	/* If larger than max, replace */
			max_result = data->num_nodes;
			max_cluster_idx = i;	/* Save cluster index (relative to list) */
		}

		current=current->next;
	}

	*spanning = span_result;	/* Save span result to variable */
	*max_num_nodes = max_result;	/* Save max num nodes to variable */

	return max_cluster_idx;
}

/* Helper function for linkedlist_to_array below */
void linkedlist_malloc(int num_elements, int N, int* num_nodes, int* top_row_idx, int* bottom_row_idx, int**top_bounds, int** bottom_bounds, int** cols_spanned, int** rows_spanned)
{
        num_nodes = malloc(num_elements * sizeof(int));
        top_row_idx = malloc(num_elements * sizeof(int));
        bottom_row_idx = malloc(num_elements * sizeof(int));
	if (num_nodes == NULL || top_row_idx == NULL || bottom_row_idx == NULL)	{
		fprintf(stderr,"MALLOC CALL UNSUCCESSFUL.\n");
	}

	int* topb_data = malloc(num_elements * N * sizeof(int));	
	int* bottomb_data = malloc(num_elements * N * sizeof(int));
	int* cols_data = malloc(num_elements * N * sizeof(int));
	int* rows_data = malloc(num_elements * N * sizeof(int));
	if (topb_data == NULL || bottomb_data == NULL || cols_data == NULL || rows_data == NULL)	{
		fprintf(stderr,"MALLOC CALL UNSUCCESSFUL.\n");
		exit(EXIT_FAILURE);
	}

	top_bounds = malloc(num_elements * sizeof(int*));
	bottom_bounds = malloc(num_elements * sizeof(int*));
	cols_spanned = malloc(num_elements * sizeof(int*));
	rows_spanned = malloc(num_elements * sizeof(int*));

	int i;
	for (i = 0; i < num_elements; i++)	{
		top_bounds[i] = &(topb_data[N*i]);
		bottom_bounds[i] = &(bottomb_data[N*i]);
		cols_spanned[i] = &(cols_data[N*i]);
		rows_spanned[i] = &(rows_data[N*i]);
	}	
}

/* Converts a linked list into an array (in order to transfer over mpi) */
void linkedlist_to_array(NODE* head, int num_elements, int N, int* num_nodes, int* top_row_idx, int* bottom_row_idx, int** top_bounds, int** bottom_bounds, int** cols_spanned, int** rows_spanned)
{
	linkedlist_malloc(num_elements, N, num_nodes, top_row_idx, bottom_row_idx, top_bounds, bottom_bounds, cols_spanned, rows_spanned);

	int k;
	int i = 0;	/* Count place in cluster array */
	/* Allocate a contiguous array */
	while( head != NULL)	{
		CLUSTER* current = head->data;

/*		num_nodes[i] = current->num_nodes;
		printf(" - num_nodes transferred\n");
		top_row_idx[i] = current->top_row_idx;
		printf(" - top_row_idx transferred\n");
		bottom_row_idx[i] = current->bottom_row_idx;
		printf(" - bottom_row_idx transferred\n");

		printf("hello\n");

		for (k = 0; k < N; k++)	{
			top_bounds[i][k] = current->top_bounds[k];
			bottom_bounds[i][k] = current->bottom_bounds[k];
			cols_spanned[i][k] = current->cols_reached[k];
			rows_spanned[i][k] = current->rows_reached[k];
		}
*/
		head = head->next;
		i++;
	}
}	
