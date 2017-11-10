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
