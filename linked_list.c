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
void display_list(struct NODE* head, int num_elements)
{
	NODE* current;
	CLUSTER* data;
	current = head;

	int i=num_elements;
	
	int max_print_clusters = 30;
	
	if (current == NULL)	{
		printf("No clusters in the list.\n");
		return;
	}

	printf("CLUSTERS -->\n");
	while ( current != NULL )	{
		data = current->data;
		printf("  Cluster %i with %i nodes\n",i,data->num_nodes);
		i--;
		current=current->next;
		if (i < (num_elements - max_print_clusters))	{
			printf("Warning! Will only display first %i clusters.\n\n",max_print_clusters);
			return;
		}
	}
	printf("\n");
}
