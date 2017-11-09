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

