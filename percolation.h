/* Header file for percolation project in C for SHPC4002 
 *	
 *	Emily Hackett, 21489688
 *	Computational Physics Honours
 */

/* External libraries used */
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <omp.h>

typedef struct CLUSTER	{
	int num_nodes;	
	int* cols_reached;
	int* rows_reached;
} CLUSTER;

typedef struct NODE	{
	CLUSTER* data;
	struct NODE* next;
} NODE;

/* Command line argument variables */
extern float occ_prob;
extern char perc_type;
extern int span_type;
extern int N;
extern int NUM_THREADS;

/* Globals */
extern int** sites;	/* Lattice sites */
extern int** hbonds;	/* Horizontal bonds */
extern int** vbonds;	/* Vertical bonds */

extern int* rows_reached;
extern int* cols_reached;

/* Miscellaneous helper functions */
extern int 	wrap(int val, int max);

/* Functions associated with initialising the 2D lattice and bonds */
extern int**	allocate_lattice(int M);
extern void	initialise_lattice(int** sites, int N, float occ_prob);
extern void	display_lattice(int** sites, int** hbonds, int** vbonds, int N);

/* Functions associated with the depth first search */
extern CLUSTER*	initialise_cluster(int N, int row, int col);
extern CLUSTER*	depth_first_search(int** sites, int** hbonds, int** vbonds, int N, int* chunk, int row, int col, CLUSTER* tmp);

/* Check for spanning clusters */
extern int	check_spanning(CLUSTER* tmp, int N, int span_type);

/* Functions associated with linked list */
extern NODE* 	push(struct NODE* head, struct CLUSTER* data);
extern NODE*	pop(struct NODE* head, struct CLUSTER* data);
extern void	display_list(struct NODE* head, int num_elements);
extern int	traverse_list(struct NODE* head, int N, int span_type, int* spanning, int* max_nu_nodes);
