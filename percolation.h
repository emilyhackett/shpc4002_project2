/* Header file for percolation project in C for SHPC4002 
 *	
 *	Emily Hackett, 21489688
 *	Computational Physics Honours
 */

/* External libraries used */
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <string.h>

#include <omp.h>
#include <mpi.h>

typedef struct CLUSTER	{
	int num_nodes;	
	int* cols_reached;
	int* rows_reached;
	int top_row_idx;
	int bottom_row_idx;
	int* top_bounds;
	int* bottom_bounds;
	int merged;
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
extern int**	allocate_lattice(int N, int chunk);
extern void	initialise_lattice(int** sites, int N, float occ_prob);
extern void	display_lattice(int** sites, int** hbonds, int** vbonds, int N, int chunk);

/* Functions associated with the depth first search */
extern CLUSTER*	initialise_cluster(int N, int row, int col);
extern CLUSTER*	depth_first_search(int** sites, int** hbonds, int** vbonds, int N, int* chunk, int row, int col, CLUSTER* tmp, int rank, int id);

/* Combination of split clusters */
extern NODE*	merge_cluster_lists(NODE* bottom_seg, NODE* top_seg,int N, int end_idx_A, int* num_clusters_in_list);
extern int	check_bounds_crossover(int* bounds_1, int* bounds_2, int N);
extern CLUSTER*	merge(CLUSTER* current, NODE* head, CLUSTER* new_cluster, int N);

/* Check for spanning clusters */
extern int	check_spanning(CLUSTER* tmp, int N, int span_type);

/* Functions associated with linked list */
extern NODE* 	push(struct NODE* head, struct CLUSTER* data);
extern NODE*	pop(struct NODE* head, struct CLUSTER* data);
extern void	display_list(struct NODE* heads);
extern int	traverse_list(struct NODE* head, int N, int span_type, int* spanning, int* max_nu_nodes);
extern void	linkedlist_to_array(NODE* head, int num_elements, int N, int* num_nodes, int* top_row_idx, int* bottom_row_idx, int** top_bounds, int** bottom_bounds, int** cols_spanned, int** rows_spanned);
extern void	linkedlist_allocate(int num_elements, int N, int* num_nodes, int* top_row_idx, int* bottom_row_idx, int** top_bounds, int** bottom_bounds, int** cols_spanned, int** rows_spanned);

