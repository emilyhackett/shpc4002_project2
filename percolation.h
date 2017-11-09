/* Header file for percolation project in C for SHPC4002 
 *	
 *	Emily Hackett, 21489688
 *	Computational Physics Honours
 */

/* External libraries used */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/* Command line argument variables */
extern float occ_prob;
extern char perc_type;
extern int span_type;
extern int N;

/* Globals */
extern int** sites;	/* Lattice sites */
extern int** hbonds;	/* Horizontal bonds */
extern int** vbonds;	/* Vertical bonds */

extern int* rows_reached;
extern int* cols_reached;

/* Extra functions */
extern int 	wrap(int val, int max);

/* Functions associated with initialising the 2D lattice and bonds */
extern int**	allocate_lattice(int M);
extern void	initialise_lattice(int** sites, int N, float occ_prob);
extern void	display_lattice(int** sites, int** hbonds, int** vbonds, int N);

/* Functions associated with the depth first search */
extern int	depth_first_search(int** sites, int** hbonds, int** vbonds, int N, int row, int col);

/* Check for spanning clusters */
extern int	check_spanning(int* rows_reached, int* cols_reached, int N, int span_type);
