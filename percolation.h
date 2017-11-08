/* Header file for percolation project in C for SHPC4002 */

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

/* Functions associated with initialising the 2D lattice and bonds */
extern int**	allocate_lattice(int M);
extern void	initialise_lattice(int** lattice, int N, float occ_prob);
extern void	display_lattice(int** lattice, int** hbonds, int** vbonds, int N);
