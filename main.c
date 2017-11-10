#include "percolation.h"

int main(int argc, char *argv[])
{
	printf("SHPC4002, PROJECT 2: Emily Hackett, 21489688\n\n");

	int N = 16;	/* Define default lattice size */
	int NUM_THREADS = 2;	/* Define default num threads */
	
	struct timeval start, end;	/* Allocate start and end time vals */
	gettimeofday(&start, NULL);	/* Begin timing */

	/* CHECK COMMAND LINE ARGUMENTS:
	 * If 3 command line arguments read, assume they are:
	 * 	argv[1]: Occupancy probability, btwn 0 and 1	 
	 * 	argv[2]: Percolation type, 's' or 'b'
	 * 	argv[3]: Spanning type, 0 (col), 1 (row), 2 (both)
	 * If 4 command line arguments read, assume:
	 * 	argv[4]: Lattice size (must be power of 2)
	 */
	float occ_prob;
	char perc_type;
	int span_type;
	
	if (argc > 3)	{
		occ_prob=atof(argv[1]);	/* Occupancy probability */
		if (occ_prob<0||occ_prob>1)	{
			fprintf(stderr,"ERROR: Occupancy probability must lie within [0,1]\n");
			exit(EXIT_FAILURE);
		}
		printf("* Occupancy probability = %5.4f\n",occ_prob);
		
		perc_type=*argv[2];	/* Percolation type, site or bond */
		if (perc_type=='s'||perc_type=='b')	{
			printf("* Percolation type = %c\n",perc_type);
		}
		else	{
			fprintf(stderr,"ERROR: Percolation type must be 's' (site) or 'b' (bond)\n");
			exit(EXIT_FAILURE);
		}
		
		span_type=atoi(argv[3]);	/* Spanning type */
		if (span_type==0||span_type==1||span_type==2)	{
			printf("* Spanning type = %i\n",span_type);
		}
		else	{
			fprintf(stderr,"ERROR: Spanning type must be 0 (col), 1 (row) or 2 (both)\n");
			exit(EXIT_FAILURE);
		}
		
		if (argc > 4)	{
			N = atoi(argv[4]);	/* Lattice size */
			if (( N & (N - 1)) != 0 )	{	/* Check if power of 2 */
				fprintf(stderr,"ERROR: Lattice size must be a power of 2\n");
				exit(EXIT_FAILURE);
			}
		}
		printf("* Lattice size = %i\n",N);
	}
	else	{
		fprintf(stderr, "Usage:	./main occ_prob perc_type span_type [size]\n");
		exit(EXIT_FAILURE);
	}
	
	int** sites = allocate_lattice(N);	/* Allocate the site lattice */
	int** hbonds = allocate_lattice(N);	/* Allocate horizontal bonds */
	int** vbonds = allocate_lattice(N);	/* Allocate vertical bonds */
	
	/* Initialise the lattice and bonds based on percolation type */
	if (perc_type == 's') {
		initialise_lattice(sites,N,occ_prob);	/* Fill lattice according to occ_prob */
		initialise_lattice(hbonds,N,-1);	/* Set all horizontal bonds equal to 1 */
		initialise_lattice(vbonds,N,-1);	/* Set all vertical bonds equal to 1 */
	}
	if (perc_type == 'b') {
		initialise_lattice(sites,N,-1);		/* Set all sites to occupied */
		initialise_lattice(hbonds,N,occ_prob);	/* Fill bonds according to occ_prob */
		initialise_lattice(vbonds,N,occ_prob);
	}

	/* Print the lattice */
	display_lattice(sites,hbonds,vbonds,N);
	
	/* Conduct the depth-first search */
	int i,j;
	
	int max_num_nodes = 0;	/* Tracks the maximum number of nodes */
	int spanning = 0;	/* Checks whether spanning or not */
	
	/* Have an array that holds the head of each cluster list for each thread.
	 * NOTE! This needs to be in the shared memory, hence why declared here */
	NODE** head_list = malloc(NUM_THREADS * sizeof(NODE*));	/* Initialise the start of cluster linked list */
	for (i = 0; i < NUM_THREADS; i++)	{
		head_list[i]= NULL;
	}
	
	int chunk_size = N/NUM_THREADS;	/* Split the cluster depending on num_threads */
	/* NOTE! Can check for not evenly divided, and possibly give the last thread/master thread 
	 * more rows to handle later
	 */

	int num_threads_running = NUM_THREADS;

	/* *** START OF PARALLEL *** */
	printf("* Beginning parallel region with %i threads ...\n",NUM_THREADS);
	#pragma omp parallel num_threads(NUM_THREADS) shared(head_list,num_threads_running,N) private(i,j)
	{
		int num_threads = omp_get_num_threads();
		int id = omp_get_thread_num();

		int* chunk = malloc(2 * sizeof(int));
		chunk[0] = id * chunk_size;	/* Start row */
		chunk[1] = chunk[0] + chunk_size - 1;	/* End row */

		printf("Thread %i of %i taking lattice rows %i -> %i\n",id,num_threads,chunk[0],chunk[1]);

		/* Loop over all lattice points and conduct a depth first search */		
		for (i = chunk[0]; i < chunk[1]; i++)	{	/* Check the rows relevant to this thread */
			for (j = 0; j < N; j++)	{
			
				if (sites[i][j] == 1) 	{	
	
					/* Allocate new cluster */
					CLUSTER* tmp = initialise_cluster(N,i,j);

					/* If the site is occupied, conduct depth_first_search */
					tmp = depth_first_search(sites,hbonds,vbonds,N,chunk,i,j,tmp);
				
					head_list[id] = push(head_list[id], tmp);	/* Push this cluster onto the list */
					printf("- Thread %i found cluster at [%i,%i] with %i nodes\n",id,i,j,tmp->num_nodes);
						
				}
			}
		}

		/* BEGIN TO COMBINE CLUSTERS */
		int divisor = 2;

		#pragma omp wait
		while (num_threads_running > 1)	{	/* Consolidate into master thread */
			if (id % divisor == 0)	{	/* Take half of the running threads */
				#pragma omp master
					num_threads_running = num_threads_running/2;

				int neighbour_id = id + divisor/2;
				printf("CHECK -> num_threads_running = %i, id = %i, neighbour_id = %i\n",num_threads_running,id,neighbour_id);				

				//head_list[id] = merge_cluster_lists(head_list[id],head_list[neighbour_id],N);
			}
		}

	}
	/* *** END OF PARALLEL *** */
	printf("* ... end of parallel region.\n");	
	
	printf("RESULTS:\n");	/* Determine if spanning, max nodes */

	display_list(head_list[0]);	/* Display the found cluster information */

	traverse_list(head_list[0],N,span_type,&spanning,&max_num_nodes);	/* Search clusters for span/max */

	printf("Maximum number of nodes in a cluster is %i.\n",max_num_nodes);
	
	if (spanning == 1)	{
		printf("Spanning cluster of type %i exists.\n",span_type);
	}
	else	{
		printf("No spanning cluster of type %i exists.\n",span_type);
	}

	gettimeofday(&end, NULL);	/* End the timer */
	double time_spent = ((end.tv_sec - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;	
	printf("\nTOTAL TIME: %10.6f\n",time_spent);
	
	return 0;
	
}
