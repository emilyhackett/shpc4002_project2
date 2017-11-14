#include "percolation.h"

int main(int argc, char *argv[])
{
	int i,j;
	
	float occ_prob;
	char perc_type;
	int span_type;

	int rank,size;

	MPI_Init( &argc, &argv);
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );
        MPI_Comm_size( MPI_COMM_WORLD, &size );
//	printf("I am %i of %i\n",rank,size);

        MPI_Status stats[3];
        MPI_Request reqs[3];

	int N = 4;	/* Define default lattice size */
	int NUM_THREADS = 1;	/* Define default num threads */
	
	struct timeval start,end;	/* Allocate start and end time vals */	
	
	if (rank == 0)	{
		printf("SHPC4002, PROJECT 2: Emily Hackett, 21489688\n\n");
		printf("	Running on %i nodes in the cluster.\n",size);
		gettimeofday(&start, NULL);	/* Begin timing */
	}
		
	/* CHECK COMMAND LINE ARGUMENTS:
	 * If 3 command line arguments read, assume they are:
	 * 	argv[1]: Occupancy probability, btwn 0 and 1	 
	 * 	argv[2]: Percolation type, 's' or 'b'
	 * 	argv[3]: Spanning type, 0 (col), 1 (row), 2 (both)
	 * If 4 command line arguments read, assume:
	 * 	argv[4]: Lattice size (must be power of 2)
	 */

	if (argc > 3)	{
		occ_prob=atof(argv[1]);	/* Occupancy probability */
		if (occ_prob<0||occ_prob>1)	{
			fprintf(stderr,"ERROR: Occupancy probability must lie within [0,1]\n");
			exit(EXIT_FAILURE);
		}
		
		perc_type=*argv[2];	/* Percolation type, site or bond */
		if (perc_type!='s' && perc_type!='b')	{
			fprintf(stderr,"ERROR: Percolation type must be 's' (site) or 'b' (bond)\n");
			exit(EXIT_FAILURE);
		}
		
		span_type=atoi(argv[3]);	/* Spanning type */
		if (span_type!=0 && span_type!=1 && span_type!=2)	{
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
	}
	else	{
		fprintf(stderr, "Usage:	./main occ_prob perc_type span_type [size]\n");
		exit(EXIT_FAILURE);
	}

	int MPI_chunk_size = N/size;
	if (N % size != 0)	{
		fprintf(stderr,"ERROR: Number of MPI processes must divide into lattice.\n");
		exit(EXIT_FAILURE);
	}
	int OMP_chunk_size = MPI_chunk_size/NUM_THREADS;
	if (MPI_chunk_size % NUM_THREADS != 0)	{
		fprintf(stderr,"ERROR: Number of OMP threads must divided into MPI chunk.\n");
		exit(EXIT_FAILURE);
	}

	/* Print these and command line argument results */
	if (rank == 0)	{
		printf("	Lattice split into MPI chunks of size %i.\n",MPI_chunk_size);
		printf("	Lattice split into OMP chunks of size %i.\n\n",OMP_chunk_size);
                
		printf("* Occupancy probability = %5.4f\n",occ_prob);
                printf("* Percolation type = %c\n",perc_type);
                printf("* Spanning type = %i\n",span_type);
                printf("* Percolation type = %c\n",perc_type);
                printf("* Lattice size = %i\n\n",N);
	}

	MPI_Barrier(MPI_COMM_WORLD);	/* Wait to check command line arguments */

	int** sites = allocate_lattice(N,MPI_chunk_size+2);	/* Pointer for the site lattice */
	int** hbonds = allocate_lattice(N,MPI_chunk_size+2);	/* Pointer for horizontal bonds */
	int** vbonds = allocate_lattice(N,MPI_chunk_size+2);	/* Pointer for vertical bonds */

	if (rank == 0)	{	
		
		int** main_sites = allocate_lattice(N,N);	/* Allocate the site lattice */
		int** main_hbonds = allocate_lattice(N,N);	/* Allocate horizontal bonds */
		int** main_vbonds = allocate_lattice(N,N);	/* Allocate vertical bonds */
	
		/* Initialise the lattice and bonds based on percolation type */
		if (perc_type == 's') {
			initialise_lattice(main_sites,N,occ_prob);	/* Fill lattice according to occ_prob */
			initialise_lattice(main_hbonds,N,-1);	/* Set all horizontal bonds equal to 1 */
			initialise_lattice(main_vbonds,N,-1);	/* Set all vertical bonds equal to 1 */
		}
		if (perc_type == 'b') {
			initialise_lattice(main_sites,N,-1);		/* Set all sites to occupied */
			initialise_lattice(main_hbonds,N,occ_prob);	/* Fill bonds according to occ_prob */
			initialise_lattice(main_vbonds,N,occ_prob);
		}
		
		/* Print the lattice */
		display_lattice(main_sites,main_hbonds,main_vbonds,N,N);

		/* Loop over all other processes and send parts of the lattice to them */
		for (i = 1; i < size; i++ )	{
			int start_idx = MPI_chunk_size*i - 1;

			/* If not at the end of the lattice, send through normally (taking extra top/bottom row */
			if (i != size-1)	{
                        	printf("%i: Sending rows [%i -> %i] to node %i.\n",rank,start_idx,start_idx+MPI_chunk_size+1,i);

				MPI_Isend(&main_sites[start_idx][0], N*(MPI_chunk_size+2), MPI_INT, i, 1, MPI_COMM_WORLD, &reqs[0]);
				MPI_Isend(&main_hbonds[start_idx][0], N*(MPI_chunk_size+2), MPI_INT, i, 2, MPI_COMM_WORLD, &reqs[1]);
				MPI_Isend(&main_vbonds[start_idx][0], N*(MPI_chunk_size+2), MPI_INT, i, 3, MPI_COMM_WORLD, &reqs[2]);
			}
			else	{
				int** end_sites = allocate_lattice(N,MPI_chunk_size+2);
				int** end_hbonds = allocate_lattice(N,MPI_chunk_size+2);
				int** end_vbonds = allocate_lattice(N,MPI_chunk_size+2);
			
				for (j = 0; j < MPI_chunk_size+2; j++)	{
					/* Copy across data from main to end sites */
					int k;
					if (j != MPI_chunk_size+1) 	{
						for (k = 0; k < N; k++)	{
							end_sites[j][k] = main_sites[start_idx+j][k];
							end_hbonds[j][k] = main_hbonds[start_idx+j][k];
							end_vbonds[j][k] = main_vbonds[start_idx+j][k];
						}
					}	
					else {
						for (k = 0; k < N; k++)	{
							end_sites[j][k] = main_sites[0][k];
							end_hbonds[j][k] = main_hbonds[0][k];
							end_vbonds[j][k] = main_vbonds[0][k];
						}
					}
				}

	                        printf("%i: Sending rows [%i -> %i] to node %i.\n",rank,start_idx,0,i);

				MPI_Isend(&end_sites[0][0], N*(MPI_chunk_size+2), MPI_INT, i, 1, MPI_COMM_WORLD, &reqs[0]);
				MPI_Isend(&end_hbonds[0][0], N*(MPI_chunk_size+2), MPI_INT, i, 2, MPI_COMM_WORLD, &reqs[1]);
				MPI_Isend(&end_vbonds[0][0], N*(MPI_chunk_size+2), MPI_INT, i, 3, MPI_COMM_WORLD, &reqs[2]);
			}

		}

		/* Copy over from master lattice data (main_lattice) */
		memcpy(&sites[1],&main_sites[0],N*(MPI_chunk_size+1));
		memcpy(&sites[0],&main_sites[N-1],N);
		memcpy(&hbonds[1],&main_hbonds[0],N*(MPI_chunk_size+1));
		memcpy(&hbonds[0],&main_hbonds[N-1],N);
		memcpy(&vbonds[1],&main_vbonds[0],N*(MPI_chunk_size+1));
		memcpy(&vbonds[0],&main_vbonds[N-1],N);

	}
	else	{
		
		printf("%i: Called MPI_Irecv for lattice ...\n",rank);
		MPI_Irecv(&sites[0][0], N*(MPI_chunk_size+2), MPI_INT, 0, 1, MPI_COMM_WORLD, &reqs[0]);
		MPI_Irecv(&hbonds[0][0], N*(MPI_chunk_size+2), MPI_INT, 0, 2, MPI_COMM_WORLD, &reqs[1]);
		MPI_Irecv(&vbonds[0][0], N*(MPI_chunk_size+2), MPI_INT, 0, 3, MPI_COMM_WORLD, &reqs[2]);
	}

	MPI_Waitall(3,reqs,stats);
	printf("%i: Lattice chunk recieved, beginning DFS.\n",rank);

	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Finalize();

	return 0;

	/* Conduct the depth-first search */	
	int max_num_nodes = 0;	/* Tracks the maximum number of nodes */
	int spanning = 0;	/* Checks whether spanning or not */
	
	/* Have an array that holds the head of each cluster list for each thread.
	 * NOTE! This needs to be in the shared memory, hence why declared here */
	NODE** head_list = malloc(NUM_THREADS * sizeof(NODE*));	/* Initialise the start of cluster linked list */
	for (i = 0; i < NUM_THREADS; i++)	{
		head_list[i]= NULL;
	}
	
	int num_threads_running = NUM_THREADS;

	/* *** START OF PARALLEL *** */
	printf("--- BEGIN PARALLEL REGION ---\n");
	#pragma omp parallel num_threads(NUM_THREADS) shared(head_list,num_threads_running,N) private(i,j)
	{
		int num_threads = omp_get_num_threads();
		int id = omp_get_thread_num();

		int* chunk = malloc(2 * sizeof(int));
		chunk[0] = id * OMP_chunk_size;	/* Start row */
		chunk[1] = chunk[0] + OMP_chunk_size - 1;	/* End row */

		printf("Node %i: Thread %i of %i taking lattice rows %i -> %i\n",rank,id,num_threads,chunk[0],chunk[1]);

		/* Loop over all lattice points and conduct a depth first search */		
		for (i = chunk[0]; i <= chunk[1]; i++)	{	/* Check the rows relevant to this thread */
			for (j = 0; j < N; j++)	{
			
				if (sites[i][j] == 1) 	{	
	
					/* Allocate new cluster */
					CLUSTER* tmp = initialise_cluster(N,i,j);

					/* If the site is occupied, conduct depth_first_search */
					tmp = depth_first_search(sites,hbonds,vbonds,N,chunk,i,j,tmp);
					tmp->top_row_idx = chunk[0];
					tmp->bottom_row_idx = chunk[1];
				
					head_list[id] = push(head_list[id], tmp);	/* Push this cluster onto the list */
				}
			}
		}

		/* BEGIN TO COMBINE CLUSTERS */
		int divisor = 2;
		int bound_idx = chunk[1];

		#pragma omp barrier
		while (num_threads_running > 1)	{	/* Consolidate into master thread */
			if (id % divisor == 0)	{	/* Take half of the running threads */
				#pragma omp master
					num_threads_running = num_threads_running/2;

				int neighbour_id = id + divisor/2;

				head_list[id] = merge_cluster_lists(head_list[id],head_list[neighbour_id], N, bound_idx);
				
				bound_idx = bound_idx + OMP_chunk_size;

				#pragma omp master
					divisor = divisor*2;
			}
		}

	}
	/* *** END OF PARALLEL *** */
	if (rank == 0)	{
		printf("--- END OF PARALLEL REGION ---\n\n");	
	
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
	}	

	MPI_Finalize();

	return 0;
}
