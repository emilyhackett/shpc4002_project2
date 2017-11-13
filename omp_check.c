#include <stdio.h>
#include <omp.h>

int main()
{
	#pragma omp parallel num_threads(8)
	{
		int id = omp_get_thread_num();
		int num_threads = omp_get_num_threads();
		
		printf("Hello from thread %i of %i\n",id,num_threads);
	}
	
	return 0;
}
