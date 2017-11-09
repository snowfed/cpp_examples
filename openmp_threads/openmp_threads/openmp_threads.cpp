/*
* Simple OpenMP Test
* Compilation:
* g++ -std=c++11 -fopenmp -O3 openmp_threads.cpp -o openmp_threads
*/

#include <cstdlib>
#include <iostream>
#include <string>
#include <omp.h>

const int nthreads = 10;

int main()
{
	omp_set_num_threads(nthreads);

	const int N = (nthreads << 1);
#pragma omp parallel for
	for (int i = 0; i < N; ++i) {
		std::string s = "Iteration " + std::to_string(i) + ": Thread " + std::to_string(omp_get_thread_num()) + ".\n";
		std::cout << s;
	}

	return EXIT_SUCCESS;
}

