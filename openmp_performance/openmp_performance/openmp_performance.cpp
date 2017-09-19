/*
 * Simple OpenMP Test
 * Compilation:
 * g++ -std=c++11 -fopenmp -O3 openmp_performance.cpp -o openmp_performance
 */
#include <cstddef>
#include <cstdlib>
#include <chrono>
#include <iostream>
#include <limits>
#include <algorithm>
#include <omp.h>

using std::cerr;
using std::cout;
using std::cin;
using std::endl;
using std::size_t;

typedef double real_t;

const size_t nops0 = 50000000; // initial number of operations
const size_t nms_stop = 10000; // max execution time in ms

real_t arithm_op_bulk (real_t input, size_t N, bool use_openmp);

int main ()
{
	const size_t multiplier = 2;
	const size_t maxN = std::numeric_limits<size_t>::max() / multiplier;
	double nms1 = 0, nms2 = 1;
	size_t N = 0;
	real_t val1 = 0, val2 = 0;
	for (N = nops0; N < maxN; N *= 2) {
		auto start = std::chrono::high_resolution_clock::now();
		val1 = arithm_op_bulk(val1, N, false);
		std::chrono::duration<double, std::milli> sequential_time =
			std::chrono::high_resolution_clock::now() - start;

		start = std::chrono::high_resolution_clock::now();
		val2 = arithm_op_bulk(val2, N, true);
		std::chrono::duration<double, std::milli> parallel_time = 
			std::chrono::high_resolution_clock::now() - start;

		if (std::abs((val2 - val1) / std::max(val1, val2)) > 0.0001) {
			cerr << "ERROR (Parallel != Sequential)" << endl;
			return EXIT_FAILURE;
		}

		nms1 = sequential_time.count();
		nms2 = parallel_time.count();
		if (nms1 > nms_stop || nms2 > nms_stop)
			break;
		(cout << ".").flush();
	}
	cout << "\rFLOPS (sequential): " << 1000 * real_t(N) / nms1 << "." << endl;
	cout << "FLOPS (parallel):   " << 1000 * real_t(N) / nms2 << "." << endl;
	cout << "OPENMP SPEEDUP:     " << (nms2 == 0 ? 1 : real_t(nms1) / nms2) << "." << endl;
	cout << "omp_get_num_procs() returns " << omp_get_num_procs() << "." << endl;
	return EXIT_SUCCESS;
}

real_t arithm_op (real_t input)
{
	return input * 0.5 + 1;
}

real_t arithm_op_bulk_helper (real_t value, size_t N)
{
	for (size_t i = 0; i < N; ++i)
		value = arithm_op(value);
	return value;
}

real_t arithm_op_bulk (real_t input, size_t N, bool use_openmp)
{
	int nprocs = omp_get_num_procs();
	size_t K = N / nprocs;
	real_t total = 0;
	if (use_openmp) {
#pragma omp parallel for reduction(+:total)
		for (int ithread = 0; ithread < nprocs; ++ithread)
			total += arithm_op_bulk_helper(input+ithread, K);
	} else {
		for (int ithread = 0; ithread < nprocs; ++ithread)
			total += arithm_op_bulk_helper(input+ithread, K);
	}
	return total;
}

