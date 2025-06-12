/*
*  This file is part of Christian's OpenMP software lab 
*
*  Copyright (C) 2016 by Christian Terboven <terboven@itc.rwth-aachen.de>
*  Copyright (C) 2016 by Jonas Hahnfeld <hahnfeld@itc.rwth-aachen.de>
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with this program; if not, write to the Free Software
*  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
*
*/

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <sys/time.h>

#include <iostream>
#include <algorithm>

#include <cstdlib>
#include <cstdio>

#include <cmath>
#include <ctime>
#include <cstring>

#include <omp.h>

// sum of the length of the 2 arrays to merge below which the merge algorithm is serial
const int MERGE_THRESHOLD_L = 100000;

/**
  * helper routine: check if array is sorted correctly
  */
bool isSorted(int ref[], int data[], const size_t size){
	std::sort(ref, ref + size);
	for (size_t idx = 0; idx < size; ++idx){
		if (ref[idx] != data[idx]) {
			return false;
		}
	}
	return true;
}


/**
  * sequential merge step (straight-forward implementation)
  */
void MsMergeSequential(int *out, int *in, long begin1, long end1, long begin2, long end2, long outBegin) {
	long left = begin1;
	long right = begin2;

	long idx = outBegin;

	while (left < end1 && right < end2) {
		if (in[left] <= in[right]) {
			out[idx] = in[left];
			left++;
		} else {
			out[idx] = in[right];
			right++;
		}
		idx++;
	}

	while (left < end1) {
		out[idx] = in[left];
		left++, idx++;
	}

	while (right < end2) {
		out[idx] = in[right];
		right++, idx++;
	}
}


/**
  * sequential MergeSort
  */
// TODO: remember one additional parameter (depth)
// TODO: recursive calls could be taskyfied
// TODO: task synchronization also is required
void MsSequential(int *array, int *tmp, bool inplace, long begin, long end) {
	if (begin < (end - 1)) {
		const long half = (begin + end) / 2;
		MsSequential(array, tmp, !inplace, begin, half);
		MsSequential(array, tmp, !inplace, half, end);
		if (inplace) {
			MsMergeSequential(array, tmp, begin, half, half, end, begin);
		} else {
			MsMergeSequential(tmp, array, begin, half, half, end, begin);
		}
	} else if (!inplace) {
		tmp[begin] = array[begin];
	}
}


/**
  * Serial MergeSort
  */
// TODO: this function should create the parallel region
// TODO: good point to compute a good depth level (cut-off)
void MsSerial(int *array, int *tmp, const size_t size) {

   // TODO: parallel version of MsSequential will receive one more parameter: 'depth' (used as cut-off)
	MsSequential(array, tmp, true, 0, size);
}

/*
 * @param array   increasingly sorted array where value must be searched
 * @param begin   search from begin index
 * @param end 	  search before end index
 * @param value   value to search
 * @returns the index where the value is located into the array. However, there are 3 special cases. Please, read below.
 * If array contains value, returns the index of the value (or the index of one of its appearances).
 * If array[begin] > value returns begin - 1
 * If array[end-1] < value returns end
 */
inline int binarySearchIndex(int* array, int begin, int end, int value) {
	int currLeft = 0;
	int currRight = end - begin;

	int index = (end - begin) / 2;
	int currValue = array[begin+index];

	if ( array[begin] > value ) {
		return begin-1;
	}
	if ( array[end-1] < value ) {
		return end;
	}

	while ( currValue != value ) {
		if ( currValue == value ) {
			break;
		}
		else if ( currValue < value ) {
			// se non esiste currMedian, trovo l'intervallo di 2 valori in cui sarebbe compreso
			if ( value < array[begin + index + 1] ) {
				break;
			}
			// looking after
			currLeft = index+1;
			index = index + ((currRight - index) / 2);
		} else {
			// se non esiste currMedian, trovo l'intervallo di 2 valori in cui sarebbe compreso
			if ( value > array[begin + index - 1] ) {
				index = index-1;
				break;
			}
			// looking before
			currRight = index;
			index /= 2;
		}
		currValue = array[begin + index];
	}
	return begin + index;
}

/*
 @brief merges the 2 given array parallelizing the merge in numThreads threads through binary search splitting
 
 The merge is parallelized in the following way:
 The array `[begin1:end1]` is diveded in numThreads subarrays of equal dimension. 
 They are separated by the median[i] values (the name median is meaningful only when numThreads=2) 
 The second array is separated by searching for this medians. 
 
 Some medians will be before `[begin2:end2]`, because their value is less (or equal) than in[begin2]. 
 Likewise, some medians will be after `[begin2:end2]` because their value is greater (or equal) than in[begin2].
 The other medians will fall inside `[begin2:end2]` partitioning it in k subarrays (numThreads >= k).

 k threads will merge the now coupled subarrays (i.e. `[median[i]:median[i+1]]` will be merged with `[binarySearch(median[i]):beinarySearch(median[i+1])]`).
 2 other threads will directly copy the subarrays of '[begin1:end1]` whose median falls outside the array `[begin1:end1]`
 */
void MsMergeParallel(int *out, int *in, long begin1, long end1, long begin2, long end2, long outBegin, int numThreads) {

	if ( (end1 - begin1 + 1) + (end2 - begin2 + 1) <= MERGE_THRESHOLD_L  ) {
		MsMergeSequential(out, in, begin1, end1, begin2, end2, outBegin);
		return;
	}

	#pragma omp parallel for 
	for ( int i=0; i<numThreads; i++ ) {

	}

	int medianIndex, currMedian, index;
	int* begins1 	= 	(int*) malloc(numThreads * sizeof(int));
	int* ends1 		= 	(int*) malloc(numThreads * sizeof(int));
	int* begins2 	= 	(int*) malloc(numThreads * sizeof(int));
	int* ends2 		= 	(int*) malloc(numThreads * sizeof(int));
	int* outBegins 	= 	(int*) malloc(numThreads * sizeof(int));

	int numSectionsBefore = 0, numSectionsAfter = 0;
	begins1[0] = begin1;
	ends1[numThreads-1] = end1;
	begins2[0] = begin2;
	ends2[numThreads-1] = end2;
	outBegins[0] = outBegin;

	/* dividing array [begin1:end1] in constant sections and [begin2:end2] with binary search */
	float array1SectionSize = ((float) end1-begin1)/numThreads;
	#pragma omp parallel for reduction(+: numSectionsBefore, numSectionsAfter)
	for ( int i=0; i<numThreads-1; i++ ) {

		medianIndex = begin1 + (i+1)*array1SectionSize;
		currMedian = in[medianIndex];
		index = binarySearchIndex(in, begin2, end2, currMedian);

		// dividing first array
		ends1[i] = medianIndex;
		begins1[i+1] = medianIndex;

		// handling edge case in which currMedian > in[end2-1]
		if ( index == end2 ) {
			numSectionsAfter += 1;
			index == end2 - 1;
		}

		// dividing second array
		ends2[i] = index+1;
		begins2[i+1] = index+1;

		numSectionsBefore += (index == begin2);
	}

	/* calculating the output indexes */
	for ( int i=1; i<numThreads; i++ ) {
		outBegins[i] = outBegins[i-1] + (ends2[i-1] - begins2[i-1]) + (ends1[i-1] - begins1[i-1]);
	}

	int firstSectionSize, secondSectionSize, finalSectionSize;
	#pragma omp parallel num_threads(numThreads)
	{
		#pragma omp single
		{
			/* if array [begin1:end1] has some sections with values less than in[begin2], those sections are copied directly in the output array*/
			if ( numSectionsBefore > 0 ) {
				firstSectionSize = begins1[numSectionsBefore-1] - begins1[0];
				#pragma omp task firstprivate(firstSectionSize)
				{
					std::memcpy(&(out[outBegins[0]]), &(in[begins1[0]]), firstSectionSize);
				}
			} else {
				firstSectionSize = 0;
			}

			/* parallel merge of couples of sections */
			for ( int i=numSectionsBefore; i<numThreads-numSectionsAfter; i++ ) {
				#pragma omp task 
				{
					MsMergeSequential(out, in, 
									begins1[i], ends1[i], 
									begins2[i], ends2[i], 
									outBegins[i]);
				}
			}

			/* if array [begin1:end1] has some sections with values greater than in[end2-1], those sections are copied directly in the output array*/
			if ( numSectionsAfter > 0 ) {
				#pragma omp task 
				{
					finalSectionSize = ends1[numThreads-1]-begins1[numThreads-numSectionsAfter];
					std::memcpy(&(out[outBegins[numThreads-numSectionsAfter]]), 
								&(in[begins1[numThreads-numSectionsAfter]]),
								finalSectionSize);
				}
			}
		}
			
	}

	free(begins1);
	free(begins2);
	free(ends1);
	free(ends2);
	free(outBegins);

}

/*
 @brief merges the 2 given array parallelizing the merge in 2 threads. Those thread will recursively call this funciton until a certain threshold on the arrays length

 Merges the 2 given array parallelizing the merge in 2 threads. 
 Those thread will recursively call this funciton until a certain threshold on the arrays length.
 Tries balancing the array division by dividing the first array in 2 equal halves and using the binary search to divide the other one.
 Does the same inverting the 2 arrays. Then, the most balanced division is choosen.
*/
void MsMergeParallel2(int *out, int *in, long begin1, long end1, long begin2, long end2, long outBegin) {
	if ( (end1 - begin1 + 1) + (end2 - begin2 + 1) <= MERGE_THRESHOLD_L  ) {
		MsMergeSequential(out, in, begin1, end1, begin2, end2, outBegin);
		return;
	}

	// division1: median of array1, searching for it in array 2
	int medianIndex1 = begin1 + (end1 - begin1)/2;
	int medianValue1 = in[medianIndex1];
	int foundIndex1 = binarySearchIndex(in, begin2, end2, medianValue1);

	// division2: median in array2, searching for it in array 1
	int medianIndex2 = begin2 + (end2 - begin2)/2;
	int medianValue2 = in[medianIndex2];
	int foundIndex2 = binarySearchIndex(in, begin1, end1, medianValue2);

	// choosing the most balanced division
	// following ratio express how far the index is from the center of the array
	int division1EquityRatio = abs(0.5 - (foundIndex1-begin2)/(end2-begin2));
	int division2EquityRation = abs(0.5 - (foundIndex2-begin1)/(end2-begin2));
	if ( division1EquityRatio < division2EquityRation ) {
		// choosing the first division

		#pragma omp parallel num_threads(2) 
		{
			#pragma omp single 
			{
				if ( foundIndex1 > begin2 - 1 ) {
					#pragma omp task
					{
						MsMergeParallel2(out, in, begin1, medianIndex1, begin2, foundIndex1+1, outBegin);
					}
				} else {
					medianIndex1 = begin1;
				}
				#pragma omp task 
				{
					MsMergeParallel2(out, in, medianIndex1, end1, foundIndex1+1, end2, outBegin + (medianIndex1-begin1)+(foundIndex1+1-begin2));
				}
			}
		}
	} else {
		// choosing the second division
		#pragma omp parallel num_threads(2)
		{
			#pragma omp single
			{
				if ( foundIndex2 > begin1 - 1 ) {
					#pragma omp task
					{
						MsMergeParallel2(out, in, begin1, foundIndex2+1, begin2, medianIndex2, outBegin);
					}
				} else {
					medianIndex2 = begin2;
				}
				#pragma omp task 
				{
					MsMergeParallel2(out, in, foundIndex2+1, end1, medianIndex2, end2, outBegin + (medianIndex2-begin2)+(foundIndex2+1-begin1));
				}
			}
		}
	}


}

void MsParallel(int *array, int *tmp, bool inplace, long begin, long end, int depth) {

	if (begin < (end - 1)) {
		const long half = (begin + end) / 2;
		if ( depth > 0 ) {
			#pragma omp parallel num_threads(2) shared(array, tmp) firstprivate(inplace, begin, half, depth)
			{
				#pragma omp single
				{
					#pragma omp task 
					{
						MsParallel(array, tmp, !inplace, begin, half, depth-1);
					}
					#pragma omp task 
					{
						MsParallel(array, tmp, !inplace, half, end, depth-1);
					}
				}
			}

			if (inplace) {
				MsMergeParallel2(array, tmp, begin, half, half, end, begin);
				//MsMergeParallel(array, tmp, begin, half, half, end, begin, pow(2, depth));
			} else {
				MsMergeParallel2(tmp, array, begin, half, half, end, begin);
				//MsMergeParallel(tmp, array, begin, half, half, end, begin, pow(2, depth));
			}

		} else {
			MsSequential(array, tmp, !inplace, begin, half);
			MsSequential(array, tmp, !inplace, half, end);

			if (inplace) {
				MsMergeSequential(array, tmp, begin, half, half, end, begin);
			} else {
				MsMergeSequential(tmp, array, begin, half, half, end, begin);
			}

		}


	} else if (!inplace) {
		tmp[begin] = array[begin];
	}
}

void MsParallelAlgorithm(int *array, int *tmp, const size_t size) {
	/* calculating the maximum depth */
	int numThreads = omp_get_num_procs();
	int depth = floor(log2(numThreads));
	if ( size < 17654 ) {
		MsSerial(array, tmp, size);
	}
	else {
		omp_set_max_active_levels(256);
		if ( size < 100000000 ) {
			omp_set_dynamic(true);
		}

		MsParallel(array, tmp, true, 0, size, depth);
	}
}

/** 
  * @brief program entry point
  */
int main(int argc, char* argv[]) {
	// variables to measure the elapsed time
	struct timeval t1, t2;
	double etime;
	double serialTime, parallelTime;

	// expect one command line arguments: array size
	if (argc != 2) {
		printf("Usage: MergeSort.exe <array size> \n");
		printf("\n");
		return EXIT_FAILURE;
	}
	else {
		const size_t stSize = strtol(argv[1], NULL, 10);
		int *data = (int*) malloc(stSize * sizeof(int));
		int *data1 = (int*) malloc(stSize * sizeof(int));
		int *tmp = (int*) malloc(stSize * sizeof(int));
		int *ref = (int*) malloc(stSize * sizeof(int));

		printf("Initialization...\n");

		//srand(95);
		int elem;
		for (size_t idx = 0; idx < stSize; ++idx){
			elem = (int) (stSize * (double(rand()) / RAND_MAX));
			data[idx] = elem;
			data1[idx] = elem;
		}
		std::copy(data, data + stSize, ref);

		double dSize = (stSize * sizeof(int)) / 1024 / 1024;
		printf("Sorting %zu elements of type int (%f MiB)...\n", stSize, dSize);
		
		/* -------------------- SERIAL ALGORITHM -------------------------*/
		/*
		printf("Serial algorithm:\n");
		gettimeofday(&t1, NULL);
		MsSerial(data, tmp, stSize);
		gettimeofday(&t2, NULL);

		serialTime = (t2.tv_sec - t1.tv_sec) * 1000 + (t2.tv_usec - t1.tv_usec) / 1000;
		serialTime = serialTime / 1000;

		printf("done, took %f sec. Verification...", serialTime);
		if (isSorted(ref, data, stSize)) {
			printf(" successful.\n");
		}
		else {
			printf(" FAILED.\n");
		}
		*/

		/* -------------------- PARALLEL ALGORITHM -------------------------*/
		printf("Parallel algorithm: \n");
		gettimeofday(&t1, NULL);
		MsParallelAlgorithm(data1, tmp, stSize);
		gettimeofday(&t2, NULL);

		parallelTime = (t2.tv_sec - t1.tv_sec) * 1000 + (t2.tv_usec - t1.tv_usec) / 1000;
		parallelTime = parallelTime / 1000;

		printf("done, took %f sec. Verification...", parallelTime);
		if (isSorted(ref, data1, stSize)) {
			printf(" successful.\n");
		}
		else {
			printf(" FAILED.\n");
			for ( int i=0; i<stSize; i++) {
				if ( data1[i] != ref[i] ) {
					printf(">> %d\n", i);
					break;
				}
			}
		}

		//printf("speedup is %lf\n", serialTime/parallelTime);

		free(data);
		free(data1);
		free(tmp);
		free(ref);
	}



	return EXIT_SUCCESS;
}
