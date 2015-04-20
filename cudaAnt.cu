/**
 *  Project TACO: Parallel ACO algorithm for TSP
 *  15-418 Parallel Algorithms - Final Project
 *  Ivan Wang, Carl Lin
 */

#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime.h>

#include "ants.h"

#define MAX_THREADS 1024

__device__ double aantProduct(int from, int to) {
  //return (pow(phero[from][to], ALPHA) * pow((1.0 / dist[from][to]), BETA));
  return 0;
}

// Randomly select a city based off an array of probabilities (return the index)
__device__ int selectCity(double *start, double *end) {
  return 0;
}

__global__ void constructAntTour(double *tourResults, int *pathResults) {
    __shared__ int tabu[MAX_CITIES]; //TODO: put in register wtf is that
    __shared__ int current_city;
    __shared__ int path[MAX_CITIES];
    __shared__ int num_visited;
    __shared__ double tour_length;
    __shared__ int bestCities[MAX_THREADS];
    __shared__ double cityProb[MAX_THREADS];

 
    if (threadIdx.x == 0) {
      current_city = 0; //TODO: random make it random random
      num_visited = 0;
    }

    const int citiesPerThread = (MAX_CITIES + MAX_THREADS - 1) / MAX_THREADS;
    int startCityIndex = threadIdx.x * citiesPerThread;
    int antId = blockIdx.x;

    double localCityProb[citiesPerThread];

    //initiailize tabu list to zero
    for (int i = startCityIndex; i < startCityIndex + citiesPerThread; i++) {
      tabu[i] = 0;
    }

    //check if we have finished the tour (can you run into dead ends?)
    while (num_visited != MAX_CITIES) {
      //pick next city
      for (int i = startCityIndex; i < startCityIndex + citiesPerThread; i++) {
        if (tabu[i] != 0) {
          localCityProb[i] = 0.0;
          continue;
        }

        localCityProb[i] = aantProduct(current_city, i);
      }

      bestCities[threadIdx.x] = 
         selectCity(&localCityProb[startCityIndex], 
                    &localCityProb[startCityIndex + citiesPerThread]);

      cityProb[threadIdx.x] = localCityProb[bestCities[threadIdx.x]];

      __syncthreads();

      //reduce over cityProb and randomly pick a city based on 
      //those probabilities

      if(threadIdx.x == 0) {
        int next_city = selectCity(&cityProb[0], &cityProb[MAX_THREADS]);
        tour_length += 1; //TODO: find some way to access global graph edge weights
        path[num_visited++] = next_city;
        current_city = next_city;
      }

      __syncthreads();
    }

    //extract best ant tour length and write the paths out to global memory
    if (threadIdx.x == 0) {
      tourResults[antId] = tour_length;
      memcpy(pathResults + antId * MAX_CITIES, 
             path,
             MAX_CITIES * sizeof(int));
    }
}


double cuda_ACO(cityType *cities) {
  int best_index = -1;
  double best = (double) MAX_TOUR;
  int* bestPath[MAX_CITIES]; 
  dim3 numBlocks(MAX_ANTS);
  dim3 threadsPerBlock(MAX_THREADS);

  double *copiedTourResults = new double[MAX_ANTS];

  double *tourResults;
  int* pathResults;
  cudaMalloc((void**)&pathResults, sizeof(int) * MAX_ANTS * MAX_CITIES);
  cudaMalloc((void**)&tourResults, sizeof(double) * MAX_ANTS);

  for (int i = 0; i < MAX_TIME; i++) {
    constructAntTour<<<numBlocks, threadsPerBlock>>>(tourResults, pathResults);
    cudaThreadSynchronize();

    cudaMemcpy(copiedTourResults, tourResults, MAX_ANTS * sizeof(double),
               cudaMemcpyDeviceToHost);

    //find the best tour result from all the ants
    for (int j = 0; j < MAX_ANTS; j++) {
      if (copiedTourResults[j] < best) {
        best = copiedTourResults[j];
        best_index = j;
      }
    }

    //copy the corresponding tour for the best ant
    if (best_index != -1) {
      cudaMemcpy(bestPath, 
                 &pathResults[MAX_CITIES * best_index], 
                 MAX_CITIES * sizeof(int),
                 cudaMemcpyDeviceToHost);
    }

    //TODO: pheromone update
  }

  cudaFree(pathResults);
  cudaFree(tourResults);
  delete copiedTourResults;
  return best;
}
