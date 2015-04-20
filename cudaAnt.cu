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

__device__ int selectCity(int *start, int *end) {
  return 0;
}

__global__ void constructAntTour(double *tourResults, int *pathResults) {
    __shared__ int tabu[MAX_CITIES]; //TODO: put in register wtf is that
    __shared__ int current_city;
    __shared__ int path[MAX_CITIES];
    __shared__ int num_visited;
    __shared__ double tour_length;
    __shared__ int bestCities[MAX_THREADS];
    __shared__ int cityProb[MAX_THREADS];
    
 
    if (threadIdx.x == 0) {
      current_city = 0; //TODO: random make it random random
      num_visited = 0;
    }
    
    const int citiesPerThread = (MAX_CITIES + MAX_THREADS - 1) / MAX_THREADS;
    int startCityIndex = threadIdx.x * citiesPerThread;
    int antId = blockIdx.x;
    
    int localCityProb[citiesPerThread];
    
    //initiailize tabu list to zero
    for (int i = startCityIndex; i < startCityIndex + citiesPerThread; i++) {
      tabu[i] = 0;
    }
    
    //check if we have finished the tour (can you run into dead ends?)
    while (num_visited != MAX_CITIES) {
      //pick next city
      for (int i = startCityIndex; i < startCityIndex + citiesPerThread; i++) {
        if (tabu[i] != 0) {
          localCityProb[i] = 0;
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
    
    if(threadIdx.x == 0) {
      //extract best ant tour length and write the paths out to global memory
      tourResults[blockIdx.x] = tour_length;
      memcpy(pathResults + blockIdx.x * MAX_CITIES, 
             path, 
             MAX_CITIES * sizeof(int));
    }
}


void cuda_ACO(cityType *cities) {
  double best = (double) MAX_TOUR;
  dim3 numBlocks(MAX_ANTS);
  dim3 threadsPerBlock(MAX_THREADS);
  
  double *tourResults;
  int* pathResults;
  cudaMalloc((void**)&pathResults, sizeof(int) * MAX_ANTS * MAX_CITIES);
  cudaMalloc((void**)&tourResults, sizeof(double) * MAX_ANTS);
  
  for (int i = 0; i < MAX_TIME; i++) {
    constructAntTour<<<numBlocks, threadsPerBlock>>>(tourResults, pathResults);
    cudaThreadSynchronize();
    //TODO: pheromone update
  }
  
  cudaFree(pathResults);
  cudaFree(tourResults);
  return;
}
