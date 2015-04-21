/**
 *  Project TACO: Parallel ACO algorithm for TSP
 *  15-418 Parallel Algorithms - Final Project
 *  Ivan Wang, Carl Lin
 */

#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <math.h>

#include "ants.h"

#define MAX_THREADS 100

__device__ static inline int toIndex(int i, int j) {
  return i * MAX_CITIES + j;
}

__device__ double cudaAntProduct(double *edges, double *phero, int from, int to) {
  return (pow(phero[toIndex(from, to)], ALPHA) * pow((1.0 / edges[toIndex(from, to)]), BETA));
}

__global__ void init_rand(curandState *state) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  curand_init(15418, idx, 0, &state[idx]);
}

__device__ void make_rand(curandState *state, double *randArray) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  randArray[idx] = curand_uniform_double(&state[idx]);
}

// Randomly select a city based off an array of values (return the index)
__device__ int selectCity(curandState *state, double *randArray, double *start, int length) {
  double sum = 0;
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  
  for (int i = 0; i < length; i++) {
    if (start[i] > 0) {
      //printf("%1.15f\n", start[i]);
    }
    sum += start[i];
  }
  
  if (sum == 0) {
    return 0;
  }
  
  make_rand(state, randArray);
  double luckyNumber = (double)randArray[idx] / RAND_MAX;
  double acc = 0;
  
  for (int i = 0; i < length; i++) {
    acc += start[i]/sum;
    if (acc >= luckyNumber) {
      return i;
    }
  }
  
  return 0;
}

__global__ void initPhero(double *phero) {
  for(int i = 0; i < MAX_CITIES * MAX_CITIES; i++){
    phero[i] = INIT_PHER;
  }
}

__global__ void constructAntTour(double *edges, double *phero, 
                                 curandState *state, double *randArray, 
                                 double *tourResults, int *pathResults) {
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
      tour_length = 0.0;
    }

    const int citiesPerThread = (MAX_CITIES + MAX_THREADS - 1) / MAX_THREADS;
    int startCityIndex = threadIdx.x * citiesPerThread;
    int antId = blockIdx.x;

    double localCityProb[citiesPerThread];

    //initiailize tabu list to zero
    for (int i = 0; i < citiesPerThread; i++) {
      int city = i + startCityIndex;
      if (city >= MAX_CITIES) {
        break;
      }
      tabu[city] = 0;
    }

    //check if we have finished the tour (can you run into dead ends?)
    while (num_visited != MAX_CITIES) {
      //pick next (unvisited) city
      for (int i = 0; i < citiesPerThread; i++) {
        int city = i + startCityIndex;

        if (city >= MAX_CITIES || tabu[city] != 0) {
          localCityProb[i] = 0.0;
        } else {
          localCityProb[i] = cudaAntProduct(edges, phero, current_city, city);
          //printf("city prob: %1.15f\n", localCityProb[i]);
        }
      }

      //for each thread, look through cities and stochastically select one
      int localCity = selectCity(state, randArray, localCityProb, citiesPerThread);
      cityProb[threadIdx.x] = localCityProb[localCity];
      bestCities[threadIdx.x] = localCity + startCityIndex;
      //printf("BESTCITIES: %d\n", localCity + startCityIndex);
      __syncthreads();

      //reduce over cityProb and randomly pick a city based on
      //those probabilities
      if (threadIdx.x == 0) {
        int nextIndex = selectCity(state, randArray, cityProb, MAX_THREADS);
        //printf("next city index: %d\n", nextIndex);
        int next_city = bestCities[nextIndex];
        tour_length += edges[toIndex(current_city, next_city)]; 
        path[num_visited++] = next_city;
        current_city = next_city;
      }

      __syncthreads();
    }

    //extract best ant tour length and write the paths out to global memory
    if (threadIdx.x == 0) {
      tourResults[antId] = tour_length;
      memcpy(pathResults + antId * MAX_CITIES, path, MAX_CITIES * sizeof(int));
    }
}

__global__ void updateTrails(double *phero, int *paths, double *tourLengths)
{
  int antId = blockIdx.x;
  int from, to;
  
  // Pheromone Evaporation
  if (antId == 0) {
    for (from = 0; from < MAX_CITIES; from++) {
      for (to = 0; to < MAX_CITIES; to++) {
        if (from != to) {
          phero[toIndex(from, to)] *= 1.0 - RHO;

          if (phero[toIndex(from, to)] < 0.0) {
            phero[toIndex(from, to)] = INIT_PHER;
          }
        }
      }
    }
  }
  
  __syncthreads();
  //Add new pheromone to the trails
  for (int i = 0; i < MAX_CITIES; i++) {
    if (i < MAX_CITIES - 1) {
      from = paths[toIndex(antId, i)];
      to = paths[toIndex(antId, i+1)];
    } else {
      from = paths[toIndex(antId, i)];
      to = paths[toIndex(antId, 0)];
    }

    phero[toIndex(from, to)] += (QVAL / tourLengths[antId]);
    phero[toIndex(to, from)] = phero[toIndex(from, to)];
  }
  

  for (from = 0; from < MAX_CITIES; from++) {
    for (to = 0; to < MAX_CITIES; to++) {
      phero[toIndex(from, to)] *= RHO;
    }
  }
}

double cuda_ACO(EdgeMatrix *dist, int *bestPath) {
  int best_index = -1;
  double best = (double) MAX_TOUR;
  dim3 numBlocks(MAX_ANTS);
  dim3 threadsPerBlock(MAX_THREADS);
  dim3 single(1);
  
  double *copiedTourResults = new double[MAX_ANTS];

  // malloc device memory
  double *tourResults;
  int* pathResults;
  double *deviceEdges;
  double *phero;
  double *randArray;
  curandState *randState;
  
  cudaMalloc((void**)&pathResults, sizeof(int) * MAX_ANTS * MAX_CITIES);
  cudaMalloc((void**)&tourResults, sizeof(double) * MAX_ANTS);
  cudaMalloc((void**)&deviceEdges, sizeof(double) * MAX_CITIES * MAX_CITIES);
  cudaMalloc((void**)&phero, sizeof(double) * MAX_CITIES * MAX_CITIES);
  cudaMalloc(&randState, sizeof(curandState) * MAX_ANTS * MAX_THREADS);
  cudaMalloc((void**)&randArray, sizeof(double) * MAX_ANTS * MAX_THREADS);
  init_rand<<<numBlocks, threadsPerBlock>>>(randState);
  
  cudaMemcpy(deviceEdges, dist->get_array(), sizeof(double) * MAX_ANTS,
            cudaMemcpyHostToDevice);
  
  initPhero<<<single, single>>>(phero);
  
  for (int i = 0; i < MAX_TIME; i++) {
    constructAntTour<<<numBlocks, threadsPerBlock>>>(deviceEdges, phero, randState, randArray, tourResults, pathResults);
    cudaThreadSynchronize();

    cudaMemcpy(copiedTourResults, tourResults, sizeof(double) * MAX_ANTS,
               cudaMemcpyDeviceToHost);

    //find the best tour result from all the ants
    for (int j = 0; j < MAX_ANTS; j++) {
      if (copiedTourResults[j] < best) {
        printf("new best: %f\n", best);
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
    updateTrails<<<numBlocks, single>>>(phero, pathResults, tourResults); //TODO: change single for optimization
    cudaThreadSynchronize();
  }

  cudaFree(pathResults);
  cudaFree(tourResults);
  cudaFree(deviceEdges);
  cudaFree(randArray);
  cudaFree(randState);
  delete copiedTourResults;
  return best;
}
