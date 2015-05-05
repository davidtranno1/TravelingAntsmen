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

#include <math_functions.h>

#include "ants.h"

#define MAX_THREADS 256

__device__ static inline int toIndex(int i, int j) {
  return i * MAX_CITIES + j;
}

__device__ double cudaAntProduct(double *edges, double *phero, int from, int to) {
  // TODO: delete this when we're sure it's fixed
  if (isinf(pow(1.0 / edges[toIndex(from, to)], BETA))) {
    printf("OH NO INFINITY: dist = %1.15f\n", edges[toIndex(from, to)]);
  }
  /*if (pow(phero[toIndex(from, to)], ALPHA) * pow(1.0 / edges[toIndex(from, to)], BETA) == 0) {
    printf("I'M ZERO\n");
  }*/

  return (pow(phero[toIndex(from, to)], ALPHA) * pow(1.0 / edges[toIndex(from, to)], BETA));
}

__global__ void init_rand(curandState *state) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  curand_init(1230, idx, 0, &state[idx]);
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
    /*if (start[i] > 0) {
      printf("%1.15f\n", start[i]);
    }*/
    sum += start[i];
  }

  if (sum == 0.0) {
    return 0;
  }

  make_rand(state, randArray);
  double luckyNumber = (double)randArray[idx];
  double acc = 0;

  for (int i = 0; i < length; i++) {
    acc += start[i] / sum;
    if (acc >= luckyNumber) {
      /*if (idx == 0) {
        printf("SUM: %1.15f, ACC: %1.15f, LUCKYNUM: %1.15f, i: %d, length: %d\n", sum, acc, luckyNumber, i, length);
      }*/
      return i;
    }
  }

  return 0;
}

__device__ int calculateFrom(int i){
  //find least triangle number less than i
  int row = (-1 + (sqrt((float)(1 + 8 * i)))) / 2;
  int tnum = (row * (row + 1)) / 2;
  int remain = i - tnum;
  return remain;
}

__device__ int calculateTo(int i){
  //find least triangle number less than i
  int row = (-1 + (sqrt((float)(1 + 8 * i)))) / 2;
  int tnum = (row * (row + 1)) / 2;
  int remain = i - tnum;
  return MAX_CITIES - (row - remain) - 1;
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
    __shared__ int path[MAX_CITIES];
    __shared__ int current_city;
    __shared__ int num_visited;
    __shared__ double tour_length;
    __shared__ int bestCities[MAX_THREADS];
    __shared__ double cityProb[MAX_THREADS];

    const int citiesPerThread = (MAX_CITIES + MAX_THREADS - 1) / MAX_THREADS;
    int startCityIndex = threadIdx.x * citiesPerThread;
    int antId = blockIdx.x;

    double localCityProb[citiesPerThread];

    if (threadIdx.x == 0) {
      current_city = 0; //TODO: random make it random random
      num_visited = 1;
      tour_length = 0.0;
      tabu[current_city] = 1;
      path[0] = current_city;
    }
    __syncthreads();

    //initiailize tabu list to zero
    for (int i = 0; i < citiesPerThread; i++) {
      int city = i + startCityIndex;
      if (city >= MAX_CITIES) {
        break;
      }
      if (city != current_city) {
        tabu[city] = 0;
      }
    }

    __syncthreads();

    //check if we have finished the tour
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
      /*if (antId == 0 && threadIdx.x == 0) {
        for (int i = 0; i < citiesPerThread; i++) {
          printf("localCityProb[%d] = %1.15f\n", i, localCityProb[i]);
        }
      }*/
      bestCities[threadIdx.x] = localCity + startCityIndex;
      //printf("BESTCITIES: %d\n", localCity + startCityIndex);
      __syncthreads();

      //reduce over cityProb and randomly pick a city based on
      //those probabilities
      if (threadIdx.x == 0) {
        int nextIndex = selectCity(state, randArray, cityProb, MAX_THREADS);
        /*if (antId == 0) {
          printf("next city index: %d\n", nextIndex);
          printf("next city prob: %1.15f\n", cityProb[nextIndex]);
        }*/
        int next_city = bestCities[nextIndex];
        tour_length += edges[toIndex(current_city, next_city)];
        path[num_visited++] = next_city;
        current_city = next_city;
        tabu[current_city] = 1;
      }

      __syncthreads();
    }

    //extract best ant tour length and write the paths out to global memory
    if (threadIdx.x == 0) {
      tour_length += edges[toIndex(current_city, path[0])];
      tourResults[antId] = tour_length;
      memcpy(pathResults + antId * MAX_CITIES, path, MAX_CITIES * sizeof(int));
    }
}

// Evaporate pheromones along each edge
__global__ void evaporatePheromones(double *phero) {
  int cityId = blockDim.x * blockIdx.x + threadIdx.x;

  for (int to = 0; to < MAX_CITIES; to++) {
    if (cityId == to) {
      continue;
    }

    int idx = toIndex(cityId, to);
    phero[idx] *= 1.0 - RHO;

    if (phero[idx] < 0.0) {
      phero[idx] = INIT_PHER;
    }
  }
}


// Add new pheromone to the trails
__global__ void updateTrails(double *phero, int *paths, double *tourLengths)
{
  //int antId = threadIdx.x;
  int from, to;
 
  __shared__ int numPhero;
  __shared__ int blockStartPhero;
  
  if (threadIdx.x == 0) {
    numPhero = (((MAX_CITIES * MAX_CITIES) / 2) + (MAX_THREADS * MAX_ANTS - 1))
               / (MAX_THREADS * MAX_ANTS);
    blockStartPhero = numPhero * MAX_THREADS * blockIdx.x;
  }
  
  __syncthreads();
  
  int cur_phero;
  // TODO: need to be atomic
  for (int i = 0; i < MAX_ANTS; i++) {
    for (int j = 0; j < numPhero; j++) {
      cur_phero = blockStartPhero + j + numPhero * threadIdx.x;
      
      if (cur_phero >= (MAX_CITIES * MAX_CITIES) / 2) {
        break;
      }
      
      from = calculateFrom(cur_phero); //triangle number thing
      to = calculateTo(cur_phero);
      bool touched = false;
      int checkTo;
      int checkFrom;
      for (int k = 0; k < MAX_CITIES; k++) {
        checkFrom = paths[toIndex(i, k)];
        if (k < MAX_CITIES - 1) {
          checkTo = paths[toIndex(i, k + 1)];
        } else {
          checkTo = paths[toIndex(i, 0)];
        }
        
        //printf("NEW VALUE: to: %d, from: %d", to, from);
        if ((checkFrom == from && checkTo == to) ||
            (checkFrom == to && checkTo == from)) 
        {
          touched = true;
          break;      
        }
      }
      
      if (touched) {
        phero[toIndex(from, to)] += (QVAL / tourLengths[i]);
        phero[toIndex(from, to)] *= RHO;
        phero[toIndex(to, from)] = phero[toIndex(from, to)];
        /*if (i == 0) {
          printf("NEW VALUE: to: %d, from: %d, value: %f\n", to, from, phero[toIndex(to, from)]);
        }*/
      }
    }
  }
  
  /*for (int i = 0; i < MAX_CITIES; i++) {
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
  }*/
}

double cuda_ACO(EdgeMatrix *dist, int *bestPath) { 
  dim3 numAntBlocks(MAX_ANTS);
  dim3 numCityBlocks((MAX_CITIES + MAX_THREADS - 1) / MAX_THREADS);
  dim3 threadsPerBlock(MAX_THREADS);
  dim3 single(1);

  int best_index = -1;
  double best = (double) MAX_TOUR;

  // allocate host memory
  double *copiedTourResults = new double[MAX_ANTS];

  // allocate device memory
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
  init_rand<<<numAntBlocks, threadsPerBlock>>>(randState);

  cudaMemcpy(deviceEdges, dist->get_array(), sizeof(double) * MAX_CITIES * MAX_CITIES,
             cudaMemcpyHostToDevice);

  initPhero<<<single, single>>>(phero);

  for (int i = 0; i < MAX_TOURS; i++) {
    best_index = -1;
    constructAntTour<<<numAntBlocks, threadsPerBlock>>>(deviceEdges, phero, randState, randArray, tourResults, pathResults);
    cudaThreadSynchronize();

    cudaMemcpy(copiedTourResults, tourResults, sizeof(double) * MAX_ANTS,
               cudaMemcpyDeviceToHost);

    //find the best tour result from all the ants
    for (int j = 0; j < MAX_ANTS; j++) {
      if (copiedTourResults[j] < best) {
        best = copiedTourResults[j];
        printf("new best: %f\n", best);
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

    //evaporate pheromones in parallel
    evaporatePheromones<<<numCityBlocks, threadsPerBlock>>>(phero);
    cudaThreadSynchronize();

    //TODO: pheromone update
    updateTrails<<<numAntBlocks, threadsPerBlock>>>(phero, pathResults, tourResults); 
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
