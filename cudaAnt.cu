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

#include "CycleTimer.h"
#include "ants.h"

#define MAX_THREADS 128

__device__ static inline int toIndex(int i, int j) {
  return i * MAX_CITIES + j;
}

__device__ static inline float cudaAntProduct(float *edges, float *phero, int city) {
  // TODO: delete this when we're sure it's fixed
  /*if (isinf(pow(1.0 / edges[toIndex(from, to)], BETA))) {
    printf("OH NO INFINITY: dist = %1.15f\n", edges[toIndex(from, to)]);
  }
  if (pow(phero[toIndex(from, to)], ALPHA) * pow(1.0 / edges[toIndex(from, to)], BETA) == 0) {
    printf("I'M ZERO\n");
  }
  if (isnan(powf(1.0 / edges[city], BETA))) {
    printf("IS NAN: city %d\n", city);
    return 0;
  }*/
  return (powf(phero[city], ALPHA) * powf(1.0 / edges[city], BETA));
}

__global__ void init_rand(curandState *state) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  curand_init(418, idx, 0, &state[idx]);
}

__device__ static inline void make_rand(curandState *state, float *randArray) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  randArray[idx] = curand_uniform(&state[idx]);
}

// Randomly select a city based off an array of values (return the index)
__device__ int selectCity(curandState *state, float *randArray, float *start, int length) {
  float sum = 0;
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
  /*if (isnan(sum)) {
    printf("error; value is nan!\n");
    return 0;
  }*/

  make_rand(state, randArray);
  float luckyNumber = (float)randArray[idx];
  float acc = 0;

  int lastBestIndex = 0;
  for (int i = 0; i < length; i++) {
    float value = start[i] / sum;
    if (value > 0) {
      acc += value;
      lastBestIndex = i;

      if (acc >= luckyNumber) {
        /*if (idx == 0) {
          printf("SUM: %1.15f, ACC: %1.15f, LUCKYNUM: %1.15f, i: %d, length: %d\n", sum, acc, luckyNumber, i, length);
        }*/
        return i;
      }

    }
  }

  //printf("warning: acc did not reach luckyNumber in selectNextCity\n");
  //printf("sum: %1.15f, acc: %1.15f, luckyNumber: %1.15f\n", sum, acc, luckyNumber);
  return lastBestIndex;
}

__device__ static inline int calculateFrom(int i) {
  //find least triangle number less than i
  int row = (int)(-1 + (sqrt((float)(1 + 8 * i)))) >> 1;
  int tnum = (row * (row + 1)) >> 1;
  int remain = i - tnum;
  return MAX_CITIES - 1 - remain;
}

__device__ static inline int calculateTo(int i) {
  //find least triangle number less than i
  int row = (int)(-1 + (sqrt((float)(1 + 8 * i)))) >> 1;
  int tnum = (row * (row + 1)) >> 1;
  int remain = i - tnum;
  return row - remain;
}

__global__ void initPhero(float *phero) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx >= MAX_CITIES * MAX_CITIES) {
    return;
  }
  phero[idx] = INIT_PHER;
}

__global__ void copyBestPath(int i, int *bestPathResult, int *pathResults) {
  memcpy(bestPathResult, &pathResults[i * MAX_ANTS], MAX_CITIES * sizeof(int));
}

__global__ void constructAntTour(float *edges, float *phero,
                                 curandState *state, float *randArray,
                                 float *tourResults, int *pathResults) {
    __shared__ bool tabu[MAX_CITIES]; //TODO: put in register wtf is that
    //__shared__ int path[MAX_CITIES];
    __shared__ int current_city;
    __shared__ int num_visited;
    __shared__ int bestCities[MAX_THREADS];
    __shared__ float cityProb[MAX_THREADS];

    __shared__ float localEdges[MAX_CITIES];
    __shared__ float localPhero[MAX_CITIES];

    const int citiesPerThread = (MAX_CITIES + MAX_THREADS - 1) / MAX_THREADS;
    const int startCityIndex = threadIdx.x * citiesPerThread;
    const int antId = blockIdx.x;
    float tour_length;

    if (startCityIndex >= MAX_CITIES) {
      //printf("i'm over. startCityIndex: %d, MAX_CITIES: %d\n", startCityIndex, MAX_CITIES);
      cityProb[threadIdx.x] = 0;
      return;
    }

    float localCityProb[citiesPerThread];

    if (threadIdx.x == 0) {
      make_rand(state, randArray);
      current_city = randArray[antId * blockDim.x + threadIdx.x] * MAX_CITIES;
      num_visited = 1;
      tour_length = 0.0;
      tabu[current_city] = true;
      pathResults[antId * MAX_CITIES] = current_city;
    }
    __syncthreads();

    //initiailize tabu list to zero
    for (int i = 0; i < citiesPerThread; i++) {
      int city = i + startCityIndex;
      if (city >= MAX_CITIES) {
        break;
      }
      if (city != current_city) {
        tabu[city] = false;
      }
    }

    __syncthreads();

    //check if we have finished the tour
    while (num_visited != MAX_CITIES) {

      int tile;
      if (startCityIndex + citiesPerThread >= MAX_CITIES) {
        tile = MAX_CITIES - startCityIndex;
      } else {
        tile = citiesPerThread;
      }
      memcpy(&localEdges[startCityIndex], &edges[current_city * MAX_CITIES + startCityIndex], tile * sizeof(float));
      memcpy(&localPhero[startCityIndex], &phero[current_city * MAX_CITIES + startCityIndex], tile * sizeof(float));

      __syncthreads();

      //pick next (unvisited) city
      for (int i = 0; i < citiesPerThread; i++) {
        int city = i + startCityIndex;

        if (city >= MAX_CITIES || tabu[city]) {
          localCityProb[i] = 0.0;
        } else {
          localCityProb[i] = cudaAntProduct(localEdges, localPhero, city);
          //printf("city prob: %1.15f\n", localCityProb[i]);
        }
      }
      //if (threadIdx.x == 0)
      //  printf("cuda ant product done\n");

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

      //reduce over bestCities and pick city with best absolute heuristic
      if (threadIdx.x == 0) {
        //int nextIndex = selectCity(state, randArray, cityProb, MAX_THREADS);
        //int next_city = bestCities[nextIndex];
        //printf("best cities done in block %d\n", blockIdx.x);
        float best_distance = MAX_DIST * 2;
        int next_city = -1;
        for (int i = 0; i < MAX_THREADS && i < MAX_CITIES; i++) {
          if (cityProb[i] == 0) {
            continue;
          }
          //printf("best city[%d]: %d\n", i, bestCities[i]);
          float distance = localEdges[bestCities[i]];
          if (distance < best_distance) {
            best_distance = distance;
            next_city = bestCities[i];
          }
        }
        /*if (next_city == -1) {
          printf("OH NO\n");
        }*/
        /*if (antId == 0) {
          printf("current: %d, next: %d\n", current_city, next_city);
        }*/
        tour_length += localEdges[next_city];
        pathResults[antId * MAX_CITIES + num_visited] = next_city;
        num_visited++;
        current_city = next_city;
        tabu[current_city] = true;
      }

      __syncthreads(); //TODO: move this syncthreads?
    }

    //extract best ant tour length and write the paths out to global memory
    if (threadIdx.x == 0) {
      tour_length += edges[toIndex(current_city, pathResults[antId * MAX_CITIES])];
      tourResults[antId] = tour_length;
    }
}

// Evaporate pheromones along each edge
__global__ void evaporatePheromones(float *phero) {
  int current_phero = blockIdx.x * blockDim.x + threadIdx.x;
  if (current_phero >= NUM_EDGES) {
    return;
  }
  int from = calculateFrom(current_phero); //triangle number thing
  int to = calculateTo(current_phero);

  int idx = toIndex(from, to);
  phero[idx] *= 1.0 - RHO;

  if (phero[idx] < 0.0) {
    phero[idx] = INIT_PHER;
  }
  phero[toIndex(to, from)] = phero[idx];
}


// Add new pheromone to the trails
__global__ void updateTrailsAtomic(float *phero, int *paths, float *tourLengths)
{
  int antId = blockIdx.x;
  int from, to;

  for (int i = 0; i < MAX_CITIES; i++) {
    from = paths[toIndex(antId, i)];
    if (i < MAX_CITIES - 1) {
      to = paths[toIndex(antId, i+1)];
    } else {
      to = paths[toIndex(antId, 0)];
    }

    if (from < to) {
      int tmp = from;
      from = to;
      to = tmp;
    }
    atomicAdd(&phero[toIndex(from, to)], QVAL / tourLengths[antId]);
    //atomicExch(&phero[toIndex(to, from)], phero[toIndex(from, to)]);
  }
}

__global__ void updateSymmetricPhero(float *phero) {
  for (int i = 0; i < MAX_CITIES; i++) {
    for (int j = 0; j < i; j++) {
      //phero[toIndex(i, j)] *= RHO;
      phero[toIndex(j, i)] = phero[toIndex(i, j)];
    }
  }
}

__global__ void updateTrails(float *phero, int *paths, float *tourLengths)
{
  //int antId = threadIdx.x;
  //__shared__ float localPaths[MAX_CITIES];

  int numPhero = (NUM_EDGES + (blockDim.x * (MAX_ANTS * 2) - 1)) /
                 (blockDim.x * (MAX_ANTS * 2));
  int blockStartPhero = numPhero * blockDim.x * blockIdx.x;
  int from, to;

  int cur_phero;
  for (int i = 0; i < MAX_ANTS; i++) {
    // For each ant, cache paths in shared memory
    /*int tile;
    if (startCityIndex + citiesPerThread >= MAX_CITIES) {
      tile = MAX_CITIES - startCityIndex;
    } else {
      tile = citiesPerThread;
    }
    memcpy(&localPaths[startCityIndex], &paths[i * MAX_CITIES + startCityIndex], tile * sizeof(float));
    */
    // TODO: figure out tiling
    /*if (threadIdx.x == 0) {
      memcpy(&localPaths, &paths[i * MAX_CITIES], MAX_CITIES * sizeof(float));
    }

    __syncthreads();
    */

    for (int j = 0; j < numPhero; j++) {
      cur_phero = blockStartPhero + j + numPhero * threadIdx.x;

      if (cur_phero >= NUM_EDGES) {
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
        int idx = toIndex(from, to);
        phero[idx] += (QVAL / tourLengths[i]);
        //phero[idx] *= RHO;
        phero[toIndex(to, from)] = phero[idx];
        /*if (i == 0) {
          printf("NEW VALUE: to: %d, from: %d, value: %f\n", to, from, phero[toIndex(to, from)]);
        }*/
      }
    }
    //__syncthreads();
  }
}


__global__ void checkPhero(float *pheroSeq, float *phero) {
  for (int i = 0; i < MAX_CITIES; i++) {
    for (int j = 0; j < MAX_CITIES; j++) {
      if (i == j) continue;
      int idx = toIndex(i, j);
      if (fabsf(pheroSeq[idx] - phero[idx]) > 0.001) {
        printf("PHERO IS BROKEN at (%d, %d); expected: %1.15f, actual: %1.15f\n", i, j, pheroSeq[idx], phero[idx]);
      }
    }
  }
}
__global__ void seqPheroUpdate(float *phero, float *pheroReal, int *paths, float *tourLengths) {
  memcpy(phero, pheroReal, sizeof(float) * MAX_CITIES * MAX_CITIES);

  int from, to;
  // evaporate
  for (from = 0; from < MAX_CITIES; from++) {
    for (to = 0; to < from; to++) {
      phero[toIndex(from, to)] *= 1.0 - RHO;

      if (phero[toIndex(from, to)] < 0.0) {
        phero[toIndex(from, to)] = INIT_PHER;
      }
      phero[toIndex(to, from)] = phero[toIndex(from, to)];
    }
  }

  //Add new pheromone to the trails
  for (int ant = 0; ant < MAX_ANTS; ant++) {
    for (int i = 0; i < MAX_CITIES; i++) {
      from = paths[toIndex(ant, i)];
      if (i < MAX_CITIES - 1) {
        to = paths[toIndex(ant, i+1)];
      } else {
        to = paths[toIndex(ant, 0)];
      }

      phero[toIndex(from, to)] += (QVAL / tourLengths[ant]);
      phero[toIndex(to, from)] = phero[toIndex(from, to)];
    }
  }

}

float cuda_ACO(EdgeMatrix *dist, int *bestPath) {
  dim3 numAntBlocks(MAX_ANTS);
  dim3 numTwoAntBlocks(MAX_ANTS * 2);
  dim3 numCityBlocks((MAX_CITIES + MAX_THREADS - 1) / MAX_THREADS);
  dim3 numEdgesBlocks((MAX_CITIES * MAX_CITIES + MAX_THREADS - 1) / MAX_THREADS);
  dim3 numPheroBlocks((NUM_EDGES + MAX_THREADS - 1) / MAX_THREADS);
  dim3 threadsPerBlock(MAX_THREADS);
  dim3 single(1);

  int best_index;
  float best = (float) MAX_TOUR;

  // allocate host memory
  float *copiedTourResults = new float[MAX_ANTS];

  // allocate device memory
  float *tourResults;
  int *pathResults;
  int *bestPathResult;
  float *deviceEdges;
  float *phero;
  float *testPhero;
  float *randArray;
  curandState *randState;

  cudaMalloc((void**)&pathResults, sizeof(int) * MAX_ANTS * MAX_CITIES);
  cudaMalloc((void**)&tourResults, sizeof(float) * MAX_ANTS);
  cudaMalloc((void**)&deviceEdges, sizeof(float) * MAX_CITIES * MAX_CITIES);
  cudaMalloc((void**)&phero, sizeof(float) * MAX_CITIES * MAX_CITIES);
  cudaMalloc((void**)&testPhero, sizeof(float) * MAX_CITIES * MAX_CITIES);
  cudaMalloc(&randState, sizeof(curandState) * MAX_ANTS * MAX_THREADS);
  cudaMalloc((void**)&randArray, sizeof(float) * MAX_ANTS * MAX_THREADS);
  cudaMalloc((void**)&bestPathResult, sizeof(int) * MAX_CITIES);
  init_rand<<<numAntBlocks, threadsPerBlock>>>(randState);

  cudaMemcpy(deviceEdges, dist->get_array(), sizeof(float) * MAX_CITIES * MAX_CITIES,
             cudaMemcpyHostToDevice);

  initPhero<<<numEdgesBlocks, threadsPerBlock>>>(phero);
  cudaThreadSynchronize();

  float pathTime = 0;
  float pheroTime = 0;
  float sBegin;
  float sEnd;
  for (int i = 0; i < MAX_TOURS; i++) {
    best_index = -1;

    sBegin = CycleTimer::currentSeconds();
    constructAntTour<<<numAntBlocks, threadsPerBlock>>>(deviceEdges, phero, randState, randArray, tourResults, pathResults);
    cudaThreadSynchronize();
    sEnd = CycleTimer::currentSeconds();

    pathTime += (sEnd - sBegin);

    cudaMemcpy(copiedTourResults, tourResults, sizeof(float) * MAX_ANTS,
               cudaMemcpyDeviceToHost);

    //find the best tour result from all the ants
    for (int j = 0; j < MAX_ANTS; j++) {
      if (copiedTourResults[j] < best) {
        best = copiedTourResults[j];
        printf("new best: %1.f\n", best);
        best_index = j;
      }
    }

    //copy the corresponding tour for the best ant
    if (best_index != -1) {
      copyBestPath<<<single, single>>>(best_index, bestPathResult, pathResults);
    }

    //seqPheroUpdate<<<single, single>>>(testPhero, phero, pathResults, tourResults);
    //cudaThreadSynchronize();

    //evaporate pheromones in parallel
    sBegin = CycleTimer::currentSeconds();
    evaporatePheromones<<<numPheroBlocks, threadsPerBlock>>>(phero);
    cudaThreadSynchronize();

    //pheromone update
    updateTrailsAtomic<<<numAntBlocks, single>>>(phero, pathResults, tourResults);
    cudaThreadSynchronize();
    updateSymmetricPhero<<<single, single>>>(phero);
    cudaThreadSynchronize();
    sEnd = CycleTimer::currentSeconds();
    pheroTime += (sEnd - sBegin);

    //checkPhero<<<single, single>>>(testPhero, phero);
    //cudaThreadSynchronize();
  }

  printf("PATHTIME: %f, PHEROTIME: %f\n", pathTime, pheroTime);

  cudaMemcpy(bestPath, bestPathResult, MAX_CITIES * sizeof(int), cudaMemcpyDeviceToHost);

  cudaFree(bestPathResult);
  cudaFree(pathResults);
  cudaFree(tourResults);
  cudaFree(deviceEdges);
  cudaFree(randArray);
  cudaFree(randState);
  delete copiedTourResults;
  return best;
}
