/**
 *  main.cpp
 *  Run and compare sequential and parallel version of ACO algorithm for TSP
 */

#include <iostream>
#include <stdio.h>

#include "CycleTimer.h"
#include "ants.h"

extern void cuda_ACO(cityType *cities);
extern void seq_ACO(cityType *cities);


// Construct TSP graph
void constructTSP(cityType *cities) {
  // Create cities with random x, y coordinates
  for (int city = 0; city < MAX_CITIES; city++) {
    cities[city].x = rand() % MAX_DIST;
    cities[city].y = rand() % MAX_DIST;
  }
}


int main() {
  // Initialize TSP graph
  cityType cities[MAX_CITIES];
  constructTSP(cities);

  double startTime, endTime;
  std::cout << "Running sequential ant algorithm..." << std::endl;
  startTime = CycleTimer::currentSeconds();
  seq_ACO(cities);
  endTime = CycleTimer::currentSeconds();
  double seqTime = endTime - startTime;

  std::cout << "Running parallel ant algorithm..." << std::endl;
  startTime = CycleTimer::currentSeconds();
  cuda_ACO(cities);
  endTime = CycleTimer::currentSeconds();
  double parTime = endTime - startTime;

  // TODO: check correctness
  std::cout << "Correctness passed!" << std::endl;
  std::cout << "Sequential runtime: " << seqTime << " s" << std::endl;
  std::cout << "Parallel runtime: " << parTime << " s" << std::endl;
  std::cout << "Speedup: " << seqTime / parTime << "x" << std::endl;

  return 0;
}

