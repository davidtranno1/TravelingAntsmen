/**
 *  main.cpp
 *  Run and compare sequential and parallel version of ACO algorithm for TSP
 */

#include <iostream>
#include <stdio.h>

#include "CycleTimer.h"
#include "ants.h"

extern double cuda_ACO(cityType *cities);
extern double seq_ACO(cityType *cities);


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
  double seqTourLength = seq_ACO(cities);
  endTime = CycleTimer::currentSeconds();
  double seqTime = endTime - startTime;
  std::cout << "Found tour of length " << seqTourLength << std::endl;

  std::cout << "Running parallel ant algorithm..." << std::endl;
  startTime = CycleTimer::currentSeconds();
  double parTourLength = cuda_ACO(cities);
  endTime = CycleTimer::currentSeconds();
  double parTime = endTime - startTime;
  std::cout << "Found tour of length " << parTourLength << std::endl;

  // check correctness
  // TODO: check that tours are actually equal
  if (seqTourLength == parTourLength) {
    std::cout << "Correctness passed!" << std::endl;
  } else {
    std::cout << "Uh oh! Found two different tours..." << std::endl;
  }
  std::cout << std::endl;
  std::cout << "Sequential runtime: " << seqTime << " s" << std::endl;
  std::cout << "Parallel runtime: " << parTime << " s" << std::endl;
  std::cout << "Speedup: " << seqTime / parTime << "x" << std::endl;

  return 0;
}

