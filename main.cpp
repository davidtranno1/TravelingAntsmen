/**
 *  main.cpp
 *  Run and compare sequential and parallel version of ACO algorithm for TSP
 */

#include <iostream>
#include <stdio.h>
#include <math.h>

#include "CycleTimer.h"
#include "ants.h"

extern double cuda_ACO(cityType *cities, EdgeMatrix *dist);
extern double seq_ACO(cityType *cities, EdgeMatrix *dist);


// Construct TSP graph
void constructTSP(cityType *cities, EdgeMatrix *dist) {
  // Create cities with random x, y coordinates
  for (int city = 0; city < MAX_CITIES; city++) {
    cities[city].x = rand() % MAX_DIST;
    cities[city].y = rand() % MAX_DIST;
  }

  // Compute distances between cities (edge weights)
  for (int from = 0; from < MAX_CITIES; from++) {
    (*dist)[from][from] = 0.0;

    for (int to = from + 1; to < MAX_CITIES; to++) {
      int xd = pow(abs(cities[from].x - cities[to].x), 2);
      int yd = pow(abs(cities[from].y - cities[to].y), 2);

      (*dist)[from][to] = sqrt(xd + yd);
      (*dist)[to][from] = (*dist)[from][to];
    }
  }
}


int main() {
  // Initialize TSP graph
  cityType cities[MAX_CITIES];
  EdgeMatrix *dist = new EdgeMatrix();

  std::cout << "Constructing graph..." << std::endl;
  constructTSP(cities, dist);

  double startTime, endTime;
  std::cout << "Running sequential ant algorithm..." << std::endl;
  startTime = CycleTimer::currentSeconds();
  double seqTourLength = seq_ACO(cities, dist);
  endTime = CycleTimer::currentSeconds();
  double seqTime = endTime - startTime;
  std::cout << "Found tour of length " << seqTourLength << std::endl;

  std::cout << "Running parallel ant algorithm..." << std::endl;
  startTime = CycleTimer::currentSeconds();
  double parTourLength = cuda_ACO(cities, dist);
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

  delete dist;
  return 0;
}

