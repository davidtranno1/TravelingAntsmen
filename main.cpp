/**
 *  main.cpp
 *  Run and compare sequential and parallel version of ACO algorithm for TSP
 */

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <unistd.h>

#include "CycleTimer.h"
#include "ants.h"

extern float cuda_ACO(EdgeMatrix *dist, int *bestPath);
extern float seq_ACO(EdgeMatrix *dist, int *bestPath);


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

      // if both cities lie on top of each other, manually set edge weight to 1
      float edge_dist = sqrt(xd + yd);
      if (sqrt(xd + yd) == 0) {
        edge_dist = 1;
      }
      (*dist)[from][to] = edge_dist;
      (*dist)[to][from] = edge_dist;
    }
  }
}

void savePathDataFile(int *path, char *filename)
{
  std::ofstream f1;
  f1.open(filename);
  for (int i = 0; i < MAX_CITIES; i++) {
    f1 << path[i] << " ";
  }
  f1 << std::endl;
  f1.close();
}

void saveCityDataFile(cityType *cities)
{
  std::ofstream f1;
  f1.open("city_data.txt");
  for (int i = 0; i < MAX_CITIES; i++) {
    f1 << cities[i].x << " " << cities[i].y << std::endl;
  }
  f1.close();
}

// Check if two paths match
bool matchPaths(int *path1, int *path2) {
  for (int i = 0; i < MAX_CITIES; i++) {
    if (path1[i] != path2[i])
      return false;
  }
  return true;
}

// Check if path is a legitimate tour (all unique vertices)
bool checkTourUnique(int *path) {
  bool cities[MAX_CITIES];
  memset(&cities, 0, sizeof(bool) * MAX_CITIES);

  for (int i = 0; i < MAX_CITIES; i++) {
    if (cities[path[i]]) {
      return false;
    }
    cities[path[i]] = true;
  }
  return true;
}

// Check if path is a legitimate tour (correct distance)
bool checkTourLength(int *path, EdgeMatrix *dist, float length) {
  float distance = 0.0;

  for (int i = 0; i < MAX_CITIES; i++) {
    distance += (*dist)[path[i]][path[(i+1) % MAX_CITIES]];
  }

  return distance == length;
}


int main(int argc, char *argv[]) {
  int opt;
  bool par_only = false;
  while ((opt = getopt(argc, argv, "p")) != -1) {
    switch (opt) {
      case 'p':
        par_only = true;
        break;
      default:
        fprintf(stderr, "Usage: %s [-p]\n", argv[0]);
        exit(EXIT_FAILURE);
    }
  }


  // Initialize TSP graph
  cityType cities[MAX_CITIES];
  EdgeMatrix *dist = new EdgeMatrix();
  int seqPath[MAX_CITIES];
  int parPath[MAX_CITIES];

  std::cout.precision(12);
  std::cout << "Constructing graph..." << std::endl;
  constructTSP(cities, dist);
  saveCityDataFile(cities);

  // Sequential algorithm
  float startTime, endTime;
  float seqTourLength, seqTime;
  if (!par_only) {
    std::cout << "Running sequential ant algorithm..." << std::endl;
    startTime = CycleTimer::currentSeconds();
    seqTourLength = seq_ACO(dist, seqPath);
    endTime = CycleTimer::currentSeconds();
    seqTime = endTime - startTime;
    std::cout << "Found tour of length " << seqTourLength << std::endl;
    if (!checkTourUnique(seqPath)) {
      std::cout << "Error: invalid tour (repeated cities!)" << std::endl;
    }
    if (!checkTourLength(seqPath, dist, seqTourLength)) {
      std::cout << "Error: invalid tour (length mismatch!)" << std::endl;
    }
    savePathDataFile(seqPath, (char *)"path_seq.txt");
  }

  // Parallel algorithm
  std::cout << "Running parallel ant algorithm..." << std::endl;
  startTime = CycleTimer::currentSeconds();
  float parTourLength = cuda_ACO(dist, parPath);
  endTime = CycleTimer::currentSeconds();
  float parTime = endTime - startTime;
  std::cout << "Found tour of length " << parTourLength << std::endl;
  if (!checkTourUnique(parPath)) {
    std::cout << "Error: invalid tour (repeated cities!)" << std::endl;
  }
  if (!checkTourLength(parPath, dist, parTourLength)) {
    std::cout << "Error: invalid tour (length mismatch!)" << std::endl;
  }
  savePathDataFile(parPath, (char *)"path_par.txt");

  // check correctness and print data
  if (!par_only) {
    if (seqTourLength == parTourLength && matchPaths(seqPath, parPath)) {
      std::cout << "Correctness passed!" << std::endl;
    } else {
      std::cout << "Uh oh! Found two different tours..." << std::endl;
    }
  }
  std::cout << std::endl;
  if (!par_only) {
    std::cout << "Sequential runtime: " << seqTime << " s" << std::endl;
    std::cout << "Parallel runtime: " << parTime << " s" << std::endl;
    std::cout << "Speedup: " << seqTime / parTime << "x" << std::endl;
  } else {
    std::cout << "Parallel runtime: " << parTime << " s" << std::endl;
  }

  delete dist;
  return 0;
}

