/**
 *  main.cpp
 *  Run and compare sequential and parallel version of ACO algorithm for TSP
 */

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <unistd.h>

#include "CycleTimer.h"
#include "ants.h"

extern float cuda_ACO(EdgeMatrix *dist, int *bestPath);
extern float seq_ACO(EdgeMatrix *dist, int *bestPath);

float euclideanDistance(int x1, int x2, int y1, int y2) {
  int xd = x1 - x2;
  int yd = y1 - y2;
  return (int) (sqrt(xd * xd + yd * yd) + 0.5);
}

float pseudoEuclideanDistance(int x1, int x2, int y1, int y2) {
  int xd = x1 - x2;
  int yd = y1 - y2;
  float rij = sqrt((xd * xd + yd * yd) / 10.0);
  return ceil(rij);
}

// Construct TSP graph
void constructTSP(std::string graph, cityType *cities, EdgeMatrix *dist) {
  // Load cities from file
  std::ifstream infile(("graphs/" + graph + ".tsp").c_str());
  std::string line;
  bool euclidean = true; // whether to run EUC_2D or ATT distance metric

  int city, x, y;
  bool reading_nodes = false;
  while (std::getline(infile, line)) {
    std::istringstream iss(line);
    std::string word;
    if (!reading_nodes) {
      iss >> word;
      if (word.compare("EDGE_WEIGHT_TYPE") == 0) {
        iss >> word >> word;
        std::cout << "edge type: " << word << std::endl;
        euclidean = !word.compare("EUC_2D");
      } else if (word.compare("NODE_COORD_SECTION") == 0) {
        reading_nodes = true;
      }
    } else if (iss >> city >> x >> y) {
      cities[city-1].x = x;
      cities[city-1].y = y;
    }
  }
  infile.close();

  // Compute distances between cities (edge weights)
  for (int from = 0; from < MAX_CITIES; from++) {
    (*dist)[from][from] = 0.0;

    for (int to = from + 1; to < MAX_CITIES; to++) {
      float edge_dist;
      if (euclidean) {
        edge_dist = euclideanDistance(cities[from].x, cities[to].x, cities[from].y, cities[to].y);
      } else {
        edge_dist = pseudoEuclideanDistance(cities[from].x, cities[to].x, cities[from].y, cities[to].y);
      }
      if (edge_dist == 0) {
        edge_dist = 1.0;
      }
      //printf("edge[%d][%d] = %f\n", from, to, edge_dist);
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
    f1 << path[i] + 1 << " ";
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

float checkTourSolution(std::string graph, EdgeMatrix *dist) {
  // Load cities from file
  std::ifstream infile(("graphs/" + graph + ".opt.tour").c_str());
  if (!infile.good()) {
    printf("No solution to verify.\n");
    return 0.0;
  }
  std::string line;

  float distance = 0.0;
  int last_city = -1;
  int first_city = -1;
  int city;
  while (std::getline(infile, line)) {
    std::istringstream iss(line);
    if (iss >> city) {
      city--;
      if (last_city == -1) {
        first_city = last_city = city;
        continue;
      }
      if (city == -1) {
        city = first_city;
      }
      distance += (*dist)[last_city][city];
      last_city = city;
    }
  }
  infile.close();

  return distance;
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
  std::string graph = "att48";
  while ((opt = getopt(argc, argv, "pg:")) != -1) {
    switch (opt) {
      case 'p':
        par_only = true;
        break;
      case 'g':
        graph = optarg;
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
  std::cout << "Constructing graph " << graph << "..." << std::endl;
  constructTSP(graph, cities, dist);
  float optimal_length = checkTourSolution(graph, dist);

  // Sequential algorithm
  float startTime, endTime, parTime, seqTime;
  float approximation;
  float seqTourLength = -1;
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
    if (optimal_length > 0) {
      approximation = seqTourLength / optimal_length;
      std::cout << "Observed " << approximation << "x optimal path" << std::endl;
    }
    savePathDataFile(seqPath, (char *)"path_seq.txt");
  }

  // Parallel algorithm
  std::cout << "Running parallel ant algorithm..." << std::endl;
  startTime = CycleTimer::currentSeconds();
  float parTourLength = cuda_ACO(dist, parPath);
  endTime = CycleTimer::currentSeconds();
  parTime = endTime - startTime;
  std::cout << "Found tour of length " << parTourLength << std::endl;
  if (!checkTourUnique(parPath)) {
    std::cout << "Error: invalid tour (repeated cities!)" << std::endl;
  }
  if (!checkTourLength(parPath, dist, parTourLength)) {
    std::cout << "Error: invalid tour (length mismatch!)" << std::endl;
  }
  if (optimal_length > 0) {
    approximation = parTourLength / optimal_length;
    std::cout << "Observed " << approximation << "x optimal path" << std::endl;
  }
  savePathDataFile(parPath, (char *)"path_par.txt");

  // check correctness and print data
  if (!par_only) {
    approximation = parTourLength / seqTourLength;
    std::cout << "Parallel path is " << approximation << "x sequential path" << std::endl;
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

