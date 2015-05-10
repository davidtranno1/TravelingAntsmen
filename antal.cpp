// Sequential ant algorithm for Traveling Salesman Problem

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "CycleTimer.h"
#include "ants.h"

// stores ant data (visited cities and current path)
struct antType {
  int curCity, nextCity, pathIndex;
  int tabu[MAX_CITIES];
  int path[MAX_CITIES];
  float tourLength;
};

// runtime structures and global variables
antType ants[MAX_ANTS];

EdgeMatrix *dist;
float phero[MAX_CITIES][MAX_CITIES];

float best = (float) MAX_TOUR;
int bestIndex;

// initializes the entire graph
void init() {
  int from, to, ant;

  for (from = 0; from < MAX_CITIES; from++) {
    for (to = 0; to < MAX_CITIES; to++) {
      phero[from][to] = INIT_PHER;
    }
  }

  // initializing the ants
  to = 0;
  for (ant = 0; ant < MAX_ANTS; ant++) {
    if (to == MAX_CITIES) {
      to = 0;
    }

    ants[ant].curCity = to++;
    for (from = 0; from < MAX_CITIES; from++) {
      ants[ant].tabu[from] = 0;
      ants[ant].path[from] = -1;
    }

    ants[ant].pathIndex = 1;
    ants[ant].path[0] = ants[ant].curCity;
    ants[ant].nextCity = -1;
    ants[ant].tourLength = 0;

    // load first city into tabu list
    ants[ant].tabu[ants[ant].curCity] = 1;
  }
}

// reinitialize all ants and redistribute them
void restartAnts() {
  for (int ant = 0; ant < MAX_ANTS; ant++) {
    ants[ant].nextCity = -1;
    ants[ant].tourLength = 0.0;

    for (int i = 0; i < MAX_CITIES; i++) {
      ants[ant].tabu[i] = 0;
      //ants[ant].path[i] = -1;
    }

    ants[ant].curCity = rand() % MAX_CITIES;
    ants[ant].pathIndex = 1;
    ants[ant].path[0] = ants[ant].curCity;
    ants[ant].tabu[ants[ant].curCity] = 1;
  }
}

float antProduct(int from, int to) {
  if (isnan((*dist)[from][to])) {
    printf("NAN (%d, %d)\n", from, to);
  }
  return (pow(phero[from][to], ALPHA) * pow((1.0 / (*dist)[from][to]), BETA));
}

int selectNextCity(int ant) {
  int from = ants[ant].curCity;
  float sum = 0.0;

  for (int to = 0; to < MAX_CITIES; to++) {
    if (ants[ant].tabu[to] == 0) {
      sum += antProduct(from, to);
    }
  }

  if (sum == 0) {
    printf("warning: zero sum in selectNextCity\n");
    return 0;
  }

  int lastBestIndex = 0;
  float acc = 0;
  float luckyNumber = (float)rand() / RAND_MAX;

  for (int to = 0; to < MAX_CITIES; to++) {
    if (ants[ant].tabu[to] == 0) {
      float product = antProduct(from, to) / sum;
      if (product > 0) {
        acc += product;
        lastBestIndex = to;

        if (acc >= luckyNumber) {
          return to;
        }
      }
    }
  }

  //if we get here (floating point errors), return last best city
  printf("warning: acc did not reach luckyNumber in selectNextCity\n");
  printf("sum: %1.15f, acc: %1.15f, luckyNumber: %1.15f\n", sum, acc, luckyNumber);
  return lastBestIndex;
}

void simulateAnts(int k) {
  //for (int k = 0; k < MAX_ANTS; k++) {
    // check if there are any more cities to visit
    
  while(ants[k].pathIndex < MAX_CITIES) {
    ants[k].nextCity = selectNextCity(k);
    ants[k].tabu[ants[k].nextCity] = 1;
    ants[k].path[ants[k].pathIndex++] = ants[k].nextCity;

    ants[k].tourLength += (*dist)[ants[k].curCity][ants[k].nextCity];

    //handle last case->last city to first
    if (ants[k].pathIndex == MAX_CITIES) {
      ants[k].tourLength += (*dist)[ants[k].path[MAX_CITIES -1]][ants[k].path[0]];
    }

    ants[k].curCity = ants[k].nextCity;
  }
}

// Updating trails
void updateTrails()
{
  int from, to, i, ant;

  // Pheromone Evaporation
  for (from = 0; from < MAX_CITIES; from++) {
    for (to = 0; to < from; to++) {
      phero[from][to] *= 1.0 - RHO;

      if (phero[from][to] < 0.0) {
        phero[from][to] = INIT_PHER;
      }
      phero[to][from] = phero[from][to];
    }
  }

  //Add new pheromone to the trails
  for (ant = 0; ant < MAX_ANTS; ant++) {
    for (i = 0; i < MAX_CITIES; i++) {
      from = ants[ant].path[i];
      if (i < MAX_CITIES - 1) {
        to = ants[ant].path[i+1];
      } else {
        to = ants[ant].path[0];
      }

      phero[from][to] += (QVAL / ants[ant].tourLength);
      phero[to][from] = phero[from][to];
    }
  }

  /*for (from = 0; from < MAX_CITIES; from++) {
    for (to = 0; to < from; to++) {
      phero[from][to] *= RHO;
      phero[to][from] *= RHO;
    }
  }*/
}

float seq_ACO(EdgeMatrix *d, int *bestPath) {
  dist = d;
  int curTime = 0;
  float pathTime = 0;
  float pheroTime = 0;
  float sBegin, sEnd;

  srand(time(NULL));

  init();

  while (curTime++ < MAX_TOURS) {
    bestIndex = -1;
    sBegin = CycleTimer::currentSeconds();
   
    for (int i = 0; i < MAX_ANTS; i++) {
      simulateAnts(i);
    }
    
    sEnd = CycleTimer::currentSeconds();
    pathTime += sEnd - sBegin;

    for (int ant = 0; ant < MAX_ANTS; ant++) {
      if (ants[ant].tourLength < best) {
        best = ants[ant].tourLength;
        printf("new best: %1.f\n", best);
        bestIndex = ant;
      }
    }

    if (bestIndex != -1) {
      memcpy(bestPath, ants[bestIndex].path, sizeof(int) * MAX_CITIES);
    }

    sBegin = CycleTimer::currentSeconds();
    updateTrails();
    sEnd = CycleTimer::currentSeconds();
    pheroTime += sEnd - sBegin;

    if (curTime != MAX_TIME) {
      restartAnts();
    }

    //cout << "\nTime is " << curTime << "(" << best << ")";
  }

  printf("PATHTIME: %f, PHEROTIME: %f\n", pathTime, pheroTime);
  return best;
}
