// Traveling Antsmen constants and struct declarations
#include <algorithm>

#define MAX_CITIES 300 // eventually test on 250-800
#define MAX_DIST 100
#define MAX_TOUR (MAX_CITIES * MAX_DIST)
#define MAX_ANTS MAX_CITIES
#define NUM_EDGES ((MAX_CITIES * MAX_CITIES - MAX_CITIES) / 2)


struct cityType {
  int x, y;
};

// allows edge queries as a 2D array
class EdgeMatrix {
  float *dist;
public:
  EdgeMatrix() {
    dist = new float[NUM_EDGES];
  }
  ~EdgeMatrix() {
    delete dist;
  }

  /* Edges are stored as an array of size (n^2 + n) / 2, with indices:
   * 0
   * 1 2
   * 3 4 5
   * 6 7 8 9
   */
  float edge(unsigned int x, unsigned int y) {
    unsigned int i = std::max(x, y) - 1;
    unsigned int j = std::min(x, y);
    return dist[(i * (i + 1) / 2) + j];
  }
  void set_edge(unsigned int x, unsigned int y, float value) {
    unsigned int i = std::max(x, y) - 1;
    unsigned int j = std::min(x, y);
    dist[(i * (i + 1) / 2) + j] = value;
  }

  float *get_array(){
    return dist;
  }
};

//Ant algorithm problem parameters
#define ALPHA 1.0
#define BETA 5.0 // this parameter raises the weight of distance over pheromone
#define RHO 0.5 // evaporation rate
#define QVAL 100
#define MAX_TOURS 50
#define MAX_TIME (MAX_TOURS * MAX_CITIES)
#define INIT_PHER (1.0 / MAX_CITIES)
