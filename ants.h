// Traveling Antsmen constants and struct declarations

#define MAX_CITIES 100 // eventually test on 250-800
#define MAX_DIST 100
#define MAX_TOUR (MAX_CITIES * MAX_DIST)
#define MAX_ANTS MAX_CITIES

struct cityType {
  int x, y;
};

// allows edge queries as a 2D array
class EdgeMatrix {
  double *dist;
public:
  EdgeMatrix() {
    dist = new double[MAX_CITIES * MAX_CITIES];
  }
  ~EdgeMatrix() {
    delete dist;
  }
  double* operator[](unsigned int i) {
    return &dist[MAX_CITIES * i];
  }
  
  double *get_array(){
    return dist;
  }
};

//Ant algorithm problem parameters
#define ALPHA 1.0
#define BETA 5.0 // this parameter raises the weight of distance over pheromone
#define RHO 0.5 // evaporation rate
#define QVAL 100
#define MAX_TOURS 20
#define MAX_TIME (MAX_TOURS * MAX_CITIES)
#define INIT_PHER (1.0 / MAX_CITIES)
