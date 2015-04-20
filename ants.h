// Traveling Antsmen constants and struct declarations

#define MAX_CITIES 100 // eventually test on 250-400
#define MAX_DIST 100
#define MAX_TOUR (MAX_CITIES * MAX_DIST)
#define MAX_ANTS MAX_CITIES

struct cityType {
  int x, y;
};

struct antType {
  int curCity, nextCity, pathIndex;
  int tabu[MAX_CITIES];
  int path[MAX_CITIES];
  double tourLength;
};

//Ant algorithm problem parameters
#define ALPHA 1.0
#define BETA 5.0 // this parameter raises the weight of distance over pheromone
#define RHO 0.5 // evaporation rate
#define QVAL 100
#define MAX_TOURS 20
#define MAX_TIME (MAX_TOURS * MAX_CITIES)
#define INIT_PHER (1.0 / MAX_CITIES)