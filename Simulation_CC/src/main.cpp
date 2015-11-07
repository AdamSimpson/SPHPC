#include "simulation.h"
#include "dimension.h"

int main(int argc, char *argv[]) {
  Simulation<float, int, two_dimensional> simulation("params.ini");

  return 0;
}
