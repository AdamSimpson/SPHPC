#include "simulation.h"
#include "dimension.h"
#include "exception"

int main(int argc, char *argv[]) {
  try {
    Simulation<float, two_dimensional> simulation("params.ini");

  } catch(std::exception const& exception) {
      std::cout << "Aborting: " << exception.what() << std::endl;
      return 1;
  }

  return 0;
}
