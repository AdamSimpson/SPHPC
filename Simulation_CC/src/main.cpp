#include "simulation.h"
#include "dimension.h"
#include "exception"

int main(int argc, char *argv[]) {
  try {
    Simulation<float, two_dimensional> simulation("params.ini");

    simulation.WriteParticles();

    for(std::size_t i=0; i<simulation.GetTimeStepCount(); ++i) {
      simulation.AdvanceParticles();

      simulation.WriteParticles();
    }

  } catch(std::exception const& exception) {
      std::cout << "Aborting: " << exception.what() << std::endl;
      return 1;
  }

  return 0;
}
