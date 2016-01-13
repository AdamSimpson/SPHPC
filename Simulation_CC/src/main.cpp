#include <exception>
#include <iostream>
#include "simulation.h"

int main(int argc, char *argv[]) {
  try {
    Simulation<float, three_dimensional> simulation("params.ini");

    simulation.write_particles();

    for(std::size_t i=0; i<simulation.time_step_count(); ++i) {
      simulation.advance_particles();

      if(i%4 == 0)
        simulation.write_particles();
    }

  } catch(std::exception const& exception) {
      std::cout << "Aborting: " << exception.what() << std::endl;
      return 1;
  } catch(...) {
    std::cout << "Aborting: unknown exception" <<std::endl;
    return 1;
  }

  return 0;
}
