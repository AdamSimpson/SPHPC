#include "particles.h"
#include "parameters.h"
#include <iostream>
#include "hemi/array.h"

typedef double real;

int main(int argc, char *argv[]) {
  Parameters<real> parameters;

  Particles<real, three_dimensional> particles(parameters, 100000);

  particles.ApplyGravity();

  return 0;
}
