#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdbool.h>
#include "safe_alloc.h"
#include "print_macros.h"
#include "geometry.h"
#include "obstacle_cuda.h"
#include "obstacle.h"

// SDF file assumed to be ASCII format produced by the following command line
// tool: https://github.com/christopherbatty/SDFGen
void AllocInitObstacle(struct Obstacle *obstacle) {
  // Open SDF file
  FILE *const file = fopen(obstacle->SDF_file_name, "r");
  if(!file)
    EXIT_PRINT("Can't open %s\n", obstacle->SDF_file_name);

  // Extract grid dimensions
  char *line = NULL;
  size_t length = 0;
  ssize_t char_count = getline(&line, &length, file);
  if(char_count == -1)
      EXIT_PRINT("Error parsing SDF file: %s\n", obstacle->SDF_file_name);
  int num_params = sscanf(line, "%d %d %d",
                                &(obstacle->SDF_dim_x),
                                &(obstacle->SDF_dim_y),
                                &(obstacle->SDF_dim_z));
  if(num_params == EOF)
    EXIT_PRINT("Error parsing SDF file: %s\n", obstacle->SDF_file_name);

  // Extract origin coordinate
  char_count = getline(&line, &length, file);
  if(char_count == -1)
      EXIT_PRINT("Error parsing SDF file: %s\n", obstacle->SDF_file_name);
  num_params = sscanf(line, "%lf %lf %lf",
                                &(obstacle->origin_x),
                                &(obstacle->origin_y),
                                &(obstacle->origin_z));
  if(num_params == EOF)
    EXIT_PRINT("Error parsing SDF file: %s\n", obstacle->SDF_file_name);

  // Extract spacing
  char_count = getline(&line, &length, file);
  if(char_count == -1)
      EXIT_PRINT("Error parsing SDF file: %s\n", obstacle->SDF_file_name);
  num_params = sscanf(line, "%lf", &(obstacle->SDF_spacing));
  if(num_params == EOF)
    EXIT_PRINT("Error parsing SDF file: %s\n", obstacle->SDF_file_name);

  // Allocate space for SDF data
  int SDF_grid_count = obstacle->SDF_dim_x *
                       obstacle->SDF_dim_y *
                       obstacle->SDF_dim_z;
  obstacle->SDF_buffer = SAFE_ALLOC(SDF_grid_count ,sizeof(float));

  // Read SDF distance values
  for(int i=0; i<SDF_grid_count; i++) {
    size_t length = 0;
    const ssize_t char_count = getline(&line, &length, file);
    if(char_count == -1)
      EXIT_PRINT("Error parsing SDF file: %s\n", obstacle->SDF_file_name);
    const int num_params = sscanf(line, "%f", &(obstacle->SDF_buffer[i]));
    if(num_params == EOF)
      EXIT_PRINT("Error parsing SDF file: %s\n", obstacle->SDF_file_name);
  }

  free(line);
  fclose(file);

  // Set obstacle bounds to the "correct" aspect ratio
  // Note this may be off by a little as we're using the SDF grid dimensions
  double obstacle_length = obstacle->world_bounds.max_x - obstacle->world_bounds.min_x;
  obstacle->world_bounds.max_y = obstacle_length*obstacle->SDF_dim_y/(float)obstacle->SDF_dim_x;
  obstacle->world_bounds.max_z = obstacle_length*obstacle->SDF_dim_z/(float)obstacle->SDF_dim_x;

  // Initialize CUDA obstacle implimentation
  AllocInitObstacle_CUDA(obstacle);
}

void FinalizeObstacle(struct Obstacle *obstacle) {
  FinalizeObstacle_CUDA(&(obstacle->obstacle_cuda));
  free(obstacle->SDF_buffer);
}
