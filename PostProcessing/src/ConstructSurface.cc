#include <iostream>
#include <vector>
#include "openvdb/openvdb.h"
#include "openvdb/tools/ParticlesToLevelSet.h"
#include "openvdb/tools/LevelSetFilter.h"
#include "openvdb/Types.h"

class Particles {
  private:
    std::vector<openvdb::Vec3R> particles;
    openvdb::Real radius;
  public:
    Particles(const openvdb::Real &radius) : radius(radius) {};
    void add(const openvdb::Vec3R &coord) {this->particles.push_back(coord); }
    void getPos(size_t n, openvdb::Vec3R &pos) const { pos = particles[n]; }
    size_t size() const { return particles.size(); }
};

int main(int argc, char **argv) {
  openvdb::initialize();

  const openvdb::Real radius = 1.0;
  Particles particles(radius);

  // Read data from file
  const int particle_count = 241776;
  std::ifstream coordinate_file;
  coordinate_file.open("./sim-0.bin", std::ios::binary);
  for(int i=0; i<particle_count; i++) {
    char coord[3*sizeof(double)];
    coordinate_file.read(coord, sizeof(coord));
    particles.add(openvdb::Vec3R((double*)coord));
  }

  const openvdb::Real particle_radius = 1.0f;

  // Create empty level set grid
  const float voxelSize = particle_radius/2.0f; // World units
  const float halfWidth = 2.0f;                 // Voxel units
  openvdb::FloatGrid::Ptr level_set = openvdb::createLevelSet<openvdb::FloatGrid>
                                        (voxelSize, halfWidth);

  // Create level set from rasterized particles
  openvdb::tools::ParticlesToLevelSet<openvdb::FloatGrid> raster(*level_set);
  raster.setGrainSize(1); //a value of zero disables threading
  raster.rasterizeSpheres(particles, particle_radius);
  raster.finalize();

  // Create filter for level set
  openvdb::tools::LevelSetFilter<openvdb::FloatGrid> filter(*level_set);

  // Dilate level set
  filter.offset(-particle_radius);

  // Smooth level set
  filter.gaussian();
//  filter.laplacian();
  openvdb::FloatGrid::Ptr grid = level_set;

  // Erode level set
  filter.offset(particle_radius);

  // Save mesh
  grid->setName("ParticlesToLevelSet");
  openvdb::GridPtrVec grids;
  grids.push_back(grid);
  std::string fileName("sph");
  openvdb::io::File file(fileName + ".vdb");
  file.write(grids);
  file.close();

  openvdb::uninitialize();

  return 0;
}
