#include "particle_renderer.h"

#include "parameters.h"
#include "adios_reader.h"

#include "boost/mpi.hpp"
#include "boost/filesystem.hpp"
#include <vector>
#include <algorithm>
#include <regex>

#include <vtkSmartPointer.h>
#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkOBJReader.h>
#include <vtkFloatArray.h>
#include <vtkPolyData.h>
#include <vtkOpenGLSphereMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkWindowToImageFilter.h>
#include <vtkCamera.h>
#include <vtkOggTheoraWriter.h>
#include <vtkTransform.h>
#include <vtkLookupTable.h>

// Needed to use OpenGL2 sphere imposters
#include <map>
#include <vtkShader.h>
#include <vtkOpenGLHelper.h>
#include <vtkOpenGLSphereMapper.h>
//

#include "../../../src/vec.h"

namespace mpi = boost::mpi;

void VtkParticleRenderer::PrepareBuffers(std::size_t step, const mpi::communicator& comm) {
  AdiosReader reader(params.bp_file_name, comm);

  std::vector< Vec<float,2> > positions = reader.FetchValue< Vec<float,2> >("positions", step);

  points->Reset();
  colors->Reset();
  radii->Reset();

  const int dim = 2;
  for(int i=0; i<positions.size(); ++i) {
    points->InsertNextPoint(static_cast<float>(positions[i].x),
                            static_cast<float>(0.0),
                            static_cast<float>(-1.0*positions[i].y + 100.0));
  }

  colors->SetNumberOfComponents(3);
  for(int i=0; i<points->GetNumberOfPoints(); i++) {
    colors->InsertNextTuple3(static_cast<unsigned char>(0),
                             static_cast<unsigned char>(0),
                             static_cast<unsigned char>(255) );
  }

  double radius = 0.5;
  for(int i=0; i<points->GetNumberOfPoints(); i++) {
    radii->InsertNextValue(static_cast<float>(radius));
  }

}

void VtkParticleRenderer::PreparePolyData() {
  // Create polydata object and add the points to it.
  fluid_polydata->SetPoints(points);

  // Attach radii and colors to points
  fluid_polydata->GetPointData()->AddArray(radii);
  fluid_polydata->GetPointData()->AddArray(colors);
}

void VtkParticleRenderer::PrepareScene() {
  // Create fluid particle actor
  fluid_mapper->SetInputData(fluid_polydata);
  fluid_mapper->SetScalarModeToUsePointFieldData();
  colors->SetName("colors");
  fluid_mapper->SelectColorArray("colors");
  radii->SetName("radii");
  fluid_mapper->SetScaleArray("radii");
  fluid_actor->SetMapper(fluid_mapper);

  // Enable offscreen rendering
  render_window->SetOffScreenRendering(1);

  // Add actor to new render window
  render_window->AddRenderer(renderer);
  renderer->AddActor(fluid_actor);
  renderer->AddActor(obstacle_actor);

  renderer->GradientBackgroundOn();
  renderer->SetBackground(82/255.0,82/255.0,110/255.0);
  renderer->SetBackground2(0,0,42/255.0);

  // Setup Camera
  double3 view_up = params.camera.view_up;
  camera->SetViewUp(view_up.x, view_up.y, view_up.z); // Set z as up to match simulation
  double3 position = params.camera.position;
  camera->SetPosition(position.x, position.y, position.z);
  double3 focal_point = params.camera.focal_point;
  camera->SetFocalPoint(focal_point.x, focal_point.y, focal_point.z);
  renderer->SetActiveCamera(camera);

  window_to_image_filter->SetInput(render_window);
  window_to_image_filter->SetMagnification(3); //set the resolution of the output image (3 times the current resolution of vtk render window)
  window_to_image_filter->ReadFrontBufferOff(); // read from the back buffer
}

void VtkParticleRenderer::RenderScene() {
  render_window->Render();
  window_to_image_filter->Modified();
  window_to_image_filter->Update();
  ogg_writer->Write();
}

void VtkParticleRenderer::BeginAnimation() {
  ogg_writer->SetInputConnection(window_to_image_filter->GetOutputPort());
  ogg_writer->SetFileName("cool.ogg");
  ogg_writer->SetRate(30);
  ogg_writer->Start();
}

void VtkParticleRenderer::EndAnimation() {
  ogg_writer->End();
}

void VtkParticleRenderer::ReadParameters() {
  params.ReadParameters();
}

int main (int argc, char **argv) {

  mpi::environment env;
  mpi::communicator world;

  // Create communicator for each rank
  mpi::communicator my_comm = world.split(world.rank());

  VtkParticleRenderer renderer;
  renderer.ReadParameters();

  renderer.BeginAnimation();

  renderer.PreparePolyData();
  renderer.PrepareScene();

  const int num_steps = 100;
  for(int step=0; step<num_steps; step++) {
    std::cout<<"Processing step "<<step<<std::endl;
    renderer.PrepareBuffers(step, my_comm);
    renderer.RenderScene();
  }

  renderer.EndAnimation();

  return 0;
}
