#include "particle_renderer.h"
#include "parameters.h"
#include "adios_reader.h"

#include <vector>
#include "boost/mpi.hpp"

#include <vtkSmartPointer.h>
#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkPolyData.h>
#include <vtkOpenGLSphereMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkCamera.h>
//#include <vtkWindowToImageFilter.h>
//#include <vtkPNGWriter.h>
//#include <vtkCamera.h>

// Needed to use OpenGL2 sphere imposters
#include <map>
#include <vtkShader.h>
#include <vtkOpenGLHelper.h>
#include <vtkOpenGLSphereMapper.h>
//

namespace mpi = boost::mpi;

void VtkParticleRenderer::PrepareBuffers(const std::string& file_name, const mpi::communicator& comm) {
  AdiosReader reader(file_name, comm);

  std::vector<double> positions = reader.FetchValue<double>("positions");
  for(int i=0; i<positions.size(); i+=3) {
    points->InsertNextPoint(static_cast<float>(positions[i]),
                                 static_cast<float>(positions[i+1]),
                                 static_cast<float>(positions[i+2]) );
  }

  colors->SetNumberOfComponents(3);
  std::vector<double> densities = reader.FetchValue<double>("densities");
  for(const double& density : densities) {
    colors->InsertNextTuple3(255, 0, 0);
  }

  double radius = reader.FetchValue<double>("particle_radius")[0];
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
  // Create OpenGL sphere imposter mapper
  fluid_mapper->SetInputData(fluid_polydata);
  fluid_mapper->SetScalarModeToUsePointFieldData();
  colors->SetName("colors");
  fluid_mapper->SelectColorArray("colors");
  radii->SetName("radii");
  fluid_mapper->SetScaleArray("radii");
  fluid_actor->SetMapper(fluid_mapper);

  // Add actor to new render window
  render_window->AddRenderer(renderer);
  renderer->AddActor(fluid_actor);
  renderer->SetBackground(.3, .6, .3);

  // Setup Camera
  camera = vtkSmartPointer<vtkCamera>::New();
  camera->SetViewUp(0,0,1); // Set z as up to match simulation
  camera->SetPosition(30, -150, 15);
  camera->SetFocalPoint(30, 15, 30);
  renderer->SetActiveCamera(camera);
}

void VtkParticleRenderer::RenderScene() {
  render_window->Render();
}

/*
void BeginAnimation() {

}

void EndEnimation() {

}
*/

int main (int argc, char **argv) {

  mpi::environment env;
  mpi::communicator world;

  // Create communicator for each rank
  mpi::communicator my_comm = world.split(world.rank());

  Parameters params("./post-config.ini");
  params.ReadParameters();

  VtkParticleRenderer renderer;
  renderer.PrepareBuffers(params.input_file_path, my_comm);
  renderer.PreparePolyData();
  renderer.PrepareScene();
  renderer.RenderScene();

/*
  // Save image
  vtkSmartPointer<vtkWindowToImageFilter> windowToImageFilter =
    vtkSmartPointer<vtkWindowToImageFilter>::New();
  windowToImageFilter->SetInput(renderWindow);
  windowToImageFilter->SetMagnification(3); //set the resolution of the output image (3 times the current resolution of vtk render window)
  windowToImageFilter->SetInputBufferTypeToRGBA(); //also record the alpha (transparency) channel
  windowToImageFilter->ReadFrontBufferOff(); // read from the back buffer
  windowToImageFilter->Update();

  vtkSmartPointer<vtkPNGWriter> writer = vtkSmartPointer<vtkPNGWriter>::New();
  writer->SetFileName("screenshot2.png");
  writer->SetInputConnection(windowToImageFilter->GetOutputPort());
  writer->Write();
*/

  return 0;
}
