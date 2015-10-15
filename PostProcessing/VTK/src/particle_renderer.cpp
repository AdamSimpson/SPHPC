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


namespace mpi = boost::mpi;

void VtkParticleRenderer::PrepareBuffers(const std::string& file_name, const mpi::communicator& comm) {
  AdiosReader reader(file_name, comm);

  std::vector<double> x = reader.FetchValue<double>("x");
  std::vector<double> y = reader.FetchValue<double>("y");
  std::vector<double> z = reader.FetchValue<double>("z");

  points->Reset();
  colors->Reset();
  radii->Reset();

  for(int i=0; i<x.size(); i++) {
    points->InsertNextPoint(static_cast<float>(x[i]),
                            static_cast<float>(y[i]),
                            static_cast<float>(z[i]));
 }

  colors->SetNumberOfComponents(3);
  std::vector<double> densities = reader.FetchValue<double>("density");
  vtkSmartPointer<vtkLookupTable> lookup_table = vtkSmartPointer<vtkLookupTable>::New();
  lookup_table->SetHueRange(0.6, 0.6);
  lookup_table->SetSaturationRange(0.0, 1.0);
  lookup_table->SetTableRange(500, 800);
  lookup_table->Build();
  for(int i=0; i<points->GetNumberOfPoints(); i++) {
    double color[3];
    lookup_table->GetColor(densities[i], color);
    colors->InsertNextTuple3(static_cast<unsigned char>(color[0]*255),
                             static_cast<unsigned char>(color[1]*255),
                             static_cast<unsigned char>(color[2]*255) );
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

  // Setup obstacle poly data
  obstacle_reader->SetFileName(params.obstacle.OBJ_file_name.c_str());
  obstacle_reader->Update();
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

  // Create obstacle actor
  obstacle_mapper->SetInputConnection(obstacle_reader->GetOutputPort());
  obstacle_actor->SetMapper(obstacle_mapper);

  double obstacle_length = params.obstacle.max_x - params.obstacle.min_coord.x;
  double *x_range = obstacle_actor->GetXRange();
  double original_length = x_range[1] - x_range[0];
  double obstacle_scale = obstacle_length/original_length;
  obstacle_actor->SetScale(obstacle_scale);

  // In general OBJ model coordinates do not have (0,0,0) origin so we subtract minimum value
  x_range = obstacle_actor->GetXRange();
  double *y_range = obstacle_actor->GetYRange();
  double *z_range = obstacle_actor->GetZRange();
  obstacle_actor->SetPosition(params.obstacle.min_coord.x - x_range[0],
                              params.obstacle.min_coord.y - y_range[0],
                              params.obstacle.min_coord.z - z_range[0]);

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

std::vector<std::string> VtkParticleRenderer::GetFiles() {
  boost::filesystem::path path(params.input_dir);

  std::vector<std::string> file_names;

  for(const auto& file_name : boost::filesystem::directory_iterator(path)) {
    file_names.push_back(file_name.path().string());
  }

  std::sort(file_names.begin(), file_names.end(), [](std::string a, std::string b) {
    std::regex regex_int("\\d+"); // Match integer in filename
    std::smatch match;

    std::regex_search(a, match, regex_int);
    int a_i = std::stoi(match.str());
    std::regex_search(b, match, regex_int);
    int b_i = std::stoi(match.str());

    return a_i < b_i;
  });

  return file_names;
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

  std::vector<std::string> files = renderer.GetFiles();
  renderer.BeginAnimation();

  renderer.PreparePolyData();
  renderer.PrepareScene();

  for(const std::string& file : files) {
    std::cout<<"Processing: "<<file<<std::endl;
    renderer.PrepareBuffers(file, my_comm);
    renderer.RenderScene();
  }

  renderer.EndAnimation();

  return 0;
}
