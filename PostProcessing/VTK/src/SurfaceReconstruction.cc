#include <vtkSmartPointer.h>
#include <vtkNew.h>
#include <vtkVersion.h>
#include <vtkCellArray.h>
#include <vtkFloatArray.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include "vtkLookupTable.h"
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkWindowToImageFilter.h>
#include <vtkPNGWriter.h>
#include <vtkCamera.h>

// Needed to use OpenGL2 sphere imposters
#include <map>
#include <vtkShader.h>
#include <vtkOpenGLHelper.h>
#include <vtkOpenGLSphereMapper.h>
//

#include <iostream>
#include <fstream>

int main (int argc, char **argv) {
  //Create points array
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

  // Create radii array
  vtkSmartPointer<vtkFloatArray> radii = vtkSmartPointer<vtkFloatArray>::New();
  radii->SetName("radii");

  // Create color array
  vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
  colors->SetName("colors");
  colors->SetNumberOfComponents(3);

  // Setup file to process
  // Place cursor at end of file to easily get size
  ifstream file ("sim-1.bin", ios::in|ios::binary|ios::ate);

  std::streampos size;
  size = file.tellg();
  int num_components = 3;
  int num_particles = size/(num_components*sizeof(double));

  std::cout<<"number of particles: "<<num_particles<<std::endl;

  // Allocate almost 1 GB for buffer
  // Needs to be divisible by num_components*8 bytes so xyz double coords fits
  unsigned int buffer_bytes = 1073741808;

  // Number of bytes read in
  std::streamsize bytes_read;
  // Number of particles read in
  unsigned int num_particles_read;
  // Total particles read in so far
  unsigned int total_particles_read = 0;

  char *buffer = new char[buffer_bytes];
  if(file.is_open()) {
    // Go to begining of file
    file.seekg (0, ios::beg);

    // Read in up to buffer_bytes and write to points
    while(!file.eof()) {
      file.read(buffer, buffer_bytes);
      std::streamsize bytes_read = file.gcount();
      num_particles_read = bytes_read/(num_components*sizeof(double));
      total_particles_read += num_particles_read;

      for(int j=0; j<num_particles_read; j++) {
        double* buffer_pos = (double*)(buffer+(num_components*j*sizeof(double)));

        // Convert to float
        float float_buffer[3];
        float_buffer[0] = (float)buffer_pos[0];
        float_buffer[1] = (float)buffer_pos[1];
        float_buffer[2] = (float)buffer_pos[2];

        points->InsertNextPoint((float_buffer));
        radii->InsertNextValue(0.35);
        colors->InsertNextTuple3(255, 0, 0);
      }
    }
  }
  file.close();
  delete[] buffer;

  // Create polydata object and add the points to it.
  vtkSmartPointer<vtkPolyData> fluid_polydata = vtkSmartPointer<vtkPolyData>::New();
  fluid_polydata->SetPoints(points);

  // Attach radii and colors to points
  fluid_polydata->GetPointData()->AddArray(radii);
  fluid_polydata->GetPointData()->AddArray(colors);

  // Create OpenGL sphere imposter mapper
  vtkSmartPointer<vtkOpenGLSphereMapper> fluid_mapper = vtkSmartPointer<vtkOpenGLSphereMapper>::New();
  fluid_mapper->SetInputData(fluid_polydata);

  // Set to use colors array for coloring and radii array for radius
  fluid_mapper->SetScalarModeToUsePointFieldData();
  fluid_mapper->SelectColorArray("colors");
  fluid_mapper->SetScaleArray("radii");

  // Create fluid particles actor
  vtkSmartPointer<vtkActor> fluid_actor = vtkSmartPointer<vtkActor>::New();
  fluid_actor->SetMapper(fluid_mapper);

  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);
 
  renderer->AddActor(fluid_actor);
  renderer->SetBackground(.3, .6, .3); // Background color green

  // Setup Camera
  vtkSmartPointer<vtkCamera> camera = vtkSmartPointer<vtkCamera>::New();
  camera->SetViewUp(0,0,1); // Set z as up to match simulation
  camera->SetPosition(30, -150, 15);
  camera->SetFocalPoint(30, 15, 30);
  renderer->SetActiveCamera(camera);

  renderWindow->Render();

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

  return 0;
}
