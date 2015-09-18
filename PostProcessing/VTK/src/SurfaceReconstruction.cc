#include <vtkVersion.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <iostream>
#include <fstream>

int main ( int, char *[] )
{
  //Create points array
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

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

      for(int j=0; j<num_particles_read; j++)
        points->InsertNextPoint(((double*)buffer)+(num_components*j));
    }
  }
  file.close();
  delete[] buffer;

  // Create a polydata object and add the points to it.
  vtkSmartPointer<vtkPolyData> fluid_polydata = vtkSmartPointer<vtkPolyData>::New();
  fluid_polydata->SetPoints(points);

  // Create OpenGL sphere imposter mapper
  vtkSmartPointer<vtkOpenGLSphereMapper> fluid_mapper = vtkSmartPointer<vtkOpenGLSphereMapper>::New();
  fluid_mapper->SetInputData(fluid_polydata);

  // Create fluid particles actor
  vtkSmartPointer<vtkActor> fluid_actor = vtkSmartPointer<vtkActor>::New();
  fluid_actor->SetMapper(fluid_mapper);

  // Add renderer and render window
  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
//  renderWindow->SetOffScreenRendering( 1 );
  renderWindow->AddRenderer(renderer);

  // Render fluid actor
  renderer->AddActor(fluid_actor);
  renderer->SetBackground(1,1,1);

  renderWindow->Render();

/*
  // Write the VTP file
  vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName("sim-1.vtp");
  writer->SetInputData(polydata);
  writer->SetDataModeToBinary();
  writer->Write();
*/

  return 0;
}
