#ifndef VTK_SRC_PARTICLE_RENDERER_H_
#define VTK_SRC_PARTICLE_RENDERER_H_

#include <vtkSmartPointer.h>
#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkFloatArray.h>
#include <vtkPolyData.h>
#include <vtkOpenGLSphereMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkCamera.h>


#include "boost/mpi.hpp"


class VtkParticleRenderer {
  public:
    VtkParticleRenderer() : points{vtkSmartPointer<vtkPoints>::New()},
                            radii{vtkSmartPointer<vtkFloatArray>::New()},
                            colors{vtkSmartPointer<vtkUnsignedCharArray>::New()},
                            fluid_polydata{vtkSmartPointer<vtkPolyData>::New()},
                            fluid_mapper{vtkSmartPointer<vtkOpenGLSphereMapper>::New()},
                            fluid_actor{vtkSmartPointer<vtkActor>::New()},
                            renderer{vtkSmartPointer<vtkRenderer>::New()},
                            render_window{vtkSmartPointer<vtkRenderWindow>::New()},
                            camera{vtkSmartPointer<vtkCamera>::New()} {};
    void PrepareBuffers(const std::string& file_name, const boost::mpi::communicator& comm);
    void PreparePolyData();
    void PrepareScene();
    void RenderScene();

  private:
    vtkSmartPointer<vtkPoints> points;
    vtkSmartPointer<vtkFloatArray> radii;
    vtkSmartPointer<vtkUnsignedCharArray> colors;
    vtkSmartPointer<vtkPolyData> fluid_polydata;
    vtkSmartPointer<vtkOpenGLSphereMapper> fluid_mapper;
    vtkSmartPointer<vtkActor> fluid_actor;
    vtkSmartPointer<vtkRenderer> renderer;
    vtkSmartPointer<vtkRenderWindow> render_window;
    vtkSmartPointer<vtkCamera> camera;
};

#endif
