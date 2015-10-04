#ifndef VTK_SRC_PARAMETERS_H_
#define VTK_SRC_PARAMETERS_H_

#include <string>

struct double3 {
  double x;
  double y;
  double z;
};
double3 ToDouble3(const std::string input_string);

struct Camera {
  double3 view_up;
  double3 position;
  double3 focal_point;
};

struct Boundary {
  double3 min_coord;
  double3 max_coord;
};

struct Obstacle {
  std::string OBJ_file_name;
  std::string MTL_file_name;
  double3 min_coord;
  double max_x;
};

class Parameters {
  public:
    Parameters(const std::string ini_name) : ini_file_name{ini_name} {};
    void ReadParameters();
    Camera camera;
    Boundary boundary;
    Obstacle obstacle;
    std::string input_file_path;

  private:
    std::string ini_file_name;
};

#endif
