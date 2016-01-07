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

class Parameters {
  public:
    Parameters(const std::string ini_name) : ini_file_name{ini_name} {};
    void ReadParameters();
    Camera camera;
    Boundary boundary;
    std::string bp_file_name;

  private:
    std::string ini_file_name;
};

#endif
