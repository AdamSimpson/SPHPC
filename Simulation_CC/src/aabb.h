#include "dimension.h"
#include "vec.h"

template<typename Real, Dimension Dim>
struct AABB {
  Vec<Real,Dim> min;
  Vec<Real,Dim> max;
};
