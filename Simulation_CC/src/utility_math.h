#pragma once

namespace Utility {

template<typename Util_T>
Util_T Clamp(const Util_T n, const Util_T lower, const Util_T upper) {
  return n < lower ? lower : n > upper ? upper : n;
}

template<typename Util_T>
void ClampInPlace(Util_T& n, const Util_T lower, const Util_T upper) {
  if(n < lower)
    n = lower;
  else if(n > upper)
    n = upper;
}

}
