#ifndef SPH_SRC_UTILITY_MATH_H_
#define SPH_SRC_UTILITY_MATH_H_

template<typename TUtil>
HEMI_DEV_CALLABLE_INLINE
TUtil Clamp(TUtil n, TUtil lower, TUtil upper) {
  return n <= lower ? lower : n >= upper ? upper : n;
}

#endif
