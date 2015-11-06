#include "utility_math.h"

/**
template <typename TReal>
**/
public:
HEMI_DEV_CALLABLE_INLINE_MEMBER
const TReal GetGravity() const { return gravity; }

HEMI_DEV_CALLABLE_INLINE_MEMBER
void SetGravity(const TReal new_gravity) {
  const TReal min = static_cast<TReal>(-9.8);
  const TReal max = static_cast<TReal>(9.8);
  gravity = Clamp<TReal>(new_gravity, min, max);
}

HEMI_DEV_CALLABLE_INLINE_MEMBER
const TReal GetDt() const { return dt; }

HEMI_DEV_CALLABLE_INLINE_MEMBER
void setDt(const TReal new_dt) { dt = new_dt; }
