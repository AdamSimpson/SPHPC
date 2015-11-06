#ifndef SPH_SRC_PARAMETERS_H_
#define SPH_SRC_PARAMETERS_H_

template <typename TReal>
class Parameters {
  public:
    HEMI_DEV_CALLABLE_INLINE_MEMBER
    Parameters() : gravity{static_cast<TReal>(9.8)},
                   dt{static_cast<TReal>(0.01)} {};

  private:
    TReal gravity;
    TReal dt;

  #include "parameters.hxx"
};

#endif
