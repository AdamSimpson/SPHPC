/**
template<typename TReal, int TDim>
**/

public:
Particles(const Parameters<TReal>& params,
          const int max_local) : parameters{params},
                                 max_count_local{max_local},
                                 pos{hemi::Array< Vector<TReal, TDim> >(max_count_local)},
                                 v{hemi::Array< Vector<TReal, TDim> >(max_count_local)},
                                 density{hemi::Array<TReal>(max_count_local)} {

  // Initialize particle values
  for (int i=0; i<max_count_local; ++i) {
    v.writeOnlyPtr(hemi::host)[i] = Vector<TReal, TDim>(0.0);
    pos.writeOnlyPtr(hemi::host)[i] = Vector<TReal, TDim>(7.0);
    density.writeOnlyPtr(hemi::host)[i] = 0.0;
  }
}

void ApplyGravity() {
  const TReal g = parameters.GetGravity();
  const TReal dt = parameters.GetDt();
  auto *const v = this->v.ptr();

  hemi::parallel_for(0, max_count_local, [=] HEMI_LAMBDA (int i) {
      v[i].y += g*dt;
  });

  checkCudaErrors();

  std::cout<<this->v.readOnlyPtr(hemi::host)[100].y<<std::endl;

}
