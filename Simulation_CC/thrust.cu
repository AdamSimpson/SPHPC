#include <thrust/for_each.h>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/iterator/counting_iterator.h>
#include <cstdio>

/*
typedef thrust::tuple<double&, double&> Particle;
__host__ __device__
double position(const Particle& p) {
  return thrust::get<0>(p);
}

__host__ __device__
double velocity(const Particle& p){
  return thrust::get<1>(p);
}

template<typename Iterator>
void func(Iterator begin, Iterator end) {
  auto lam = [=] __device__(Particle p) {
    double p_position, p_vel;
    thrust::tie(p_position, p_vel) = p;
    p_position = 7.0;
    printf("%f, %f\n", p_position, p_vel); 
  };

  thrust::for_each(begin, end, lam);
}
*/

int main() {
  thrust::host_vector<double> position(3);
  thrust::host_vector<double> velocity(3);
  thrust::device_vector<double> d_position(3);
  thrust::device_vector<double> d_velocity(3);
 
  position[0] = 0; position[1] = 1; position[2] = 2;
  velocity[0] = 3.0; velocity[1] = 4.0; velocity[2] = 5.0;

  d_position = position;
  d_velocity = velocity;

  if(thrust::host == thrust::device)
    std::cout<<"stuff"<<std::endl;

/*
  auto particle = thrust::make_zip_iterator(thrust::make_tuple(position.begin(), velocity.begin()));
  auto d_particle = thrust::make_zip_iterator(thrust::make_tuple(d_position.begin(), d_velocity.begin()));

  func(d_particle, d_particle + 3);
*/
  return 0;
}
