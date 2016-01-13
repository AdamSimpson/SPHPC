#include <iostream>
#include "cuda_runtime.h"

template <typename T>
__global__ void kernel(T lam) { lam(); }

int main(int argc, char **argv) {
  int *arr;
  cudaMallocManaged(&arr, 1*sizeof(int));

  auto body = void { arr[0] = 7; };
  auto lambda = [=] __device__ { arr[0] = 7; };

  kernel<<<1,1>>>(lambda);
  lambda();

  cudaDeviceSynchronize();
  std::cout<<arr[0]<<std::endl;

  return 0;
}
