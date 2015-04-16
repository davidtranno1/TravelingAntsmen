#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime.h>

__global__ void kernelSmallComputePixels() {
    __shared__ uint circleDoesIntersect[2];
    __syncthreads();
}

void cuda_ACO(void) {
  printf("HIII\n");
  return;
}
