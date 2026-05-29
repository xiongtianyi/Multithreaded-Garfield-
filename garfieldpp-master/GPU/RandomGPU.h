#ifndef G_RANDOMGPU_H
#define G_RANDOMGPU_H

#include "Garfield/GPUInterface.hh"

namespace Garfield {

  __device__ GPUFLOAT RndmUniformGPU();
  __device__ GPUFLOAT RndmUniformPosGPU();
  __device__ void RndmDirectionGPU(GPUFLOAT& dx, GPUFLOAT& dy, GPUFLOAT& dz, const GPUFLOAT length = 1.);

}

#endif
