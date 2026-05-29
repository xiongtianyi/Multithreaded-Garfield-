#ifndef G_RANDOM_ENGINE_GPU_H
#define G_RANDOM_ENGINE_GPU_H

#include "Garfield/GPUInterface.hh"
#include <curand.h>
#include <curand_kernel.h>

namespace Garfield {

class RandomEngineGPU {
 public:
  /// Constructor
  RandomEngineGPU();
  /// Destructor
  ~RandomEngineGPU();
  /// Call the random number generator.
  __device__ GPUFLOAT Draw();
  void setRandomEngineOnDevice();

  // cuRAND specifics
  curandState* d_curand_states{nullptr};
  double initCURandStates(const unsigned int seed);
};

}

#endif
