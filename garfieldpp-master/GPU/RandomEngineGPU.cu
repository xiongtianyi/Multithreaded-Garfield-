#include "RandomEngineGPU.h"
#include "GPUFunctions.h"
#include <iostream>

namespace Garfield {

    static constexpr double Pi = 3.1415926535897932384626433832795;
    static constexpr double TwoPi = 2. * Pi;

    __device__ RandomEngineGPU *randomEngineGPU{nullptr};

    __global__ void setRandomEngine_d(RandomEngineGPU *engine)
    {
        randomEngineGPU = engine;
    }

    __device__ GPUFLOAT RndmUniformGPU()
    {
        return randomEngineGPU->Draw();
    }

    __device__ GPUFLOAT RndmUniformPosGPU() {
        GPUFLOAT r = RndmUniformGPU();
        while (r <= 0.) r = RndmUniformGPU();
        return r;
    }

    __device__ void RndmDirectionGPU(GPUFLOAT& dx, GPUFLOAT& dy, GPUFLOAT& dz,
        const GPUFLOAT length) {
        const GPUFLOAT phi = TwoPi * RndmUniformGPU();
        const GPUFLOAT ctheta = 2 * RndmUniformGPU() - 1.;
        const GPUFLOAT stheta = sqrt(1. - ctheta * ctheta);
        dx = length * cos(phi) * stheta;
        dy = length * sin(phi) * stheta;
        dz = length * ctheta;
    }

    __global__ void initCURandStates_d(curandState *state, const unsigned int seed) {
        int tid = threadIdx.x + blockIdx.x * blockDim.x;
        if (tid < MAXSTACKSIZE) curand_init(seed, tid, 0, &state[tid]);
    }

    RandomEngineGPU::RandomEngineGPU(){
    }

    RandomEngineGPU::~RandomEngineGPU()
    {
    }

    double RandomEngineGPU::initCURandStates(const unsigned int seed)
    {
        checkCudaErrors( cudaMalloc( &d_curand_states, MAXSTACKSIZE * sizeof(curandState) ) );
        initCURandStates_d<<< 1 + MAXSTACKSIZE/256, 256 >>> (d_curand_states, (seed == 0 ? time(NULL) : seed));
        return MAXSTACKSIZE * sizeof(curandState);
    }

    void RandomEngineGPU::setRandomEngineOnDevice()
    {
        setRandomEngine_d<<<1,1>>>(this);
    }

    __device__ GPUFLOAT RandomEngineGPU::Draw()
    {
        unsigned int tid = (threadIdx.x + blockIdx.x * blockDim.x);
        return curand_uniform(&(d_curand_states[tid]));
    }

}
