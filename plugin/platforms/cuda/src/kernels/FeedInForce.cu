
extern "C" __global__ void computeForceKernel(
unsigned long long* __restrict__ forceBuffer,
const double* __restrict__ forces_in,
real* __restrict__ energyBuffer
) {
    int threadIndex = threadIdx.x;
    for (int index=blockIdx.x * blockDim.x + threadIndex; index<NUM_ATOMS; index+=blockDim.x * gridDim.x) {
        atomicAdd(&forceBuffer[index], static_cast<unsigned long long>((long long)(forces_in[index*3]*0x100000000)));                         // X component of the force
        atomicAdd(&forceBuffer[index+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long)(forces_in[index*3+1]*0x100000000)));      // Y component of the force
        atomicAdd(&forceBuffer[index+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long)(forces_in[index*3+2]*0x100000000)));    // Z component of the force
    }
    if (threadIndex == 0) {
        energyBuffer[0] += 0;
    }
}
