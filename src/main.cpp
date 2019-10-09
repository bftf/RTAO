#include "RayGenerator.h"
#include "BVHManager.h"
#include "RayTraceManager.h"

#include <iostream>

int main()
{
  RayGenerator rg = RayGenerator();
  rg.loadModelOBJ();
  rg.generateObjectRays();
  rg.uploadRaysToGPU();
  rg.fillWithListOfRays(); // fill custom list of rays for debugging
  // rg.printRaysForVisualization(); // prints list of rays for visualization

  BVHManager bvh_manager = BVHManager();
  // bvh_manager.buildBVH2(rg);
  bvh_manager.buildCWBVH(rg);

  RayTraceManager rt_manager = RayTraceManager(rg);
  // rt_manager.traceOptiXPrime(rg);
  // rt_manager.traceAila(rg); 
  rt_manager.traceCWBVH(rg);

  // rt_manager.evaluateAndPrintForPLYVisualization(rg); // prints outcome from kernels for visualization
  // rt_manager.debugging(rg); // for debugging - prints output of kernels into command line
  
  printf("Done, traced %i rays \n", rg.getRayCount());
  return 0;
}

// /home/francois/cuda_10_0/lib64:/home/francois/Documents/RayTracing/NVIDIA-OptiX-SDK-5.1.0-linux64/lib64/
// /home/francois/Documents/gpgpu-sim_distribution/lib/gcc-5.4.0/cuda-10000/debug:/home/francois/Documents/RayTracing/NVIDIA-OptiX-SDK-5.1.0-linux64/lib64/