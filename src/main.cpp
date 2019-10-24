#include "RayGenerator.h"
#include "BVHManager.h"
#include "RayTraceManager.h"

#include <iostream>
#include <string>

int main()
{
  // const std::string model_path = "/home/francois/Documents/RayTracing/models/teapot/teapot.obj";
  // const std::string out_ply_path = "/home/francois/Documents/RayTracing/CWBVH/build_make/ply_files/teapot.ply";

  // const std::string model_path = "/home/francois/Documents/RayTracing/models/Sponza/models/sponza.obj";
  // const std::string out_ply_path = "/home/francois/Documents/RayTracing/CWBVH/build_make/ply_files/sponza.ply";

  const std::string model_path = "/home/francois/Documents/RayTracing/models/Dragon/dragon.obj";
  const std::string out_ply_path = "/home/francois/Documents/RayTracing/CWBVH/build_make/ply_files/dragon.ply";

  RayGenerator rg = RayGenerator();
  rg.loadModelOBJ(model_path);
  printf("Done loading obj\n");
  
  rg.generateObjectRays();
  rg.downsizeRayVector(5000);
  printf("Done generating rays\n");
  
  // rg.fillWithListOfRays(); // fill custom list of rays for debugging
  // rg.printRaysForVisualization(); // prints list of rays for visualization
  rg.uploadRaysToGPU();

  BVHManager bvh_manager = BVHManager();
  RayTraceManager rt_manager = RayTraceManager(rg);

  /* OptiX */
  rt_manager.traceOptiXPrime(rg);

  /* Aila */
  // bvh_manager.buildBVH2(rg);
  // rt_manager.traceAila(rg); 
  
  /* CWBVH */
  bvh_manager.buildCWBVH(rg);
  printf("Start tracing %i rays \n", rg.getRayCount());
  rt_manager.traceCWBVH(rg);


  // rt_manager.evaluateAndPrintForPLYVisualization(rg, out_ply_path); // prints outcome from kernels for visualization
  // rt_manager.debugging(rg); // for debugging - prints output of kernels into command line
  
  printf("Done, traced %i rays \n", rg.getRayCount());
  return 0;
}

// /home/francois/cuda_10_0/lib64:/home/francois/Documents/RayTracing/NVIDIA-OptiX-SDK-5.1.0-linux64/lib64/
// /home/francois/Documents/gpgpu-sim_distribution/lib/gcc-5.4.0/cuda-10000/debug:/home/francois/Documents/RayTracing/NVIDIA-OptiX-SDK-5.1.0-linux64/lib64/
// /home/francois/Documents/gpgpu-sim_distribution/lib/gcc-5.4.0/cuda-10000/release:/home/francois/Documents/RayTracing/NVIDIA-OptiX-SDK-5.1.0-linux64/lib64/