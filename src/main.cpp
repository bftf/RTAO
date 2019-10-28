#include "RayGenerator.h"
#include "BVHManager.h"
#include "RayTraceManager.h"

#include <iostream>
#include <string>

int main()
{

  const std::string ply_base_path = "/home/francois/Documents/RayTracing/CWBVH/build_make/ply_files/";
  const std::string model_base_path = "/home/francois/Documents/RayTracing/models/";
  const std::string ray_base_path = "/home/francois/Documents/RayTracing/CWBVH/build_make/ray_files/"; 

  // const std::string model_name = "teapot";
  // const std::string obj_path_addition = "/teapot/teapot.obj";

  // const std::string model_name = "sponza";
  // const std::string obj_path_addition = "sponza.obj";

  // const std::string model_name = "dragon";
  // const std::string obj_path_addition = "Dragon/dragon.obj";

  const std::string model_name = "san-miguel";
  const std::string obj_path_addition = "San_Miguel/san-miguel.obj";

  const std::string out_ply_path = ply_base_path + model_name + ".ply";
  const std::string model_path = model_base_path + obj_path_addition;

  RayGenerator rg = RayGenerator(4, 4, 0.1, 10); /*spp, spt, t_min, t_max*/
  rg.loadModelOBJ(model_path);
  printf("Done loading obj\n");
  
  rg.generateObjectRays(500000);
  printf("Done generating rays\n");

  rg.saveRaysToFile(ray_base_path, model_name);
  //rg.readRaysFromFile(ray_base_path + "teapot_100_28_9.ray_file", 100);
  
  // rg.fillWithListOfRays(); // fill custom list of rays for debugging
  
  rg.uploadRaysToGPU();

  BVHManager bvh_manager = BVHManager();
  RayTraceManager rt_manager = RayTraceManager(rg);

  /* OptiX */
  // rt_manager.traceOptiXPrime(rg);

  /* Aila */
  // bvh_manager.buildBVH2(rg);
  // rt_manager.traceAila(rg); 
  
  /* CWBVH */
  bvh_manager.buildCWBVH(rg);
  printf("Start tracing %i rays \n", rg.getRayCount());
  rt_manager.traceCWBVH(rg);


  rt_manager.evaluateAndPrintForPLYVisualization(rg, out_ply_path); // prints outcome from kernels for visualization
  // rt_manager.debugging(rg); // for debugging - prints output of kernels into command line
  
  printf("Done, traced %i rays \n", rg.getRayCount());
  return 0;
}

// aamodt-pc12
// /home/francois/cuda_10_0/lib64:/home/francois/Documents/RayTracing/NVIDIA-OptiX-SDK-5.1.0-linux64/lib64/
// /home/francois/Documents/gpgpu-sim_distribution/lib/gcc-5.4.0/cuda-10000/debug:/home/francois/Documents/RayTracing/NVIDIA-OptiX-SDK-5.1.0-linux64/lib64/
// /home/francois/Documents/gpgpu-sim_distribution/lib/gcc-5.4.0/cuda-10000/release:/home/francois/Documents/RayTracing/NVIDIA-OptiX-SDK-5.1.0-linux64/lib64/

// aamodt-pc16
// /home/francois/Documents/gpgpu-sim_cuda10/lib/gcc-7.4.0/cuda-10000/debug:/home/francois/Documents/embree/build/
// /home/common/modules/cuda/10.0.130/lib64:/home/francois/Documents/embree/build/