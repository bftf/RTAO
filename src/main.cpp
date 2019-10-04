#include "RayGenerator.h"
#include "BVHManager.h"
#include "RayTraceManager.h"

#include <iostream>

int main()
{
  RayGenerator rg = RayGenerator();
  rg.loadModelOBJ();
  rg.generateObjectRays();
  // rg.debugging();

  BVHManager bvh_manager = BVHManager();
  bvh_manager.buildBVH2(rg);
  bvh_manager.buildCWBVH(rg);

  RayTraceManager rt_manager = RayTraceManager(rg);
  rt_manager.traceOptiXPrime(rg);
  rt_manager.traceAila(rg); 
  rt_manager.traceCWBVH(rg);
  // rt_manager.evaluate(rg);
 
  return 0;
}
