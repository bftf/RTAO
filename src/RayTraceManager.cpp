
#include "RayTraceManager.h"
#include "RayGenerator.h"
#include "CUDAAssert.h"
#include "TraversalKernelBVH2.h"
#include "TraversalKernelCWBVH.h"

#include <cuda_runtime.h>
#include <cuda_profiler_api.h>
#include <optix_prime/optix_primepp.h>
#include <optix.h>

RayTraceManager::RayTraceManager(RayGenerator& rg)
{
  numRays = rg.ray_helper_vec.size();
  cudaCheck(cudaMalloc(&cudaHits, sizeof(Hit) * numRays));
};

void RayTraceManager::traceOptiXPrime(RayGenerator& rg)
{
  optix::prime::Context OptiXContext = optix::prime::Context::create(RTP_CONTEXT_TYPE_CUDA);
  optix::prime::Model SceneModel = OptiXContext->createModel();
  SceneModel->setTriangles(rg.IndexBuffer.size(), RTP_BUFFER_TYPE_HOST, rg.IndexBuffer.data(), rg.VertexBuffer.size(), RTP_BUFFER_TYPE_HOST, rg.VertexBuffer.data());
  SceneModel->update(RTP_MODEL_HINT_NONE);
  SceneModel->finish();

  optix::prime::Query query = SceneModel->createQuery(RTP_QUERY_TYPE_CLOSEST);
  query->setRays(numRays, RTP_BUFFER_FORMAT_RAY_ORIGIN_TMIN_DIRECTION_TMAX, RTP_BUFFER_TYPE_CUDA_LINEAR, rg.cudaRays);
  query->setHits(numRays, RTP_BUFFER_FORMAT_HIT_T_TRIID_U_V, RTP_BUFFER_TYPE_CUDA_LINEAR, cudaHits);

  cudaProfilerStart();
  {
    float elapsedTime;
    cudaEvent_t startEvent, stopEvent;
    cudaCheck(cudaEventCreate(&startEvent));
    cudaCheck(cudaEventCreate(&stopEvent));
    cudaCheck(cudaEventRecord(startEvent, 0));
    query->execute(0);
    cudaCheck(cudaEventRecord(stopEvent, 0));
    cudaCheck(cudaEventSynchronize(stopEvent));
    cudaCheck(cudaEventElapsedTime(&elapsedTime, startEvent, stopEvent));

    Log("%.3fMS, %.2fMRays/s (OptiXPrime)", elapsedTime, (float)numRays / 1000000.0f / (elapsedTime / 1000.0f));
  }
  cudaProfilerStop();
}

void RayTraceManager::traceAila(RayGenerator& rg)
{
  rtTraceBVH2(rg.cudaRays, cudaHits, numRays);
}

void RayTraceManager::traceCWBVH(RayGenerator& rg)
{
  rtTraceCWBVH(rg.cudaRays, cudaHits, numRays); 
}

void RayTraceManager::traceOptiX(RayGenerator& rg)
{
  RTcontext context = 0;
  // RT_CHECK_ERROR(rtContextCreate(&context)); 
}

void RayTraceManager::evaluate(RayGenerator& rg)
{
  std::vector<Hit> hostHits(numRays);
  cudaCheck(cudaMemcpy(hostHits.data(), cudaHits, sizeof(Hit) * numRays, cudaMemcpyDeviceToHost));
  assert(numRays % rg.spp == 0);

  for (int i = 0; i < numRays/rg.spp; i++)
  {
    float sum = 0.f;
    for (int j = 0; j < rg.spp; j++) 
    { 
      float cur_t = hostHits[i * rg.spp + j].t_triId_u_v.x;
      if (cur_t > rg.t_min && cur_t < rg.t_max)
      {
        assert(cur_t >= rg.t_min);
        assert(cur_t <= rg.t_max);

        sum += 1;  
        
      }
    }

    if (sum > 0) 
    {
      uint color = (uint)((255.f / (float)rg.spp) * sum);
      printf("%f %f %f %u %u %u\n", 
        rg.ray_helper_vec[i * rg.spp].origin_tmin.x, rg.ray_helper_vec[i * rg.spp].origin_tmin.y, rg.ray_helper_vec[i * rg.spp].origin_tmin.z, 
        color, color, color);
    }
    
  }
  printf("done!");

  // Print out the first 10 results to validate by eye
  // for (int rayIndex = 0; rayIndex < 10; rayIndex++)
  // printf("%.2f %d\t", hostHits[rayIndex].t_triId_u_v.x, *(int*)&hostHits[rayIndex].t_triId_u_v.y);
}

void RayTraceManager::debugging(RayGenerator& rg)
{
  std::vector<Hit> hostHits(numRays);
  cudaCheck(cudaMemcpy(hostHits.data(), cudaHits, sizeof(Hit) * numRays, cudaMemcpyDeviceToHost));
  assert(numRays % rg.spp == 0);

  printf("numrays traced: %i \n", (int)numRays); 

  for( int i = 0; i < numRays; i++)
  {
    printf("%.2f %d\t", hostHits[i].t_triId_u_v.x, *(int*)&hostHits[i].t_triId_u_v.y);

    printf("%f %f %f %f %f %f \n", 
    rg.ray_helper_vec[i * rg.spp].origin_tmin.x, rg.ray_helper_vec[i * rg.spp].origin_tmin.y, rg.ray_helper_vec[i * rg.spp].origin_tmin.z,
    rg.ray_helper_vec[i * rg.spp].dir_tmax.x, rg.ray_helper_vec[i * rg.spp].dir_tmax.y, rg.ray_helper_vec[i * rg.spp].dir_tmax.z );
  }
  
}
