#include <vector>
#include <memory>

#include <cuda_runtime.h>
#include <cuda_profiler_api.h>


#include "CUDAAssert.h"



#include <optix_prime/optix_primepp.h>


#include "RayGenerator.h"
#include "BVHManager.h"

int main()
{
 
  // init the Ray Generator
  RayGenerator rg;
  rg.loadModelOBJ();
  rg.generateObjectRays();
  rg.debugging();

  BVHManager bvh_manager = BVHManager();


  // generate cudaRays
  {
    // numRays = ray_helper_vec.size();
    // cudaRays = GenerateRaysFromFile(ray_helper_vec, numRays);
    
  }

  {
    // BVH2  
  }

  {
    // CWBVH
  }

  /* TODO! REFACTOR THIS



  Hit* cudaHits;
  cudaCheck(cudaMalloc(&cudaHits, sizeof(Hit) *  numRays));

  printf("Launching kernels\n");

  rtTraceBVH2(cudaRays, cudaHits, numRays);
  rtTraceCWBVH(cudaRays, cudaHits, numRays);

#if 1 // Optional OptiX Prime comparison
  {
    optix::prime::Context OptiXContext = optix::prime::Context::create(RTP_CONTEXT_TYPE_CUDA);
    optix::prime::Model SceneModel = OptiXContext->createModel();
    SceneModel->setTriangles(IndexBuffer.size(), RTP_BUFFER_TYPE_HOST, IndexBuffer.data(), VertexBuffer.size(), RTP_BUFFER_TYPE_HOST, VertexBuffer.data());
    SceneModel->update(RTP_MODEL_HINT_NONE);
    SceneModel->finish();

    optix::prime::Query query = SceneModel->createQuery(RTP_QUERY_TYPE_CLOSEST);
    query->setRays(numRays, RTP_BUFFER_FORMAT_RAY_ORIGIN_TMIN_DIRECTION_TMAX, RTP_BUFFER_TYPE_CUDA_LINEAR, cudaRays);
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
#endif

  std::vector<Hit> hostHits(numRays);
  cudaCheck(cudaMemcpy(hostHits.data(), cudaHits, sizeof(Hit) * numRays, cudaMemcpyDeviceToHost));

  // fdemoullin
  for (int i = 0; i < numRays/spp; i++)
  {
    float sum = 0.f;
    for (int j = 0; j < spp; j++) 
    {
      if (hostHits[i*spp+j].t_triId_u_v.x > 0)
      {
        assert( hostHits[i * spp + j].t_triId_u_v.x >= t_min);
        assert(hostHits[i * spp + j].t_triId_u_v.x <= t_max);

        // sum += hostHits[i*spp+j].t_triId_u_v.x;
        sum++;
      }
    }
    if (sum > 0) 
    {
      std::cout << i << " " << reference_vec[i].x << " " << reference_vec[i].y << " " << (int)(sum) << std::endl;
    }
    
  }

  // Print out the first 10 results to validate by eye
  // for (int rayIndex = 0; rayIndex < 10; rayIndex++)
    // printf("%.2f %d\t", hostHits[rayIndex].t_triId_u_v.x, *(int*)&hostHits[rayIndex].t_triId_u_v.y);
  
  */
  

  return 0;
}
