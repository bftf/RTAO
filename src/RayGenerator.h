#include "helper_math.h"
#include <vector>

class BVHManager;

class RayGenerator
{
  public:
  
  RayGenerator();
  ~RayGenerator();

  // populates VertexBuffer and IndexBuffer
  void loadModelOBJ();

  // generate ray_helper_vec and cudaRays - chose either one!
  void generateAORaysFromFile();
  void generateObjectRays();

  void debugging();

  private:

  const float t_min = 0.1;
  const float t_max = 10;
  const uint spp = 16;
  const uint samples_per_triangle = 64;

  std::vector<Ray> ray_helper_vec;
  Ray* cudaRays = nullptr;

  std::vector<float3> VertexBuffer;
  std::vector<float3> NormalBuffer;
  std::vector<int3> IndexBuffer;

  // give rays the proper ray direction
  float3 CosineSampleHemisphere();
  void GetTangentBasis(const float3 &n, float3 &b1, float3& b2);
  void generateSPPRaysFromWorldPosAndDir(float4 world_pos, float4 world_normal, std::vector<Ray>& out_rays);
  void generatePointInsideTriangle(
    const float3 a, const float3 b, const float3 c,
    const float3 n_a, const float3 n_b, const float3 n_c, 
    float3& out_point, float3& out_normal);

  enum ray_sorting { no_sort, random_shuffle, direction, origin };
  ray_sorting m_ray_sorting_strategy = direction;

  void raySorting(std::vector<Ray>& v);
  
  friend class BVHManager;
  friend class RayTraceManager;
};