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
  const float t_max = 100;
  const uint spp = 1;
  const uint samples_per_triangle = 4;

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

  friend class BVHManager;
};