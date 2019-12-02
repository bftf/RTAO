#include "helper_math.h"
#include <vector>
#include <string>

class BVHManager;

class RayGenerator
{
  public:
  
  RayGenerator(uint p_spp=1, uint p_spt=1, float p_t_min=0.1, float p_t_max=10);
  ~RayGenerator();

  // populates VertexBuffer and IndexBuffer
  int loadModelOBJ(const std::string& model_path);

  // generate ray_helper_vec and cudaRays - chose either one!
  void generateAORaysFromFile();
  void generatePointsAndNormals(int number_of_rays=-1);

  /* Generates points and normals */
  // Generate 1 point for each traingle on model, repeat spt times
  int generatePointsAndNormals_model(int number_of_rays=-1);
  // Generate spt points for each triangle
  int generatePointsAndNormals_spt(int number_of_rays=-1);

  void RandomizeAndDownsizeRays();
  void uploadRaysToGPU();
  void uploadSingleRayToGPU(int ray_id);
  
  int getRayCount() { return ray_helper_vec.size(); }

  int generate_entangled_spp(int entangle_size);
  int generate_detangled_spp();

  // ray file IO - this is used to make rays deterministic
  void saveRaysToFile(const std::string& file_path, const std::string& model_name);
  void readRaysFromFile(const std::string& file_path, const uint number_of_rays);

  enum ray_sorting { no_sort, random_shuffle, direction, origin, origin_chunk};
  void raySorting(const ray_sorting sorting_strategy = m_ray_sorting_strategy);

  void clear_rays() { ray_helper_vec.clear(); }

  // debugging
  void fillWithListOfRays();

  // Hack! Invert the normal for the Tepot model because of the OBJ file
  // const int invertNormal = -1; // teapot
  const int invertNormal = 1; // sponza + dragon + san-mig

  private:

  float t_min;
  float t_max;
  uint spp;
  uint samples_per_triangle;

  std::vector<Ray> ray_helper_vec;
  Ray* cudaRays = nullptr;

  std::vector<float3> VertexBuffer;
  std::vector<float3> NormalBuffer;
  std::vector<int3> IndexBuffer;

  std::vector<float3> points_vec;
  std::vector<float3> normal_vec;

  // give rays the proper ray direction
  float3 CosineSampleHemisphere();
  void GetTangentBasis(const float3 &n, float3 &b1, float3& b2);
  void generateSPPRaysFromWorldPosAndDir(float4 world_pos, float4 world_normal, std::vector<Ray>& out_rays, int my_spp=0);
  void generatePointInsideTriangle(
    const float3 a, const float3 b, const float3 c,
    const float3 n_a, const float3 n_b, const float3 n_c, 
    float3& out_point, float3& out_normal);

  static const ray_sorting m_ray_sorting_strategy = direction;
  
  friend class BVHManager;
  friend class RayTraceManager;
};