#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <math.h>
#include <random>
#include <ctime>
#include <assert.h>

#include "RayGenerator.h"
#include "OBJ_Loader.h"
#include "ValidationKernels.h"


RayGenerator::RayGenerator()
{
  srand(time(0));
  return;
}

RayGenerator::~RayGenerator()
{
  return;
}

void RayGenerator::loadModelOBJ()
{

  objl::Loader loader;
  assert(loader.LoadFile("/home/francois/Documents/RayTracing/models/teapot/teapot.obj"));

  // loader.LoadFile("C:/Users/f.demoullin/Documents/graphics_assets/teapot.obj");
  // loader.LoadFile("C:/Users/f.demoullin/Documents/graphics_assets/conferenceRoomModel/conference.obj");
  // loader.LoadFile("C:/Users/f.demoullin/Documents/graphics_assets/AmazonLumberjard/exterior.obj");

  // generate Vertex and Index Buffers from obj loader:
  for (int i = 0; i < loader.LoadedIndices.size(); i += 3) 
  {
    int3 cur_index;
    cur_index.x = loader.LoadedIndices[i];
    cur_index.y = loader.LoadedIndices[i+1];
    cur_index.z = loader.LoadedIndices[i+2];
    IndexBuffer.push_back(cur_index);
  }
  for (int i = 0; i < loader.LoadedVertices.size(); i++) 
  {
    float3 cur_vertex;
    cur_vertex.x = loader.LoadedVertices[i].Position.X;
    cur_vertex.y = loader.LoadedVertices[i].Position.Y;
    cur_vertex.z = loader.LoadedVertices[i].Position.Z;
    VertexBuffer.push_back(cur_vertex);

    float3 cur_normal;
    cur_normal.x = loader.LoadedVertices[i].Normal.X;
    cur_normal.y = loader.LoadedVertices[i].Normal.Y;
    cur_normal.z = loader.LoadedVertices[i].Normal.Z;
    float length_cur_normal = sqrt(cur_normal.x*cur_normal.x + cur_normal.y*cur_normal.y + cur_normal.z*cur_normal.z);
    NormalBuffer.push_back(-cur_normal/length_cur_normal); // HACK HACK HACK
  }
}

// see pixar paper!
// Building an orthonormal basis revisited
void RayGenerator::GetTangentBasis(const float3 &n, float3 &b1, float3& b2)
{
  float sign = copysignf(1.0f, n.z);
  const float a = -1.0f / (sign + n.z);
  const float b = n.x * n.y * a;

  b1.x = 1.f + sign * n.x * n.x * a;
  b1.y = sign * b;
  b1.z = -sign * n.x;
  
  b2.x = b;
  b2.y = sign + n.y * n.y * a;
  b2.z = -n.y;
}

float3 RayGenerator::CosineSampleHemisphere()
{
  // generate random sample
  float2 E;
  
  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<float> dist(0.0, 1.0);

  E.x = dist(mt);
  E.y = dist(mt);
  
  // do the cosing hemisphere sampling!
  float Phi = 2 * M_PI * E.x;
  float CosTheta = sqrt(E.y);
  float SinTheta = sqrt(1 - CosTheta * CosTheta);

  float3 H = make_float3(SinTheta * cos(Phi), SinTheta * sin(Phi), CosTheta);
  return H;
}

void RayGenerator::generateSPPRaysFromWorldPosAndDir(float4 world_pos, float4 world_normal, std::vector<Ray>& out_rays)
{
  float3 n, b1, b2;

  const float world_normal_magnitude = sqrt(world_normal.x * world_normal.x + world_normal.y * world_normal.y + world_normal.z * world_normal.z);
  n.x = world_normal.x / world_normal_magnitude;
  n.y = world_normal.y / world_normal_magnitude;
  n.z = world_normal.z / world_normal_magnitude;

  // get orthonormal basis of the normal
  GetTangentBasis(n, b1, b2);
  
  for(int i=0; i < spp; i++)
  {
    // sample the hemisphere, get a tangent
    float3 tangent = CosineSampleHemisphere();

    // translate tangent into world_normal space
    // tangent * float3x3(row0 = b1, row1 = b2, row2 = n)
    float4 ray_direction;
    ray_direction.x = tangent.x * b1.x + tangent.y * b2.x + tangent.z * n.x;
    ray_direction.y = tangent.x * b1.y + tangent.y * b2.y + tangent.z * n.y;
    ray_direction.z = tangent.x * b1.z + tangent.y * b2.z + tangent.z * n.z;

    // normalize ray_direction - should not be required?
    ray_direction = ray_direction / sqrt(ray_direction.x * ray_direction.x + ray_direction.y * ray_direction.y + ray_direction.z * ray_direction.z);

    // encode for Ray struct
    ray_direction.w = t_max;
    world_pos.w = t_min;

    Ray cur_ray;
    cur_ray.origin_tmin = world_pos;
    cur_ray.dir_tmax = ray_direction;

    out_rays.push_back(cur_ray); 
  }
}

void RayGenerator::generateAORaysFromFile()
{
  std::string line;
  std::ifstream infile("C:/Users/f.demoullin/Documents/graphics_assets/ray_file_teapot.txt");
  // std::ifstream infile("C:/Users/f.demoullin/Documents/graphics_assets/ray_file_conference.txt");

  while (std::getline(infile, line))
  {
    std::stringstream iss(line);
    uint2 cur_coords;
    float4 world_pos;
    float4 world_normal;

    iss >> cur_coords.x;
    iss >> cur_coords.y;

    iss >> world_pos.x;
    iss >> world_pos.y;
    iss >> world_pos.z;

    iss >> world_normal.x;
    iss >> world_normal.y;
    iss >> world_normal.z;

    std::vector<Ray> temp;
    generateSPPRaysFromWorldPosAndDir(world_pos, world_normal, temp);
    ray_helper_vec.insert(std::end(ray_helper_vec), std::begin(temp), std::end(temp));
  }

  cudaRays = GenerateRaysFromFile(ray_helper_vec, ray_helper_vec.size());
}

/*
* See: https://www.cs.princeton.edu/~funk/tog02.pdf, Chapter 4.2
* For the proof, start here: https://math.stackexchange.com/questions/18686/uniform-random-point-in-triangle
*/
void RayGenerator::generatePointInsideTriangle(
  const float3 a, const float3 b, const float3 c,
  const float3 n_a, const float3 n_b, const float3 n_c, 
  float3& out_point, float3& out_normal)
  {     
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<float> dist(0.0, 1.0);

    float r1 = dist(mt);
    float r2 = dist(mt);

    out_point = (1 - sqrt(r1)) * a + (sqrt(r1) * (1 - r2)) * b + (sqrt(r1) * r2) * c;
    // assumption: the same linear combination can be used for the normal at out_point
    out_normal = (1 - sqrt(r1)) * n_a + (sqrt(r1) * (1 - r2)) * n_b + (sqrt(r1) * r2) * n_c;
  }

void RayGenerator::generateObjectRays()
{
  for (auto cur_index : IndexBuffer)
  {
    for (int cur_samples_per_triangle = 0; cur_samples_per_triangle < samples_per_triangle; cur_samples_per_triangle++)
    {
      float3 cur_point, cur_normal;
      generatePointInsideTriangle(
        VertexBuffer[cur_index.x], VertexBuffer[cur_index.y], VertexBuffer[cur_index.z],
        NormalBuffer[cur_index.x], NormalBuffer[cur_index.y], NormalBuffer[cur_index.z],
        cur_point, cur_normal);

      std::vector<Ray> temp;
      generateSPPRaysFromWorldPosAndDir(make_float4(cur_point), make_float4(cur_normal), temp);
      ray_helper_vec.insert(std::end(ray_helper_vec), std::begin(temp), std::end(temp));
    }
  }

  cudaRays = GenerateRaysFromFile(ray_helper_vec, ray_helper_vec.size());
}

void RayGenerator::debugging()
{
  std::cout << ray_helper_vec.size() << std::endl;

  /*
  for (int i = 0; i < ray_helper_vec.size(); i += spp)
  {
    float3 direction_point = VertexBuffer[i] + 1 * NormalBuffer[i]; 
    printf("%f %f %f %u %u %u\n", direction_point.x, direction_point.y, direction_point.z, 0, 255, 0);
  }
  */
  
  for (int i = 0; i < ray_helper_vec.size(); i += spp)
  {
    Ray cur_ray = ray_helper_vec[i];
    // printf("%f %f %f %u %u %u\n", cur_ray.origin_tmin.x, cur_ray.origin_tmin.y, cur_ray.origin_tmin.z, 255, 0, 0);
    
    for (int j = 0; j < spp; j++)
    {
      Ray cur_ray = ray_helper_vec[i + j];
      float4 direction_point = cur_ray.origin_tmin + 1 * cur_ray.dir_tmax; 
      printf("%f %f %f %u %u %u\n", direction_point.x, direction_point.y, direction_point.z, 0, 255, 0);
    }  
  }

}