#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <math.h>
#include <random>
#include <ctime>
#include <assert.h>
#include <algorithm>

#include "RayGenerator.h"
#include "OBJ_Loader.h"
#include "ValidationKernels.h"

static bool compareByDirection(const Ray &a, const Ray &b)
{
  float4 a_d = a.dir_tmax;
  float4 b_d = b.dir_tmax;

  double a_d_d = a_d.x * 1/sqrt(3) + a_d.y * 1/sqrt(3) + a_d.z * 1/sqrt(3);
  a_d_d = (a_d_d < -1.0 ? -1.0 : (a_d_d > 1.0 ? 1.0 : a_d_d));
  
  double b_d_d = b_d.x * 1/sqrt(3) + b_d.y * 1/sqrt(3) + b_d.z * 1/sqrt(3);
  b_d_d = (b_d_d < -1.0 ? -1.0 : (b_d_d > 1.0 ? 1.0 : b_d_d));

  double a_angle = acos(a_d_d);
  double b_angle = acos(b_d_d);

  return a_angle < b_angle;
}

static bool compareByOrigin(const Ray &a, const Ray &b)
{
  float4 a_o = a.origin_tmin;
  float4 b_o = b.origin_tmin;

  float a_o_d = sqrt(a_o.x * a_o.x + a_o.y * a_o.y + a_o.z * a_o.z);
  float b_o_d = sqrt(b_o.x * b_o.x + b_o.y * b_o.y + b_o.z * b_o.z);

  return a_o_d < b_o_d;
}

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

void RayGenerator::raySorting(std::vector<Ray>& v)
{
  switch(m_ray_sorting_strategy)
  {
    case no_sort: { break; }
    case random_shuffle: 
    {
      // sorting goes here!
      std::random_device rd;
      std::mt19937 g(rd());
      std::shuffle(v.begin(), v.end(), g);
      break;
    }
    case direction:
    {
      sort(v.begin(), v.end(), compareByDirection);
      break;
    }
    case origin:
    {
      sort(v.begin(), v.end(), compareByOrigin);
      break;
    }
  }
}

void RayGenerator::generateObjectRays()
{
  // /* Commenting out for testbench
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

  raySorting(ray_helper_vec);

  // */

  /*
  Ray debug_ray; 
  // this ray misses!
  debug_ray.origin_tmin = make_float4(1.382953, 2.411897, -0.156976, t_min);
  debug_ray.dir_tmax = make_float4(-0.417863, 0.406558, 0.812466, t_max);

  // this ray hits
  debug_ray.origin_tmin = make_float4(1.379896, 2.418430, -0.162858, t_min);
  debug_ray.dir_tmax = make_float4(-0.743145, -0.584695, -0.325373, t_max);
  
  ray_helper_vec.push_back(debug_ray); 
  */
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

void RayGenerator::uploadRays()
{
  cudaRays = GenerateRaysFromFile(ray_helper_vec, ray_helper_vec.size());
}

void RayGenerator::fillWithListOfRays()
{
  ray_helper_vec.resize(1);

  int i = 0; 
  // doesn't finish!! 120 rays
  ray_helper_vec[i].make_ray(1.382365, 2.400508, -0.222089, -0.688209, 0.090686, 0.719823, 0.100000, 10.000000); i++;

/*
ray_helper_vec[i].make_ray(1.371162, 2.435336, -0.188570, -0.796378, 0.598703, 0.085652, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(1.339913, 2.402648, -0.399668, -0.815969, -0.577133, 0.033343, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(1.376911, 2.410918, -0.232221, 0.238469, -0.957137, 0.164379, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(1.253219, 2.421222, -0.608350, 0.057156, -0.704196, 0.707702, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(1.283489, 2.422934, -0.533878, -0.959833, 0.254433, 0.118254, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(1.174238, 2.413205, -0.754423, -0.676332, 0.717534, 0.166492, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(1.138589, 2.430622, -0.799165, -0.034084, 0.313031, 0.949131, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(1.078339, 2.405096, -0.892272, -0.445994, -0.470326, 0.761500, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(1.039050, 2.425531, -0.925709, -0.699509, -0.524711, -0.485146, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(0.931053, 2.404622, -1.045419, -0.127005, -0.545635, 0.828343, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(0.920085, 2.426615, -1.043293, 0.066757, 0.319319, 0.945293, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(0.695562, 2.416121, -1.209086, 0.109543, -0.977200, -0.181883, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(0.691388, 2.427629, -1.206285, -0.419841, -0.895149, 0.149808, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(0.499359, 2.403853, -1.306064, -0.081860, -0.189182, 0.978524, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(0.580183, 2.419757, -1.265576, -0.179151, 0.694778, 0.696555, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(0.320691, 2.416475, -1.353315, -0.518742, -0.488586, 0.701563, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(0.423454, 2.416458, -1.328547, 0.146167, -0.073696, 0.986511, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(0.086985, 2.411746, -1.388507, 0.677337, -0.583175, 0.448465, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(0.160039, 2.425794, -1.377192, 0.413330, -0.899393, -0.142307, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(1.364147, 2.456148, -0.224541, -0.830631, -0.493349, 0.258183, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(1.375882, 2.457595, -0.071359, -0.656876, -0.699146, 0.282326, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(1.335821, 2.435544, -0.361090, -0.975234, -0.194924, -0.104513, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(1.359038, 2.450637, -0.251908, -0.489339, -0.036897, 0.871313, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(1.249422, 2.444684, -0.597969, -0.489667, -0.044997, 0.870748, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(1.270999, 2.450718, -0.542998, -0.906015, 0.035522, 0.421753, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(1.122965, 2.454505, -0.813614, -0.059985, -0.631531, 0.773027, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(1.218555, 2.462325, -0.655226, -0.573796, -0.249331, 0.780123, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(1.049258, 2.444765, -0.904834, -0.642764, -0.612675, 0.459874, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(1.065757, 2.456927, -0.881788, 0.057473, -0.566358, 0.822153, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(0.864820, 2.446386, -1.083053, -0.568192, -0.166257, 0.805926, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(0.880485, 2.458230, -1.066528, -0.318442, 0.079607, 0.944594, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(0.722468, 2.436638, -1.183240, 0.039048, -0.675986, 0.735879, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(0.804025, 2.438656, -1.132595, -0.025968, 0.291573, 0.956196, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(0.538525, 2.439370, -1.275324, 0.811958, 0.228253, 0.537238, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(0.566208, 2.456754, -1.260027, -0.178294, -0.678749, 0.712399, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(0.373423, 2.438268, -1.332287, -0.556184, -0.660611, 0.504235, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(0.401095, 2.442338, -1.324778, 0.185710, -0.394163, 0.900082, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(0.008792, 2.446363, -1.383045, -0.162731, 0.585481, 0.794186, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(0.169131, 2.451705, -1.369387, -0.180841, 0.520501, 0.834491, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(1.365512, 2.465453, -0.192558, -0.688415, 0.075266, -0.721402, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(1.377684, 2.469195, -0.041914, -0.825403, -0.550056, -0.127075, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(1.336182, 2.472784, -0.340019, -0.419911, 0.557017, 0.716524, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(1.349698, 2.472884, -0.283996, -0.394168, -0.903887, -0.166194, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(1.237200, 2.463438, -0.617883, -0.944089, -0.325421, -0.052881, 0.100000, 10.000000); i++;

ray_helper_vec[i].make_ray(1.280202, 2.474146, -0.517180, -0.342843, -0.078731, 0.936088, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(1.151700, 2.470037, -0.764895, -0.982272, -0.005419, 0.187382, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(1.208346, 2.470816, -0.672948, -0.003655, 0.787580, 0.616202, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(1.002231, 2.479640, -0.956565, -0.907686, -0.409110, 0.093463, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(1.115499, 2.481941, -0.824459, 0.303308, 0.142749, 0.942140, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(0.864919, 2.470056, -1.079443, -0.975839, -0.162202, 0.146383, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(0.945800, 2.475055, -1.010878, -0.699795, 0.373258, 0.609070, 0.100000, 10.000000); i++;

ray_helper_vec[i].make_ray(0.722108, 2.465477, -1.177524, -0.893769, -0.331318, -0.302334, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(0.775222, 2.473571, -1.145732, -0.560550, -0.770175, 0.304326, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(0.601899, 2.463311, -1.243842, -0.594625, -0.715747, 0.366234, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(0.621572, 2.480246, -1.237351, -0.518732, -0.269353, 0.811397, 0.100000, 10.000000); i++;

// THIS IS THE FUCKER!! WHY DO YOU SUCK SO MUCH?
ray_helper_vec[i].make_ray(0.372081, 2.463659, -1.327584, 0.020548, 0.182130, 0.983060, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(0.396929, 2.476619, -1.322826, -0.306163, 0.199582, 0.930823, 0.100000, 10.000000); i++;

ray_helper_vec[i].make_ray(0.126553, 2.464505, -1.370605, -0.941568, -0.141493, 0.305661, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(0.208290, 2.470536, -1.364748, -0.904755, 0.298510, 0.303827, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(1.380467, 2.484768, -0.040303, -0.784370, -0.264529, 0.561060, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(1.383408, 2.490822, -0.054886, -0.899118, -0.166731, 0.404707, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(1.331774, 2.489940, -0.383060, -0.573806, 0.755571, -0.316004, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(1.359466, 2.490663, -0.270277, -0.213397, -0.009348, 0.976921, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(1.258764, 2.489914, -0.583417, -0.592507, 0.730646, 0.339253, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(1.311077, 2.485118, -0.449349, 0.190031, 0.832224, 0.520856, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(1.140074, 2.489575, -0.794913, -0.514756, 0.358657, -0.778712, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(1.210014, 2.487439, -0.678464, -0.024960, -0.281267, 0.959305, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(0.997742, 2.492829, -0.972794, -0.995247, -0.073582, 0.063782, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(1.094083, 2.488571, -0.855723, -0.998993, -0.043421, -0.011273, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(0.847196, 2.482900, -1.096333, -0.565081, 0.823242, 0.054365, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(0.874623, 2.493861, -1.082623, 0.263231, 0.908025, 0.325882, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(0.661722, 2.484344, -1.217859, -0.694097, -0.191643, 0.693903, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(0.714963, 2.493906, -1.192690, -0.440307, 0.771007, 0.460085, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(0.501336, 2.487133, -1.290904, -0.376027, -0.429361, 0.821129, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(0.591139, 2.489311, -1.255108, -0.943175, 0.143159, 0.299876, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(0.333236, 2.485739, -1.340875, 0.136001, 0.385712, 0.912541, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(0.384558, 2.489928, -1.331405, -0.344317, 0.874083, -0.342674, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(0.219015, 2.482983, -1.365233, 0.041543, 0.888577, 0.456842, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(0.104639, 2.494024, -1.381665, -0.931253, 0.329087, 0.156427, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(1.385211, 2.497266, -0.179202, -0.162394, 0.986723, -0.002326, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(1.391847, 2.497941, -0.122920, -0.479718, 0.713756, 0.510316, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(1.363333, 2.494731, -0.268500, -0.936592, 0.260244, 0.234666, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(1.371831, 2.496788, -0.262002, -0.548817, 0.817585, -0.174229, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(1.281657, 2.494660, -0.537847, -0.081307, 0.652274, -0.753610, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(1.320699, 2.496408, -0.458994, -0.214505, 0.950564, 0.224534, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(1.134184, 2.495100, -0.814595, -0.023993, 0.822633, -0.568067, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(1.160269, 2.497526, -0.787398, 0.421235, 0.876771, -0.232023, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(1.076956, 2.495727, -0.888115, 0.712873, 0.451198, 0.536873, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(1.038512, 2.497703, -0.943055, 0.153454, 0.708451, 0.688875, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(0.860647, 2.496916, -1.105588, -0.904409, 0.424549, 0.042457, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(0.986412, 2.496607, -0.996691, -0.805381, 0.589447, 0.062556, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(0.796394, 2.494724, -1.143931, -0.592169, 0.761600, 0.263253, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(0.769034, 2.496223, -1.166543, -0.129882, 0.847135, 0.515260, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(0.480209, 2.494820, -1.306221, -0.890209, 0.431911, 0.144849, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(0.599574, 2.497161, -1.264844, -0.153242, 0.527624, 0.835542, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(0.317985, 2.496695, -1.358021, -0.419866, 0.667025, 0.615459, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(0.394298, 2.496449, -1.338794, -0.132340, 0.937681, 0.321311, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(0.098900, 2.496718, -1.389716, 0.000690, 0.999331, 0.036556, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(0.179557, 2.497617, -1.386334, -0.263937, 0.761865, -0.591522, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(1.393938, 2.497982, -0.141648, 0.546686, 0.772309, 0.323532, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(1.406069, 2.494613, -0.168849, -0.066986, 0.898117, -0.434625, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(1.349841, 2.496479, -0.411606, -0.145883, 0.984283, -0.099522, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(1.378829, 2.495545, -0.308200, -0.640871, 0.764785, 0.066237, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(1.296466, 2.497827, -0.541246, -0.077202, 0.895130, -0.439071, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(1.312988, 2.495625, -0.525806, 0.447233, 0.526566, -0.722988, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(1.183028, 2.498130, -0.758609, 0.604682, 0.513864, 0.608526, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(1.248736, 2.495857, -0.670171, 0.745680, 0.594470, -0.300945, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(1.068156, 2.497882, -0.915730, -0.169043, 0.948512, 0.267862, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(1.097606, 2.494823, -0.901253, 0.834846, 0.512218, 0.201655, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(0.853036, 2.495761, -1.133613, -0.716585, 0.660479, -0.224217, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(0.982569, 2.496132, -1.020757, -0.134374, 0.596684, -0.791146, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(0.676570, 2.495806, -1.245052, -0.267629, 0.943132, -0.197173, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(0.813200, 2.495568, -1.162185, -0.231892, 0.270603, -0.934345, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(0.509116, 2.496655, -1.315217, 0.362312, 0.753689, -0.548345, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(0.521995, 2.494901, -1.317895, -0.175084, 0.975671, -0.131957, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(0.265862, 2.497072, -1.382402, -0.545838, 0.663295, -0.511957, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(0.282140, 2.494875, -1.388023, -0.248369, 0.618754, -0.745289, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(0.008503, 2.495338, -1.415581, -0.884264, 0.428765, -0.185034, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(0.193555, 2.495647, -1.399748, -0.673487, 0.669298, -0.313777, 0.100000, 10.000000); i++;
*/

  /*
  if (i != ray_helper_vec.size())
  {
    printf("%i, %i\n", i, (int)ray_helper_vec.size());
    assert(i == ray_helper_vec.size());  
  }

    for( int i = 0; i < ray_helper_vec.size(); i++)
    {
      printf("ray_helper_vec[i].make_ray(%f, %f, %f, %f, %f, %f, %f, %f); i++;\n", 
        ray_helper_vec[i].origin_tmin.x, ray_helper_vec[i].origin_tmin.y, ray_helper_vec[i].origin_tmin.z,
        ray_helper_vec[i].dir_tmax.x, ray_helper_vec[i].dir_tmax.y, ray_helper_vec[i].dir_tmax.z, 
        ray_helper_vec[i].origin_tmin.w, ray_helper_vec[i].dir_tmax.w
      );
    }
  }
  */
}

/*
// 240 - finishes 6 out of 6 (105 secs)
ray_helper_vec[0].make_ray(1.392292, 2.406880, -0.063378, -0.641424, -0.352417, 0.681453, 0.100000, 10.000000);
ray_helper_vec[1].make_ray(1.375892, 2.430883, -0.150853, -0.575911, -0.810474, -0.107049, 0.100000, 10.000000);
ray_helper_vec[2].make_ray(1.330851, 2.405409, -0.432583, -0.889826, -0.402806, 0.214376, 0.100000, 10.000000);
ray_helper_vec[3].make_ray(1.337846, 2.434403, -0.354529, -0.893245, -0.449197, 0.018283, 0.100000, 10.000000);
ray_helper_vec[4].make_ray(1.264574, 2.420126, -0.582208, -0.669584, -0.633028, 0.388501, 0.100000, 10.000000);
ray_helper_vec[5].make_ray(1.267083, 2.428584, -0.567446, -0.509429, -0.219582, 0.832025, 0.100000, 10.000000);
ray_helper_vec[6].make_ray(1.212835, 2.407542, -0.695974, -0.869675, -0.303751, -0.389104, 0.100000, 10.000000);
ray_helper_vec[7].make_ray(1.215293, 2.416025, -0.685545, -0.835568, 0.004158, 0.549371, 0.100000, 10.000000);
ray_helper_vec[8].make_ray(1.015043, 2.404240, -0.966802, -0.861588, -0.213038, 0.460739, 0.100000, 10.000000);
ray_helper_vec[9].make_ray(1.083520, 2.424719, -0.874210, -0.822016, -0.446118, -0.353932, 0.100000, 10.000000);
ray_helper_vec[10].make_ray(0.859410, 2.402907, -1.107590, -0.743499, -0.600204, -0.294896, 0.100000, 10.000000);
ray_helper_vec[11].make_ray(0.957435, 2.427222, -1.011032, -0.914424, -0.230013, 0.333051, 0.100000, 10.000000);
ray_helper_vec[12].make_ray(0.706357, 2.401536, -1.209250, -0.019657, -0.598681, 0.800746, 0.100000, 10.000000);
ray_helper_vec[13].make_ray(0.759488, 2.430387, -1.163107, 0.020499, 0.675263, 0.737292, 0.100000, 10.000000);
ray_helper_vec[14].make_ray(0.542886, 2.409120, -1.285678, 0.363726, 0.535899, 0.761915, 0.100000, 10.000000);
ray_helper_vec[15].make_ray(0.587496, 2.414554, -1.264768, 0.299796, -0.700031, 0.648135, 0.100000, 10.000000);
ray_helper_vec[16].make_ray(0.302609, 2.406819, -1.361612, -0.477909, -0.665262, 0.573611, 0.100000, 10.000000);
ray_helper_vec[17].make_ray(0.438831, 2.418477, -1.324016, -0.214190, 0.000917, 0.976792, 0.100000, 10.000000);
ray_helper_vec[18].make_ray(0.096856, 2.406073, -1.389985, 0.180378, -0.651526, 0.736871, 0.100000, 10.000000);
ray_helper_vec[19].make_ray(0.132040, 2.417585, -1.382651, 0.025912, -0.985231, 0.169256, 0.100000, 10.000000);
ray_helper_vec[20].make_ray(1.382188, 2.437813, -0.041538, -0.906852, 0.057575, -0.417498, 0.100000, 10.000000);
ray_helper_vec[21].make_ray(1.371029, 2.459855, -0.127406, -0.446238, -0.810724, -0.378943, 0.100000, 10.000000);
ray_helper_vec[22].make_ray(1.352844, 2.441952, -0.285012, -0.514179, -0.248048, -0.821031, 0.100000, 10.000000);
ray_helper_vec[23].make_ray(1.328143, 2.459620, -0.372385, -0.365427, -0.879244, 0.305604, 0.100000, 10.000000);
ray_helper_vec[24].make_ray(1.296394, 2.441381, -0.486871, -0.852907, -0.175259, 0.491766, 0.100000, 10.000000);
ray_helper_vec[25].make_ray(1.257410, 2.461162, -0.570195, -0.859777, 0.474335, -0.189183, 0.100000, 10.000000);
ray_helper_vec[26].make_ray(1.208971, 2.437404, -0.680346, -0.609118, -0.784954, -0.113239, 0.100000, 10.000000);
ray_helper_vec[27].make_ray(1.178054, 2.460121, -0.721909, -0.844891, -0.527867, 0.086696, 0.100000, 10.000000);
ray_helper_vec[28].make_ray(1.028068, 2.437071, -0.931986, -0.961760, -0.178514, -0.207727, 0.100000, 10.000000);
ray_helper_vec[29].make_ray(1.063950, 2.460931, -0.882665, -0.454549, -0.241321, 0.857409, 0.100000, 10.000000);
ray_helper_vec[30].make_ray(0.929062, 2.441089, -1.029507, -0.884585, 0.389048, -0.257199, 0.100000, 10.000000);
ray_helper_vec[31].make_ray(0.970652, 2.450878, -0.991352, 0.260295, -0.438691, 0.860114, 0.100000, 10.000000);
ray_helper_vec[32].make_ray(0.669968, 2.443246, -1.213980, -0.417654, 0.130088, 0.899245, 0.100000, 10.000000);
ray_helper_vec[33].make_ray(0.756334, 2.458468, -1.157267, -0.552767, -0.696863, 0.456980, 0.100000, 10.000000);
ray_helper_vec[34].make_ray(0.534371, 2.447715, -1.275243, -0.045690, -0.649792, 0.758738, 0.100000, 10.000000);
ray_helper_vec[35].make_ray(0.625972, 2.456219, -1.235260, 0.197230, -0.979955, 0.028088, 0.100000, 10.000000);
ray_helper_vec[36].make_ray(0.343260, 2.441013, -1.338994, 0.378462, -0.735719, 0.561679, 0.100000, 10.000000);
ray_helper_vec[37].make_ray(0.346265, 2.455473, -1.335294, -0.025096, -0.586984, 0.809209, 0.100000, 10.000000);
ray_helper_vec[38].make_ray(0.145044, 2.443914, -1.372840, 0.472540, -0.706862, 0.526357, 0.100000, 10.000000);
ray_helper_vec[39].make_ray(0.205649, 2.458117, -1.365235, 0.281047, 0.722535, 0.631630, 0.100000, 10.000000);
ray_helper_vec[40].make_ray(1.371671, 2.471542, -0.121285, -0.788661, 0.575374, 0.216700, 0.100000, 10.000000);
ray_helper_vec[41].make_ray(1.377863, 2.476803, -0.048615, -0.864242, -0.493562, 0.097383, 0.100000, 10.000000);
ray_helper_vec[42].make_ray(1.341712, 2.463960, -0.313598, -0.567863, 0.186239, 0.801777, 0.100000, 10.000000);
ray_helper_vec[43].make_ray(1.358384, 2.467112, -0.245688, -0.503415, 0.787405, -0.355763, 0.100000, 10.000000);
ray_helper_vec[44].make_ray(1.246343, 2.463216, -0.595870, -0.617979, -0.397439, 0.678339, 0.100000, 10.000000);
ray_helper_vec[45].make_ray(1.310967, 2.468275, -0.441872, -0.792183, -0.243903, 0.559426, 0.100000, 10.000000);
ray_helper_vec[46].make_ray(1.157707, 2.469567, -0.755047, -0.606396, -0.759642, 0.235008, 0.100000, 10.000000);
ray_helper_vec[47].make_ray(1.167095, 2.482026, -0.741990, 0.068254, -0.291419, 0.954157, 0.100000, 10.000000);
ray_helper_vec[48].make_ray(1.031651, 2.471097, -0.920947, 0.010529, 0.152242, 0.988287, 0.100000, 10.000000);
ray_helper_vec[49].make_ray(1.073173, 2.471776, -0.872496, -0.713158, 0.079450, 0.696486, 0.100000, 10.000000);
ray_helper_vec[50].make_ray(0.886163, 2.473118, -1.061647, -0.430068, 0.427306, 0.795268, 0.100000, 10.000000);
ray_helper_vec[51].make_ray(0.844121, 2.480847, -1.098548, -0.356027, 0.826389, 0.436263, 0.100000, 10.000000);
ray_helper_vec[52].make_ray(0.705698, 2.463287, -1.187380, -0.255099, 0.396584, 0.881842, 0.100000, 10.000000);
ray_helper_vec[53].make_ray(0.807865, 2.478608, -1.126200, -0.078894, 0.987282, 0.138021, 0.100000, 10.000000);
ray_helper_vec[54].make_ray(0.462249, 2.474161, -1.303073, -0.018424, 0.924636, 0.380405, 0.100000, 10.000000);
ray_helper_vec[55].make_ray(0.543106, 2.479522, -1.269947, 0.220949, 0.036251, 0.974611, 0.100000, 10.000000);
ray_helper_vec[56].make_ray(0.294246, 2.470012, -1.346954, 0.742198, 0.315986, 0.591012, 0.100000, 10.000000);
ray_helper_vec[57].make_ray(0.283935, 2.478410, -1.350239, -0.671659, 0.405223, 0.620217, 0.100000, 10.000000);
ray_helper_vec[58].make_ray(0.126889, 2.470006, -1.371089, 0.192754, 0.933363, 0.302787, 0.100000, 10.000000);
ray_helper_vec[59].make_ray(0.182643, 2.478538, -1.367503, -0.034233, 0.765739, 0.642240, 0.100000, 10.000000);
ray_helper_vec[60].make_ray(1.380018, 2.485112, -0.048965, -0.220701, 0.899924, -0.376068, 0.100000, 10.000000);
ray_helper_vec[61].make_ray(1.375973, 2.491078, -0.151805, 0.211478, 0.653196, -0.727058, 0.100000, 10.000000);
ray_helper_vec[62].make_ray(1.349111, 2.483034, -0.291301, -0.034464, 0.096714, 0.994715, 0.100000, 10.000000);
ray_helper_vec[63].make_ray(1.364611, 2.489255, -0.244892, -0.040405, 0.790668, 0.610910, 0.100000, 10.000000);
ray_helper_vec[64].make_ray(1.248300, 2.489001, -0.606949, -0.885067, 0.314769, 0.342895, 0.100000, 10.000000);
ray_helper_vec[65].make_ray(1.304689, 2.486138, -0.466482, -0.842632, 0.341482, -0.416367, 0.100000, 10.000000);
ray_helper_vec[66].make_ray(1.127481, 2.484725, -0.809132, 0.243172, 0.842423, 0.480823, 0.100000, 10.000000);
ray_helper_vec[67].make_ray(1.202785, 2.492563, -0.696821, 0.203217, 0.965374, 0.163573, 0.100000, 10.000000);
ray_helper_vec[68].make_ray(1.062570, 2.482835, -0.886604, -0.860361, 0.494356, 0.124059, 0.100000, 10.000000);
ray_helper_vec[69].make_ray(1.060781, 2.490168, -0.896321, -0.310583, 0.935670, 0.167513, 0.100000, 10.000000);
ray_helper_vec[70].make_ray(0.917813, 2.483467, -1.036441, 0.204139, 0.974820, 0.089744, 0.100000, 10.000000);
ray_helper_vec[71].make_ray(0.972034, 2.487555, -0.993702, -0.916429, 0.374401, -0.141360, 0.100000, 10.000000);
ray_helper_vec[72].make_ray(0.651707, 2.490490, -1.228894, -0.384847, -0.617438, 0.686049, 0.100000, 10.000000);
ray_helper_vec[73].make_ray(0.797068, 2.492015, -1.140683, -0.087846, 0.898606, 0.429872, 0.100000, 10.000000);
ray_helper_vec[74].make_ray(0.542252, 2.488451, -1.274833, 0.382534, -0.556152, 0.737809, 0.100000, 10.000000);
ray_helper_vec[75].make_ray(0.579456, 2.489432, -1.260060, -0.937874, 0.013995, 0.346694, 0.100000, 10.000000);
ray_helper_vec[76].make_ray(0.228225, 2.486313, -1.366590, 0.095984, 0.356418, 0.929384, 0.100000, 10.000000);
ray_helper_vec[77].make_ray(0.369185, 2.489439, -1.334772, -0.613698, 0.671932, 0.414586, 0.100000, 10.000000);
ray_helper_vec[78].make_ray(0.183972, 2.483694, -1.368463, -0.187164, -0.193178, 0.963147, 0.100000, 10.000000);
ray_helper_vec[79].make_ray(0.208855, 2.485208, -1.367532, -0.317205, -0.527531, 0.788094, 0.100000, 10.000000);
ray_helper_vec[80].make_ray(1.377277, 2.495839, -0.220625, -0.334351, 0.936996, 0.101236, 0.100000, 10.000000);
ray_helper_vec[81].make_ray(1.393668, 2.497807, -0.094088, 0.480358, 0.762963, 0.432601, 0.100000, 10.000000);
ray_helper_vec[82].make_ray(1.335210, 2.498160, -0.433082, 0.527606, 0.601598, 0.599760, 0.100000, 10.000000);
ray_helper_vec[83].make_ray(1.350583, 2.497788, -0.364122, -0.653715, 0.744602, 0.134998, 0.100000, 10.000000);
ray_helper_vec[84].make_ray(1.292860, 2.494987, -0.513733, -0.338343, 0.940854, -0.017839, 0.100000, 10.000000);
ray_helper_vec[85].make_ray(1.324345, 2.495937, -0.446215, -0.584069, 0.751729, -0.306214, 0.100000, 10.000000);
ray_helper_vec[86].make_ray(1.152184, 2.495265, -0.786371, 0.680488, 0.726684, -0.094159, 0.100000, 10.000000);
ray_helper_vec[87].make_ray(1.223908, 2.495914, -0.673846, -0.085170, 0.991273, 0.100617, 0.100000, 10.000000);
ray_helper_vec[88].make_ray(1.022041, 2.494758, -0.947425, -0.611477, 0.497245, 0.615502, 0.100000, 10.000000);
ray_helper_vec[89].make_ray(1.097948, 2.498325, -0.876705, -0.317200, 0.863081, 0.393034, 0.100000, 10.000000);
ray_helper_vec[90].make_ray(0.841284, 2.497874, -1.126290, 0.897914, 0.439261, 0.028271, 0.100000, 10.000000);
ray_helper_vec[91].make_ray(0.937168, 2.497470, -1.042539, -0.102715, 0.847229, -0.521203, 0.100000, 10.000000);
ray_helper_vec[92].make_ray(0.722526, 2.496311, -1.195491, -0.283216, 0.933406, 0.220321, 0.100000, 10.000000);
ray_helper_vec[93].make_ray(0.819016, 2.498011, -1.142689, -0.039089, 0.485326, 0.873459, 0.100000, 10.000000);
ray_helper_vec[94].make_ray(0.460286, 2.496364, -1.320002, -0.498023, 0.768744, 0.401256, 0.100000, 10.000000);
ray_helper_vec[95].make_ray(0.561499, 2.496422, -1.278069, 0.860379, 0.505319, 0.066344, 0.100000, 10.000000);
ray_helper_vec[96].make_ray(0.330476, 2.495541, -1.351122, 0.741613, 0.372096, 0.558171, 0.100000, 10.000000);
ray_helper_vec[97].make_ray(0.370211, 2.497694, -1.348798, -0.306476, -0.106618, 0.945889, 0.100000, 10.000000);
ray_helper_vec[98].make_ray(0.030874, 2.496973, -1.395895, -0.259933, 0.941118, -0.216176, 0.100000, 10.000000);
ray_helper_vec[99].make_ray(0.078709, 2.497484, -1.393815, 0.743893, 0.668106, 0.016040, 0.100000, 10.000000);
ray_helper_vec[100].make_ray(1.389637, 2.498231, -0.182966, 0.291508, 0.628142, -0.721430, 0.100000, 10.000000);
ray_helper_vec[101].make_ray(1.407923, 2.496617, -0.037087, 0.306742, 0.878371, 0.366570, 0.100000, 10.000000);
ray_helper_vec[102].make_ray(1.350553, 2.495141, -0.432767, -0.437070, 0.256984, 0.861933, 0.100000, 10.000000);
ray_helper_vec[103].make_ray(1.374077, 2.495482, -0.329057, 0.976833, 0.212803, -0.022642, 0.100000, 10.000000);
ray_helper_vec[104].make_ray(1.279408, 2.496475, -0.597102, -0.552490, 0.741073, 0.381530, 0.100000, 10.000000);
ray_helper_vec[105].make_ray(1.334499, 2.494647, -0.484911, 0.065493, 0.957576, 0.280639, 0.100000, 10.000000);
ray_helper_vec[106].make_ray(1.205203, 2.497835, -0.724941, 0.321350, 0.393066, -0.861530, 0.100000, 10.000000);
ray_helper_vec[107].make_ray(1.213744, 2.495152, -0.732755, -0.487747, 0.832042, 0.264214, 0.100000, 10.000000);
ray_helper_vec[108].make_ray(1.088106, 2.497151, -0.897178, 0.631448, 0.159791, -0.758776, 0.100000, 10.000000);
ray_helper_vec[109].make_ray(1.117482, 2.496978, -0.863959, 0.759390, 0.486867, 0.431610, 0.100000, 10.000000);
ray_helper_vec[110].make_ray(0.957617, 2.498112, -1.031052, -0.630111, 0.542692, -0.555378, 0.100000, 10.000000);
ray_helper_vec[111].make_ray(0.963412, 2.495492, -1.040709, 0.143284, 0.961392, 0.234936, 0.100000, 10.000000);
ray_helper_vec[112].make_ray(0.702112, 2.496653, -1.225126, -0.046348, 0.518090, -0.854070, 0.100000, 10.000000);
ray_helper_vec[113].make_ray(0.671739, 2.494551, -1.254267, 0.521797, 0.205009, -0.828070, 0.100000, 10.000000);
ray_helper_vec[114].make_ray(0.550733, 2.498400, -1.289891, -0.184121, 0.398183, 0.898638, 0.100000, 10.000000);
ray_helper_vec[115].make_ray(0.563605, 2.495706, -1.296879, -0.544225, 0.783405, 0.300158, 0.100000, 10.000000);
ray_helper_vec[116].make_ray(0.271127, 2.495705, -1.387073, -0.115337, 0.742379, -0.659978, 0.100000, 10.000000);
ray_helper_vec[117].make_ray(0.319871, 2.495412, -1.376595, -0.264170, 0.580593, -0.770147, 0.100000, 10.000000);
ray_helper_vec[118].make_ray(0.082511, 2.496985, -1.402798, 0.061583, 0.996348, 0.059143, 0.100000, 10.000000);
ray_helper_vec[119].make_ray(0.120508, 2.495697, -1.405272, 0.623816, 0.778405, 0.070289, 0.100000, 10.000000);
ray_helper_vec[120].make_ray(1.419548, 2.494145, -0.010557, -0.135449, 0.304631, -0.942790, 0.100000, 10.000000);
ray_helper_vec[121].make_ray(1.422235, 2.483898, -0.188780, 0.950847, -0.285924, -0.118904, 0.100000, 10.000000);
ray_helper_vec[122].make_ray(1.373723, 2.493747, -0.353447, 0.732360, 0.249652, -0.633501, 0.100000, 10.000000);
ray_helper_vec[123].make_ray(1.393683, 2.488341, -0.308120, 0.635669, 0.565814, 0.525147, 0.100000, 10.000000);
ray_helper_vec[124].make_ray(1.319263, 2.492214, -0.532796, 0.673715, 0.702009, -0.230848, 0.100000, 10.000000);
ray_helper_vec[125].make_ray(1.349388, 2.486781, -0.483429, -0.493646, 0.648766, -0.579152, 0.100000, 10.000000);
ray_helper_vec[126].make_ray(1.239760, 2.493308, -0.699437, 0.602142, 0.047661, -0.796965, 0.100000, 10.000000);
ray_helper_vec[127].make_ray(1.253273, 2.486987, -0.697111, -0.121360, 0.701837, 0.701923, 0.100000, 10.000000);
ray_helper_vec[128].make_ray(1.044914, 2.490230, -0.975668, 0.252263, 0.958325, 0.134079, 0.100000, 10.000000);
ray_helper_vec[129].make_ray(1.063528, 2.483907, -0.969755, -0.366889, 0.927729, 0.068641, 0.100000, 10.000000);
ray_helper_vec[130].make_ray(0.893717, 2.492407, -1.110337, 0.852452, 0.479208, -0.209010, 0.100000, 10.000000);
ray_helper_vec[131].make_ray(0.926265, 2.484474, -1.099506, 0.528661, 0.409777, -0.743371, 0.100000, 10.000000);
ray_helper_vec[132].make_ray(0.809478, 2.492574, -1.173471, 0.742706, 0.553597, 0.376721, 0.100000, 10.000000);
ray_helper_vec[133].make_ray(0.726268, 2.482974, -1.243009, 0.416412, 0.753715, -0.508444, 0.100000, 10.000000);
ray_helper_vec[134].make_ray(0.536873, 2.490142, -1.321216, -0.331112, 0.342631, -0.879187, 0.100000, 10.000000);
ray_helper_vec[135].make_ray(0.587117, 2.488577, -1.303055, 0.094370, 0.986701, 0.132344, 0.100000, 10.000000);
ray_helper_vec[136].make_ray(0.239431, 2.486864, -1.412711, 0.946983, 0.214858, -0.238873, 0.100000, 10.000000);
ray_helper_vec[137].make_ray(0.337335, 2.486137, -1.390321, -0.115398, 0.309050, -0.944019, 0.100000, 10.000000);
ray_helper_vec[138].make_ray(0.045275, 2.490427, -1.422874, -0.003551, 0.879041, -0.476732, 0.100000, 10.000000);
ray_helper_vec[139].make_ray(0.078473, 2.483393, -1.431716, 0.632476, 0.195796, -0.749425, 0.100000, 10.000000);
ray_helper_vec[140].make_ray(1.424597, 2.478899, -0.234223, 0.539899, 0.821233, -0.184623, 0.100000, 10.000000);
ray_helper_vec[141].make_ray(1.442555, 2.464552, -0.196476, 0.850011, -0.384822, 0.359712, 0.100000, 10.000000);
ray_helper_vec[142].make_ray(1.397087, 2.471287, -0.383920, -0.014410, 0.605137, 0.795991, 0.100000, 10.000000);
ray_helper_vec[143].make_ray(1.416144, 2.480077, -0.265760, 0.843418, 0.461838, 0.274502, 0.100000, 10.000000);
ray_helper_vec[144].make_ray(1.348927, 2.480955, -0.506571, 0.095227, 0.973100, -0.209779, 0.100000, 10.000000);
ray_helper_vec[145].make_ray(1.361048, 2.465232, -0.520184, 0.278929, 0.955370, -0.097296, 0.100000, 10.000000);
ray_helper_vec[146].make_ray(1.222714, 2.478408, -0.768691, 0.094663, 0.989658, 0.107779, 0.100000, 10.000000);
ray_helper_vec[147].make_ray(1.268401, 2.478305, -0.694629, 0.641592, 0.464691, -0.610264, 0.100000, 10.000000);
ray_helper_vec[148].make_ray(1.075523, 2.475763, -0.969931, 0.821407, 0.484401, -0.301075, 0.100000, 10.000000);
ray_helper_vec[149].make_ray(1.159623, 2.470558, -0.879978, -0.043256, 0.999002, -0.011180, 0.100000, 10.000000);
ray_helper_vec[150].make_ray(0.917621, 2.476925, -1.118663, -0.785228, 0.299206, -0.542119, 0.100000, 10.000000);
ray_helper_vec[151].make_ray(1.010612, 2.470565, -1.047886, -0.205584, 0.877045, 0.434197, 0.100000, 10.000000);
ray_helper_vec[152].make_ray(0.682924, 2.474653, -1.280084, 0.448440, -0.095903, -0.888653, 0.100000, 10.000000);
ray_helper_vec[153].make_ray(0.799554, 2.473957, -1.209193, 0.444545, 0.848380, -0.287457, 0.100000, 10.000000);
ray_helper_vec[154].make_ray(0.480471, 2.471260, -1.370762, -0.428292, 0.123379, -0.895178, 0.100000, 10.000000);
ray_helper_vec[155].make_ray(0.615711, 2.469254, -1.316725, -0.155417, -0.227108, -0.961388, 0.100000, 10.000000);
ray_helper_vec[156].make_ray(0.350267, 2.473876, -1.402423, 0.379992, 0.772358, -0.508988, 0.100000, 10.000000);
ray_helper_vec[157].make_ray(0.266496, 2.464126, -1.433081, 0.430264, 0.902444, -0.021617, 0.100000, 10.000000);
ray_helper_vec[158].make_ray(0.027307, 2.472506, -1.447521, -0.087110, 0.978023, 0.189428, 0.100000, 10.000000);
ray_helper_vec[159].make_ray(0.118244, 2.467414, -1.445705, -0.189151, 0.015463, -0.981826, 0.100000, 10.000000);
ray_helper_vec[160].make_ray(1.456477, 2.451214, -0.152684, -0.284511, 0.955668, -0.075846, 0.100000, 10.000000);
ray_helper_vec[161].make_ray(1.478017, 2.437023, -0.014216, 0.854063, 0.295056, 0.428390, 0.100000, 10.000000);
ray_helper_vec[162].make_ray(1.404747, 2.445600, -0.444747, 0.362028, 0.477450, -0.800611, 0.100000, 10.000000);
ray_helper_vec[163].make_ray(1.452761, 2.443533, -0.252215, 0.879524, 0.114128, -0.461965, 0.100000, 10.000000);
ray_helper_vec[164].make_ray(1.339111, 2.446979, -0.610264, 0.715233, 0.675882, -0.177839, 0.100000, 10.000000);
ray_helper_vec[165].make_ray(1.362350, 2.436308, -0.575310, 0.571130, 0.818314, 0.064594, 0.100000, 10.000000);
ray_helper_vec[166].make_ray(1.241895, 2.459467, -0.773331, 0.384185, 0.755395, -0.530829, 0.100000, 10.000000);
ray_helper_vec[167].make_ray(1.296306, 2.446749, -0.703140, 0.945939, 0.286902, -0.151286, 0.100000, 10.000000);
ray_helper_vec[168].make_ray(1.104706, 2.451327, -0.969871, -0.524855, 0.792172, -0.311433, 0.100000, 10.000000);
ray_helper_vec[169].make_ray(1.133207, 2.447800, -0.940632, 0.398737, 0.001329, -0.917064, 0.100000, 10.000000);
ray_helper_vec[170].make_ray(0.880642, 2.456492, -1.175905, 0.100864, 0.947814, -0.302449, 0.100000, 10.000000);
ray_helper_vec[171].make_ray(1.015152, 2.458698, -1.058678, 0.211479, 0.532367, -0.819672, 0.100000, 10.000000);
ray_helper_vec[172].make_ray(0.821253, 2.457471, -1.214180, 0.613688, -0.301231, -0.729826, 0.100000, 10.000000);
ray_helper_vec[173].make_ray(0.816212, 2.437855, -1.234605, 0.729172, 0.549202, -0.408271, 0.100000, 10.000000);
ray_helper_vec[174].make_ray(0.530850, 2.458932, -1.362444, 0.192239, 0.981195, -0.017308, 0.100000, 10.000000);
ray_helper_vec[175].make_ray(0.637866, 2.442541, -1.331231, 0.560147, -0.317806, -0.765007, 0.100000, 10.000000);
ray_helper_vec[176].make_ray(0.370324, 2.456438, -1.414323, -0.735730, 0.659356, -0.154758, 0.100000, 10.000000);
ray_helper_vec[177].make_ray(0.429430, 2.441175, -1.411855, -0.637937, -0.039385, -0.769081, 0.100000, 10.000000);
ray_helper_vec[178].make_ray(0.040708, 2.459337, -1.459159, 0.844893, 0.207555, -0.493027, 0.100000, 10.000000);
ray_helper_vec[179].make_ray(0.130963, 2.445404, -1.462551, 0.734987, 0.330305, -0.592193, 0.100000, 10.000000);
ray_helper_vec[180].make_ray(1.471978, 2.434966, -0.109661, 0.468735, 0.568971, 0.675692, 0.100000, 10.000000);
ray_helper_vec[181].make_ray(1.493566, 2.408638, -0.020869, 0.494409, -0.168669, -0.852708, 0.100000, 10.000000);
ray_helper_vec[182].make_ray(1.422364, 2.418210, -0.444945, 0.142539, 0.527392, 0.837580, 0.100000, 10.000000);
ray_helper_vec[183].make_ray(1.444876, 2.402413, -0.388922, 0.499289, 0.676342, -0.541546, 0.100000, 10.000000);
ray_helper_vec[184].make_ray(1.354838, 2.419467, -0.618110, 0.931311, 0.249661, 0.265195, 0.100000, 10.000000);
ray_helper_vec[185].make_ray(1.397680, 2.409196, -0.530037, 0.094789, 0.951252, -0.293486, 0.100000, 10.000000);
ray_helper_vec[186].make_ray(1.252441, 2.429477, -0.797009, 0.934642, 0.344788, -0.086977, 0.100000, 10.000000);
ray_helper_vec[187].make_ray(1.269262, 2.416852, -0.783048, 0.782634, -0.012938, -0.622347, 0.100000, 10.000000);
ray_helper_vec[188].make_ray(1.115018, 2.417997, -0.991129, -0.472600, 0.495252, -0.728955, 0.100000, 10.000000);
ray_helper_vec[189].make_ray(1.102796, 2.405776, -1.015869, 0.868085, 0.464404, -0.175378, 0.100000, 10.000000);
ray_helper_vec[190].make_ray(0.941887, 2.415640, -1.158857, 0.145006, 0.903801, -0.402639, 0.100000, 10.000000);
ray_helper_vec[191].make_ray(1.035078, 2.402729, -1.088596, 0.539098, 0.074614, -0.838932, 0.100000, 10.000000);
ray_helper_vec[192].make_ray(0.704047, 2.418029, -1.317092, 0.354264, -0.156758, -0.921913, 0.100000, 10.000000);
ray_helper_vec[193].make_ray(0.753060, 2.407225, -1.293989, 0.983486, 0.038450, 0.176855, 0.100000, 10.000000);
ray_helper_vec[194].make_ray(0.561610, 2.420023, -1.378028, 0.463521, -0.261676, -0.846566, 0.100000, 10.000000);
ray_helper_vec[195].make_ray(0.680606, 2.405461, -1.337237, 0.498553, 0.535128, -0.681970, 0.100000, 10.000000);
ray_helper_vec[196].make_ray(0.276402, 2.422688, -1.460446, -0.306577, 0.731260, -0.609319, 0.100000, 10.000000);
ray_helper_vec[197].make_ray(0.389666, 2.419404, -1.435011, 0.773831, 0.252364, -0.580946, 0.100000, 10.000000);
ray_helper_vec[198].make_ray(0.090710, 2.423551, -1.479803, -0.454611, 0.835152, -0.309596, 0.100000, 10.000000);
ray_helper_vec[199].make_ray(0.232971, 2.402086, -1.480554, 0.559620, 0.453501, -0.693659, 0.100000, 10.000000);
ray_helper_vec[200].make_ray(-0.134246, 2.410408, -1.385328, -0.297511, -0.950003, 0.094768, 0.100000, 10.000000);
ray_helper_vec[201].make_ray(-0.045810, 2.420457, -1.388279, 0.048892, -0.547323, 0.835492, 0.100000, 10.000000);
ray_helper_vec[202].make_ray(-0.344609, 2.406566, -1.351590, -0.391912, -0.918047, 0.059959, 0.100000, 10.000000);
ray_helper_vec[203].make_ray(-0.260126, 2.421040, -1.366056, -0.431592, -0.800584, 0.415683, 0.100000, 10.000000);
ray_helper_vec[204].make_ray(-0.576613, 2.403958, -1.273855, 0.418521, 0.515239, 0.747910, 0.100000, 10.000000);
ray_helper_vec[205].make_ray(-0.474065, 2.431474, -1.304719, 0.443555, -0.309531, 0.841100, 0.100000, 10.000000);
ray_helper_vec[206].make_ray(-0.781897, 2.425588, -1.151560, -0.063162, 0.456468, 0.887495, 0.100000, 10.000000);
ray_helper_vec[207].make_ray(-0.691185, 2.427330, -1.206549, -0.394186, -0.282794, 0.874440, 0.100000, 10.000000);
ray_helper_vec[208].make_ray(-0.904999, 2.410619, -1.064565, 0.716448, -0.677515, -0.166363, 0.100000, 10.000000);
ray_helper_vec[209].make_ray(-0.842814, 2.430272, -1.107467, 0.667364, -0.652175, 0.359573, 0.100000, 10.000000);
ray_helper_vec[210].make_ray(-1.060995, 2.415567, -0.906146, 0.025499, -0.834067, 0.551074, 0.100000, 10.000000);
ray_helper_vec[211].make_ray(-1.065577, 2.431562, -0.891004, 0.821417, 0.501850, -0.270963, 0.100000, 10.000000);
ray_helper_vec[212].make_ray(-1.234715, 2.424866, -0.647267, -0.130090, -0.263379, 0.955881, 0.100000, 10.000000);
ray_helper_vec[213].make_ray(-1.168000, 2.424964, -0.755646, 0.393687, -0.882698, 0.256622, 0.100000, 10.000000);
ray_helper_vec[214].make_ray(-1.310971, 2.400304, -0.491240, 0.573223, -0.744893, 0.341394, 0.100000, 10.000000);
ray_helper_vec[215].make_ray(-1.257030, 2.428040, -0.592155, 0.758164, -0.643481, 0.105448, 0.100000, 10.000000);
ray_helper_vec[216].make_ray(-1.362515, 2.423745, -0.270236, 0.673934, -0.098869, 0.732146, 0.100000, 10.000000);
ray_helper_vec[217].make_ray(-1.326944, 2.409286, -0.442231, 0.550306, -0.561010, 0.618410, 0.100000, 10.000000);
ray_helper_vec[218].make_ray(-1.384918, 2.421718, -0.082240, 0.648321, 0.654407, 0.389141, 0.100000, 10.000000);
ray_helper_vec[219].make_ray(-1.378744, 2.423402, -0.152380, 0.842665, -0.179927, -0.507486, 0.100000, 10.000000);
ray_helper_vec[220].make_ray(-0.136576, 2.450971, -1.372090, -0.077415, -0.885134, 0.458852, 0.100000, 10.000000);
ray_helper_vec[221].make_ray(-0.060629, 2.457107, -1.376822, -0.767855, -0.248971, 0.590265, 0.100000, 10.000000);
ray_helper_vec[222].make_ray(-0.340977, 2.441322, -1.339481, 0.197120, -0.166908, 0.966067, 0.100000, 10.000000);
ray_helper_vec[223].make_ray(-0.270395, 2.459050, -1.352849, -0.061596, -0.744617, 0.664644, 0.100000, 10.000000);
ray_helper_vec[224].make_ray(-0.610606, 2.436899, -1.245851, -0.089958, 0.329884, 0.939726, 0.100000, 10.000000);
ray_helper_vec[225].make_ray(-0.573078, 2.460800, -1.256288, 0.130534, 0.531777, 0.836764, 0.100000, 10.000000);
ray_helper_vec[226].make_ray(-0.764596, 2.441509, -1.156178, 0.390150, -0.259726, 0.883361, 0.100000, 10.000000);
ray_helper_vec[227].make_ray(-0.752067, 2.459263, -1.159705, -0.157984, -0.833826, 0.528937, 0.100000, 10.000000);
ray_helper_vec[228].make_ray(-0.911251, 2.449453, -1.042532, 0.244028, -0.729491, 0.638978, 0.100000, 10.000000);
ray_helper_vec[229].make_ray(-0.847704, 2.461054, -1.093819, 0.813343, -0.018754, 0.581482, 0.100000, 10.000000);
ray_helper_vec[230].make_ray(-1.033293, 2.436273, -0.926122, 0.255715, -0.777318, 0.574792, 0.100000, 10.000000);
ray_helper_vec[231].make_ray(-1.041997, 2.458355, -0.909128, 0.167285, -0.011984, 0.985836, 0.100000, 10.000000);
ray_helper_vec[232].make_ray(-1.197917, 2.436435, -0.698687, 0.924554, -0.380872, -0.011703, 0.100000, 10.000000);
ray_helper_vec[233].make_ray(-1.197855, 2.462609, -0.688769, -0.063678, -0.493763, 0.867262, 0.100000, 10.000000);
ray_helper_vec[234].make_ray(-1.312239, 2.449799, -0.444425, 0.861186, -0.378050, -0.339760, 0.100000, 10.000000);
ray_helper_vec[235].make_ray(-1.244292, 2.455151, -0.604834, 0.593733, -0.328040, 0.734760, 0.100000, 10.000000);
ray_helper_vec[236].make_ray(-1.332231, 2.443689, -0.369027, 0.442339, -0.843032, 0.305995, 0.100000, 10.000000);
ray_helper_vec[237].make_ray(-1.344812, 2.456559, -0.305856, 0.623447, -0.779304, 0.063240, 0.100000, 10.000000);
ray_helper_vec[238].make_ray(-1.378955, 2.449168, -0.053729, 0.981087, 0.185502, 0.055295, 0.100000, 10.000000);
ray_helper_vec[239].make_ray(-1.374021, 2.449745, -0.115104, 0.493733, -0.812873, -0.308975, 0.100000, 10.000000);
*/