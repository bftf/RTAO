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

void RayGenerator::loadModelOBJ(const std::string& model_path)
{

  objl::Loader loader;
  assert(loader.LoadFile(model_path));

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
}

void RayGenerator::printRaysForVisualization()
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

void RayGenerator::uploadRaysToGPU()
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