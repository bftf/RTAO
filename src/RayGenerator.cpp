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
    NormalBuffer.push_back(invertNormal * cur_normal/length_cur_normal); // HACK HACK HACK
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

void RayGenerator::downsizeRayVector(const uint64_t number_of_rays)
{
  assert(m_ray_sorting_strategy == random_shuffle);
  ray_helper_vec.resize(number_of_rays);
}

void RayGenerator::uploadRaysToGPU()
{
  cudaRays = GenerateRaysFromFile(ray_helper_vec, ray_helper_vec.size());
}

void RayGenerator::fillWithListOfRays()
{
  ray_helper_vec.resize(50);
  int i = 0; 

ray_helper_vec[i].make_ray(539.231323, 198.971649, 204.659210, 0.901625, 0.370240, -0.223594, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(806.181396, 158.500748, -276.240814, -0.039514, -0.682485, -0.729831, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(-666.666199, 290.715454, -280.687988, -0.478790, 0.821487, -0.309708, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(-1039.442627, 674.629150, 182.131470, -0.814641, 0.476629, -0.330431, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(-639.919983, 123.229416, -229.858459, -0.524756, -0.298174, -0.797322, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(674.080200, 106.372047, 201.033081, 0.900200, 0.159620, -0.405168, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(-345.748535, 393.215942, -238.874512, 0.592021, 0.286174, -0.753402, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(518.723755, 200.837112, -645.093018, 0.540266, 0.739437, -0.401678, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(1152.068970, 131.655334, 347.817749, -0.386070, 0.821806, -0.419028, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(-167.892609, 512.504211, 182.253250, 0.405867, 0.904507, 0.130915, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(-464.455444, 511.735046, 218.740768, 0.995961, -0.069651, 0.056660, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(-220.943512, 172.881821, 180.896973, -0.397235, -0.699856, -0.593638, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(489.539673, 107.036201, -220.235626, 0.214076, -0.887672, -0.407689, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(174.476562, 450.102295, 162.743042, 0.357323, -0.896706, -0.261226, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(-681.083191, 486.827454, 161.382095, -0.283187, -0.957095, -0.061443, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(112.707443, 122.987320, 234.094879, -0.065004, -0.771857, 0.632465, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(1160.653442, 376.500793, 18.280207, -0.049902, -0.975628, -0.213681, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(510.528442, 133.237457, -218.204681, 0.420336, 0.498902, 0.757902, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(88.151955, 126.069489, 220.765656, -0.952481, 0.033372, -0.302763, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(-450.642517, 676.236450, -293.027832, -0.975175, 0.210094, -0.069965, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(-208.003204, 172.822861, 206.457458, 0.633983, -0.647134, 0.423419, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(492.537170, 846.500488, 159.650848, -0.238495, -0.224538, 0.944830, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(-329.631653, 265.009338, -279.456543, 0.005927, 0.678398, -0.734671, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(-755.765625, 512.504211, 210.231171, -0.438008, 0.890666, 0.121912, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(-777.207581, 236.728821, -271.497681, -0.257632, -0.770708, 0.582782, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(60.983471, 461.603210, 161.231445, 0.488689, -0.855586, -0.170751, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(396.922729, 440.277832, -238.567780, 0.139870, -0.026229, -0.989822, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(485.435303, 101.288239, 144.192032, -0.293868, -0.795055, 0.530593, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(482.536133, 168.317368, 227.481415, -0.343009, -0.891052, -0.297272, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(-1289.165649, 1397.564331, -699.827881, 0.240359, -0.970048, -0.035129, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(-196.120483, 161.241623, 199.261261, -0.028170, -0.765883, -0.642363, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(-273.029419, 155.501968, 180.590988, -0.799474, -0.398845, -0.449181, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(112.559555, 25.668892, -225.175323, -0.370606, -0.517494, -0.771266, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(-593.292542, 384.034637, -238.763672, 0.513614, -0.793446, -0.326566, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(-257.280701, 1257.654297, -255.967117, 0.143681, -0.756651, 0.637836, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(-250.407913, 354.360046, -234.496445, 0.111284, 0.641038, 0.759398, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(812.876221, 287.619354, 195.726837, -0.391578, 0.595090, -0.701808, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(841.026001, 39.495056, -238.539337, 0.812093, 0.580554, 0.058834, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(-584.725891, 127.543472, -221.312897, 0.942757, 0.296819, -0.152013, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(-215.163849, 170.872620, -303.035370, 0.093228, 0.835610, -0.541355, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(856.716736, 48.432487, -213.840744, 0.775054, -0.304437, 0.553724, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(-238.895996, 55.726707, 170.768982, 0.267282, 0.660908, 0.701257, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(-1231.386108, 83.354019, 472.906921, 0.016037, 0.762957, -0.646250, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(-298.555420, 431.608826, 162.162003, -0.639714, -0.635952, -0.431662, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(365.803589, 63.834991, 194.711639, -0.114968, 0.265883, -0.957125, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(1308.769775, 687.891296, -26.383528, 0.619478, -0.699526, -0.356245, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(-839.464600, 254.645599, -270.074127, -0.305232, 0.347513, 0.886605, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(-1423.480225, 250.247543, 150.899155, 0.347356, -0.465243, -0.814182, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(-262.780518, 63.595718, 180.296509, -0.198311, 0.455250, -0.867998, 0.100000, 10.000000); i++;
ray_helper_vec[i].make_ray(-1339.224854, 1398.982300, -698.798035, 0.870580, 0.245951, -0.426143, 0.100000, 10.000000); i++;

/*
  10.00 -1  539.231323 198.971649 204.659210 0.901625 0.370240 -0.223594 
10.00 -1  806.181396 158.500748 -276.240814 -0.039514 -0.682485 -0.729831 
10.00 -1  -666.666199 290.715454 -280.687988 -0.478790 0.821487 -0.309708 
10.00 -1  -1039.442627 674.629150 182.131470 -0.814641 0.476629 -0.330431 
10.00 -1  -639.919983 123.229416 -229.858459 -0.524756 -0.298174 -0.797322 
10.00 -1  674.080200 106.372047 201.033081 0.900200 0.159620 -0.405168 
10.00 -1  -345.748535 393.215942 -238.874512 0.592021 0.286174 -0.753402 
10.00 -1  518.723755 200.837112 -645.093018 0.540266 0.739437 -0.401678 
10.00 -1  1152.068970 131.655334 347.817749 -0.386070 0.821806 -0.419028 
10.00 -1  -167.892609 512.504211 182.253250 0.405867 0.904507 0.130915 
6.47 899346 -464.455444 511.735046 218.740768 0.995961 -0.069651 0.056660 
7.72 845226 -220.943512 172.881821 180.896973 -0.397235 -0.699856 -0.593638 
0.30 460851 489.539673 107.036201 -220.235626 0.214076 -0.887672 -0.407689 
10.00 -1  174.476562 450.102295 162.743042 0.357323 -0.896706 -0.261226 
10.00 -1  -681.083191 486.827454 161.382095 -0.283187 -0.957095 -0.061443 
10.00 -1  112.707443 122.987320 234.094879 -0.065004 -0.771857 0.632465 
0.00 0  1160.653442 376.500793 18.280207 -0.049902 -0.975628 -0.213681 
10.00 -1  510.528442 133.237457 -218.204681 0.420336 0.498902 0.757902 
10.00 -1  88.151955 126.069489 220.765656 -0.952481 0.033372 -0.302763 
10.00 -1  -450.642517 676.236450 -293.027832 -0.975175 0.210094 -0.069965 
3.94 840672 -208.003204 172.822861 206.457458 0.633983 -0.647134 0.423419 
10.00 -1  492.537170 846.500488 159.650848 -0.238495 -0.224538 0.944830 
10.00 -1  -329.631653 265.009338 -279.456543 0.005927 0.678398 -0.734671 
10.00 -1  -755.765625 512.504211 210.231171 -0.438008 0.890666 0.121912 
10.00 -1  -777.207581 236.728821 -271.497681 -0.257632 -0.770708 0.582782 
10.00 -1  60.983471 461.603210 161.231445 0.488689 -0.855586 -0.170751 
10.00 -1  396.922729 440.277832 -238.567780 0.139870 -0.026229 -0.989822 
10.00 -1  485.435303 101.288239 144.192032 -0.293868 -0.795055 0.530593 
10.00 -1  482.536133 168.317368 227.481415 -0.343009 -0.891052 -0.297272 
4.05 1088988  -1289.165649 1397.564331 -699.827881 0.240359 -0.970048 -0.035129 
10.00 -1  -196.120483 161.241623 199.261261 -0.028170 -0.765883 -0.642363 
10.00 -1  -273.029419 155.501968 180.590988 -0.799474 -0.398845 -0.449181 
10.00 -1  112.559555 25.668892 -225.175323 -0.370606 -0.517494 -0.771266 
10.00 -1  -593.292542 384.034637 -238.763672 0.513614 -0.793446 -0.326566 
10.00 -1  -257.280701 1257.654297 -255.967117 0.143681 -0.756651 0.637836 
10.00 -1  -250.407913 354.360046 -234.496445 0.111284 0.641038 0.759398 
10.00 -1  812.876221 287.619354 195.726837 -0.391578 0.595090 -0.701808 
10.00 -1  841.026001 39.495056 -238.539337 0.812093 0.580554 0.058834 
10.00 -1  -584.725891 127.543472 -221.312897 0.942757 0.296819 -0.152013 
3.85 666885 -215.163849 170.872620 -303.035370 0.093228 0.835610 -0.541355 
10.00 -1  856.716736 48.432487 -213.840744 0.775054 -0.304437 0.553724 
0.48 872007 -238.895996 55.726707 170.768982 0.267282 0.660908 0.701257 
10.00 -1  -1231.386108 83.354019 472.906921 0.016037 0.762957 -0.646250 
10.00 -1  -298.555420 431.608826 162.162003 -0.639714 -0.635952 -0.431662 
10.00 -1  365.803589 63.834991 194.711639 -0.114968 0.265883 -0.957125 
9.92 57354  1308.769775 687.891296 -26.383528 0.619478 -0.699526 -0.356245 
10.00 -1  -839.464600 254.645599 -270.074127 -0.305232 0.347513 0.886605 
0.00 0  -1423.480225 250.247543 150.899155 0.347356 -0.465243 -0.814182 
10.00 -1  -262.780518 63.595718 180.296509 -0.198311 0.455250 -0.867998 
10.00 -1  -1339.224854 1398.982300 -698.798035 0.870580 0.245951 -0.426143
  */
  assert(i == ray_helper_vec.size()); 

/*
  for( int i = 0; i < ray_helper_vec.size(); i++)
  {
    printf("ray_helper_vec[i].make_ray(%f, %f, %f, %f, %f, %f, %f, %f); i++;\n",
      ray_helper_vec[i].origin_tmin.x, ray_helper_vec[i].origin_tmin.y, ray_helper_vec[i].origin_tmin.z,
      ray_helper_vec[i].dir_tmax.x, ray_helper_vec[i].dir_tmax.y, ray_helper_vec[i].dir_tmax.z,
      ray_helper_vec[i].origin_tmin.w, ray_helper_vec[i].dir_tmax.w
    );
  }
*/

}