#include "helper_math.h"
#include <string>

class RayGenerator;
class BVHManager;

class RayTraceManager {
  public:
    RayTraceManager(RayGenerator& rg);
    ~RayTraceManager() {};

    void traceOptiXPrime(RayGenerator& rg);
    void traceAila(RayGenerator& rg);
    void traceCWBVH(RayGenerator& rg);
    void traceCWBVHSingleRay(RayGenerator& rg);
    void traceOptiX(RayGenerator& rg);

    void evaluateAndPrintForPLYVisualization(RayGenerator& rg, const std::string& out_ply_path);
    void debugging(RayGenerator& rg);

  private:

    uint numRays;
    Hit* cudaHits;
};