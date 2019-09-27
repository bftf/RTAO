#include "helper_math.h"

class RayGenerator;
class BVHManager;

class RayTraceManager {
  public:
    RayTraceManager(RayGenerator& rg);
    RayTraceManager() {};

    void traceOptiXPrime(RayGenerator& rg);
    void traceAila(RayGenerator& rg);
    void traceCWBVH(RayGenerator& rg);
    void evaluate(RayGenerator& rg);

  private:

    uint numRays;
    Hit* cudaHits;
};