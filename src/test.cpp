#include <math.h>
#include <iostream>
#include <string>
#include "particle_filter.h"

// for convenience
using std::string;
using std::vector;

int main{
  // Create particle filter
  ParticleFilter pf;
  double std[] = {2,2,0.05};
  pf.init(4983, 5029, 1.201, std[]);
  return 0;
}