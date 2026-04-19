#ifndef PHYSICS_CONSTANTS_H
#define PHYSICS_CONSTANTS_H

#include <cmath>

namespace PhysConst {

const double EPS0 = 8.8541878128e-12;
const double MU0 = 1.25663706212e-6;
const double QE = 1.602176634e-19;
const double ME = 9.1093837015e-31;
const double KB = 1.380649e-23;
const double EV = 1.602176634e-19;
const double PI = 3.14159265358979323846;
const double C = 1.0 / std::sqrt(EPS0 * MU0);

}

#endif
