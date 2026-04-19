#ifndef MAXWELL_SOLVER_H
#define MAXWELL_SOLVER_H

#include <vector>

namespace MaxwellSolver {

void update_magnetic_field(int nx,
                           double dt,
                           double dz,
                           const std::vector<double>& Ex,
                           std::vector<double>& By);

void update_transverse_electric_field(int nx,
                                      double dt,
                                      double dz,
                                      double laser_ex_t,
                                      const std::vector<double>& jx,
                                      const std::vector<double>& By,
                                      std::vector<double>& Ex);

void update_longitudinal_electric_field(int nx,
                                        double dt,
                                        const std::vector<double>& jz,
                                        std::vector<double>& Ez);

}

#endif
