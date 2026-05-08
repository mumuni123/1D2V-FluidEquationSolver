#ifndef MAXWELL_SOLVER_H
#define MAXWELL_SOLVER_H

#include <vector>

namespace MaxwellSolver {

void update_magnetic_field(int nx,
                           double dt,
                           double dz,
                           const std::vector<double>& Ex,
                           bool drive_left_boundary,
                           double left_laser_ex_t,
                           bool drive_right_boundary,
                           double right_laser_ex_t,
                           std::vector<double>& By);

void update_transverse_electric_field(int nx,
                                      double dt,
                                      double dz,
                                      bool drive_left_boundary,
                                      double left_laser_ex_t,
                                      bool drive_right_boundary,
                                      double right_laser_ex_t,
                                      const std::vector<double>& jx,
                                      const std::vector<double>& By,
                                      std::vector<double>& Ex);

void update_longitudinal_electric_field(int nx,
                                        double dt,
                                        const std::vector<double>& jz,
                                        std::vector<double>& Ez);

}

#endif
