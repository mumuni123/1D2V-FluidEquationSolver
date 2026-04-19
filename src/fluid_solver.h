#ifndef FLUID_SOLVER_H
#define FLUID_SOLVER_H

#include <vector>

namespace FluidSolver {

void update_derived_from_n_vx_qz(int nx,
                                      const std::vector<double>& n,
                                      const std::vector<double>& vx,
                                      const std::vector<double>& qz,
                                      std::vector<double>& gamma,
                                      std::vector<double>& inv_n,
                                      std::vector<double>& vz,
                                      std::vector<double>& jx,
                                      std::vector<double>& jz);

void update_n(int nx,
              double dt,
              double dz,
              double n_bath,
              double vth_tilde,
              const std::vector<double>& n_old,
              const std::vector<double>& vz,
              const std::vector<double>& jz,
              std::vector<double>& n_new);

void update_vx(int nx,
               double dt,
               double dz,
               const std::vector<double>& Ex,
               const std::vector<double>& By,
               const std::vector<double>& vz,
               std::vector<double>& vx);

void update_qz(int nx,
               double dt,
               double dz,
               double beta,
               const std::vector<double>& Ez,
               const std::vector<double>& By,
               const std::vector<double>& n,
               const std::vector<double>& inv_n,
               const std::vector<double>& vx,
               const std::vector<double>& vz,
               std::vector<double>& qz);

void apply_thermal_both_sides(int nx,
                              double vth_tilde,
                              std::vector<double>& n,
                              std::vector<double>& vx,
                              std::vector<double>& qz);

void enforce_continuity_constraint_from_jz(int nx,
                                           double dt,
                                           double dz,
                                           double n_bath,
                                           double vth_tilde,
                                           const std::vector<double>& n_old,
                                           const std::vector<double>& vz,
                                           const std::vector<double>& jz,
                                           std::vector<double>& n_new);

}

#endif
