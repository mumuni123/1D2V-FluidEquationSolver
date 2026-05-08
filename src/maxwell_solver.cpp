#include "maxwell_solver.h"
namespace {

void apply_zero_gradient(std::vector<double>& a) {
    if (a.size() >= 2U) {
        a[0] = a[1];
        a[a.size() - 1] = a[a.size() - 2];
    }
}

}

namespace MaxwellSolver {

void update_magnetic_field(int nx,
                           double dt,
                           double dz,
                           const std::vector<double>& Ex,
                           bool drive_left_boundary,
                           double left_laser_ex_t,
                           bool drive_right_boundary,
                           double right_laser_ex_t,
                           std::vector<double>& By) {
    for (int i = 1; i < nx; ++i) {
        By[i] -= (dt / dz) * (Ex[i] - Ex[i - 1]);
    }

    By[0] = drive_left_boundary ? left_laser_ex_t : Ex[0];
    By[nx] = drive_right_boundary ? -right_laser_ex_t : By[nx - 1];
}

void update_transverse_electric_field(int nx,
                                      double dt,
                                      double dz,
                                      bool drive_left_boundary,
                                      double left_laser_ex_t,
                                      bool drive_right_boundary,
                                      double right_laser_ex_t,
                                      const std::vector<double>& jx,
                                      const std::vector<double>& By,
                                      std::vector<double>& Ex) {

    for (int i = 0; i < nx; ++i) {
        const double curl_by = (By[i + 1] - By[i]) / dz;
        Ex[i] -= dt * curl_by - dt * jx[i];
    }

    // simple_flow style: optional left drive, copy-outflow at right.
    if (drive_left_boundary) {
        Ex[0] = left_laser_ex_t;
    } else {
        Ex[0] = Ex[1];
    }
    if (drive_right_boundary) {
        Ex[nx - 1] = right_laser_ex_t;
    } else {
        Ex[nx - 1] = Ex[nx - 2];
    }
}

void update_longitudinal_electric_field(int nx,
                                        double dt,
                                        const std::vector<double>& jz,
                                        std::vector<double>& Ez) {
    for (int i = 0; i < nx; ++i) {
        Ez[i] += dt * jz[i];
    }
    apply_zero_gradient(Ez);
}

}
