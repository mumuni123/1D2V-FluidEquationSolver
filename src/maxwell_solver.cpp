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
                           std::vector<double>& By) {
    for (int i = 1; i < nx; ++i) {
        By[i] -= (dt / dz) * (Ex[i] - Ex[i - 1]);
    }

    By[0] = Ex[0];
    By[nx] = By[nx - 1];
}

void update_transverse_electric_field(int nx,
                                      double dt,
                                      double dz,
                                      double laser_ex_t,
                                      bool drive_left_boundary,
                                      const std::vector<double>& jx,
                                      const std::vector<double>& By,
                                      std::vector<double>& Ex) {

    for (int i = 0; i < nx; ++i) {
        const double curl_by = (By[i + 1] - By[i]) / dz;
        Ex[i] -= dt * curl_by - dt * jx[i];
    }

    // simple_flow style: optional left drive, copy-outflow at right.
    if (drive_left_boundary) {
        Ex[0] = laser_ex_t;
    } else {
        Ex[0] = Ex[1];
    }
    Ex[nx - 1] = Ex[nx - 2];
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
