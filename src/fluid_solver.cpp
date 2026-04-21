#include "fluid_solver.h"
#include <algorithm>
#include <cmath>

namespace {

constexpr double kInvSqrtTwoPi = 0.39894228040143267794;
constexpr double kSqrtHalfPi = 1.25331413731550025121;

void apply_zero_gradient(std::vector<double>& a) {
    if (a.size() >= 2U) {
        a[0] = a[1];
        a[a.size() - 1] = a[a.size() - 2];
    }
}

void central_grad(const std::vector<double>& phi, double dz, std::vector<double>& g) {
    const int n = static_cast<int>(phi.size());
    g.resize(n);
    if (n < 2) {
        return;
    }

    for (int i = 1; i < n - 1; ++i) {
        g[i] = (phi[i + 1] - phi[i - 1]) / (2.0 * dz);
    }
    g[0] = (phi[1] - phi[0]) / dz;
    g[n - 1] = (phi[n - 1] - phi[n - 2]) / dz;
}

void upwind_grad(const std::vector<double>& phi,
                 const std::vector<double>& v,
                 double dz,
                 std::vector<double>& g) {
    const int n = static_cast<int>(phi.size());
    g.resize(n);
    if (n < 2) {
        return;
    }

    for (int i = 1; i < n - 1; ++i) {
        const double back = (phi[i] - phi[i - 1]) / dz;
        const double fwd = (phi[i + 1] - phi[i]) / dz;
        g[i] = (v[i] >= 0.0) ? back : fwd;
    }
    g[0] = (phi[1] - phi[0]) / dz;
    g[n - 1] = (phi[n - 1] - phi[n - 2]) / dz;
}

}

namespace FluidSolver {

void update_derived_from_n_vx_qz(int nx,
                                      const std::vector<double>& n,
                                      const std::vector<double>& vx,
                                      const std::vector<double>& qz,
                                      std::vector<double>& gamma,
                                      std::vector<double>& inv_n,
                                      std::vector<double>& vz,
                                      std::vector<double>& jx,
                                      std::vector<double>& jz) {
    if (static_cast<int>(gamma.size()) != nx) gamma.assign(nx, 1.0);
    if (static_cast<int>(inv_n.size()) != nx) inv_n.assign(nx, 1.0);

    vz.resize(nx);
    jx.resize(nx);
    jz.resize(nx);

    for (int i = 0; i < nx; ++i) {
        gamma[i] = 1.0;
        if (n[i] <= 0.0) {
            vz[i] = 0.0;
            inv_n[i] = 0.0;
            jx[i] = 0.0;
            jz[i] = 0.0;
        } else {
            vz[i] = qz[i];
            inv_n[i] = 1.0 / n[i];
            jx[i] = n[i] * vx[i];
            jz[i] = n[i] * vz[i];
        }
    }

    apply_zero_gradient(gamma);
    apply_zero_gradient(inv_n);
    apply_zero_gradient(vz);
    apply_zero_gradient(jx);
    apply_zero_gradient(jz);
}

void update_n(int nx,
              double dt,
              double dz,
              double n_bath,
              double vth_tilde,
              const std::vector<double>& n_old,
              const std::vector<double>& vz,
              const std::vector<double>& jz,
              std::vector<double>& n_new) {
    std::vector<double> flux(nx + 1, 0.0);

    for (int i = 0; i < nx - 1; ++i) {
        const double UL = n_old[i];
        const double UR = n_old[i + 1];
        const double FL = jz[i];
        const double FR = jz[i + 1];
        const double a = std::max(std::fabs(vz[i]), std::fabs(vz[i + 1]));
        flux[i + 1] = 0.5 * (FL + FR) - 0.5 * a * (UR - UL);
    }
    const double thermal_in = n_bath * vth_tilde * kInvSqrtTwoPi;
    const double outgoing_left = std::min(jz[0], 0.0);
    flux[0] = thermal_in + outgoing_left;
    const double outgoing_right = std::max(jz[nx - 1], 0.0);
    flux[nx] = outgoing_right - thermal_in;

    n_new.resize(nx);
    for (int i = 0; i < nx; ++i) {
        n_new[i] = n_old[i] - (dt / dz) * (flux[i + 1] - flux[i]);
    }
    apply_zero_gradient(n_new);
}

void enforce_continuity_constraint_from_jz(int nx,
                                           double dt,
                                           double dz,
                                           double n_bath,
                                           double vth_tilde,
                                           const std::vector<double>& n_old,
                                           const std::vector<double>& vz,
                                           const std::vector<double>& jz,
                                           std::vector<double>& n_new) {
    std::vector<double> flux(nx + 1, 0.0);
    for (int i = 0; i < nx - 1; ++i) {
        const double a = std::max(std::fabs(vz[i]), std::fabs(vz[i + 1]));
        const double dn = n_old[i + 1] - n_old[i];
        flux[i + 1] = 0.5 * (jz[i] + jz[i + 1]) - 0.5 * a * dn;
    }

    const double thermal_in = n_bath * vth_tilde * kInvSqrtTwoPi;
    const double outgoing_left = std::min(jz[0], 0.0);
    flux[0] = thermal_in + outgoing_left;
    const double outgoing_right = std::max(jz[nx - 1], 0.0);
    flux[nx] = outgoing_right - thermal_in;

    std::vector<double> n_proj(nx, 0.0);
    for (int i = 0; i < nx; ++i) {
        n_proj[i] = n_old[i] - (dt / dz) * (flux[i + 1] - flux[i]);
    }
    apply_zero_gradient(n_proj);

    // Convex relaxation keeps the constraint correction from over-shooting density.
    double theta = 1.0;
    const double n_eps = 1.0e-10;
    for (int i = 0; i < nx; ++i) {
        const double dn = n_proj[i] - n_new[i];
        if (dn < 0.0) {
            const double denom = n_new[i] - n_proj[i];
            if (denom > 0.0) {
                const double ti = (n_new[i] - n_eps) / denom;
                theta = std::min(theta, std::max(0.0, std::min(1.0, ti)));
            }
        }
    }

    for (int i = 0; i < nx; ++i) {
        n_new[i] += theta * (n_proj[i] - n_new[i]);
    }
    apply_zero_gradient(n_new);
}

void update_vx(int nx,
               double dt,
               double dz,
               const std::vector<double>& Ex,
               const std::vector<double>& By,
               const std::vector<double>& vz,
               std::vector<double>& vx) {
    std::vector<double> by_center(nx, 0.0);
    std::vector<double> dvx_dz(nx, 0.0);

    for (int i = 0; i < nx; ++i) {
        by_center[i] = 0.5 * (By[i] + By[i + 1]);
    }

    upwind_grad(vx, vz, dz, dvx_dz);

    for (int i = 1; i < nx - 1; ++i) {
        const double src = -(Ex[i] - vz[i] * by_center[i]);
        vx[i] += dt * (src - vz[i] * dvx_dz[i]);
    }

    apply_zero_gradient(vx);
}

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
               std::vector<double>& qz) {
    std::vector<double> by_center(nx, 0.0);
    std::vector<double> dqz_dz(nx, 0.0);
    std::vector<double> dn_dz(nx, 0.0);

    for (int i = 0; i < nx; ++i) {
        by_center[i] = 0.5 * (By[i] + By[i + 1]);
    }

    upwind_grad(qz, vz, dz, dqz_dz);
    central_grad(n, dz, dn_dz);

    for (int i = 1; i < nx - 1; ++i) {
        const double src = -(Ez[i] + vx[i] * by_center[i]) - beta * inv_n[i] * dn_dz[i];
        qz[i] += dt * (src - vz[i] * dqz_dz[i]);
    }

    apply_zero_gradient(qz);
}

void apply_thermal_both_sides(int nx,
                              double vth_tilde,
                              std::vector<double>& n,
                              std::vector<double>& vx,
                              std::vector<double>& qz) {
    if (nx < 2) {
        return;
    }

    // Half-space Maxwellian bath at left boundary: <vx>=0, <vz>_incoming = sqrt(pi/2) * vth.
    n[0] = n[1];
    vx[0] = 0.0;
    qz[0] = kSqrtHalfPi * vth_tilde;

    // Symmetric thermal bath at right boundary: incoming mean velocity points to -z.
    n[nx - 1] = n[nx - 2];
    vx[nx - 1] = 0.0;
    qz[nx - 1] = -kSqrtHalfPi * vth_tilde;
}

}
