#include "fluid_solver.h"
#include <algorithm>
#include <cmath>

namespace {

constexpr double kInvSqrtTwoPi = 0.39894228040143267794;
constexpr double kSqrtHalfPi = 1.25331413731550025121;

inline bool in_plasma(int i, const std::vector<char>* plasma_mask) {
    return (plasma_mask == nullptr) || ((*plasma_mask)[i] != 0);
}

void apply_zero_gradient(std::vector<double>& a) {
    if (a.size() >= 2U) {
        a[0] = a[1];
        a[a.size() - 1] = a[a.size() - 2];
    }
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
                                      std::vector<double>& jz,
                                      const std::vector<char>* plasma_mask) {
    update_kinematics_from_n_qz(nx, n, qz, gamma, inv_n, vz, plasma_mask);
    update_current_from_n_v(nx, n, vx, vz, jx, jz, plasma_mask);
}

void update_kinematics_from_n_qz(int nx,
                                 const std::vector<double>& n,
                                 const std::vector<double>& qz,
                                 std::vector<double>& gamma,
                                 std::vector<double>& inv_n,
                                 std::vector<double>& vz,
                                 const std::vector<char>* plasma_mask) {
    if (static_cast<int>(gamma.size()) != nx) gamma.assign(nx, 1.0);
    if (static_cast<int>(inv_n.size()) != nx) inv_n.assign(nx, 1.0);

    vz.resize(nx);

    #pragma omp parallel for
    for (int i = 0; i < nx; ++i) {
        if (!in_plasma(i, plasma_mask)) {
            gamma[i] = 1.0;
            vz[i] = 0.0;
            inv_n[i] = 0.0;
            continue;
        }

        gamma[i] = 1.0;
        if (n[i] <= 0.0) {
            vz[i] = 0.0;
            inv_n[i] = 0.0;
        } else {
            vz[i] = qz[i];
            inv_n[i] = 1.0 / n[i];
        }
    }

    if (plasma_mask == nullptr)
    {
        apply_zero_gradient(gamma);
        apply_zero_gradient(inv_n);
        apply_zero_gradient(vz);
    }
}

void update_current_from_n_v(int nx,
                             const std::vector<double>& n,
                             const std::vector<double>& vx,
                             const std::vector<double>& vz,
                             std::vector<double>& jx,
                             std::vector<double>& jz,
                             const std::vector<char>* plasma_mask) {
    jx.resize(nx);
    jz.resize(nx);

    #pragma omp parallel for
    for (int i = 0; i < nx; ++i) {
        if (!in_plasma(i, plasma_mask) || n[i] <= 0.0) {
            jx[i] = 0.0;
            jz[i] = 0.0;
        } else {
            jx[i] = n[i] * vx[i];
            jz[i] = n[i] * vz[i];
        }
    }

    if (plasma_mask == nullptr)
    {
        apply_zero_gradient(jx);
        apply_zero_gradient(jz);
    }
}

void update_n(int nx,
              double dt,
              double dz,
              double n_bath,
              double vth_tilde,
              const std::vector<double>& n_old,
              const std::vector<double>& vz,
              const std::vector<double>& jz,
              std::vector<double>& n_new,
              const std::vector<char>* plasma_mask,
              bool closed_boundary,
              std::vector<double>* face_flux_workspace) {
    const double thermal_in = n_bath * vth_tilde * kInvSqrtTwoPi;

    std::vector<double> local_face_flux;
    std::vector<double>& face_flux = face_flux_workspace ? *face_flux_workspace : local_face_flux;
    face_flux.resize(nx + 1);

    face_flux[0] = closed_boundary ? 0.0 : thermal_in + std::min(jz[0], 0.0);
    face_flux[nx] = closed_boundary ? 0.0 : std::max(jz[nx - 1], 0.0) - thermal_in;

    #pragma omp parallel for
    for (int face = 1; face < nx; ++face) {
        const int left = face - 1;
        const int right = face;
        if (!in_plasma(left, plasma_mask) || !in_plasma(right, plasma_mask)) {
            face_flux[face] = 0.0;
            continue;
        }

        const double a = std::max(std::fabs(vz[left]), std::fabs(vz[right]));
        face_flux[face] = 0.5 * (jz[left] + jz[right]) - 0.5 * a * (n_old[right] - n_old[left]);
    }

    n_new.resize(nx);
    #pragma omp parallel for
    for (int i = 0; i < nx; ++i) {
        if (!in_plasma(i, plasma_mask)) {
            n_new[i] = 0.0;
        } else {
            n_new[i] = n_old[i] - (dt / dz) * (face_flux[i + 1] - face_flux[i]);
        }
    }

    if (plasma_mask == nullptr) {
        apply_zero_gradient(n_new);
    }
}

void enforce_continuity_constraint_from_jz(int nx,
                                           double dt,
                                           double dz,
                                           double n_bath,
                                           double vth_tilde,
                                           const std::vector<double>& n_old,
                                           const std::vector<double>& vz,
                                           const std::vector<double>& jz,
                                           std::vector<double>& n_new,
                                           const std::vector<char>* plasma_mask,
                                           bool closed_boundary,
                                           std::vector<double>* n_proj_workspace) {
    const double thermal_in = n_bath * vth_tilde * kInvSqrtTwoPi;
    auto face_flux = [&](int face) -> double {
        if (face == 0) {
            return closed_boundary ? 0.0 : thermal_in + std::min(jz[0], 0.0);
        }
        if (face == nx) {
            return closed_boundary ? 0.0 : std::max(jz[nx - 1], 0.0) - thermal_in;
        }
        const int left = face - 1;
        const int right = face;
        if (!in_plasma(left, plasma_mask) || !in_plasma(right, plasma_mask)) {
            return 0.0;
        }

        const double a = std::max(std::fabs(vz[left]), std::fabs(vz[right]));
        return 0.5 * (jz[left] + jz[right]) - 0.5 * a * (n_old[right] - n_old[left]);
    };

    std::vector<double> local_n_proj;
    std::vector<double>& n_proj = n_proj_workspace ? *n_proj_workspace : local_n_proj;
    n_proj.resize(nx);
    double flux_left = face_flux(0);
    for (int i = 0; i < nx; ++i) {
        const double flux_right = face_flux(i + 1);
        if (!in_plasma(i, plasma_mask)) {
            n_proj[i] = 0.0;
        } else {
            n_proj[i] = n_old[i] - (dt / dz) * (flux_right - flux_left);
        }
        flux_left = flux_right;
    }
    if (plasma_mask == nullptr) {
        apply_zero_gradient(n_proj);
    }

    // Convex relaxation keeps the constraint correction from over-shooting density.
    double theta = 1.0;
    const double n_eps = 1.0e-10;
    #pragma omp parallel for reduction(min:theta)
    for (int i = 0; i < nx; ++i) {
        if (!in_plasma(i, plasma_mask)) {
            n_new[i] = 0.0;
            continue;
        }
        const double dn = n_proj[i] - n_new[i];
        if (dn < 0.0) {
            const double denom = n_new[i] - n_proj[i];
            if (denom > 0.0) {
                const double ti = (n_new[i] - n_eps) / denom;
                theta = std::min(theta, std::max(0.0, std::min(1.0, ti)));
            }
        }
    }

    #pragma omp parallel for
    for (int i = 0; i < nx; ++i) {
        if (in_plasma(i, plasma_mask)) {
            n_new[i] += theta * (n_proj[i] - n_new[i]);
        } else {
            n_new[i] = 0.0;
        }
    }

    if (plasma_mask == nullptr) {
        apply_zero_gradient(n_new);
    }
}

void update_vx(int nx,
               double dt,
               double dz,
               const std::vector<double>& Ex,
               const std::vector<double>& By,
               const std::vector<double>& vz,
               std::vector<double>& vx,
               const std::vector<char>* plasma_mask,
               std::vector<double>* dvx_dz_workspace) {
    std::vector<double> local_dvx_dz;
    std::vector<double>& dvx_dz = dvx_dz_workspace ? *dvx_dz_workspace : local_dvx_dz;
    dvx_dz.resize(nx);

    if (nx >= 2) {
        dvx_dz[0] = (vx[1] - vx[0]) / dz;
        dvx_dz[nx - 1] = (vx[nx - 1] - vx[nx - 2]) / dz;
    }
    #pragma omp parallel for
    for (int i = 1; i < nx - 1; ++i) {
        const double back = (vx[i] - vx[i - 1]) / dz;
        const double fwd = (vx[i + 1] - vx[i]) / dz;
        dvx_dz[i] = (vz[i] >= 0.0) ? back : fwd;
    }

    #pragma omp parallel for
    for (int i = 1; i < nx - 1; ++i) {
        if (!in_plasma(i, plasma_mask)) {
            vx[i] = 0.0;
            continue;
        }
        const double by_center = 0.5 * (By[i] + By[i + 1]);
        const double src = -(Ex[i] - vz[i] * by_center);
        vx[i] += dt * (src - vz[i] * dvx_dz[i]);
    }

    if (plasma_mask == nullptr) {
        apply_zero_gradient(vx);
    } else {
        #pragma omp parallel for
        for (int i = 0; i < nx; ++i) {
            if (!in_plasma(i, plasma_mask)) {
                vx[i] = 0.0;
            }
        }
    }
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
               std::vector<double>& qz,
               const std::vector<char>* plasma_mask,
               std::vector<double>* dqz_dz_workspace,
               std::vector<double>* dn_dz_workspace) {
    std::vector<double> local_dqz_dz;
    std::vector<double> local_dn_dz;
    std::vector<double>& dqz_dz = dqz_dz_workspace ? *dqz_dz_workspace : local_dqz_dz;
    std::vector<double>& dn_dz = dn_dz_workspace ? *dn_dz_workspace : local_dn_dz;
    dqz_dz.resize(nx);
    dn_dz.resize(nx);

    if (nx >= 2) {
        dqz_dz[0] = (qz[1] - qz[0]) / dz;
        dqz_dz[nx - 1] = (qz[nx - 1] - qz[nx - 2]) / dz;
        dn_dz[0] = (n[1] - n[0]) / dz;
        dn_dz[nx - 1] = (n[nx - 1] - n[nx - 2]) / dz;
    }
    #pragma omp parallel for
    for (int i = 1; i < nx - 1; ++i) {
        const double back = (qz[i] - qz[i - 1]) / dz;
        const double fwd = (qz[i + 1] - qz[i]) / dz;
        dqz_dz[i] = (vz[i] >= 0.0) ? back : fwd;
        dn_dz[i] = (n[i + 1] - n[i - 1]) / (2.0 * dz);
    }

    #pragma omp parallel for
    for (int i = 1; i < nx - 1; ++i) {
        if (!in_plasma(i, plasma_mask)) {
            qz[i] = 0.0;
            continue;
        }
        const double by_center = 0.5 * (By[i] + By[i + 1]);
        const double src = -(Ez[i] + vx[i] * by_center) - beta * inv_n[i] * dn_dz[i];
        qz[i] += dt * (src - vz[i] * dqz_dz[i]);
    }

    if (plasma_mask == nullptr) {
        apply_zero_gradient(qz);
    } else {
        #pragma omp parallel for
        for (int i = 0; i < nx; ++i) {
            if (!in_plasma(i, plasma_mask)) {
                qz[i] = 0.0;
            }
        }
    }
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
