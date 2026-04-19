#include "fluid_maxwell_1d.h"
#include "fluid_solver.h"
#include "maxwell_solver.h"
#include "physics_constants.h"
#include "result_output.h"

#include <cmath>
#include <iostream>
#include <random>

FluidMaxwell1D::FluidMaxwell1D(const Config& cfg)
    : cfg_(cfg),
      nx_(static_cast<int>(std::floor(cfg.length / cfg.dz + 0.5))),
      z_(nx_, 0.0),
      Ex_(nx_, 0.0),
      Ez_(nx_, 0.0),
            n_(nx_, 1.0),
        vx_(nx_, 0.0),
        qz_(nx_, 0.0),
        gamma_(nx_, 1.0),
            inv_n_(nx_, 1.0),
            vz_(nx_, 0.0),
            jx_(nx_, 0.0),
            jz_(nx_, 0.0),
      By_(nx_ + 1, 0.0),
      dt_(0.0),
      dt_tilde_(0.0),
      dz_tilde_(0.0),
                strict_cfl_(0.05),
      n_steps_(0),
      snapshot_every_(1),
      omega_pe_(0.0),
      omega0_(0.0),
      omega0_tilde_(0.0),
      E_scale_(0.0),
      B_scale_(0.0),
      ne0_(0.0),
      beta_(0.0),
      vth_tilde_(0.0),
      rng_(std::random_device{}()),
    time_tilde_(0.0) {
    for (int i = 0; i < nx_; ++i) {
        z_[i] = i * cfg_.dz;
    }

    omega_pe_ = std::sqrt(cfg_.electron_density0 * PhysConst::QE * PhysConst::QE /
                          (PhysConst::EPS0 * PhysConst::ME));
    dt_ = strict_cfl_ * cfg_.dt_multiplier * cfg_.dz / PhysConst::C;
    dt_tilde_ = dt_ * omega_pe_;
    dz_tilde_ = cfg_.dz * omega_pe_ / PhysConst::C;
    n_steps_ = static_cast<int>(std::ceil(cfg_.t_end / dt_));
    snapshot_every_ = static_cast<int>(std::floor(cfg_.snapshot_dt / dt_ + 0.5));
    if (snapshot_every_ < 1) {
        snapshot_every_ = 1;
    }

    omega0_ = 2.0 * PhysConst::PI * PhysConst::C / cfg_.lambda0;
    omega0_tilde_ = omega0_ / omega_pe_;
    {
        const double intensity_w_m2 = cfg_.intensity_w_cm2 * 1.0e4;
        E0_ = std::sqrt(2.0 * intensity_w_m2 / (PhysConst::C * PhysConst::EPS0));
    }

    E_scale_ = PhysConst::ME * PhysConst::C * omega_pe_ / PhysConst::QE;
    B_scale_ = PhysConst::ME * omega_pe_ / PhysConst::QE;
    ne0_ = cfg_.electron_density0;

    {
        const double Te_J = cfg_.electron_temperature_ev * PhysConst::EV;
        beta_ = cfg_.gamma_e * Te_J / (PhysConst::ME * PhysConst::C * PhysConst::C);
        if (cfg_.gamma_e > 0.0) {
            vth_tilde_ = std::sqrt(beta_ / cfg_.gamma_e);
        } else {
            vth_tilde_ = 0.0;
        }
    }
}

void FluidMaxwell1D::run() {
    if (!ResultOutput::ensure_output_dir(cfg_.output_dir)) {
        std::cerr << "Failed to create output directory: " << cfg_.output_dir << std::endl;
        return;
    }

    std::ofstream diag;
    if (!ResultOutput::open_diagnostics(cfg_.output_dir, diag)) {
        std::cerr << "Failed to open diagnostics.csv" << std::endl;
        return;
    }

    auto write_snapshot = [&](int step) {
        std::vector<double> Ex_phys(nx_, 0.0);
        std::vector<double> Ez_phys(nx_, 0.0);
        std::vector<double> ne_phys(nx_, 0.0);
        std::vector<double> vx_phys(nx_, 0.0);
        std::vector<double> vz_phys(nx_, 0.0);
        std::vector<double> Pe_phys(nx_, 0.0);
        std::vector<double> By_phys(nx_ + 1, 0.0);

        for (int i = 0; i < nx_; ++i) {
            Ex_phys[i] = Ex_[i] * E_scale_;
            Ez_phys[i] = Ez_[i] * E_scale_;
            ne_phys[i] = n_[i] * ne0_;
            const double vx_tilde = vx_[i];
            const double vz_tilde = vz_[i];
            vx_phys[i] = vx_tilde * PhysConst::C;
            vz_phys[i] = vz_tilde * PhysConst::C;
            Pe_phys[i] = (beta_ * PhysConst::ME * PhysConst::C * PhysConst::C / cfg_.gamma_e) * ne_phys[i];
        }
        for (int i = 0; i < nx_ + 1; ++i) {
            By_phys[i] = By_[i] * B_scale_;
        }

        ResultOutput::write_state_csv(cfg_.output_dir, step, nx_, z_, Ex_phys, Ez_phys, By_phys, ne_phys, vx_phys, vz_phys, Pe_phys);
        ResultOutput::append_diagnostics(diag, time_tilde_ / omega_pe_, nx_, cfg_.dz, Ex_phys, Ez_phys, By_phys, ne_phys, vx_phys, vz_phys);
    };

    write_snapshot(0);

    for (int step = 1; step <= n_steps_; ++step) {
        step_once();
        if ((step % snapshot_every_ == 0) || (step == n_steps_)) {
            write_snapshot(step);
        }
    }
}

void FluidMaxwell1D::print_summary() const {
    std::cout << "1D Maxwell-Fluid simulation" << std::endl;
    std::cout << "nx=" << nx_ << ", dz=" << cfg_.dz << " m, dt=" << dt_ << " s" << std::endl;
    std::cout << "dt_tilde=" << dt_tilde_ << ", dz_tilde=" << dz_tilde_ << std::endl;
    std::cout << "omega_pe=" << omega_pe_ << " rad/s, beta=" << beta_ << std::endl;
    std::cout << "lambda0=" << cfg_.lambda0 << " m, a0=" << (E0_ / E_scale_) << std::endl;
    std::cout << "ne0=" << cfg_.electron_density0 << " m^-3, Te="
              << cfg_.electron_temperature_ev << " eV" << std::endl;
}

double FluidMaxwell1D::laser_ex(double t_tilde) const { return (E0_ / E_scale_) * std::sin(omega0_tilde_ * t_tilde); }

void FluidMaxwell1D::step_once() {
    const double t_next = time_tilde_ + dt_tilde_;
    const std::vector<double> n_old = n_;
    const std::vector<double> Ez_old = Ez_;

    FluidSolver::update_derived_from_n_vx_qz(nx_, n_, vx_, qz_, gamma_, inv_n_, vz_, jx_, jz_);

    MaxwellSolver::update_magnetic_field(nx_, dt_tilde_, dz_tilde_, Ex_, By_);

    MaxwellSolver::update_transverse_electric_field(nx_, dt_tilde_, dz_tilde_,
                                                    laser_ex(t_next), jx_, By_, Ex_);
    MaxwellSolver::update_longitudinal_electric_field(nx_, dt_tilde_, jz_, Ez_);

    std::vector<double> jz_from_ez(nx_, 0.0);
    for (int i = 0; i < nx_; ++i) {
        jz_from_ez[i] = (Ez_[i] - Ez_old[i]) / dt_tilde_;
    }

    std::vector<double> n_new;
    FluidSolver::update_n(nx_, dt_tilde_, dz_tilde_, 1.0, vth_tilde_, n_, vz_, jz_, n_new);
    FluidSolver::enforce_continuity_constraint_from_jz(nx_, dt_tilde_, dz_tilde_,
                                                       1.0, vth_tilde_, n_old, vz_, jz_from_ez, n_new);
    n_.swap(n_new);

    FluidSolver::update_derived_from_n_vx_qz(nx_, n_, vx_, qz_, gamma_, inv_n_, vz_, jx_, jz_);
    FluidSolver::update_vx(nx_, dt_tilde_, dz_tilde_, Ex_, By_, vz_, vx_);

    FluidSolver::update_derived_from_n_vx_qz(nx_, n_, vx_, qz_, gamma_, inv_n_, vz_, jx_, jz_);
    FluidSolver::update_qz(nx_, dt_tilde_, dz_tilde_, beta_, Ez_, By_, n_, inv_n_, vx_, vz_, qz_);

    FluidSolver::apply_thermal_both_sides(nx_, vth_tilde_, n_, vx_, qz_);

    FluidSolver::update_derived_from_n_vx_qz(nx_, n_, vx_, qz_, gamma_, inv_n_, vz_, jx_, jz_);

    time_tilde_ = t_next;
}

