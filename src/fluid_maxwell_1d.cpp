#include "fluid_maxwell_1d.h"
#include "fluid_solver.h"
#include "maxwell_solver.h"
#include "physics_constants.h"
#include "result_output.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>

FluidMaxwell1D::FluidMaxwell1D(const Config &cfg)
    : cfg_(cfg),
      nx_(static_cast<int>(std::floor(cfg.length / cfg.dz + 0.5))),
      z_(nx_, 0.0),
      Ex_(nx_, 0.0),
      Ez_(nx_, 0.0),
      n_(nx_, 0.0),
      vx_(nx_, 0.0),
      qz_(nx_, 0.0),
      gamma_(nx_, 1.0),
      inv_n_(nx_, 1.0),
      vz_(nx_, 0.0),
      jx_(nx_, 0.0),
      jz_(nx_, 0.0),
      By_(nx_ + 1, 0.0),
      plasma_mask_(nx_, 0),
      interface_faces_(),
      n_old_(nx_, 0.0),
      Ez_old_(nx_, 0.0),
      jz_from_ez_(nx_, 0.0),
      n_new_(nx_, 0.0),
      n_proj_(nx_, 0.0),
      dvx_dz_(nx_, 0.0),
      dqz_dz_(nx_, 0.0),
      dn_dz_(nx_, 0.0),
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
      n_bath_tilde_(0.0),
      plasma_left_(0.0),
      plasma_right_(0.0),
      plasma_left_index_(-1),
      plasma_right_index_(-1),
      rng_(std::random_device{}()),
      time_tilde_(0.0)
{
    for (int i = 0; i < nx_; ++i)
    {
        z_[i] = i * cfg_.dz;
    }

    omega_pe_ = std::sqrt(cfg_.electron_density0 * PhysConst::QE * PhysConst::QE /
                          (PhysConst::EPS0 * PhysConst::ME));
    dt_ = strict_cfl_ * cfg_.dt_multiplier * cfg_.dz / PhysConst::C;
    dt_tilde_ = dt_ * omega_pe_;
    dz_tilde_ = cfg_.dz * omega_pe_ / PhysConst::C;
    n_steps_ = static_cast<int>(std::ceil(cfg_.t_end / dt_));
    snapshot_every_ = static_cast<int>(std::floor(cfg_.snapshot_dt / dt_ + 0.5));
    if (snapshot_every_ < 1)
    {
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
        if (cfg_.gamma_e > 0.0)
        {
            vth_tilde_ = std::sqrt(beta_ / cfg_.gamma_e);
        }
        else
        {
            vth_tilde_ = 0.0;
        }
    }

    plasma_left_ = cfg_.plasma_left;
    plasma_right_ = cfg_.plasma_right;
    if (!(plasma_left_ >= 0.0 && plasma_right_ > plasma_left_ && plasma_right_ <= cfg_.length))
    {
        plasma_left_ = 0.25 * cfg_.length;
        plasma_right_ = 0.75 * cfg_.length;
    }

    n_bath_tilde_ = std::max(cfg_.n_floor_ratio, 0.0);

    for (int i = 0; i < nx_; ++i)
    {
        const bool in_plasma = (z_[i] >= plasma_left_ && z_[i] <= plasma_right_);
        plasma_mask_[i] = in_plasma ? 1 : 0;
        n_[i] = in_plasma ? 1.0 : 0.0;
        vx_[i] = 0.0;
        qz_[i] = 0.0;
        if (in_plasma)
        {
            if (plasma_left_index_ < 0)
            {
                plasma_left_index_ = i;
            }
            plasma_right_index_ = i;
        }
    }

    for (int i = 0; i < nx_ - 1; ++i)
    {
        if (plasma_mask_[i] != plasma_mask_[i + 1])
        {
            interface_faces_.push_back(i);
        }
    }

    // Keep initial electromagnetic fields zero; laser enters from left boundary during stepping.
}

void FluidMaxwell1D::run()
{
    if (!ResultOutput::ensure_output_dir(cfg_.output_dir))
    {
        std::cerr << "Failed to create output directory: " << cfg_.output_dir << std::endl;
        return;
    }

    std::ofstream diag;
    if (!ResultOutput::open_diagnostics(cfg_.output_dir, diag))
    {
        std::cerr << "Failed to open diagnostics.csv" << std::endl;
        return;
    }

    auto write_snapshot = [&](int step)
    {
        std::vector<double> Ex_phys(nx_, 0.0);
        std::vector<double> Ez_phys(nx_, 0.0);
        std::vector<double> ne_phys(nx_, 0.0);
        std::vector<double> vx_phys(nx_, 0.0);
        std::vector<double> vz_phys(nx_, 0.0);
        std::vector<double> Pe_phys(nx_, 0.0);
        std::vector<double> By_phys(nx_ + 1, 0.0);

        for (int i = 0; i < nx_; ++i)
        {
            Ex_phys[i] = Ex_[i] * E_scale_;
            Ez_phys[i] = Ez_[i] * E_scale_;
            ne_phys[i] = n_[i] * ne0_;
            const double vx_tilde = vx_[i];
            const double vz_tilde = vz_[i];
            vx_phys[i] = vx_tilde * PhysConst::C;
            vz_phys[i] = vz_tilde * PhysConst::C;
            Pe_phys[i] = (beta_ * PhysConst::ME * PhysConst::C * PhysConst::C / cfg_.gamma_e) * ne_phys[i];
        }
        for (int i = 0; i < nx_ + 1; ++i)
        {
            By_phys[i] = By_[i] * B_scale_;
        }

        ResultOutput::write_state_csv(cfg_.output_dir, step, nx_, z_, Ex_phys, Ez_phys, By_phys, ne_phys, vx_phys, vz_phys, Pe_phys);
        ResultOutput::append_diagnostics(diag, time_tilde_ / omega_pe_, nx_, cfg_.dz, Ex_phys, Ez_phys, By_phys, ne_phys, vx_phys, vz_phys);
    };

    write_snapshot(0);

    for (int step = 1; step <= n_steps_; ++step)
    {
        step_once();
        if ((step % snapshot_every_ == 0) || (step == n_steps_))
        {
            write_snapshot(step);
        }
    }
}

void FluidMaxwell1D::print_summary() const
{
    std::cout << "1D Maxwell-Fluid simulation" << std::endl;
    std::cout << "nx=" << nx_ << ", dz=" << cfg_.dz << " m, dt=" << dt_ << " s" << std::endl;
    std::cout << "dt_tilde=" << dt_tilde_ << ", dz_tilde=" << dz_tilde_ << std::endl;
    std::cout << "omega_pe=" << omega_pe_ << " rad/s, beta=" << beta_ << std::endl;
    std::cout << "lambda0=" << cfg_.lambda0 << " m, a0=" << (E0_ / E_scale_) << std::endl;
    std::cout << "ne0=" << cfg_.electron_density0 << " m^-3, Te="
              << cfg_.electron_temperature_ev << " eV" << std::endl;
    std::cout << "domain_length=" << cfg_.length << " m" << std::endl;
    std::cout << "plasma_region=[" << plasma_left_ << ", " << plasma_right_ << "] m" << std::endl;
    std::cout << "laser_from_left_boundary=" << (cfg_.laser_from_left_boundary ? "true" : "false") << std::endl;
    std::cout << "laser_smooth_ramp=" << (cfg_.laser_smooth_ramp ? "true" : "false")
              << ", laser_ramp_cycles=" << cfg_.laser_ramp_cycles << std::endl;
    std::cout << "linear_overdense_test=" << (cfg_.linear_overdense_test ? "true" : "false") << std::endl;
}

double FluidMaxwell1D::laser_ex(double t_tilde) const
{
    double envelope = 1.0;
    if (cfg_.laser_smooth_ramp && cfg_.laser_ramp_cycles > 0.0 && omega0_tilde_ > 0.0)
    {
        const double ramp_time = cfg_.laser_ramp_cycles * 2.0 * PhysConst::PI / omega0_tilde_;
        if (t_tilde < ramp_time)
        {
            const double s = std::sin(0.5 * PhysConst::PI * t_tilde / ramp_time);
            envelope = s * s;
        }
    }

    return envelope * (E0_ / E_scale_) * std::sin(omega0_tilde_ * t_tilde);
}

void FluidMaxwell1D::step_once()
{
    const double t_next = time_tilde_ + dt_tilde_;
    n_old_ = n_;
    Ez_old_ = Ez_;

    auto is_plasma_cell = [&](int i)
    {
        return plasma_mask_[i] != 0;
    };

    auto apply_interface_jump = [&]()
    {
        if (!cfg_.enable_interface_jump_bc)
        {
            return;
        }

        for (std::vector<int>::const_iterator it = interface_faces_.begin(); it != interface_faces_.end(); ++it)
        {
            const int i = *it;
            const bool left_plasma = is_plasma_cell(i);

            const int ip = left_plasma ? i : (i + 1);
            const int iv = left_plasma ? (i + 1) : i;
            double alpha_tilde = cfg_.interface_alpha_tilde;
            double sigma_tilde = cfg_.interface_sigma_tilde;

            if (cfg_.interface_source_from_plasma)
            {
                // Local effective sheet sources from the plasma side state.
                alpha_tilde = -dt_tilde_ * jx_[ip];
                sigma_tilde = Ez_[ip];
            }

            // e_n x (E2-E1)=0 -> tangential Ex continuity.
            const double ex_avg = 0.5 * (Ex_[i] + Ex_[i + 1]);
            Ex_[i] = ex_avg;
            Ex_[i + 1] = ex_avg;

            // e_n x (H2-H1)=alpha (discrete face form in solver units).
            const int k = i + 1;
            if (k >= 1 && k <= nx_)
            {
                By_[k] = By_[k - 1] + alpha_tilde;
            }

            // e_n·(D2-D1)=sigma, with vacuum Ez pinned to 0 in this model.
            Ez_[ip] = sigma_tilde;
            Ez_[iv] = 0.0;
        }
    };

    auto apply_closed_plasma_boundary = [&]()
    {
        if (plasma_left_index_ >= 0)
        {
            qz_[plasma_left_index_] = 0.0;
            vz_[plasma_left_index_] = 0.0;
            qz_[plasma_right_index_] = 0.0;
            vz_[plasma_right_index_] = 0.0;
        }
    };

    if (cfg_.linear_overdense_test)
    {
        FluidSolver::update_derived_from_n_vx_qz(nx_, n_, vx_, qz_, gamma_, inv_n_, vz_, jx_, jz_, &plasma_mask_);

        MaxwellSolver::update_magnetic_field(nx_, dt_tilde_, dz_tilde_, Ex_, By_);
        MaxwellSolver::update_transverse_electric_field(nx_, dt_tilde_, dz_tilde_,
                                                        laser_ex(t_next), cfg_.laser_from_left_boundary, jx_, By_, Ex_);

        FluidSolver::update_vx(nx_, dt_tilde_, dz_tilde_, Ex_, By_, vz_, vx_, &plasma_mask_, &dvx_dz_);

        FluidSolver::update_derived_from_n_vx_qz(nx_, n_, vx_, qz_, gamma_, inv_n_, vz_, jx_, jz_, &plasma_mask_);
        time_tilde_ = t_next;
        return;
    }

    FluidSolver::update_derived_from_n_vx_qz(nx_, n_, vx_, qz_, gamma_, inv_n_, vz_, jx_, jz_, &plasma_mask_);

    MaxwellSolver::update_magnetic_field(nx_, dt_tilde_, dz_tilde_, Ex_, By_);

    MaxwellSolver::update_transverse_electric_field(nx_, dt_tilde_, dz_tilde_,
                                                    laser_ex(t_next), cfg_.laser_from_left_boundary, jx_, By_, Ex_);
    MaxwellSolver::update_longitudinal_electric_field(nx_, dt_tilde_, jz_, Ez_);
    apply_interface_jump();

    for (int i = 0; i < nx_; ++i)
    {
        jz_from_ez_[i] = is_plasma_cell(i) ? ((Ez_[i] - Ez_old_[i]) / dt_tilde_) : 0.0;
    }

    FluidSolver::update_n(nx_, dt_tilde_, dz_tilde_, n_bath_tilde_, vth_tilde_, n_, vz_, jz_, n_new_,
                          &plasma_mask_, true);
    FluidSolver::enforce_continuity_constraint_from_jz(nx_, dt_tilde_, dz_tilde_,
                                                       n_bath_tilde_, vth_tilde_, n_old_, vz_, jz_from_ez_, n_new_,
                                                       &plasma_mask_, true, &n_proj_);
    n_.swap(n_new_);

    FluidSolver::update_derived_from_n_vx_qz(nx_, n_, vx_, qz_, gamma_, inv_n_, vz_, jx_, jz_, &plasma_mask_);
    FluidSolver::update_vx(nx_, dt_tilde_, dz_tilde_, Ex_, By_, vz_, vx_, &plasma_mask_, &dvx_dz_);

    FluidSolver::update_derived_from_n_vx_qz(nx_, n_, vx_, qz_, gamma_, inv_n_, vz_, jx_, jz_, &plasma_mask_);
    FluidSolver::update_qz(nx_, dt_tilde_, dz_tilde_, beta_, Ez_, By_, n_, inv_n_, vx_, vz_, qz_,
                           &plasma_mask_, &dqz_dz_, &dn_dz_);

    apply_closed_plasma_boundary();

    apply_interface_jump();

    FluidSolver::update_derived_from_n_vx_qz(nx_, n_, vx_, qz_, gamma_, inv_n_, vz_, jx_, jz_, &plasma_mask_);

    time_tilde_ = t_next;
}
