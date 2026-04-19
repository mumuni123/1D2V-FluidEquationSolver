#ifndef FLUID_MAXWELL_1D_H
#define FLUID_MAXWELL_1D_H

#include <fstream>
#include <random>
#include <vector>

#include "config.h"

class FluidMaxwell1D {
public:
    explicit FluidMaxwell1D(const Config& cfg);

    void run();
    void print_summary() const;

private:
    double laser_ex(double t) const;

    void step_once();

    Config cfg_;

    int nx_;
    std::vector<double> z_;

    std::vector<double> Ex_;
    std::vector<double> Ez_;
    std::vector<double> n_;
    std::vector<double> vx_;
    std::vector<double> qz_;
    std::vector<double> gamma_;
    std::vector<double> inv_n_;
    std::vector<double> vz_;
    std::vector<double> jx_;
    std::vector<double> jz_;
    std::vector<double> By_;

    double dt_;
    double dt_tilde_;
    double dz_tilde_;
    double strict_cfl_;
    int n_steps_;
    int snapshot_every_;

    double omega_pe_;
    double omega0_;
    double omega0_tilde_;
    double E0_;
    double E_scale_;
    double B_scale_;
    double ne0_;
    double beta_;
    double vth_tilde_;

    std::mt19937 rng_;

    double time_tilde_;
};

#endif
