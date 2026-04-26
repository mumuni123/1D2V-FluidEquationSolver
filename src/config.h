#ifndef CONFIG_H
#define CONFIG_H

#include <string>

struct Config {
    // Aligned with laser.txt
    double lambda0;
    double intensity_w_cm2;
    double electron_temperature_ev;
    double electron_density0;

    double length;
    double dz;
    double t_end;
    double dt_multiplier;

    double gamma_e;

    double snapshot_dt;
    std::string output_dir;

    // Plasma occupies [plasma_left, plasma_right] in meters.
    // If invalid, solver auto-centers to [0.25*length, 0.75*length].
    double plasma_left;
    double plasma_right;

    // Laser enters from the left boundary and propagates rightward.
    bool laser_from_left_boundary;
    bool laser_smooth_ramp;
    double laser_ramp_cycles;

    // Interface jump parameters (dimensionless, solver units):
    // e_n x (H2-H1)=alpha, e_n·(D2-D1)=sigma at vacuum-plasma interface.
    bool enable_interface_jump_bc;
    bool interface_source_from_plasma;
    double interface_alpha_tilde;
    double interface_sigma_tilde;

    // Linear overdense-plasma test mode:
    // freeze n, Ez, vz and qz; only advance Ex, By and vx in the plasma slab.
    bool linear_overdense_test;

    double n_floor_ratio;

    Config()
        : lambda0(3.0e-6),
          intensity_w_cm2(1.0e17),
          electron_temperature_ev(100.0),
          electron_density0(3.0e26),
          length(25.0e-6),
          dz(0.002e-6),
          t_end(120.0e-15),
          dt_multiplier(0.1),
          gamma_e(3.0),
          snapshot_dt(0.5e-15),
          output_dir("output"),
          plasma_left(5.0e-6),
          plasma_right(20.0e-6),
          laser_from_left_boundary(true),
          laser_smooth_ramp(true),
          laser_ramp_cycles(10.0),
          enable_interface_jump_bc(false),
          interface_source_from_plasma(true),
          interface_alpha_tilde(0.0),
          interface_sigma_tilde(0.0),
          linear_overdense_test(true),
          n_floor_ratio(0.0) {}
};

#endif
