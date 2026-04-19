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

    double n_floor_ratio;

    Config()
        : lambda0(1.0e-6),
          intensity_w_cm2(1.0e17),
          electron_temperature_ev(100.0),
          electron_density0(3.0e26),
          length(8.0e-6),
          dz(0.002e-6),
          t_end(120.0e-15),
          dt_multiplier(0.1),
          gamma_e(3.0),
          snapshot_dt(0.5e-15),
          output_dir("output"),
          n_floor_ratio(0.0) {}
};

#endif
