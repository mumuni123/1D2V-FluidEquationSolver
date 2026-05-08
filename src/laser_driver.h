#ifndef LASER_DRIVER_H
#define LASER_DRIVER_H

#include "config.h"

namespace LaserDriver {

struct BoundaryDrive {
    bool drive_left;
    bool drive_right;
    double left_ex;
    double right_ex;

    BoundaryDrive()
        : drive_left(false),
          drive_right(false),
          left_ex(0.0),
          right_ex(0.0) {}
};

double beam_ex_tilde(const LaserBeamConfig& beam,
                     double t_tilde,
                     double omega_pe,
                     double e_scale);

BoundaryDrive boundary_drive(const Config& cfg,
                             double t_tilde,
                             double omega_pe,
                             double e_scale);

}

#endif
