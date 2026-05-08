#include "laser_driver.h"

#include "physics_constants.h"

#include <cmath>

namespace LaserDriver {

double beam_ex_tilde(const LaserBeamConfig& beam,
                     double t_tilde,
                     double omega_pe,
                     double e_scale) {
    if (!beam.enabled || beam.intensity_w_cm2 <= 0.0 || beam.lambda <= 0.0 ||
        omega_pe <= 0.0 || e_scale <= 0.0) {
        return 0.0;
    }

    const double omega = 2.0 * PhysConst::PI * PhysConst::C / beam.lambda;
    const double omega_tilde = omega / omega_pe;
    const double intensity_w_m2 = beam.intensity_w_cm2 * 1.0e4;
    const double e0 = std::sqrt(2.0 * intensity_w_m2 / (PhysConst::C * PhysConst::EPS0));
    return (e0 / e_scale) * std::sin(omega_tilde * t_tilde);
}

BoundaryDrive boundary_drive(const Config& cfg,
                             double t_tilde,
                             double omega_pe,
                             double e_scale) {
    BoundaryDrive drive;
    const LaserBeamConfig beams[2] = {cfg.laser_beam1, cfg.laser_beam2};

    for (int i = 0; i < 2; ++i) {
        const double ex = beam_ex_tilde(beams[i], t_tilde, omega_pe, e_scale);
        if (ex == 0.0 && !beams[i].enabled) {
            continue;
        }

        if (beams[i].from_left_boundary) {
            drive.drive_left = true;
            drive.left_ex += ex;
        } else {
            drive.drive_right = true;
            drive.right_ex += ex;
        }
    }

    return drive;
}

}
