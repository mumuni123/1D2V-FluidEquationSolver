#include "result_output.h"
#include "physics_constants.h"

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <sstream>

namespace ResultOutput {

bool ensure_output_dir(const std::string& output_dir) {
#ifdef _WIN32
    std::string cmd = "if not exist \"" + output_dir + "\" mkdir \"" + output_dir + "\"";
#else
    std::string cmd = "mkdir -p \"" + output_dir + "\"";
#endif
    return std::system(cmd.c_str()) == 0;
}

bool open_diagnostics(const std::string& output_dir, std::ofstream& diag) {
    diag.open((output_dir + "/diagnostics.csv").c_str());
    if (!diag) {
        return false;
    }
    diag << "time_fs,n_min,n_max,vx_rms,vz_rms,energy_em,energy_fluid\n";
    return true;
}

void write_state_csv(const std::string& output_dir,
                     int step,
                     int nx,
                     const std::vector<double>& z,
                     const std::vector<double>& Ex,
                     const std::vector<double>& Ez,
                     const std::vector<double>& By,
                     const std::vector<double>& ne,
                     const std::vector<double>& vx,
                     const std::vector<double>& vz,
                     const std::vector<double>& Pe) {
    std::ostringstream oss;
    oss << output_dir << "/state_" << std::setw(7) << std::setfill('0') << step << ".csv";

    std::ofstream out(oss.str().c_str());
    if (!out) {
        return;
    }

    out << "z,Ex,Ez,By,ne,vx,vz,Pe\n";
    for (int i = 0; i < nx; ++i) {
        const double by_cell = 0.5 * (By[i] + By[i + 1]);
        out << std::setprecision(10)
            << z[i] << ","
            << Ex[i] << ","
            << Ez[i] << ","
            << by_cell << ","
            << ne[i] << ","
            << vx[i] << ","
            << vz[i] << ","
            << Pe[i] << "\n";
    }
}

void append_diagnostics(std::ofstream& diag,
                        double time,
                        int nx,
                        double dz,
                        const std::vector<double>& Ex,
                        const std::vector<double>& Ez,
                        const std::vector<double>& By,
                        const std::vector<double>& ne,
                        const std::vector<double>& vx,
                        const std::vector<double>& vz) {
    if (!diag) {
        return;
    }

    double n_min = ne[0];
    double n_max = ne[0];
    double vx2_sum = 0.0;
    double vz2_sum = 0.0;
    double em_sum = 0.0;
    double fluid_sum = 0.0;

    for (int i = 0; i < nx; ++i) {
        if (ne[i] < n_min) n_min = ne[i];
        if (ne[i] > n_max) n_max = ne[i];

        vx2_sum += vx[i] * vx[i];
        vz2_sum += vz[i] * vz[i];

        const double by_cell = 0.5 * (By[i] + By[i + 1]);
        em_sum += 0.5 * PhysConst::EPS0 * (Ex[i] * Ex[i] + Ez[i] * Ez[i])
            + 0.5 * (by_cell * by_cell) / PhysConst::MU0;
        fluid_sum += 0.5 * PhysConst::ME * ne[i] * (vx[i] * vx[i] + vz[i] * vz[i]);
    }

    const double vx_rms = std::sqrt(vx2_sum / static_cast<double>(nx));
    const double vz_rms = std::sqrt(vz2_sum / static_cast<double>(nx));
    const double energy_em = em_sum * dz;
    const double energy_fluid = fluid_sum * dz;

    diag << std::setprecision(10)
         << (time * 1.0e15) << ","
         << n_min << ","
         << n_max << ","
         << vx_rms << ","
         << vz_rms << ","
         << energy_em << ","
         << energy_fluid << "\n";
}

}
