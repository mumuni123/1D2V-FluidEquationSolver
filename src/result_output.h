#ifndef RESULT_OUTPUT_H
#define RESULT_OUTPUT_H

#include <fstream>
#include <string>
#include <vector>

namespace ResultOutput {

bool ensure_output_dir(const std::string& output_dir);

bool open_diagnostics(const std::string& output_dir, std::ofstream& diag);

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
                     const std::vector<double>& Pe);

void append_diagnostics(std::ofstream& diag,
                        double time,
                        int nx,
                        double dz,
                        const std::vector<double>& Ex,
                        const std::vector<double>& Ez,
                        const std::vector<double>& By,
                        const std::vector<double>& ne,
                        const std::vector<double>& vx,
                        const std::vector<double>& vz);

}

#endif
