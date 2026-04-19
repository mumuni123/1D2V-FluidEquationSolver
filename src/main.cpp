#include <iostream>
#include <exception>

#include "config.h"
#include "fluid_maxwell_1d.h"

int main() {
    try {
        Config cfg;
        FluidMaxwell1D sim(cfg);
        sim.print_summary();
        sim.run();
        std::cout << "Finished. Output written to " << cfg.output_dir << "/" << std::endl;
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "Simulation failed: " << ex.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "Simulation failed: unknown exception" << std::endl;
        return 1;
    }
}
