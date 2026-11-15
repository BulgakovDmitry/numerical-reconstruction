#include "config.hpp"
#include "geometry.hpp"
#include "scanner.hpp"
#include "projector.hpp"
#include "regularization.hpp"
#include "io.hpp"

#include <exception>
#include <iostream>
#include <stdexcept>
#include <string_view>

int main(int argc, char* argv[]) {
    using NumReconstruction::Config;
    using NumReconstruction::Grid2D;
    using NumReconstruction::FanBeamGeometry;
    using NumReconstruction::SiddonProjector;

    if (argc == 2 && (std::string_view(argv[1]) == "--help" || std::string_view(argv[1]) == "-h")) {
        Config::print_help(std::cout);
        return 0;
    }

    try {
        // 1. Configuration
        const auto config = Config::from_command_line(argc, argv);
        if (config.is_poisson_noise_enabled()) {
            throw std::invalid_argument(
                "Poisson noise simulation is not implemented yet; use --no-poisson-noise.");
        }

        // 2. Image grid
        const auto grid = Grid2D::from_config(config);

        // 3. Scanner geometry
        const auto scanner_geometry = FanBeamGeometry::from_config(config);
        NumReconstruction::validate_geometry(grid, scanner_geometry);

        config.dump(std::cout);
        std::cout << '\n';
        grid.dump(std::cout);
        std::cout << '\n';
        scanner_geometry.dump(std::cout);

        // 4. Load the original attenuation image
        constexpr double mu_max = 1.0;
        const auto mu_true = NumReconstruction::load_image(config.get_input_path(), grid, mu_max);
        std::cout << "\nLoaded image: " << config.get_input_path()
                  << " (" << mu_true.size() << " pixels)\n";

        // 5. Build the system matrix
        SiddonProjector projector{grid, scanner_geometry};
        projector.build_system_matrix();
        std::cout << '\n';
        projector.dump(std::cout);

        // 6. Simulate ideal measurements
        const auto b = projector.project(mu_true);
        std::cout << "Ideal measurements: " << b.size() << '\n';

        // 7. Build regularization matrix
       

    } catch (const std::invalid_argument& error) {
        std::cerr << "Error: " << error.what() << '\n'
                  << "Run with --help to see available options.\n";
        return 1;
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }

    return 0;
}
