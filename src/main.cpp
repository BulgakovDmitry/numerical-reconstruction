
#include "config.hpp"
#include "geometry.hpp"
#include "scanner.hpp"
#include "projector.hpp"

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
        const auto config = Config::from_command_line(argc, argv);

        // Build the image grid and scanner from the same experiment settings.
        const auto grid = Grid2D::from_config(config);
        const auto geometry = FanBeamGeometry::from_config(config);
        NumReconstruction::validate_geometry(grid, geometry);

        config.dump(std::cout);
        std::cout << '\n';
        grid.dump(std::cout);
        std::cout << '\n';
        geometry.dump(std::cout);

        // Geometry -> discrete forward operator A. Image loading and b=A*mu
        // simulation will follow once the input module is implemented.
        SiddonProjector projector{grid, geometry};
        projector.build_system_matrix();
        std::cout << '\n';
        projector.dump(std::cout);

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
