
#include "config.hpp"
#include <iostream>
#include <stdexcept>
#include <string_view>

int main(int argc, char* argv[]) {
    using NumReconstruction::Config;

    if (argc == 2 && (std::string_view(argv[1]) == "--help" || std::string_view(argv[1]) == "-h")) {
        Config::print_help(std::cout);
        return 0;
    }

    try {
        const auto config = Config::from_command_line(argc, argv);

        config.dump(std::cout);
    } catch (const std::invalid_argument& error) {
        std::cerr << "Error: " << error.what() << '\n'
                  << "Run with --help to see available options.\n";
        return 1;
    }
}
