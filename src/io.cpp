#include "io.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace NumReconstruction {

void save_image(const std::filesystem::path& path,
                const Eigen::VectorXd& mu,
                const int width,
                const int height,
                const double mu_max) {
    if (width <= 0 || height <= 0) {
        throw std::invalid_argument("save_image: width and height must be greater than zero");
    }

    const Eigen::Index pixel_count = static_cast<Eigen::Index>(width) * height;
    if (mu.size() != pixel_count) {
        throw std::invalid_argument("save_image: vector size " + std::to_string(mu.size())
                                    + " does not match " + std::to_string(width) + "x"
                                    + std::to_string(height));
    }

    if (!std::isfinite(mu_max) || mu_max <= 0.0) {
        throw std::invalid_argument("save_image: mu_max must be finite and greater than zero");
    }

    // mu -> intensity: map [0, mu_max] linearly onto [0, 255] with clamping and rounding.
    const double scale = 255.0 / mu_max;
    std::vector<std::uint8_t> pixels(static_cast<std::size_t>(pixel_count));
    for (Eigen::Index i = 0; i < pixel_count; ++i) {
        const double value = mu[i];
        if (!std::isfinite(value)) {
            throw std::invalid_argument("save_image: mu[" + std::to_string(i) + "] is NaN or Inf");
        }
        const double level = std::clamp(value * scale, 0.0, 255.0);
        pixels[static_cast<std::size_t>(i)] = static_cast<std::uint8_t>(std::lround(level));
    }

    if (path.has_parent_path()) {
        std::filesystem::create_directories(path.parent_path());
    }

    std::ofstream file(path, std::ios::binary);
    if (!file) {
        throw std::runtime_error("save_image: cannot open " + path.string());
    }

    // PGM P5 header, then rows top to bottom - exactly the storage order of the vector.
    file << "P5\n" << width << ' ' << height << "\n255\n";
    file.write(reinterpret_cast<const char*>(pixels.data()),
               static_cast<std::streamsize>(pixels.size()));
    if (!file) {
        throw std::runtime_error("save_image: failed to write " + path.string());
    }
}

} // namespace NumReconstruction
