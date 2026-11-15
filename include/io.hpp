#ifndef IO_HPP
#define IO_HPP

#include <filesystem>

#include <Eigen/Core>

namespace NumReconstruction {

// Saves the distribution mu as an 8-bit PGM (P5) image.
// Fixed scale, inverse of loading: pixel = round(mu / mu_max * 255),
// mu <= 0 -> 0 (black), mu >= mu_max -> 255 (white).
// Pixel order matches the input data: index = row * width + col,
// row = 0 is the top row. Parent directories are created if needed.
void save_image(const std::filesystem::path& path,
                const Eigen::VectorXd& mu,
                int width,
                int height,
                double mu_max = 1.0);

} // namespace NumReconstruction

#endif // IO_HPP
