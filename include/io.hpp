#ifndef IO_HPP
#define IO_HPP

#include "linalg.hpp"

#include <filesystem>

namespace NumReconstruction {

class Grid2D;

// Loads a grayscale PGM image (P2 or P5, Maxval 1..65535).
// Dimensions must match the grid; no resizing or gamma correction is applied.
// mu[j] = pixel / Maxval * mu_max, j = row * width + col, row 0 is at the top.
[[nodiscard]] Vector load_image(const std::filesystem::path& path,
                                const Grid2D& grid,
                                double mu_max = 1.0);

// Saves the distribution mu as an 8-bit PGM (P5) image.
// Fixed scale, inverse of loading: pixel = round(mu / mu_max * 255),
// mu <= 0 -> 0 (black), mu >= mu_max -> 255 (white).
// Pixel order matches the input data: index = row * width + col,
// row = 0 is the top row. Parent directories are created if needed.
void save_image(const std::filesystem::path& path,
                const Vector& mu,
                int width,
                int height,
                double mu_max = 1.0);

// Saves a numeric CSV matrix without a header: rows are views, columns detectors.
// measurements[view * num_detectors + detector] is written without normalization.
void save_sinogram(const std::filesystem::path& path,
                   const Vector& measurements,
                   int num_angles,
                   int num_detectors);

} // namespace NumReconstruction

#endif // IO_HPP
