#ifndef IO_HPP
#define IO_HPP

#include "linalg.hpp"

#include <filesystem>
#include <string_view>

namespace NumReconstruction {

class Grid2D;

[[nodiscard]] Vector load_image(const std::filesystem::path& path,
                                const Grid2D& grid,
                                double mu_max = 1.0);

void save_image(const std::filesystem::path& path,
                const Vector& mu,
                int width,
                int height,
                double mu_max = 1.0);

void save_sinogram(const std::filesystem::path& path,
                   const Vector& measurements,
                   int num_angles,
                   int num_detectors);

// ════ Numeric results ════════════════════════════════════════════
void save_matrix(const std::filesystem::path& path,
                 const Vector& values,
                 int rows,
                 int cols);

void save_report(const std::filesystem::path& path, std::string_view report);
// ══════════════════════════════════════════════════════════════════

} // namespace NumReconstruction

#endif // IO_HPP
