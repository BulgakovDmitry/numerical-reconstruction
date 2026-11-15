#include "phantom.hpp"
#include "geometry.hpp"
#include "tools.hpp"

#include <cmath>

namespace NumReconstruction {

static bool inside_ellipse(double x, double y, double cx, double cy, double rx, double ry);


// ════ Artificial image ═══════════════════════════════════════════
Vector create_suitcase_phantom(const Grid2D& grid, const double mu_max) {
    Tools::require_positive(mu_max, "Phantom: mu_max");
    Vector mu = Vector::Zero(grid.pixel_count());

    for (int row = 0; row < grid.get_pixels_y(); ++row) {
        for (int col = 0; col < grid.get_pixels_x(); ++col) {
            const auto center = grid.pixel_center(row, col);
            const double x = center.x() / (0.5 * grid.get_width());
            const double y = center.y() / (0.5 * grid.get_height());
            double value = 0.0;

            // Outer shell and low-density contents.
            if (std::abs(x) <= 0.8 && std::abs(y) <= 0.65) {
                value = std::abs(x) > 0.74 || std::abs(y) > 0.59 ? 0.65 : 0.08;
                if (inside_ellipse(x, y, -0.38, 0.16, 0.22, 0.32)) value = 0.35;
                if (inside_ellipse(x, y, 0.35, 0.25, 0.13, 0.13)) value = 0.95;
                if (x >= 0.05 && x <= 0.60 && y >= -0.45 && y <= -0.15) value = 0.55;
                if (inside_ellipse(x, y, -0.35, -0.38, 0.18, 0.08)) value = 0.75;
            }
            // Handle, including its empty center.
            if (std::abs(x) <= 0.24 && y > 0.65 && y <= 0.86 &&
                (std::abs(x) >= 0.17 || y >= 0.79)) {
                value = 0.65;
            }
            mu[grid.index(row, col)] = value * mu_max;
        }
    }

    return mu;
}
// ══════════════════════════════════════════════════════════════════


// ════ Shape helpers ══════════════════════════════════════════════
static bool inside_ellipse(const double x, const double y, const double cx, const double cy,
                           const double rx, const double ry) {
    const double dx = (x - cx) / rx;
    const double dy = (y - cy) / ry;
    return dx * dx + dy * dy <= 1.0;
}
// ══════════════════════════════════════════════════════════════════

} // namespace NumReconstruction
