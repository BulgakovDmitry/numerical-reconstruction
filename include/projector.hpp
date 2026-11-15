#ifndef NUM_RECONSTRUCTION_PROJECTOR_HPP
#define NUM_RECONSTRUCTION_PROJECTOR_HPP

#include "geometry.hpp"
#include "scanner.hpp"

#include <iosfwd>
#include <optional>
#include <vector>

namespace NumReconstruction {

struct RayPixelIntersection {
    Index pixel_index;
    double length;
};

// Piecewise-constant forward model: A(i,j) is the path length of ray i in pixel j.
// Owns small geometry values and, after build_system_matrix(), one cached matrix.
class SiddonProjector {
public:
    SiddonProjector(Grid2D grid, FanBeamGeometry geometry);

    // Segments are returned in source-to-detector order. Misses and point contacts
    // give an empty result. Pixel ownership follows Grid2D, including edge paths.
    [[nodiscard]] std::vector<RayPixelIntersection> trace_ray(const Ray2D& ray) const;

    // Build once; subsequent calls reuse the matrix. The returned reference is
    // owned by this projector. Use a const reference to avoid a large copy.
    const SparseMatrix& build_system_matrix();
    [[nodiscard]] const SparseMatrix& get_system_matrix() const;

    // Require a built matrix and finite vectors of the corresponding size.
    // project returns logarithmic measurements b, not intensities I0*exp(-b).
    [[nodiscard]] Vector project(const Vector& mu) const;
    // Adjoint A^T*b; this is not an inverse/reconstruction.
    [[nodiscard]] Vector backproject(const Vector& measurements) const;

    void dump(std::ostream& out) const;

private:
    Grid2D grid_;
    FanBeamGeometry geometry_;
    std::optional<SparseMatrix> system_matrix_;
};

} // namespace NumReconstruction

#endif
