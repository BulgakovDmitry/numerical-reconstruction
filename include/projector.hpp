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

class SiddonProjector {
public:
    SiddonProjector(Grid2D grid, FanBeamGeometry geometry);

    [[nodiscard]] std::vector<RayPixelIntersection> trace_ray(const Ray2D& ray) const;

    const SparseMatrix& build_system_matrix();
    [[nodiscard]] const SparseMatrix& get_system_matrix() const;

    [[nodiscard]] Vector project(const Vector& mu) const;
    [[nodiscard]] Vector backproject(const Vector& measurements) const;

    void dump(std::ostream& out) const;

private:
    Grid2D grid_;
    FanBeamGeometry geometry_;
    std::optional<SparseMatrix> system_matrix_;
};

} // namespace NumReconstruction

#endif
