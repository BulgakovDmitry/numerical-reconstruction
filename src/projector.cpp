#include "projector.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <utility>

namespace NumReconstruction {
namespace {

// Extended intermediates keep very short intervals near pixel corners distinct.
// No tolerance from the reconstruction solver is used for geometric decisions.
using Parameter = long double;

int locate_axis(Parameter coordinate, int cells, const auto& edge, bool descending) {
    int first = 0;
    int last = cells;
    while (last - first > 1) {
        const int middle = first + (last - first) / 2;
        const bool beyond = descending ? coordinate <= edge(middle)
                                       : coordinate >= edge(middle);
        if (beyond) first = middle;
        else last = middle;
    }
    return first;
}

} // namespace

SiddonProjector::SiddonProjector(Grid2D grid, FanBeamGeometry geometry)
    : grid_(std::move(grid)), geometry_(std::move(geometry)) {
    validate_geometry(grid_, geometry_);
    const auto limit = std::numeric_limits<SparseMatrix::StorageIndex>::max();
    if (grid_.pixel_count() > limit || geometry_.measurement_count() > limit) {
        throw std::overflow_error("SiddonProjector: matrix dimensions exceed sparse storage index range");
    }
}

std::vector<RayPixelIntersection> SiddonProjector::trace_ray(const Ray2D& ray) const {
    const auto box = grid_.bounds();
    const auto delta = ray.displacement();
    const Parameter sx = ray.get_source().x();
    const Parameter sy = ray.get_source().y();
    const Parameter dx = delta.x();
    const Parameter dy = delta.y();

    // A segment along a grid line belongs to just one side. The right/bottom
    // external boundaries are excluded by the same convention as Grid2D::locate.
    if ((dx == 0.0L && (sx < box.xmin || sx >= box.xmax)) ||
        (dy == 0.0L && (sy <= box.ymin || sy > box.ymax))) {
        return {};
    }

    Parameter enter = 0.0L;
    Parameter exit = 1.0L;
    const auto clip_axis = [&](Parameter source, Parameter direction,
                               Parameter lower, Parameter upper) {
        if (direction == 0.0L) return;
        Parameter a = (lower - source) / direction;
        Parameter b = (upper - source) / direction;
        if (a > b) std::swap(a, b);
        enter = std::max(enter, a);
        exit = std::min(exit, b);
    };
    clip_axis(sx, dx, box.xmin, box.xmax);
    clip_axis(sy, dy, box.ymin, box.ymax);
    if (enter >= exit) return {};

    std::vector<Parameter> crossings{enter, exit};
    if (dx != 0.0L) {
        for (int col = 1; col < grid_.get_pixels_x(); ++col) {
            const Parameter t = (static_cast<Parameter>(grid_.x_edge(col)) - sx) / dx;
            if (t > enter && t < exit) crossings.push_back(t);
        }
    }
    if (dy != 0.0L) {
        for (int row = 1; row < grid_.get_pixels_y(); ++row) {
            const Parameter t = (static_cast<Parameter>(grid_.y_edge(row)) - sy) / dy;
            if (t > enter && t < exit) crossings.push_back(t);
        }
    }
    std::sort(crossings.begin(), crossings.end());
    crossings.erase(std::unique(crossings.begin(), crossings.end()), crossings.end());

    std::vector<RayPixelIntersection> result;
    result.reserve(crossings.size() - 1);
    for (std::size_t i = 1; i < crossings.size(); ++i) {
        const Parameter t0 = crossings[i - 1];
        const Parameter t1 = crossings[i];
        const Parameter middle = t0 + (t1 - t0) / 2.0L;
        const Parameter x = sx + middle * dx;
        const Parameter y = sy + middle * dy;
        const int col = locate_axis(x, grid_.get_pixels_x(),
                                    [&](int c) { return grid_.x_edge(c); }, false);
        const int row = locate_axis(y, grid_.get_pixels_y(),
                                    [&](int r) { return grid_.y_edge(r); }, true);
        const Index pixel = grid_.index(row, col);
        const double length = static_cast<double>((t1 - t0) * ray.length());
        if (length == 0.0) continue; // Below the representable double length scale.
        if (!result.empty() && result.back().pixel_index == pixel) {
            result.back().length += length;
        } else {
            result.push_back({pixel, length});
        }
    }
    return result;
}

const SparseMatrix& SiddonProjector::build_system_matrix() {
    if (system_matrix_) return *system_matrix_;

    std::vector<Triplet> entries;
    const auto limit = static_cast<std::size_t>(
        std::numeric_limits<SparseMatrix::StorageIndex>::max());
    for (int view = 0; view < geometry_.get_num_angles(); ++view) {
        for (int detector = 0; detector < geometry_.get_num_detectors(); ++detector) {
            const Index measurement = geometry_.measurement_index(view, detector);
            const auto intersections = trace_ray(geometry_.ray(view, detector));
            if (intersections.size() > limit - entries.size()) {
                throw std::overflow_error("SiddonProjector: nonzero count exceeds sparse storage index range");
            }
            for (const auto& intersection : intersections) {
                entries.emplace_back(measurement, intersection.pixel_index, intersection.length);
            }
        }
    }

    // Build locally so a failed allocation does not leave a partial cached matrix.
    SparseMatrix matrix(geometry_.measurement_count(), grid_.pixel_count());
    matrix.setFromTriplets(entries.begin(), entries.end());
    matrix.makeCompressed();
    system_matrix_.emplace(std::move(matrix));
    return *system_matrix_;
}

const SparseMatrix& SiddonProjector::get_system_matrix() const {
    if (!system_matrix_) {
        throw std::logic_error("SiddonProjector: call build_system_matrix first");
    }
    return *system_matrix_;
}

Vector SiddonProjector::project(const Vector& mu) const {
    const auto& matrix = get_system_matrix();
    if (mu.size() != matrix.cols() || !mu.allFinite()) {
        throw std::invalid_argument("SiddonProjector: mu must be finite with one value per pixel");
    }
    return matrix * mu;
}

Vector SiddonProjector::backproject(const Vector& measurements) const {
    const auto& matrix = get_system_matrix();
    if (measurements.size() != matrix.rows() || !measurements.allFinite()) {
        throw std::invalid_argument("SiddonProjector: measurements must be finite with one value per ray");
    }
    return matrix.transpose() * measurements;
}

void SiddonProjector::dump(std::ostream& out) const {
    out << "Siddon projector:\n"
        << "System matrix size: " << geometry_.measurement_count()
        << " x " << grid_.pixel_count() << '\n'
        << "System matrix built: " << (system_matrix_ ? "true" : "false") << '\n';
    if (system_matrix_) {
        out << "Nonzero elements: " << system_matrix_->nonZeros() << '\n';
    }
}

} // namespace NumReconstruction
