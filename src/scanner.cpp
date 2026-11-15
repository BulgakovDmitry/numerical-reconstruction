#include "scanner.hpp"
#include "config.hpp"
#include "tools.hpp"

#include <cmath>
#include <limits>
#include <ostream>
#include <stdexcept>

namespace NumReconstruction {
namespace {

Index checked_measurement_count(int views, int detectors) {
    Tools::require_positive(views, "FanBeamGeometry: num_angles");
    Tools::require_positive(detectors, "FanBeamGeometry: num_detectors");
    if (static_cast<Index>(views) > std::numeric_limits<Index>::max() / detectors) {
        throw std::overflow_error("FanBeamGeometry: measurement count exceeds Eigen::Index");
    }
    return static_cast<Index>(views) * detectors;
}

} // namespace

FanBeamGeometry::FanBeamGeometry(double radius, double detector_distance,
                                 int num_angles, int num_detectors, double half_fan_angle,
                                 double start_angle, double scan_angle)
    : radius_(radius), detector_distance_(detector_distance), num_angles_(num_angles),
      num_detectors_(num_detectors), half_fan_angle_(half_fan_angle),
      start_angle_(start_angle), scan_angle_(scan_angle),
      measurement_count_(checked_measurement_count(num_angles, num_detectors)) {
    Tools::require_positive(radius_, "FanBeamGeometry: radius");
    Tools::require_positive(detector_distance_, "FanBeamGeometry: distance");
    Tools::require_positive(half_fan_angle_, "FanBeamGeometry: half fan angle");
    if (half_fan_angle_ >= std::numbers::pi / 2.0) {
        throw std::invalid_argument("FanBeamGeometry: half fan angle must be less than pi/2");
    }
    Tools::require_positive(scan_angle_, "FanBeamGeometry: scan angle");
    Tools::require_finite(start_angle_, "FanBeamGeometry: start_angle");
    if (scan_angle_ > 2.0 * std::numbers::pi) {
        throw std::invalid_argument("FanBeamGeometry: invalid scan interval");
    }
    // All detector coordinates are bounded in magnitude by radius+distance.
    if (radius_ > std::numeric_limits<double>::max() - detector_distance_) {
        throw std::invalid_argument("FanBeamGeometry: coordinate scale is too large");
    }
    // Reject angular scales at which successive views cannot be represented.
    const double step = scan_angle_ / num_angles_;
    if (step == 0.0 || (num_angles_ > 1 &&
        (start_angle_ + step == start_angle_ ||
         angle(num_angles_ - 1) == angle(num_angles_ - 2)))) {
        throw std::invalid_argument("FanBeamGeometry: view angles are not distinguishable");
    }
    // In particular, a subnormal distance can round detector back to source.
    (void)ray(0, 0);
}

FanBeamGeometry FanBeamGeometry::from_config(const Config& config) {
    return {config.get_radius(), config.get_detector_distance(), config.get_num_angles(),
            config.get_num_detectors(), config.get_half_fan_angle(),
            config.get_start_angle(), config.get_scan_angle()};
}

void FanBeamGeometry::dump(std::ostream& out) const {
    out << "Fan-beam geometry:\n"
        << "Source radius: " << radius_ << " m\n"
        << "Source-detector distance: " << detector_distance_ << " m\n"
        << "Number of angles: " << num_angles_ << '\n'
        << "Number of detectors: " << num_detectors_ << '\n'
        << "Start angle: " << start_angle_ << " rad\n"
        << "Scan angle: " << scan_angle_ << " rad\n"
        << "Half fan angle: " << half_fan_angle_ << " rad\n"
        << "Ray count: " << measurement_count_ << '\n';
}

double FanBeamGeometry::angle(int view) const {
    Tools::check_index(view, num_angles_, "FanBeamGeometry: view index out of range");
    return start_angle_ + (static_cast<double>(view) / num_angles_) * scan_angle_;
}

double FanBeamGeometry::fan_angle(int detector) const {
    Tools::check_index(detector, num_detectors_, "FanBeamGeometry: detector index out of range");
    return (2.0 * ((static_cast<double>(detector) + 0.5) / num_detectors_) - 1.0)
        * half_fan_angle_;
}

Vec2 FanBeamGeometry::source_position(int view) const {
    const double theta = angle(view);
    return {radius_ * std::cos(theta), radius_ * std::sin(theta)};
}

Vec2 FanBeamGeometry::detector_position(int view, int detector) const {
    const double theta = angle(view);
    const double gamma = fan_angle(detector);
    const double c = std::cos(theta);
    const double s = std::sin(theta);
    const Vec2 source{radius_ * c, radius_ * s};
    const Vec2 u{-c, -s};
    const Vec2 v{-s, c};
    const Vec2 result = source + detector_distance_ *
        (u * std::cos(gamma) + v * std::sin(gamma));
    if (!result.allFinite()) {
        throw std::invalid_argument("FanBeamGeometry: detector coordinates are not finite");
    }
    return result;
}

Ray2D FanBeamGeometry::ray(int view, int detector) const {
    return {source_position(view), detector_position(view, detector)};
}

Index FanBeamGeometry::measurement_index(int view, int detector) const {
    Tools::check_index(view, num_angles_, "FanBeamGeometry: view index out of range");
    Tools::check_index(detector, num_detectors_, "FanBeamGeometry: detector index out of range");
    return static_cast<Index>(view) * num_detectors_ + detector;
}

MeasurementIndex FanBeamGeometry::measurement_indices(Index index) const {
    Tools::check_index(index, measurement_count_, "FanBeamGeometry: measurement index out of range");
    return {static_cast<int>(index / num_detectors_), static_cast<int>(index % num_detectors_)};
}

void validate_geometry(const Grid2D& grid, const FanBeamGeometry& scanner) {
    const double rho = std::hypot(grid.get_width() / 2.0, grid.get_height() / 2.0);
    const double radius = scanner.get_radius();
    if (radius <= rho) {
        throw std::invalid_argument("Geometry: source must lie outside the grid's circumscribed circle");
    }
    // Subtraction avoids overflowing radius+rho.
    if (scanner.get_detector_distance() <= radius ||
        scanner.get_detector_distance() - radius <= rho) {
        throw std::invalid_argument("Geometry: detector distance must exceed radius+grid radius");
    }
    if (scanner.get_half_fan_angle() < std::asin(rho / radius)) {
        throw std::invalid_argument("Geometry: fan aperture does not cover the grid's circumscribed circle");
    }
}

} // namespace NumReconstruction
