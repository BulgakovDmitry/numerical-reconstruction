#ifndef NUM_RECONSTRUCTION_SCANNER_HPP
#define NUM_RECONSTRUCTION_SCANNER_HPP

#include "geometry.hpp"

#include <numbers>

namespace NumReconstruction {

struct MeasurementIndex {
    int view;
    int detector;
    bool operator==(const MeasurementIndex&) const = default;
};

// Detector arc is centered at the source, following the explicit README formula.
// Views sample [start_angle,start_angle+scan_angle); detectors sample bin centers.
class FanBeamGeometry {
public:
    FanBeamGeometry(double radius, double detector_distance,
                    int num_angles, int num_detectors, double half_fan_angle,
                    double start_angle = 0.0,
                    double scan_angle = 2.0 * std::numbers::pi);
    [[nodiscard]] static FanBeamGeometry from_config(const Config& config);

    [[nodiscard]] double get_radius() const noexcept { return radius_; }
    [[nodiscard]] double get_detector_distance() const noexcept { return detector_distance_; }
    [[nodiscard]] int get_num_angles() const noexcept { return num_angles_; }
    [[nodiscard]] int get_num_detectors() const noexcept { return num_detectors_; }
    [[nodiscard]] double get_start_angle() const noexcept { return start_angle_; }
    [[nodiscard]] double get_scan_angle() const noexcept { return scan_angle_; }
    [[nodiscard]] double get_half_fan_angle() const noexcept { return half_fan_angle_; }
    [[nodiscard]] Index measurement_count() const noexcept { return measurement_count_; }

    [[nodiscard]] double angle(int view) const;
    [[nodiscard]] double fan_angle(int detector) const;
    [[nodiscard]] Vec2 source_position(int view) const;
    [[nodiscard]] Vec2 detector_position(int view, int detector) const;
    [[nodiscard]] Ray2D ray(int view, int detector) const;
    [[nodiscard]] Index measurement_index(int view, int detector) const;
    [[nodiscard]] MeasurementIndex measurement_indices(Index index) const;

private:
    double radius_;
    double detector_distance_;
    int num_angles_;
    int num_detectors_;
    double half_fan_angle_;
    double start_angle_;
    double scan_angle_;
    Index measurement_count_;
};

// Conservative sufficient conditions using the circumscribed circle of the grid.
// Also applied to partial scans; does not guarantee enough discrete measurements.
// Throws std::invalid_argument if source clearance, detector distance or fan
// aperture is insufficient. Constructors do not invoke this cross-object check.
void validate_geometry(const Grid2D& grid, const FanBeamGeometry& scanner);

} // namespace NumReconstruction

#endif
