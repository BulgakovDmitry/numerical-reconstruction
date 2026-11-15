#include "config.hpp"

#include <cmath>
#include <stdexcept>
#include <string>

namespace NumReconstruction {

static void require_positive(const int    value, const char* name);
static void require_positive(const double value, const char* name);
static void require_finite  (const double value, const char* name);

void Config::validate() const {
    require_positive(pixels_x_,     "pixels_x"    );
    require_positive(pixels_y_,     "pixels_y"    );
    require_positive(image_width_,  "image_width" );
    require_positive(image_height_, "image_height");

    require_positive(radius_,            "radius"           );
    require_positive(detector_distance_, "detector_distance");
    require_positive(num_angles_,        "num_angles"       );
    require_positive(num_detectors_,     "num_detectors"    );
    require_finite(start_angle_,         "start_angle"      );
    require_positive(scan_angle_,        "scan_angle"       );
    require_positive(half_fan_angle_,    "half_fan_angle"   );

    if (half_fan_angle_ >= std::numbers::pi / 2.0)
        throw std::invalid_argument("Config: half_fan_angle must be less than pi/2 radians");

    require_finite(lambda_, "lambda");
    if (lambda_ < 0.0)
        throw std::invalid_argument("Config: lambda must be nonnegative");

    switch (regularization_) {
        case RegularizationType::ZeroOrder:
        case RegularizationType::FirstOrder:
            break;
        default:
            throw std::invalid_argument("Config: regularization must be ZeroOrder or FirstOrder");
    }

    require_positive(epsilon_,               "epsilon"              );
    require_positive(max_iterations_,        "max_iterations"       );
    require_positive(incident_photon_count_, "incident_photon_count");

    if (input_path_.empty()) {
        throw std::invalid_argument("Config: input_path must not be empty");
    }
    if (output_directory_.empty()) {
        throw std::invalid_argument("Config: output_directory must not be empty");
    }
}

static void require_positive(const int value, const char* name) {
    if (value <= 0) {
        throw std::invalid_argument(std::string("Config: ") + name + " must be greater than zero");
    }
}

static void require_finite(const double value, const char* name) {
    if (!std::isfinite(value)) {
        throw std::invalid_argument(std::string("Config: ") + name + " must be finite");
    }
}

static void require_positive(const double value, const char* name) {
    require_finite(value, name);
    if (value <= 0.0) {
        throw std::invalid_argument(std::string("Config: ") + name + " must be greater than zero");
    }
}

} // namespace NumReconstruction
