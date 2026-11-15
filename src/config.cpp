#include "config.hpp"
#include "tools.hpp"

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <charconv>
#include <cmath>
#include <cstdlib>
#include <ostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <type_traits>
#include <unordered_set>

namespace NumReconstruction {

template <typename T>
static T parse_number(const std::string_view text, const std::string_view option);
static Config::RegularizationType parse_regularization(const std::string_view value);


// ════ Getters ═════════════════════════════════════════════════════
int    Config::get_pixels_x    () const noexcept { return pixels_x_;     }
int    Config::get_pixels_y    () const noexcept { return pixels_y_;     }
double Config::get_image_width () const noexcept { return image_width_;  }
double Config::get_image_height() const noexcept { return image_height_; }

double Config::get_radius           () const noexcept { return radius_;            }
double Config::get_detector_distance() const noexcept { return detector_distance_; }
int    Config::get_num_angles       () const noexcept { return num_angles_;        }
int    Config::get_num_detectors    () const noexcept { return num_detectors_;     }
double Config::get_start_angle      () const noexcept { return start_angle_;       }
double Config::get_scan_angle       () const noexcept { return scan_angle_;        }
double Config::get_half_fan_angle   () const noexcept { return half_fan_angle_;    }

double                     Config::get_lambda             () const noexcept { return lambda_;         }
Config::RegularizationType Config::get_regularization_type() const noexcept { return regularization_; }
double                     Config::get_epsilon            () const noexcept { return epsilon_;        }
int                        Config::get_max_iterations     () const noexcept { return max_iterations_; }

bool          Config::is_poisson_noise_enabled () const noexcept { return use_poisson_noise_;     }
double        Config::get_incident_photon_count() const noexcept { return incident_photon_count_; }
std::uint32_t Config::get_random_seed          () const noexcept { return random_seed_;           }

std::filesystem::path Config::get_input_path      () const { return input_path_;       }
std::filesystem::path Config::get_output_directory() const { return output_directory_; }
// ══════════════════════════════════════════════════════════════════


// ════ Builder: Image grid ═════════════════════════════════════════
Config::Builder& Config::Builder::set_pixels_x(const int pixels_x) {
    config_.pixels_x_ = pixels_x;
    return *this;
}

Config::Builder& Config::Builder::set_pixels_y(const int pixels_y) {
    config_.pixels_y_ = pixels_y;
    return *this;
}

Config::Builder& Config::Builder::set_image_width(const double image_width) {
    config_.image_width_ = image_width;
    return *this;
}

Config::Builder& Config::Builder::set_image_height(const double image_height) {
    config_.image_height_ = image_height;
    return *this;
}
// ══════════════════════════════════════════════════════════════════


// ════ Builder: Scanner geometry ═══════════════════════════════════
Config::Builder& Config::Builder::set_radius(const double radius) {
    config_.radius_ = radius;
    return *this;
}

Config::Builder& Config::Builder::set_detector_distance(const double detector_distance) {
    config_.detector_distance_ = detector_distance;
    return *this;
}

Config::Builder& Config::Builder::set_num_angles(const int num_angles) {
    config_.num_angles_ = num_angles;
    return *this;
}

Config::Builder& Config::Builder::set_num_detectors(const int num_detectors) {
    config_.num_detectors_ = num_detectors;
    return *this;
}

Config::Builder& Config::Builder::set_start_angle(const double start_angle) {
    config_.start_angle_ = start_angle;
    return *this;
}

Config::Builder& Config::Builder::set_scan_angle(const double scan_angle) {
    config_.scan_angle_ = scan_angle;
    return *this;
}

Config::Builder& Config::Builder::set_half_fan_angle(const double half_fan_angle) {
    config_.half_fan_angle_ = half_fan_angle;
    return *this;
}
// ══════════════════════════════════════════════════════════════════


// ════ Builder: Reconstruction parameters ══════════════════════════
Config::Builder& Config::Builder::set_lambda(const double lambda) {
    config_.lambda_ = lambda;
    return *this;
}

Config::Builder& Config::Builder::set_regularization(const Config::RegularizationType regularization) {
    config_.regularization_ = regularization;
    return *this;
}

Config::Builder& Config::Builder::set_epsilon(const double epsilon) {
    config_.epsilon_ = epsilon;
    return *this;
}

Config::Builder& Config::Builder::set_max_iterations(const int max_iterations) {
    config_.max_iterations_ = max_iterations;
    return *this;
}
// ══════════════════════════════════════════════════════════════════


// ════ Builder: Measurement simulation ═════════════════════════════
Config::Builder& Config::Builder::enable_poisson_noise(const bool enable) {
    config_.use_poisson_noise_ = enable;
    return *this;
}

Config::Builder& Config::Builder::set_incident_photon_count(const double incident_photon_count) {
    config_.incident_photon_count_ = incident_photon_count;
    return *this;
}

Config::Builder& Config::Builder::set_random_seed(const std::uint32_t random_seed) {
    config_.random_seed_ = random_seed;
    return *this;
}
// ══════════════════════════════════════════════════════════════════


// ════ Builder: Input / Output ═════════════════════════════════════
Config::Builder& Config::Builder::set_input_path(const std::filesystem::path& input_path) {
    config_.input_path_ = input_path;
    return *this;
}

Config::Builder& Config::Builder::set_output_directory(const std::filesystem::path& output_directory) {
    config_.output_directory_ = output_directory;
    return *this;
}
// ══════════════════════════════════════════════════════════════════


// ════ Builder: build ══════════════════════════════════════════════
Config Config::Builder::build() const {
    config_.validate();
    return config_;
}
// ══════════════════════════════════════════════════════════════════

// ════ Command line ═══════════════════════════════════════════════
Config Config::from_command_line(const int argc, const char* const argv[]) {
    Builder builder;
    std::unordered_set<std::string> seen;

    for (int i = 1; i < argc; ++i) {
        std::string option = argv[i];
        std::replace(option.begin(), option.end(), '_', '-');
        if (!option.starts_with("--")) {
            throw std::invalid_argument("Unexpected argument: '" + option + "'");
        }

        const auto key = option == "--no-poisson-noise" ? "--poisson-noise" : option;
        if (!seen.insert(key).second) {
            throw std::invalid_argument("Repeated or conflicting option: " + option);
        }

        const auto next_value = [&]() -> std::string_view {
            if (i + 1 >= argc || std::string_view(argv[i + 1]).starts_with("--")) {
                throw std::invalid_argument(option + ": missing value");
            }
            return argv[++i];
        };

        if (option == "--pixels-x") {
            builder.set_pixels_x(parse_number<int>(next_value(), option));
        } else if (option == "--pixels-y") {
            builder.set_pixels_y(parse_number<int>(next_value(), option));
        } else if (option == "--image-width") {
            builder.set_image_width(parse_number<double>(next_value(), option));
        } else if (option == "--image-height") {
            builder.set_image_height(parse_number<double>(next_value(), option));
        } else if (option == "--radius") {
            builder.set_radius(parse_number<double>(next_value(), option));
        } else if (option == "--detector-distance") {
            builder.set_detector_distance(parse_number<double>(next_value(), option));
        } else if (option == "--num-angles") {
            builder.set_num_angles(parse_number<int>(next_value(), option));
        } else if (option == "--num-detectors") {
            builder.set_num_detectors(parse_number<int>(next_value(), option));
        } else if (option == "--start-angle") {
            builder.set_start_angle(parse_number<double>(next_value(), option));
        } else if (option == "--scan-angle") {
            builder.set_scan_angle(parse_number<double>(next_value(), option));
        } else if (option == "--half-fan-angle") {
            builder.set_half_fan_angle(parse_number<double>(next_value(), option));
        } else if (option == "--lambda") {
            builder.set_lambda(parse_number<double>(next_value(), option));
        } else if (option == "--regularization") {
            builder.set_regularization(parse_regularization(next_value()));
        } else if (option == "--epsilon") {
            builder.set_epsilon(parse_number<double>(next_value(), option));
        } else if (option == "--max-iterations") {
            builder.set_max_iterations(parse_number<int>(next_value(), option));
        } else if (option == "--poisson-noise") {
            builder.enable_poisson_noise(true);
        } else if (option == "--no-poisson-noise") {
            builder.enable_poisson_noise(false);
        } else if (option == "--incident-photon-count") {
            builder.set_incident_photon_count(parse_number<double>(next_value(), option));
        } else if (option == "--random-seed") {
            builder.set_random_seed(parse_number<std::uint32_t>(next_value(), option));
        } else if (option == "--input") {
            builder.set_input_path(std::filesystem::path(next_value()));
        } else if (option == "--output") {
            builder.set_output_directory(std::filesystem::path(next_value()));
        } else {
            throw std::invalid_argument("Unknown option: " + option);
        }
    }

    return builder.build();
}

void Config::print_help(std::ostream& out) {
    const Config defaults;
    out << "Usage: numerical-reconstruction [options]\n"
        << "       numerical-reconstruction --help\n\n"
        << "Options use --name value; underscores in option names are also accepted.\n"
        << "Unspecified values keep their defaults. Repeated options are errors.\n"
        << "Distances are in meters; angles are in radians.\n\n"
        << "Image grid:\n"
        << "  --pixels-x N               Horizontal pixels (> 0); default: " << defaults.pixels_x_ << '\n'
        << "  --pixels-y N               Vertical pixels (> 0); default: " << defaults.pixels_y_ << '\n'
        << "  --image-width X            Image width (> 0); default: " << defaults.image_width_ << '\n'
        << "  --image-height X           Image height (> 0); default: " << defaults.image_height_ << '\n'
        << "\nScanner geometry:\n"
        << "  --radius X                 Source radius (> 0); default: " << defaults.radius_ << '\n'
        << "  --detector-distance X      Source-detector distance (> 0); default: " << defaults.detector_distance_ << '\n'
        << "  --num-angles N             Number of angles (> 0); default: " << defaults.num_angles_ << '\n'
        << "  --num-detectors N          Number of detectors (> 0); default: " << defaults.num_detectors_ << '\n'
        << "  --start-angle X            Start angle; default: " << defaults.start_angle_ << '\n'
        << "  --scan-angle X             Scan angle (> 0); default: " << defaults.scan_angle_ << '\n'
        << "  --half-fan-angle X         Half fan angle (0 < X < pi/2); default: " << defaults.half_fan_angle_ << '\n'
        << "\nReconstruction:\n"
        << "  --lambda X                 Regularization weight (>= 0); default: " << defaults.lambda_ << '\n'
        << "  --regularization TYPE      zero-order or first-order; default: "
        << (defaults.regularization_ == RegularizationType::ZeroOrder ? "zero-order" : "first-order") << '\n'
        << "  --epsilon X                Solver tolerance (> 0); default: " << defaults.epsilon_ << '\n'
        << "  --max-iterations N         Iteration limit (> 0); default: " << defaults.max_iterations_ << '\n'
        << "\nMeasurement simulation:\n"
        << "  --poisson-noise            Enable Poisson noise\n"
        << "  --no-poisson-noise         Disable Poisson noise; default: "
        << (defaults.use_poisson_noise_ ? "enabled" : "disabled") << '\n'
        << "  --incident-photon-count X  Incident photons (> 0); default: " << defaults.incident_photon_count_ << '\n'
        << "  --random-seed N            Seed (0..4294967295); default: " << defaults.random_seed_ << '\n'
        << "\nInput / Output:\n"
        << "  --input PATH               Input image; default: " << defaults.input_path_ << '\n'
        << "  --output PATH              Output directory; default: " << defaults.output_directory_ << '\n'
        << "\n  -h, --help                 Show this help (use on its own)\n";
}
// ══════════════════════════════════════════════════════════════════


// ════ Validation ══════════════════════════════════════════════════
void Config::validate() const {
    Tools::require_positive(pixels_x_,     "Config: pixels_x"    );
    Tools::require_positive(pixels_y_,     "Config: pixels_y"    );
    Tools::require_positive(image_width_,  "Config: image_width" );
    Tools::require_positive(image_height_, "Config: image_height");

    Tools::require_positive(radius_,            "Config: radius"           );
    Tools::require_positive(detector_distance_, "Config: detector_distance");
    Tools::require_positive(num_angles_,        "Config: num_angles"       );
    Tools::require_positive(num_detectors_,     "Config: num_detectors"    );
    Tools::require_finite  (start_angle_,       "Config: start_angle"      );
    Tools::require_positive(scan_angle_,        "Config: scan_angle"       );
    Tools::require_positive(half_fan_angle_,    "Config: half_fan_angle"   );

    if (half_fan_angle_ >= std::numbers::pi / 2.0)
        throw std::invalid_argument("Config: half_fan_angle must be less than pi/2 radians");

    Tools::require_finite(lambda_, "Config: lambda");
    if (lambda_ < 0.0)
        throw std::invalid_argument("Config: lambda must be nonnegative");

    switch (regularization_) {
        case RegularizationType::ZeroOrder:
        case RegularizationType::FirstOrder:
            break;
        default:
            throw std::invalid_argument("Config: regularization must be ZeroOrder or FirstOrder");
    }

    Tools::require_positive(epsilon_,               "Config: epsilon"              );
    Tools::require_positive(max_iterations_,        "Config: max_iterations"       );
    Tools::require_positive(incident_photon_count_, "Config: incident_photon_count");

    if (input_path_.empty()) {
        throw std::invalid_argument("Config: input_path must not be empty");
    }
    if (output_directory_.empty()) {
        throw std::invalid_argument("Config: output_directory must not be empty");
    }
}
// ══════════════════════════════════════════════════════════════════


// ════ Dump ════════════════════════════════════════════════════════
void Config::dump(std::ostream& out) const {
    out << "Configuration:\n"
        << "Pixels: " << pixels_x_ << " x " << pixels_y_ << '\n'
        << "Image size: " << image_width_ << " x " << image_height_ << " m\n"
        << "Radius: " << radius_ << " m\n"
        << "Detector distance: " << detector_distance_ << " m\n"
        << "Number of angles: " << num_angles_ << '\n'
        << "Number of detectors: " << num_detectors_ << '\n'
        << "Start angle: " << start_angle_ << " rad\n"
        << "Scan angle: " << scan_angle_ << " rad\n"
        << "Half fan angle: " << half_fan_angle_ << " rad\n"
        << "Lambda: " << lambda_ << '\n'
        << "Regularization: ";

    switch (regularization_) {
        case RegularizationType::ZeroOrder:
            out << "zero-order";
            break;
        case RegularizationType::FirstOrder:
            out << "first-order";
            break;
        default:
            out << "unknown (" << static_cast<int>(regularization_) << ')';
    }

    out << '\n'
        << "Epsilon: " << epsilon_ << '\n'
        << "Max iterations: " << max_iterations_ << '\n'
        << "Poisson noise enabled: " << (use_poisson_noise_ ? "true" : "false") << '\n'
        << "Incident photon count: " << incident_photon_count_ << '\n'
        << "Random seed: " << random_seed_ << '\n'
        << "Input path: " << input_path_ << '\n'
        << "Output directory: " << output_directory_ << '\n';
}
// ══════════════════════════════════════════════════════════════════


/// ════ Command line parsing helpers ═══════════════════════════════
template <typename T>
static T parse_number(const std::string_view text, const std::string_view option) {
    T value{};
    const char* end = text.data();
    std::errc error{};

    if constexpr (std::is_floating_point_v<T>) {
        // libc++ shipped with Apple clang (<= 15) lacks std::from_chars for double,
        // so floating-point values are parsed with strtod under the same rules:
        // no leading whitespace and no '+' sign.
        const std::string buffer(text);
        char* parse_end = nullptr;
        if (!buffer.empty() && !std::isspace(static_cast<unsigned char>(buffer.front())) && buffer.front() != '+') {
            errno = 0;
            value = static_cast<T>(std::strtod(buffer.c_str(), &parse_end));
        }
        if (parse_end == nullptr || parse_end == buffer.c_str()) {
            error = std::errc::invalid_argument;
        } else {
            end = text.data() + (parse_end - buffer.c_str());
            if (errno == ERANGE) {
                error = std::errc::result_out_of_range;
            }
        }
    } else {
        const auto result = std::from_chars(text.data(), text.data() + text.size(), value);
        end   = result.ptr;
        error = result.ec;
    }

    if (error == std::errc::result_out_of_range) {
        throw std::invalid_argument(std::string(option) + ": number out of range: '" + std::string(text) + "'");
    }
    if (error != std::errc{} || end != text.data() + text.size()) {
        const char* expected = std::is_integral_v<T> ? "an integer" : "a number";
        throw std::invalid_argument(std::string(option) + ": expected " + expected + ", got '" + std::string(text) + "'");
    }
    return value;
}

static Config::RegularizationType parse_regularization(const std::string_view value) {
    if (value == "zero-order") {
        return Config::RegularizationType::ZeroOrder;
    }
    if (value == "first-order") {
        return Config::RegularizationType::FirstOrder;
    }
    throw std::invalid_argument("--regularization: expected zero-order or first-order, got '" + std::string(value) + "'");
}
// ══════════════════════════════════════════════════════════════════

} // namespace NumReconstruction
