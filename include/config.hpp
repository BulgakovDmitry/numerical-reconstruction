#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <cstdint>
#include <filesystem>
#include <numbers>

namespace NumReconstruction {

class Config {
public:
    enum class RegularizationType {
        ZeroOrder,   
        FirstOrder,
    };
private:
    // ════ Image grid ══════════════════════════════════════════════
    int pixels_x_ = 256;
    int pixels_y_ = 256;

    double image_width_  = 2.0; //in meters
    double image_height_ = 2.0; //in meters
    // ══════════════════════════════════════════════════════════════


    // ════ scanner geometry ════════════════════════════════════════
    double radius_            = 2.5;
    double detector_distance_ = 5.5;

    int num_angles_    = 360;
    int num_detectors_ = 256;

    double start_angle_    = 0.0;
    double scan_angle_     = 2.0 * std::numbers::pi;
    double half_fan_angle_ = 35.0 * std::numbers::pi / 180.0;
    // ══════════════════════════════════════════════════════════════


    // ════ Reconstruction parameters ═══════════════════════════════
    double lambda_ = 0.01;

    RegularizationType regularization_ = RegularizationType::FirstOrder;

    double epsilon_     = 1e-8;
    int max_iterations_ = 2000;
    // ══════════════════════════════════════════════════════════════


    // ════ Measurement simulation ══════════════════════════════════
    bool use_poisson_noise_       = false;
    double incident_photon_count_ = 100000.0;
    std::uint32_t random_seed_    = 42;
    // ══════════════════════════════════════════════════════════════


    // ════ Input / Output ══════════════════════════════════════════
    std::filesystem::path input_path_       = "input/high/04_suitcase.pgm";
    std::filesystem::path output_directory_ = "output";
    // ══════════════════════════════════════════════════════════════

public:
    // ════ Builder ═════════════════════════════════════════════════
    struct Builder;
    // ══════════════════════════════════════════════════════════════

    // ════ Validation ══════════════════════════════════════════════
    void validate() const;
    // ══════════════════════════════════════════════════════════════


    // ════ Getters ═════════════════════════════════════════════════
    [[nodiscard]] int    get_pixels_x    () const noexcept { return pixels_x_;     }
    [[nodiscard]] int    get_pixels_y    () const noexcept { return pixels_y_;     }
    [[nodiscard]] double get_image_width () const noexcept { return image_width_;  }
    [[nodiscard]] double get_image_height() const noexcept { return image_height_; }

    [[nodiscard]] double get_radius           () const noexcept { return radius_;            }
    [[nodiscard]] double get_detector_distance() const noexcept { return detector_distance_; }
    [[nodiscard]] int    get_num_angles       () const noexcept { return num_angles_;        }
    [[nodiscard]] int    get_num_detectors    () const noexcept { return num_detectors_;     }
    [[nodiscard]] double get_start_angle      () const noexcept { return start_angle_;       }
    [[nodiscard]] double get_scan_angle       () const noexcept { return scan_angle_;        }
    [[nodiscard]] double get_half_fan_angle   () const noexcept { return half_fan_angle_;    }

    [[nodiscard]] double             get_lambda             () const noexcept { return lambda_;         }
    [[nodiscard]] RegularizationType get_regularization_type() const noexcept { return regularization_; }
    [[nodiscard]] double             get_epsilon            () const noexcept { return epsilon_;        }
    [[nodiscard]] int                get_max_iterations     () const noexcept { return max_iterations_; }

    [[nodiscard]] bool          is_poisson_noise_enabled () const noexcept { return use_poisson_noise_;     }
    [[nodiscard]] double        get_incident_photon_count() const noexcept { return incident_photon_count_; }
    [[nodiscard]] std::uint32_t get_random_seed          () const noexcept { return random_seed_;           }

    [[nodiscard]] std::filesystem::path get_input_path      () const { return input_path_;       }
    [[nodiscard]] std::filesystem::path get_output_directory() const { return output_directory_; }
    // ══════════════════════════════════════════════════════════════
};

struct Config::Builder {
    Config config_;

    // ════ Image grid ══════════════════════════════════════════════
    Builder& set_pixels_x(const int pixels_x) {
        config_.pixels_x_ = pixels_x;
        return *this;
    }

    Builder& set_pixels_y(const int pixels_y) {
        config_.pixels_y_ = pixels_y;
        return *this;
    }

    Builder& set_image_width(const double image_width) {
        config_.image_width_ = image_width;
        return *this;
    }

    Builder& set_image_height(const double image_height) {
        config_.image_height_ = image_height;
        return *this;
    }
    // ══════════════════════════════════════════════════════════════


    // ════ scanner geometry ════════════════════════════════════════
    Builder& set_radius(const double radius) {
        config_.radius_ = radius;
        return *this;
    }

    Builder& set_detector_distance(const double detector_distance) {
        config_.detector_distance_ = detector_distance;
        return *this;
    }

    Builder& set_num_angles(const int num_angles) {
        config_.num_angles_ = num_angles;
        return *this;
    }

    Builder& set_num_detectors(const int num_detectors) {
        config_.num_detectors_ = num_detectors;
        return *this;
    }

    Builder& set_start_angle(const double start_angle) {
        config_.start_angle_ = start_angle;
        return *this;
    }

    Builder& set_scan_angle(const double scan_angle) {
        config_.scan_angle_ = scan_angle;
        return *this;
    }

    Builder& set_half_fan_angle(const double half_fan_angle) {
        config_.half_fan_angle_ = half_fan_angle;
        return *this;
    }
    // ══════════════════════════════════════════════════════════════


    // ════ Reconstruction parameters ═══════════════════════════════
    Builder& set_lambda(const double lambda) {
        config_.lambda_ = lambda;
        return *this;
    }

    Builder& set_regularization(const Config::RegularizationType regularization) {
        config_.regularization_ = regularization;
        return *this;
    }

    Builder& set_epsilon(const double epsilon) {
        config_.epsilon_ = epsilon;
        return *this;
    }

    Builder& set_max_iterations(const int max_iterations) {
        config_.max_iterations_ = max_iterations;
        return *this;
    }
    // ══════════════════════════════════════════════════════════════


    // ════ Measurement simulation ══════════════════════════════════
    Builder& enable_poisson_noise(const bool enable) {
        config_.use_poisson_noise_ = enable;
        return *this;
    }

    Builder& set_incident_photon_count(const double incident_photon_count) {
        config_.incident_photon_count_ = incident_photon_count;
        return *this;
    }

    Builder& set_random_seed(const std::uint32_t random_seed) {
        config_.random_seed_ = random_seed;
        return *this;
    }
    // ══════════════════════════════════════════════════════════════


    // ════ Input / Output ══════════════════════════════════════════
    Builder& set_input_path(const std::filesystem::path& input_path) {
        config_.input_path_ = input_path;
        return *this;
    }

    Builder& set_output_directory(const std::filesystem::path& output_directory) {
        config_.output_directory_ = output_directory;
        return *this;
    }
    // ══════════════════════════════════════════════════════════════

    [[nodiscard]] Config build() const {
        config_.validate();
        return config_;
    }
};


} // namespace NumReconstruction

#endif // CONFIG_HPP
