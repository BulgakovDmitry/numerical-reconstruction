#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <cstdint>
#include <filesystem>
#include <iosfwd>
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

    int iteration_save_step_ = 0; // save every N-th solver iterate; 0 = off
    // ══════════════════════════════════════════════════════════════

public:
    // ════ Builder ═════════════════════════════════════════════════
    struct Builder;
    // ══════════════════════════════════════════════════════════════

    // ════ Command line ═══════════════════════════════════════════
    [[nodiscard]] static Config from_command_line(int argc, const char* const argv[]);
    static void print_help(std::ostream& out);
    // ══════════════════════════════════════════════════════════════

    // ════ Validation ══════════════════════════════════════════════
    void validate() const;
    // ══════════════════════════════════════════════════════════════

    // ════ Dump ════════════════════════════════════════════════════
    void dump(std::ostream& out) const;
    // ══════════════════════════════════════════════════════════════


    // ════ Getters ═════════════════════════════════════════════════
    [[nodiscard]] int    get_pixels_x    () const noexcept;
    [[nodiscard]] int    get_pixels_y    () const noexcept;
    [[nodiscard]] double get_image_width () const noexcept;
    [[nodiscard]] double get_image_height() const noexcept;

    [[nodiscard]] double get_radius           () const noexcept;
    [[nodiscard]] double get_detector_distance() const noexcept;
    [[nodiscard]] int    get_num_angles       () const noexcept;
    [[nodiscard]] int    get_num_detectors    () const noexcept;
    [[nodiscard]] double get_start_angle      () const noexcept;
    [[nodiscard]] double get_scan_angle       () const noexcept;
    [[nodiscard]] double get_half_fan_angle   () const noexcept;

    [[nodiscard]] double             get_lambda             () const noexcept;
    [[nodiscard]] RegularizationType get_regularization_type() const noexcept;
    [[nodiscard]] double             get_epsilon            () const noexcept;
    [[nodiscard]] int                get_max_iterations     () const noexcept;

    [[nodiscard]] bool          is_poisson_noise_enabled () const noexcept;
    [[nodiscard]] double        get_incident_photon_count() const noexcept;
    [[nodiscard]] std::uint32_t get_random_seed          () const noexcept;

    [[nodiscard]] std::filesystem::path get_input_path      () const;
    [[nodiscard]] std::filesystem::path get_output_directory() const;
    [[nodiscard]] int get_iteration_save_step() const noexcept;
    // ══════════════════════════════════════════════════════════════
};

struct Config::Builder {
    Config config_;

    // ════ Image grid ══════════════════════════════════════════════
    Builder& set_pixels_x    (const int pixels_x);
    Builder& set_pixels_y    (const int pixels_y);

    Builder& set_image_width (const double image_width );
    Builder& set_image_height(const double image_height);
    // ══════════════════════════════════════════════════════════════


    // ════ scanner geometry ════════════════════════════════════════
    Builder& set_radius           (const double radius);
    Builder& set_detector_distance(const double detector_distance);

    Builder& set_num_angles       (const int num_angles);
    Builder& set_num_detectors    (const int num_detectors);

    Builder& set_start_angle      (const double start_angle);
    Builder& set_scan_angle       (const double scan_angle);
    Builder& set_half_fan_angle   (const double half_fan_angle);
    // ══════════════════════════════════════════════════════════════


    // ════ Reconstruction parameters ═══════════════════════════════
    Builder& set_lambda        (const double                     lambda);
    Builder& set_regularization(const Config::RegularizationType regularization);
    Builder& set_epsilon       (const double                     epsilon);
    Builder& set_max_iterations(const int                        max_iterations);
    // ══════════════════════════════════════════════════════════════


    // ════ Measurement simulation ══════════════════════════════════
    Builder& enable_poisson_noise     (const bool          enable);
    Builder& set_incident_photon_count(const double        incident_photon_count);
    Builder& set_random_seed          (const std::uint32_t random_seed);
    // ══════════════════════════════════════════════════════════════


    // ════ Input / Output ══════════════════════════════════════════
    Builder& set_input_path         (const std::filesystem::path& input_path);
    Builder& set_output_directory   (const std::filesystem::path& output_directory);
    Builder& set_iteration_save_step(const int iteration_save_step);
    // ══════════════════════════════════════════════════════════════

    [[nodiscard]] Config build() const;
};


} // namespace NumReconstruction

#endif // CONFIG_HPP
