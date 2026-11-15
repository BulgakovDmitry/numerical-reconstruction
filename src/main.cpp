#include "config.hpp"
#include "geometry.hpp"
#include "scanner.hpp"
#include "projector.hpp"
#include "regularization.hpp"
#include "reconstruction.hpp"
#include "simulation.hpp"
#include "io.hpp"

#include <chrono>
#include <cmath>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

// One conjugate gradient iteration, recorded for the visualization.
struct IterationRecord {
    int iteration;
    double relative_residual;
    double rmse;
    bool frame_saved;
};

std::string iteration_frame_name(const int iteration) {
    std::ostringstream name;
    name << "iter_" << std::setw(4) << std::setfill('0') << iteration << ".pgm";
    return name.str();
}

} // namespace

int main(int argc, char* argv[]) {
    using NumReconstruction::Config;
    using NumReconstruction::Grid2D;
    using NumReconstruction::FanBeamGeometry;
    using NumReconstruction::SiddonProjector;
    using Clock = std::chrono::steady_clock;

    if (argc == 2 && (std::string_view(argv[1]) == "--help" || std::string_view(argv[1]) == "-h")) {
        Config::print_help(std::cout);
        return 0;
    }

    try {
        // 1. Configuration
        const auto config = Config::from_command_line(argc, argv);

        // 2. Image grid
        const auto grid = Grid2D::from_config(config);

        // 3. Scanner geometry
        const auto scanner_geometry = FanBeamGeometry::from_config(config);
        NumReconstruction::validate_geometry(grid, scanner_geometry);

        config.dump(std::cout);
        std::cout << '\n';
        grid.dump(std::cout);
        std::cout << '\n';
        scanner_geometry.dump(std::cout);

        // 4. Load the original attenuation image
        constexpr double mu_max = 1.0;
        const auto mu_true = NumReconstruction::load_image(config.get_input_path(), grid, mu_max);
        std::cout << "\nOriginal image: " << mu_true.size() << " pixels\n";

        // 5. Build the system matrix
        SiddonProjector projector{grid, scanner_geometry};
        std::cout << "Building system matrix..." << std::endl;
        const auto matrix_start = Clock::now();
        const auto& A = projector.build_system_matrix();
        const double matrix_seconds = std::chrono::duration<double>(Clock::now() - matrix_start).count();
        std::cout << '\n';
        projector.dump(std::cout);

        // 6. Simulate measurements
        const auto ideal_b = projector.project(mu_true);
        NumReconstruction::MeasurementResult measurements{ideal_b};
        if (config.is_poisson_noise_enabled()) {
            measurements = NumReconstruction::simulate_poisson_noise(
                ideal_b, config.get_incident_photon_count(), config.get_random_seed());
        }
        const auto& b = measurements.b;
        std::cout << "Measurements: " << b.size() << '\n';

        // 7. Build the regularization matrix
        const auto L = NumReconstruction::build_regularization_matrix(
            grid, config.get_regularization_type());
        std::cout << "\nRegularization matrix: " << L.rows() << " x " << L.cols() << '\n'
                  << "Nonzero elements: " << L.nonZeros() << '\n';

        // 8. Reconstruct the attenuation image
        const auto output = config.get_output_directory();
        const auto iterations_directory = output / "iterations";
        const int save_step = config.get_iteration_save_step();
        std::vector<IterationRecord> history;
        double frame_seconds = 0.0;

        NumReconstruction::IterationCallback on_iteration;
        if (save_step > 0) {
            on_iteration = [&](const int iteration, const NumReconstruction::Vector& mu,
                               const double relative_residual) {
                const auto frame_start = Clock::now();
                const bool save_frame = iteration % save_step == 0;
                if (save_frame) {
                    NumReconstruction::save_image(iterations_directory / iteration_frame_name(iteration),
                        mu, grid.get_pixels_x(), grid.get_pixels_y(), mu_max);
                }
                const double rmse = (mu - mu_true).stableNorm()
                                  / std::sqrt(static_cast<double>(mu_true.size()));
                history.push_back({iteration, relative_residual, rmse, save_frame});
                frame_seconds += std::chrono::duration<double>(Clock::now() - frame_start).count();
            };
        }

        std::cout << "Solving the inverse problem..." << std::endl;
        const auto solve_start = Clock::now();
        const auto result = NumReconstruction::solve_tikhonov(
            A, b, L, config.get_lambda(), config.get_epsilon(), config.get_max_iterations(),
            on_iteration);
        // Frame output is excluded so the reported solver time stays comparable.
        const double solve_seconds =
            std::chrono::duration<double>(Clock::now() - solve_start).count() - frame_seconds;

        // 9. Evaluate the reconstruction
        const NumReconstruction::Vector image_error = result.mu - mu_true;
        const NumReconstruction::Vector data_error = A * result.mu - b;
        const NumReconstruction::Vector differences = L * result.mu;
        const double rmse = image_error.stableNorm() / std::sqrt(static_cast<double>(mu_true.size()));
        const double data_residual = data_error.stableNorm();
        const double regularization_norm = differences.stableNorm();
        const double weighted_regularization = config.get_lambda() == 0.0 ? 0.0
            : std::sqrt(config.get_lambda()) * regularization_norm;
        const double objective = data_residual * data_residual
                               + weighted_regularization * weighted_regularization;
        if (!std::isfinite(rmse) || !std::isfinite(regularization_norm) || !std::isfinite(objective)) {
            throw std::runtime_error("Nonfinite reconstruction metrics; check data scaling");
        }

        std::ostringstream summary;
        summary.imbue(std::locale::classic());
        summary << std::setprecision(std::numeric_limits<double>::max_digits10);
        result.dump(summary);
        summary << "RMSE against original: " << rmse << '\n'
                << "Data residual ||A*mu - b||: " << data_residual << '\n'
                << "Regularization norm ||L*mu||: " << regularization_norm << '\n'
                << "Objective ||A*mu - b||^2 + lambda*||L*mu||^2: " << objective << '\n'
                << "Reconstructed mu range: " << result.mu.minCoeff() << " .. " << result.mu.maxCoeff() << '\n'
                << "Zero photon counts replaced by 1 before logarithm: " << measurements.zero_photon_counts << '\n'
                << "Matrix construction time (s): " << matrix_seconds << '\n'
                << "Solver time (s): " << solve_seconds << '\n';
        std::cout << '\n' << summary.str();

        // 10. Save images, sinograms and the report
        NumReconstruction::save_image(output / "original.pgm", mu_true,
            grid.get_pixels_x(), grid.get_pixels_y(), mu_max);
        NumReconstruction::save_image(output / "reconstructed.pgm", result.mu,
            grid.get_pixels_x(), grid.get_pixels_y(), mu_max);
        NumReconstruction::save_matrix(output / "original.csv", mu_true,
            grid.get_pixels_y(), grid.get_pixels_x());
        NumReconstruction::save_matrix(output / "reconstructed.csv", result.mu,
            grid.get_pixels_y(), grid.get_pixels_x());
        NumReconstruction::save_sinogram(output / "sinogram.csv", b,
            config.get_num_angles(), config.get_num_detectors());
        NumReconstruction::save_sinogram(output / "ideal_sinogram.csv", ideal_b,
            config.get_num_angles(), config.get_num_detectors());

        // Iteration frames: history.csv lists every iteration of this run and names the
        // saved frames, so stale frames from earlier runs in the same directory are ignored.
        if (!history.empty()) {
            // Always keep the final iterate, so the animation ends on the saved result.
            if (!history.back().frame_saved) {
                NumReconstruction::save_image(iterations_directory
                    / iteration_frame_name(history.back().iteration),
                    result.mu, grid.get_pixels_x(), grid.get_pixels_y(), mu_max);
                history.back().frame_saved = true;
            }

            std::ostringstream csv;
            csv.imbue(std::locale::classic());
            csv << std::setprecision(std::numeric_limits<double>::max_digits10)
                << "iteration,relative_residual,rmse,frame\n";
            for (const auto& record : history) {
                csv << record.iteration << ',' << record.relative_residual << ',' << record.rmse << ','
                    << (record.frame_saved ? iteration_frame_name(record.iteration) : "") << '\n';
            }
            NumReconstruction::save_report(iterations_directory / "history.csv", csv.str());
        }

        std::ostringstream report;
        report.imbue(std::locale::classic());
        report << std::setprecision(std::numeric_limits<double>::max_digits10);
        config.dump(report);
        report << '\n';
        projector.dump(report);
        report << "Regularization matrix: " << L.rows() << " x " << L.cols()
               << ", nonzeros: " << L.nonZeros() << "\n\n" << summary.str()
               << "\nPGM images use the fixed mu range [0, 1]; values outside it are clipped for display.\n"
               << "CSV images and metrics preserve the raw solution, without clipping.\n"
               << "Image CSV: rows from top to bottom, columns from left to right.\n"
               << "Sinogram CSV: rows are views, columns are detectors; entries are logarithmic measurements.\n";
        if (!history.empty()) {
            report << "Iteration frames: iterations/iter_NNNN.pgm (same scale as the PGM images);\n"
                   << "iterations/history.csv lists every iteration with its residual, RMSE and frame.\n";
        }
        NumReconstruction::save_report(output / "report.txt", report.str());
        std::cout << "\nResults saved to: " << output << '\n';

        if (!result.converged) {
            std::cerr << "Iteration limit reached: the saved reconstruction has not met --epsilon.\n";
            return 2;
        }
    } catch (const std::invalid_argument& error) {
        std::cerr << "Error: " << error.what() << '\n'
                  << "Run with --help to see available options.\n";
        return 1;
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }

    return 0;
}
