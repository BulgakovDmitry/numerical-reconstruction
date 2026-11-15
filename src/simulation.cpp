#include "simulation.hpp"
#include "tools.hpp"

#include <algorithm>
#include <cmath>
#include <random>
#include <stdexcept>

namespace NumReconstruction {

// ════ Measurement simulation ═════════════════════════════════════
MeasurementResult simulate_poisson_noise(const Vector& ideal_measurements,
                                         const double incident_photon_count,
                                         const std::uint32_t random_seed) {
    Tools::require_positive(incident_photon_count, "Simulation: incident_photon_count");
    // Keep integer samples and the distribution's intermediate arithmetic well within range.
    if (incident_photon_count > 1e12) {
        throw std::invalid_argument("Simulation: incident_photon_count must not exceed 1e12");
    }
    if (ideal_measurements.size() == 0 || !ideal_measurements.allFinite() ||
        (ideal_measurements.array() < 0.0).any()) {
        throw std::invalid_argument("Simulation: ideal measurements must be nonempty, finite and nonnegative");
    }

    std::mt19937_64 generator(random_seed);
    MeasurementResult result;
    result.b.resize(ideal_measurements.size());
    const double log_incident = std::log(incident_photon_count);

    for (Index i = 0; i < ideal_measurements.size(); ++i) {
        const double mean = incident_photon_count * std::exp(-ideal_measurements[i]);
        // Underflow to zero means no photons arrive; Poisson requires a positive mean.
        const auto count = mean == 0.0 ? std::int64_t{0}
                                      : std::poisson_distribution<std::int64_t>(mean)(generator);
        if (count == 0) ++result.zero_photon_counts;
        result.b[i] = log_incident - std::log(static_cast<double>(std::max(count, std::int64_t{1})));
    }

    return result;
}
// ══════════════════════════════════════════════════════════════════

} // namespace NumReconstruction
