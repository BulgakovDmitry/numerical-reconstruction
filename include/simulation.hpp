#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include "linalg.hpp"

#include <cstdint>

namespace NumReconstruction {

// ════ Measurement simulation ═════════════════════════════════════
struct MeasurementResult {
    Vector b;
    Index zero_photon_counts = 0;
};

[[nodiscard]] MeasurementResult simulate_poisson_noise(
    const Vector& ideal_measurements,
    double incident_photon_count,
    std::uint32_t random_seed);
// ══════════════════════════════════════════════════════════════════

} // namespace NumReconstruction

#endif // SIMULATION_HPP
