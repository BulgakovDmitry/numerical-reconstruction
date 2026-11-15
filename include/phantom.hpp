#ifndef PHANTOM_HPP
#define PHANTOM_HPP

#include "linalg.hpp"

namespace NumReconstruction {

class Grid2D;

// ════ Artificial image ═══════════════════════════════════════════
[[nodiscard]] Vector create_suitcase_phantom(const Grid2D& grid, double mu_max = 1.0);
// ══════════════════════════════════════════════════════════════════

} // namespace NumReconstruction

#endif // PHANTOM_HPP
