#ifndef RECONSTRUCTION_HPP
#define RECONSTRUCTION_HPP

#include "linalg.hpp"

#include <iosfwd>

namespace NumReconstruction {

// ════ Reconstruction result ═════════════════════════════════════
struct ReconstructionResult {
    Vector mu;
    int iterations = 0;
    double relative_residual = 0.0;
    bool converged = false;

    void dump(std::ostream& out) const;
};
// ══════════════════════════════════════════════════════════════════


// ════ Tikhonov reconstruction ════════════════════════════════════
[[nodiscard]] ReconstructionResult solve_tikhonov(
    const SparseMatrix& A,
    const Vector& b,
    const SparseMatrix& L,
    double lambda,
    double epsilon,
    int max_iterations);
// ══════════════════════════════════════════════════════════════════

} // namespace NumReconstruction

#endif // RECONSTRUCTION_HPP
