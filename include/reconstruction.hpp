#ifndef RECONSTRUCTION_HPP
#define RECONSTRUCTION_HPP

#include "linalg.hpp"

#include <functional>
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
// Optional progress observer, called after every conjugate gradient iteration
// with the 1-based iteration number, the current iterate and its relative
// normal-equation residual. The iterate reference is valid only during the call.
using IterationCallback =
    std::function<void(int iteration, const Vector& mu, double relative_residual)>;

[[nodiscard]] ReconstructionResult solve_tikhonov(
    const SparseMatrix& A,
    const Vector& b,
    const SparseMatrix& L,
    double lambda,
    double epsilon,
    int max_iterations,
    const IterationCallback& on_iteration = {});
// ══════════════════════════════════════════════════════════════════

} // namespace NumReconstruction

#endif // RECONSTRUCTION_HPP
