#include "reconstruction.hpp"
#include "tools.hpp"

#include <cmath>
#include <ostream>
#include <stdexcept>
#include <string>
#include <string_view>

namespace NumReconstruction {

static void validate_matrix(const SparseMatrix& matrix, std::string_view name);
static Vector build_diagonal(const SparseMatrix& A, const SparseMatrix& L, double lambda);
static Vector apply_normal_matrix(const SparseMatrix& A, const SparseMatrix& L,
                                  double lambda, const Vector& x);
static void require_positive_step(double value);


// ════ Tikhonov reconstruction ════════════════════════════════════
ReconstructionResult solve_tikhonov(const SparseMatrix& A,
                                     const Vector& b,
                                     const SparseMatrix& L,
                                     const double lambda,
                                     const double epsilon,
                                     const int max_iterations,
                                     const IterationCallback& on_iteration) {
    if (A.rows() <= 0 || A.cols() <= 0 || b.size() != A.rows() || L.cols() != A.cols()) {
        throw std::invalid_argument("Reconstruction: incompatible or empty matrix/vector dimensions");
    }
    Tools::require_finite(lambda, "Reconstruction: lambda");
    if (lambda < 0.0) {
        throw std::invalid_argument("Reconstruction: lambda must be nonnegative");
    }
    Tools::require_positive(epsilon, "Reconstruction: epsilon");
    Tools::require_positive(max_iterations, "Reconstruction: max_iterations");
    validate_matrix(A, "A");
    validate_matrix(L, "L");
    if (!b.allFinite()) {
        throw std::invalid_argument("Reconstruction: measurements must be finite");
    }

    const Vector rhs = A.transpose() * b;
    const double rhs_norm = rhs.stableNorm();
    if (!rhs.allFinite() || !std::isfinite(rhs_norm)) {
        throw std::runtime_error("Reconstruction: overflow while computing A^T*b");
    }

    ReconstructionResult result;
    result.mu = Vector::Zero(A.cols());
    result.relative_residual = rhs_norm == 0.0 ? 0.0 : 1.0;
    if (result.relative_residual <= epsilon) {
        result.converged = true;
        return result;
    }

    const Vector diagonal = build_diagonal(A, L, lambda);
    Vector residual = rhs;
    Vector preconditioned = residual.cwiseQuotient(diagonal);
    Vector direction = preconditioned;
    double rho = residual.dot(preconditioned);
    require_positive_step(rho);

    for (int iteration = 0; iteration < max_iterations; ++iteration) {
        const Vector action = apply_normal_matrix(A, L, lambda, direction);
        const double curvature = direction.dot(action);
        require_positive_step(curvature);
        const double alpha = rho / curvature;
        require_positive_step(alpha);

        result.mu += alpha * direction;
        residual -= alpha * action;
        result.iterations = iteration + 1;
        if (!result.mu.allFinite() || !residual.allFinite()) {
            throw std::runtime_error("Reconstruction: nonfinite conjugate gradient iterate");
        }

        result.relative_residual = residual.stableNorm() / rhs_norm;
        const bool last_iteration = result.iterations == max_iterations;
        if (result.relative_residual <= epsilon || last_iteration) {
            // Check the actual equation, since the recursively updated residual can drift.
            residual = rhs - apply_normal_matrix(A, L, lambda, result.mu);
            result.relative_residual = residual.stableNorm() / rhs_norm;
            if (!std::isfinite(result.relative_residual)) {
                throw std::runtime_error("Reconstruction: nonfinite final residual");
            }
            result.converged = result.relative_residual <= epsilon;
            if (on_iteration) {
                on_iteration(result.iterations, result.mu, result.relative_residual);
            }
            if (result.converged || last_iteration) {
                return result;
            }

            // Restart from the true residual if roundoff prevented convergence.
            preconditioned = residual.cwiseQuotient(diagonal);
            direction = preconditioned;
            rho = residual.dot(preconditioned);
            require_positive_step(rho);
            continue;
        }
        if (on_iteration) {
            on_iteration(result.iterations, result.mu, result.relative_residual);
        }

        preconditioned = residual.cwiseQuotient(diagonal);
        const double next_rho = residual.dot(preconditioned);
        require_positive_step(next_rho);
        const double beta = next_rho / rho;
        require_positive_step(beta);
        direction = preconditioned + beta * direction;
        rho = next_rho;
    }

    return result;
}
// ══════════════════════════════════════════════════════════════════


// ════ Dump ═══════════════════════════════════════════════════════
void ReconstructionResult::dump(std::ostream& out) const {
    out << "Reconstruction:\n"
        << "Method: conjugate gradients with diagonal preconditioning\n"
        << "Status: " << (converged ? "converged" : "iteration limit reached") << '\n'
        << "Iterations: " << iterations << '\n'
        << "Relative normal-equation residual: " << relative_residual << '\n';
}
// ══════════════════════════════════════════════════════════════════


// ════ Solver helpers ═════════════════════════════════════════════
static void validate_matrix(const SparseMatrix& matrix, const std::string_view name) {
    for (Index col = 0; col < matrix.outerSize(); ++col) {
        for (SparseMatrix::InnerIterator entry(matrix, col); entry; ++entry) {
            if (!std::isfinite(entry.value())) {
                throw std::invalid_argument("Reconstruction: " + std::string(name) + " must be finite");
            }
        }
    }
}

static Vector build_diagonal(const SparseMatrix& A, const SparseMatrix& L, const double lambda) {
    Vector diagonal = Vector::Zero(A.cols());
    const double weight = std::sqrt(lambda);
    for (Index col = 0; col < A.cols(); ++col) {
        for (SparseMatrix::InnerIterator entry(A, col); entry; ++entry) {
            diagonal[col] += entry.value() * entry.value();
        }
        if (lambda > 0.0) {
            for (SparseMatrix::InnerIterator entry(L, col); entry; ++entry) {
                const double value = weight * entry.value();
                diagonal[col] += value * value;
            }
        }
        // An unobserved, unpenalized variable has zero diagonal and zero RHS.
        if (diagonal[col] == 0.0) diagonal[col] = 1.0;
    }
    if (!diagonal.allFinite()) {
        throw std::runtime_error("Reconstruction: overflow in the diagonal preconditioner");
    }
    return diagonal;
}

static Vector apply_normal_matrix(const SparseMatrix& A, const SparseMatrix& L,
                                  const double lambda, const Vector& x) {
    const Vector projected = A * x;
    Vector result = A.transpose() * projected;
    if (lambda > 0.0 && L.rows() > 0) {
        const Vector differences = L * x;
        result += lambda * (L.transpose() * differences);
    }
    if (!result.allFinite()) {
        throw std::runtime_error("Reconstruction: overflow applying the normal-equation matrix");
    }
    return result;
}

static void require_positive_step(const double value) {
    if (!std::isfinite(value) || value <= 0.0) {
        throw std::runtime_error(
            "Reconstruction: numerical breakdown in conjugate gradients; "
            "check data scaling, regularization and solver tolerance");
    }
}
// ══════════════════════════════════════════════════════════════════

} // namespace NumReconstruction
