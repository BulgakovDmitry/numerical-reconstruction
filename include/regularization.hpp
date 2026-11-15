#ifndef NUM_RECONSTRUCTION_REGULARIZATION_HPP
#define NUM_RECONSTRUCTION_REGULARIZATION_HPP

#include "config.hpp"
#include "linalg.hpp"

namespace NumReconstruction {

class Grid2D;

[[nodiscard]] SparseMatrix build_regularization_matrix(
    const Grid2D& grid, Config::RegularizationType type);

} // namespace NumReconstruction

#endif // NUM_RECONSTRUCTION_REGULARIZATION_HPP
