#include "regularization.hpp"
#include "geometry.hpp"

#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace NumReconstruction {

static SparseMatrix build_identity(const Index pixel_count) {
    SparseMatrix matrix(pixel_count, pixel_count);
    matrix.setIdentity();
    matrix.makeCompressed();
    return matrix;
}

static SparseMatrix build_neighbor_differences(const Grid2D& grid) {
    const Index pixel_count = grid.pixel_count();
    const int nx = grid.get_pixels_x();
    const int ny = grid.get_pixels_y();
    const Index horizontal_pairs = pixel_count - ny;
    const Index vertical_pairs = pixel_count - nx;

    const Index max_pairs = std::numeric_limits<SparseMatrix::StorageIndex>::max() / 2;
    if (horizontal_pairs > max_pairs || vertical_pairs > max_pairs - horizontal_pairs) {
        throw std::overflow_error("Regularization: nonzero count exceeds sparse storage index range");
    }
    const Index pair_count = horizontal_pairs + vertical_pairs;
    std::vector<Triplet> entries;
    entries.reserve(static_cast<std::size_t>(pair_count) * 2);

    Index equation = 0;
    for (int row = 0; row < ny; ++row) {
        for (int col = 0; col < nx - 1; ++col) {
            entries.emplace_back(equation, grid.index(row, col), -1.0);
            entries.emplace_back(equation, grid.index(row, col + 1), 1.0);
            ++equation;
        }
    }
    for (int row = 0; row < ny - 1; ++row) {
        for (int col = 0; col < nx; ++col) {
            entries.emplace_back(equation, grid.index(row, col), -1.0);
            entries.emplace_back(equation, grid.index(row + 1, col), 1.0);
            ++equation;
        }
    }

    SparseMatrix matrix(pair_count, pixel_count);
    matrix.setFromTriplets(entries.begin(), entries.end());
    matrix.makeCompressed();
    return matrix;
}

SparseMatrix build_regularization_matrix(const Grid2D& grid, const Config::RegularizationType type) {
    if (grid.pixel_count() > std::numeric_limits<SparseMatrix::StorageIndex>::max()) {
        throw std::overflow_error("Regularization: pixel count exceeds sparse storage index range");
    }

    switch (type) {
        case Config::RegularizationType::ZeroOrder:
            return build_identity(grid.pixel_count());
        case Config::RegularizationType::FirstOrder:
            return build_neighbor_differences(grid);
        default:
            throw std::invalid_argument("Regularization: expected ZeroOrder or FirstOrder");
    }
}

} // namespace NumReconstruction
