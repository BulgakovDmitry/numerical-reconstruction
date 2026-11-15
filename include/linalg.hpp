#ifndef NUM_RECONSTRUCTION_LINALG_HPP
#define NUM_RECONSTRUCTION_LINALG_HPP

#include <Eigen/Core>
#include <Eigen/SparseCore>

namespace NumReconstruction {

using Vec2 = Eigen::Vector2d;
using Vector = Eigen::VectorXd;
using SparseMatrix = Eigen::SparseMatrix<double>;
using Triplet = Eigen::Triplet<double>;
using Index = Eigen::Index;

} // namespace NumReconstruction

#endif
