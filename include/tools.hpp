#ifndef NUM_RECONSTRUCTION_TOOLS_HPP
#define NUM_RECONSTRUCTION_TOOLS_HPP

#include <cstddef>
#include <string_view>

namespace NumReconstruction::Tools {

// name includes the caller's context, e.g. "Grid2D: width".
// Invalid values throw std::invalid_argument; double positivity also requires
// finiteness, so NaN and infinity are rejected.
void require_positive(int value, std::string_view name);
void require_positive(double value, std::string_view name);
void require_finite(double value, std::string_view name);

// Validate [0,count); throw std::out_of_range with the supplied message.
void check_index(std::ptrdiff_t index, std::ptrdiff_t count, std::string_view message);

} // namespace NumReconstruction::Tools

#endif
