#ifndef NUM_RECONSTRUCTION_TOOLS_HPP
#define NUM_RECONSTRUCTION_TOOLS_HPP

#include <cstddef>
#include <string_view>

namespace NumReconstruction::Tools {

void require_positive(int value, std::string_view name);
void require_positive(double value, std::string_view name);
void require_finite(double value, std::string_view name);

void check_index(std::ptrdiff_t index, std::ptrdiff_t count, std::string_view message);

} // namespace NumReconstruction::Tools

#endif
