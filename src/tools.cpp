#include "tools.hpp"

#include <cmath>
#include <stdexcept>
#include <string>

namespace NumReconstruction::Tools {

void require_positive(int value, std::string_view name) {
    if (value <= 0) {
        throw std::invalid_argument(std::string(name) + " must be greater than zero");
    }
}

void require_finite(double value, std::string_view name) {
    if (!std::isfinite(value)) {
        throw std::invalid_argument(std::string(name) + " must be finite");
    }
}

void require_positive(double value, std::string_view name) {
    require_finite(value, name);
    if (value <= 0.0) {
        throw std::invalid_argument(std::string(name) + " must be greater than zero");
    }
}

void check_index(std::ptrdiff_t index, std::ptrdiff_t count, std::string_view message) {
    if (index < 0 || index >= count) {
        throw std::out_of_range(std::string(message));
    }
}

} // namespace NumReconstruction::Tools
