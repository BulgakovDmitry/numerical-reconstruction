#ifndef NUM_RECONSTRUCTION_TEST_SUPPORT_HPP
#define NUM_RECONSTRUCTION_TEST_SUPPORT_HPP

#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>

namespace Test {

inline void check(bool condition, const char* expression, int line) {
    if (!condition) {
        throw std::runtime_error(std::string("line ") + std::to_string(line) + ": " + expression);
    }
}

inline void near(double actual, double expected, double tolerance = 1e-12) {
    if (!std::isfinite(actual) || !std::isfinite(expected) ||
        std::abs(actual - expected) > tolerance * std::max(1.0, std::abs(expected))) {
        throw std::runtime_error("expected " + std::to_string(expected) +
                                 ", got " + std::to_string(actual));
    }
}

template<class Exception, class Function>
void throws(Function&& function) {
    try {
        std::forward<Function>(function)();
    } catch (const Exception&) {
        return;
    }
    throw std::runtime_error("expected exception was not thrown");
}

class Suite {
public:
    template<class Function>
    void run(const char* name, Function&& function) {
        ++total_;
        try {
            std::forward<Function>(function)();
            std::cout << "PASS " << name << '\n';
        } catch (const std::exception& error) {
            ++failed_;
            std::cerr << "FAIL " << name << ": " << error.what() << '\n';
        }
    }

    int result() const {
        std::cout << total_ - failed_ << '/' << total_ << " tests passed\n";
        return failed_ == 0 ? 0 : 1;
    }

private:
    int total_ = 0;
    int failed_ = 0;
};

} // namespace Test

// Variadic so expressions containing braced row/col pairs remain readable.
#define CHECK(...) Test::check(static_cast<bool>((__VA_ARGS__)), #__VA_ARGS__, __LINE__)

#endif
