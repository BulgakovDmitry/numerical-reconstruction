#include "io.hpp"
#include "geometry.hpp"

#include <algorithm>
#include <charconv>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <vector>

namespace NumReconstruction {

static Index checked_size(const int rows, const int cols, const std::string_view operation) {
    if (rows <= 0 || cols <= 0) {
        throw std::invalid_argument(std::string(operation) + ": dimensions must be greater than zero");
    }
    if (static_cast<Index>(rows) > std::numeric_limits<Index>::max() / cols) {
        throw std::overflow_error(std::string(operation) + ": dimensions exceed the index range");
    }
    return static_cast<Index>(rows) * cols;
}

static void check_scale(const double mu_max, const std::string_view operation) {
    if (!std::isfinite(mu_max) || mu_max <= 0.0) {
        throw std::invalid_argument(std::string(operation) + ": mu_max must be finite and greater than zero");
    }
}

static void skip_pgm_whitespace(std::istream& file) {
    while (true) {
        const auto next = file.peek();
        if (next == std::char_traits<char>::eof()) return;
        if (std::isspace(static_cast<unsigned char>(next))) {
            file.get();
        } else if (next == '#') {
            while (true) {
                const auto ch = file.get();
                if (ch == '\n' || ch == '\r' || ch == std::char_traits<char>::eof()) break;
            }
        } else {
            return;
        }
    }
}

// Leaves the delimiter in the stream so the P5 raster separator can be consumed
// exactly once: a pixel byte may itself be whitespace or '#'.
static std::string read_pgm_token(std::istream& file, const std::string& context) {
    skip_pgm_whitespace(file);
    std::string token;
    while (true) {
        const auto next = file.peek();
        if (next == std::char_traits<char>::eof() || next == '#' ||
            std::isspace(static_cast<unsigned char>(next))) break;
        if (token.size() == 64) {
            throw std::runtime_error(context + ": invalid or excessively long PGM token");
        }
        token.push_back(static_cast<char>(file.get()));
    }
    if (token.empty()) {
        throw std::runtime_error(context + ": unexpected end of PGM data");
    }
    return token;
}

static unsigned int read_pgm_number(std::istream& file, const std::string& context) {
    const auto token = read_pgm_token(file, context);
    unsigned int value{};
    const auto [end, error] = std::from_chars(token.data(), token.data() + token.size(), value);
    if (error != std::errc{} || end != token.data() + token.size()) {
        throw std::runtime_error(context + ": expected a nonnegative integer, got '" + token + "'");
    }
    return value;
}

static std::ofstream open_output(const std::filesystem::path& path, const std::string_view operation) {
    if (path.has_parent_path()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::ofstream file(path, std::ios::binary);
    if (!file) {
        throw std::runtime_error(std::string(operation) + ": cannot open " + path.string());
    }
    file.imbue(std::locale::classic());
    return file;
}

static void close_output(std::ofstream& file, const std::filesystem::path& path,
                         const std::string_view operation) {
    file.close();
    if (!file) {
        throw std::runtime_error(std::string(operation) + ": failed to write " + path.string());
    }
}

Vector load_image(const std::filesystem::path& path, const Grid2D& grid, const double mu_max) {
    check_scale(mu_max, "load_image");
    const std::string context = "load_image: " + path.string();
    std::ifstream file(path, std::ios::binary);
    if (!file) {
        throw std::runtime_error(context + ": cannot open file");
    }

    const auto magic = read_pgm_token(file, context);
    if (magic != "P2" && magic != "P5") {
        throw std::runtime_error(context + ": expected a PGM image (P2 or P5)");
    }
    const auto width = read_pgm_number(file, context);
    const auto height = read_pgm_number(file, context);
    const auto maxval = read_pgm_number(file, context);
    if (width != static_cast<unsigned int>(grid.get_pixels_x()) ||
        height != static_cast<unsigned int>(grid.get_pixels_y())) {
        throw std::invalid_argument(context + ": image size " + std::to_string(width) + "x"
                                    + std::to_string(height) + " does not match grid "
                                    + std::to_string(grid.get_pixels_x()) + "x"
                                    + std::to_string(grid.get_pixels_y()));
    }
    if (maxval == 0 || maxval > 65535) {
        throw std::runtime_error(context + ": PGM Maxval must be in [1, 65535]");
    }

    if (magic == "P5") {
        const auto separator = file.get();
        if (separator == std::char_traits<char>::eof() ||
            !std::isspace(static_cast<unsigned char>(separator))) {
            throw std::runtime_error(context + ": missing whitespace before PGM raster");
        }
        // Check for truncation before allocating the image vector.
        const auto raster_start = file.tellg();
        file.seekg(0, std::ios::end);
        const auto raster_end = file.tellg();
        const int bytes_per_pixel = maxval < 256 ? 1 : 2;
        if (raster_start == std::streampos(-1) || raster_end == std::streampos(-1) ||
            (raster_end - raster_start) / bytes_per_pixel < grid.pixel_count()) {
            throw std::runtime_error(context + ": truncated PGM raster");
        }
        file.seekg(raster_start);
    }

    Vector mu(grid.pixel_count());
    for (Index i = 0; i < mu.size(); ++i) {
        unsigned int sample{};
        if (magic == "P2") {
            sample = read_pgm_number(file, context);
        } else {
            const auto first = file.get();
            if (first == std::char_traits<char>::eof()) {
                throw std::runtime_error(context + ": truncated PGM raster");
            }
            sample = static_cast<unsigned int>(first);
            if (maxval >= 256) {
                const auto second = file.get();
                if (second == std::char_traits<char>::eof()) {
                    throw std::runtime_error(context + ": truncated PGM raster");
                }
                sample = sample * 256 + static_cast<unsigned int>(second);
            }
        }
        if (sample > maxval) {
            throw std::runtime_error(context + ": pixel " + std::to_string(i) + " exceeds PGM Maxval");
        }
        mu[i] = (static_cast<double>(sample) / maxval) * mu_max;
    }
    return mu;
}

void save_image(const std::filesystem::path& path,
                const Vector& mu,
                const int width,
                const int height,
                const double mu_max) {
    const Index pixel_count = checked_size(height, width, "save_image");
    if (mu.size() != pixel_count) {
        throw std::invalid_argument("save_image: vector size " + std::to_string(mu.size())
                                    + " does not match " + std::to_string(width) + "x"
                                    + std::to_string(height));
    }

    check_scale(mu_max, "save_image");
    if (static_cast<std::uintmax_t>(pixel_count) >
        static_cast<std::uintmax_t>(std::numeric_limits<std::streamsize>::max())) {
        throw std::overflow_error("save_image: image exceeds the stream size range");
    }

    // Map attenuation coefficients to display levels, with a fixed physical scale.
    std::vector<std::uint8_t> pixels(static_cast<std::size_t>(pixel_count));
    for (Index i = 0; i < pixel_count; ++i) {
        const double value = mu[i];
        if (!std::isfinite(value)) {
            throw std::invalid_argument("save_image: mu[" + std::to_string(i) + "] is NaN or Inf");
        }
        const double level = (std::clamp(value, 0.0, mu_max) / mu_max) * 255.0;
        pixels[static_cast<std::size_t>(i)] = static_cast<std::uint8_t>(std::lround(level));
    }

    auto file = open_output(path, "save_image");

    // PGM P5 header, then rows top to bottom - exactly the storage order of the vector.
    file << "P5\n" << width << ' ' << height << "\n255\n";
    file.write(reinterpret_cast<const char*>(pixels.data()),
               static_cast<std::streamsize>(pixels.size()));
    close_output(file, path, "save_image");
}

void save_sinogram(const std::filesystem::path& path, const Vector& measurements,
                   const int num_angles, const int num_detectors) {
    const Index count = checked_size(num_angles, num_detectors, "save_sinogram");
    if (measurements.size() != count) {
        throw std::invalid_argument("save_sinogram: vector size does not match the scanner dimensions");
    }
    if (!measurements.allFinite()) {
        throw std::invalid_argument("save_sinogram: measurements must be finite");
    }

    auto file = open_output(path, "save_sinogram");
    file << std::setprecision(std::numeric_limits<double>::max_digits10);
    for (int view = 0; view < num_angles; ++view) {
        for (int detector = 0; detector < num_detectors; ++detector) {
            if (detector != 0) file << ',';
            file << measurements[static_cast<Index>(view) * num_detectors + detector];
        }
        file << '\n';
    }
    close_output(file, path, "save_sinogram");
}

} // namespace NumReconstruction
