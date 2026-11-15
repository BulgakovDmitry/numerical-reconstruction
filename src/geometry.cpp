#include "geometry.hpp"
#include "config.hpp"
#include "tools.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <ostream>
#include <stdexcept>

namespace NumReconstruction {
namespace {

Index checked_pixel_count(int nx, int ny) {
    Tools::require_positive(nx, "Grid2D: pixels_x");
    Tools::require_positive(ny, "Grid2D: pixels_y");
    if (static_cast<Index>(nx) > std::numeric_limits<Index>::max() / ny) {
        throw std::overflow_error("Grid2D: pixel count exceeds Eigen::Index");
    }
    return static_cast<Index>(nx) * ny;
}

} // namespace

Grid2D::Grid2D(int pixels_x, int pixels_y, double width, double height)
    : pixels_x_(pixels_x), pixels_y_(pixels_y), width_(width), height_(height),
      hx_(0.0), hy_(0.0), pixel_count_(checked_pixel_count(pixels_x, pixels_y)) {
    Tools::require_positive(width_, "Grid2D: width");
    Tools::require_positive(height_, "Grid2D: height");
    hx_ = width_ / pixels_x_;
    hy_ = height_ / pixels_y_;
    Tools::require_positive(hx_, "Grid2D: x spacing");
    Tools::require_positive(hy_, "Grid2D: y spacing");
    if (hx_ / 2.0 == 0.0 || hy_ / 2.0 == 0.0) {
        throw std::invalid_argument("Grid2D: pixel half spacing is not representable");
    }
}

Grid2D Grid2D::from_config(const Config& config) {
    return {config.get_pixels_x(), config.get_pixels_y(),
            config.get_image_width(), config.get_image_height()};
}

void Grid2D::dump(std::ostream& out) const {
    const auto box = bounds();
    out << "Grid:\n"
        << "Pixels: " << pixels_x_ << " x " << pixels_y_ << '\n'
        << "Image size: " << width_ << " x " << height_ << " m\n"
        << "Pixel count: " << pixel_count_ << '\n'
        << "Pixel spacing (x, y): " << hx_ << ", " << hy_ << " m\n"
        << "Grid bounds: x=[" << box.xmin << ", " << box.xmax
        << "], y=[" << box.ymin << ", " << box.ymax << "] m\n";
}

Bounds2D Grid2D::bounds() const noexcept {
    return {-width_ / 2.0, width_ / 2.0, -height_ / 2.0, height_ / 2.0};
}

Index Grid2D::index(int row, int col) const {
    Tools::check_index(row, pixels_y_, "Grid2D: pixel coordinates out of range");
    Tools::check_index(col, pixels_x_, "Grid2D: pixel coordinates out of range");
    return static_cast<Index>(row) * pixels_x_ + col;
}

PixelIndex Grid2D::pixel(Index value) const {
    Tools::check_index(value, pixel_count_, "Grid2D: pixel index out of range");
    return {static_cast<int>(value / pixels_x_), static_cast<int>(value % pixels_x_)};
}

Vec2 Grid2D::pixel_center(int row, int col) const {
    (void)index(row, col);
    return {((static_cast<double>(col) + 0.5) / pixels_x_ - 0.5) * width_,
            (0.5 - (static_cast<double>(row) + 0.5) / pixels_y_) * height_};
}

double Grid2D::x_edge(int col_edge) const {
    if (col_edge < 0 || col_edge > pixels_x_) {
        throw std::out_of_range("Grid2D: x edge index out of range");
    }
    return (static_cast<double>(col_edge) / pixels_x_ - 0.5) * width_;
}

double Grid2D::y_edge(int row_edge) const {
    if (row_edge < 0 || row_edge > pixels_y_) {
        throw std::out_of_range("Grid2D: y edge index out of range");
    }
    return (0.5 - static_cast<double>(row_edge) / pixels_y_) * height_;
}

std::optional<PixelIndex> Grid2D::locate(const Vec2& point) const {
    const auto box = bounds();
    if (!point.allFinite() || point.x() < box.xmin || point.x() >= box.xmax ||
        point.y() <= box.ymin || point.y() > box.ymax) {
        return std::nullopt;
    }
    // Search actual edge coordinates so an exact edge has deterministic ownership,
    // even when division-based conversion to an index would round to a neighbour.
    int left = 0;
    int right = pixels_x_;
    while (right - left > 1) {
        const int mid = left + (right - left) / 2;
        if (point.x() >= x_edge(mid)) left = mid;
        else right = mid;
    }
    const int col = left;
    left = 0;
    right = pixels_y_;
    while (right - left > 1) {
        const int mid = left + (right - left) / 2;
        if (point.y() <= y_edge(mid)) left = mid;
        else right = mid;
    }
    return PixelIndex{left, col};
}

Ray2D::Ray2D(Vec2 source, Vec2 detector)
    : source_(source), detector_(detector), length_(0.0) {
    if (!source_.allFinite() || !detector_.allFinite()) {
        throw std::invalid_argument("Ray2D: endpoints must be finite");
    }
    const Vec2 delta = displacement();
    length_ = std::hypot(delta.x(), delta.y());
    if (!delta.allFinite() || !std::isfinite(length_) || length_ <= 0.0) {
        throw std::invalid_argument("Ray2D: segment must have finite positive length");
    }
}

void Ray2D::dump(std::ostream& out) const {
    out << "Ray:\n"
        << "Source (x, y): " << source_.transpose() << " m\n"
        << "Detector (x, y): " << detector_.transpose() << " m\n"
        << "Source-detector length: " << length_ << " m\n";
}

void Ray2D::dump(std::ostream& out, const Bounds2D& box) const {
    const auto interval = clip_ray_to_box(*this, box);
    dump(out);
    if (interval) {
        out << "Entry point (x, y): " << point_at(interval->t_enter).transpose() << " m\n"
            << "Exit point (x, y): " << point_at(interval->t_exit).transpose() << " m\n"
            << "Path length inside grid: "
            << (interval->t_exit - interval->t_enter) * length_ << " m\n";
    } else {
        out << "Ray does not cross the grid.\n";
    }
}

Vec2 Ray2D::point_at(double t) const {
    if (!std::isfinite(t) || t < 0.0 || t > 1.0) {
        throw std::out_of_range("Ray2D: t must be finite and in [0,1]");
    }
    if (t == 0.0) return source_;
    if (t == 1.0) return detector_;
    return source_ + t * displacement();
}

std::optional<RayInterval> clip_ray_to_box(const Ray2D& ray, const Bounds2D& box) {
    if (!std::isfinite(box.xmin) || !std::isfinite(box.xmax) ||
        !std::isfinite(box.ymin) || !std::isfinite(box.ymax) ||
        box.xmin >= box.xmax || box.ymin >= box.ymax) {
        throw std::invalid_argument("clip_ray_to_box: invalid bounds");
    }
    // Extended intermediates avoid overflow when subtracting large coordinates
    // or dividing by a tiny nonzero direction component. No solver epsilon here.
    long double enter = 0.0L;
    long double exit = 1.0L;
    const Vec2 delta = ray.displacement();
    for (int axis = 0; axis < 2; ++axis) {
        const long double origin = ray.get_source()[axis];
        const long double direction = delta[axis];
        const long double lower = axis == 0 ? box.xmin : box.ymin;
        const long double upper = axis == 0 ? box.xmax : box.ymax;
        if (direction == 0.0L) {
            if (origin < lower || origin > upper) return std::nullopt;
            continue;
        }
        long double t0 = (lower - origin) / direction;
        long double t1 = (upper - origin) / direction;
        if (t0 > t1) std::swap(t0, t1);
        enter = std::max(enter, t0);
        exit = std::min(exit, t1);
        if (enter >= exit) return std::nullopt;
    }
    const RayInterval result{static_cast<double>(enter), static_cast<double>(exit)};
    if (result.t_enter >= result.t_exit) return std::nullopt;
    return result;
}

} // namespace NumReconstruction
