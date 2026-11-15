#ifndef NUM_RECONSTRUCTION_GEOMETRY_HPP
#define NUM_RECONSTRUCTION_GEOMETRY_HPP

#include "linalg.hpp"

#include <iosfwd>
#include <optional>

namespace NumReconstruction {

class Config;

struct PixelIndex {
    int row;
    int col;
    bool operator==(const PixelIndex&) const = default;
};

struct Bounds2D {
    double xmin;
    double xmax;
    double ymin;
    double ymax;
};

struct RayInterval {
    double t_enter;
    double t_exit;
};

// Centered rectangle; x points right, y up, image rows run from top to bottom.
class Grid2D {
public:
    Grid2D(int pixels_x, int pixels_y, double width, double height);
    [[nodiscard]] static Grid2D from_config(const Config& config);
    void dump(std::ostream& out) const;

    [[nodiscard]] int get_pixels_x() const noexcept { return pixels_x_; }
    [[nodiscard]] int get_pixels_y() const noexcept { return pixels_y_; }
    [[nodiscard]] double get_width() const noexcept { return width_; }
    [[nodiscard]] double get_height() const noexcept { return height_; }
    [[nodiscard]] double get_hx() const noexcept { return hx_; }
    [[nodiscard]] double get_hy() const noexcept { return hy_; }
    [[nodiscard]] Index pixel_count() const noexcept { return pixel_count_; }
    [[nodiscard]] Bounds2D bounds() const noexcept;

    [[nodiscard]] Index index(int row, int col) const;
    [[nodiscard]] PixelIndex pixel(Index index) const;
    [[nodiscard]] Vec2 pixel_center(int row, int col) const;
    // Edge indices include the final edge: [0,pixels_x] / [0,pixels_y].
    [[nodiscard]] double x_edge(int col_edge) const;
    [[nodiscard]] double y_edge(int row_edge) const;

    // Cells own [left,right) x (bottom,top]. Nonfinite/outside points -> nullopt.
    [[nodiscard]] std::optional<PixelIndex> locate(const Vec2& point) const;

private:
    int pixels_x_;
    int pixels_y_;
    double width_;
    double height_;
    double hx_;
    double hy_;
    Index pixel_count_;
};

// Finite source-to-detector segment: r(t)=source+t*(detector-source), t in [0,1].
class Ray2D {
public:
    Ray2D(Vec2 source, Vec2 detector);
    void dump(std::ostream& out) const;
    // Include the segment's intersection with a given closed rectangular area.
    void dump(std::ostream& out, const Bounds2D& bounds) const;
    [[nodiscard]] const Vec2& get_source() const noexcept { return source_; }
    [[nodiscard]] const Vec2& get_detector() const noexcept { return detector_; }
    [[nodiscard]] Vec2 displacement() const noexcept { return detector_ - source_; }
    [[nodiscard]] double length() const noexcept { return length_; }
    [[nodiscard]] Vec2 point_at(double t) const;

private:
    Vec2 source_;
    Vec2 detector_;
    double length_;
};

// Closed-box geometric clipping. A point contact has zero length -> nullopt.
// Segments along box edges are retained; pixel ownership belongs to Grid2D.
// Invalid/nonfinite or degenerate bounds throw std::invalid_argument.
[[nodiscard]] std::optional<RayInterval>
clip_ray_to_box(const Ray2D& ray, const Bounds2D& bounds);

} // namespace NumReconstruction

#endif
