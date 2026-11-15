#include "geometry.hpp"
#include "config.hpp"
#include "test_support.hpp"

#include <limits>
#include <vector>

using namespace NumReconstruction;
using Test::near;
using Test::throws;

namespace {

void expect_clip(const Ray2D& ray, const Bounds2D& box,
                 double enter, double exit, double length) {
    const auto interval = clip_ray_to_box(ray, box);
    CHECK(interval.has_value());
    near(interval->t_enter, enter);
    near(interval->t_exit, exit);
    near((interval->t_exit - interval->t_enter) * ray.length(), length);
}

} // namespace

int main() {
    Test::Suite suite;
    suite.run("rectangular grid dimensions and coordinates", [] {
        const Grid2D grid{3, 2, 6.0, 2.0};
        CHECK(grid.get_pixels_x() == 3 && grid.get_pixels_y() == 2);
        CHECK(grid.pixel_count() == 6);
        near(grid.get_width(), 6.0);
        near(grid.get_height(), 2.0);
        near(grid.get_hx(), 2.0);
        near(grid.get_hy(), 1.0);
        const auto box = grid.bounds();
        near(box.xmin, -3.0); near(box.xmax, 3.0);
        near(box.ymin, -1.0); near(box.ymax, 1.0);
        near(grid.x_edge(0), -3.0); near(grid.x_edge(1), -1.0);
        near(grid.x_edge(3), 3.0);
        near(grid.y_edge(0), 1.0); near(grid.y_edge(1), 0.0);
        near(grid.y_edge(2), -1.0);
        near(grid.pixel_center(0, 0).x(), -2.0);
        near(grid.pixel_center(0, 0).y(), 0.5);
        near(grid.pixel_center(1, 2).x(), 2.0);
        near(grid.pixel_center(1, 2).y(), -0.5);
    });
    suite.run("row-major indexing and center lookup round trips", [] {
        const Grid2D grid{17, 11, 3.0, 2.0};
        CHECK(grid.index(1, 0) == 17);
        CHECK(grid.pixel(186) == PixelIndex{10, 16});
        for (int row = 0; row < 11; ++row) {
            for (int col = 0; col < 17; ++col) {
                CHECK(grid.pixel(grid.index(row, col)) == PixelIndex{row, col});
                CHECK(grid.locate(grid.pixel_center(row, col)) == PixelIndex{row, col});
            }
        }
    });
    suite.run("half-open pixel boundaries and immediately adjacent points", [] {
        const Grid2D grid{3, 2, 6.0, 2.0};
        CHECK(grid.locate(Vec2{-3.0, 1.0}) == PixelIndex{0, 0});
        CHECK(grid.locate(Vec2{-1.0, 0.0}) == PixelIndex{1, 1});
        CHECK(!grid.locate(Vec2{3.0, 0.5}));
        CHECK(!grid.locate(Vec2{0.0, -1.0}));
        CHECK(!grid.locate(Vec2{-3.1, 0.0}));
        CHECK(!grid.locate(Vec2{0.0, 1.1}));
        CHECK(grid.locate(Vec2{std::nextafter(-1.0, -2.0), 0.5}) == PixelIndex{0, 0});
        CHECK(grid.locate(Vec2{std::nextafter(-1.0, 0.0), 0.5}) == PixelIndex{0, 1});
        CHECK(grid.locate(Vec2{std::nextafter(3.0, 0.0), 0.5}) == PixelIndex{0, 2});
        CHECK(grid.locate(Vec2{0.0, std::nextafter(-1.0, 0.0)}) == PixelIndex{1, 1});
        const double nan = std::numeric_limits<double>::quiet_NaN();
        CHECK(!grid.locate(Vec2{nan, 0.0}));
        CHECK(!grid.locate(Vec2{0.0, std::numeric_limits<double>::infinity()}));
        const Grid2D uneven{7, 5, 2.0, 3.0};
        for (int r = 0; r < 5; ++r)
            for (int c = 0; c < 7; ++c)
                CHECK(uneven.locate(Vec2{uneven.x_edge(c), uneven.y_edge(r)}) == PixelIndex{r,c});
    });
    suite.run("single-pixel grid", [] {
        const Grid2D grid{1, 1, 4.0, 2.0};
        CHECK(grid.pixel_count() == 1);
        CHECK(grid.pixel_center(0, 0).isZero());
        CHECK(grid.locate(Vec2{0.0, 0.0}) == PixelIndex{0, 0});
    });
    suite.run("grid rejects invalid dimensions and indexing", [] {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        const double inf = std::numeric_limits<double>::infinity();
        throws<std::invalid_argument>([] { Grid2D{0, 2, 2.0, 2.0}; });
        throws<std::invalid_argument>([] { Grid2D{2, -1, 2.0, 2.0}; });
        for (double bad : {0.0, -1.0, nan, inf}) {
            throws<std::invalid_argument>([&] { Grid2D{2, 2, bad, 2.0}; });
            throws<std::invalid_argument>([&] { Grid2D{2, 2, 2.0, bad}; });
        }
        throws<std::invalid_argument>([] {
            Grid2D{2, 2, std::numeric_limits<double>::denorm_min(), 2.0};
        });
        throws<std::invalid_argument>([] {
            Grid2D{2, 2, 2.0 * std::numeric_limits<double>::denorm_min(), 2.0};
        });
        const Grid2D grid{3, 2, 6.0, 2.0};
        throws<std::out_of_range>([&] { (void)grid.index(-1, 0); });
        throws<std::out_of_range>([&] { (void)grid.index(2, 0); });
        throws<std::out_of_range>([&] { (void)grid.index(0, 3); });
        throws<std::out_of_range>([&] { (void)grid.pixel(-1); });
        throws<std::out_of_range>([&] { (void)grid.pixel(6); });
        throws<std::out_of_range>([&] { (void)grid.pixel_center(0, -1); });
        throws<std::out_of_range>([&] { (void)grid.x_edge(-1); });
        throws<std::out_of_range>([&] { (void)grid.x_edge(4); });
        throws<std::out_of_range>([&] { (void)grid.y_edge(-1); });
        throws<std::out_of_range>([&] { (void)grid.y_edge(3); });
    });
    suite.run("grid Config adapter", [] {
        const auto config = Config::Builder{}.set_pixels_x(3).set_pixels_y(2)
            .set_image_width(6.0).set_image_height(2.0).build();
        const auto grid = Grid2D::from_config(config);
        CHECK(grid.pixel_count() == 6);
        near(grid.get_hx(), 2.0); near(grid.get_hy(), 1.0);
    });
    suite.run("segment parameterization and robust length", [] {
        const Ray2D ray{Vec2{1.0, 2.0}, Vec2{4.0, 6.0}};
        near(ray.length(), 5.0);
        CHECK(ray.displacement() == Vec2(3.0, 4.0));
        CHECK(ray.point_at(0.0) == ray.get_source());
        CHECK(ray.point_at(1.0) == ray.get_detector());
        CHECK(ray.point_at(0.5) == Vec2(2.5, 4.0));
        // Squaring the coordinates would overflow/underflow; hypot must not.
        near(Ray2D{Vec2{0.0, 0.0}, Vec2{3e200, 4e200}}.length() / 1e200, 5.0);
        near(Ray2D{Vec2{0.0, 0.0}, Vec2{3e-200, 4e-200}}.length() / 1e-200, 5.0);
    });
    suite.run("segment validation", [] {
        throws<std::invalid_argument>([] { Ray2D{Vec2{0.0,0.0}, Vec2{0.0,0.0}}; });
        throws<std::invalid_argument>([] {
            Ray2D{Vec2{std::numeric_limits<double>::infinity(),0.0}, Vec2{1.0,0.0}};
        });
        throws<std::invalid_argument>([] {
            Ray2D{Vec2{0.0,0.0}, Vec2{0.0,std::numeric_limits<double>::quiet_NaN()}};
        });
        throws<std::invalid_argument>([] { Ray2D{Vec2{-1e308,0.0}, Vec2{1e308,0.0}}; });
        const Ray2D ray{Vec2{0.0,0.0}, Vec2{1.0,0.0}};
        for (double bad : {-0.1, 1.1, std::numeric_limits<double>::quiet_NaN()})
            throws<std::out_of_range>([&] { (void)ray.point_at(bad); });
    });
    suite.run("analytical horizontal vertical diagonal clipping", [] {
        const Bounds2D box{-1.0,1.0,-1.0,1.0};
        expect_clip(Ray2D{Vec2{-2.0,0.0},Vec2{2.0,0.0}}, box, 0.25, 0.75, 2.0);
        expect_clip(Ray2D{Vec2{0.0,2.0},Vec2{0.0,-2.0}}, box, 0.25, 0.75, 2.0);
        expect_clip(Ray2D{Vec2{-2.0,-2.0},Vec2{2.0,2.0}}, box,
                    0.25, 0.75, 2.0*std::sqrt(2.0));
        expect_clip(Ray2D{Vec2{0.0,0.0},Vec2{0.5,0.0}}, box, 0.0, 1.0, 0.5);
        expect_clip(Ray2D{Vec2{0.0,0.0},Vec2{2.0,0.0}}, box, 0.0, 0.5, 1.0);
    });
    suite.run("clipping respects finite endpoints and point tangencies", [] {
        const Bounds2D box{-1.0,1.0,-1.0,1.0};
        CHECK(!clip_ray_to_box(Ray2D{Vec2{-3.0,0.0},Vec2{-2.0,0.0}}, box));
        CHECK(!clip_ray_to_box(Ray2D{Vec2{-2.0,2.0},Vec2{2.0,2.0}}, box));
        CHECK(!clip_ray_to_box(Ray2D{Vec2{2.0,-2.0},Vec2{2.0,2.0}}, box));
        CHECK(!clip_ray_to_box(Ray2D{Vec2{-2.0,0.0},Vec2{0.0,2.0}}, box));
        CHECK(!clip_ray_to_box(Ray2D{Vec2{1.0,0.0},Vec2{2.0,0.0}}, box));
        expect_clip(Ray2D{Vec2{-2.0,1.0},Vec2{2.0,1.0}}, box, 0.25, 0.75, 2.0);
        expect_clip(Ray2D{Vec2{1.0,-2.0},Vec2{1.0,2.0}}, box, 0.25, 0.75, 2.0);
    });
    suite.run("clipping reversal and scale invariance", [] {
        for (double scale : {1e-150, 1.0, 1e150}) {
            const Bounds2D box{-scale,scale,-scale,scale};
            const Ray2D forward{Vec2{-3.0*scale,0.0},Vec2{5.0*scale,0.0}};
            const Ray2D reverse{forward.get_detector(),forward.get_source()};
            const auto a = clip_ray_to_box(forward,box);
            const auto b = clip_ray_to_box(reverse,box);
            CHECK(a && b);
            near(a->t_enter, 0.25); near(a->t_exit, 0.5);
            near(b->t_enter, 1.0-a->t_exit); near(b->t_exit, 1.0-a->t_enter);
        }
        // Tiny nonzero direction components must not be treated as zero.
        const auto tiny = clip_ray_to_box(
            Ray2D{Vec2{-2.0,-1e-16},Vec2{2.0,1e-16}}, {-1.0,1.0,0.0,1.0});
        CHECK(tiny);
        near(tiny->t_enter, 0.5); near(tiny->t_exit, 0.75);
    });
    suite.run("clipping validates bounds", [] {
        const Ray2D ray{Vec2{-2.0,0.0},Vec2{2.0,0.0}};
        for (const auto& box : std::vector<Bounds2D>{
                 {1.0,-1.0,-1.0,1.0}, {-1.0,1.0,0.0,0.0},
                 {-1.0,std::numeric_limits<double>::infinity(),-1.0,1.0},
                 {-1.0,1.0,std::numeric_limits<double>::quiet_NaN(),1.0}}) {
            throws<std::invalid_argument>([&] { (void)clip_ray_to_box(ray, box); });
        }
    });
    suite.run("Eigen vector arithmetic and sparse multiplication", [] {
        const Vec2 a{3.0,4.0};
        near(a.norm(), 5.0); near(a.dot(Vec2{2.0,-1.0}), 2.0);
        CHECK((2.0*a-Vec2{1.0,2.0}) == Vec2(5.0,6.0));
        SparseMatrix matrix(3,2);
        const std::vector<Triplet> values{{0,0,2.0},{0,1,1.0},{1,1,3.0},{2,0,-1.0}};
        matrix.setFromTriplets(values.begin(),values.end());
        Vector x(2); x << 4.0,5.0;
        const Vector b = matrix*x;
        CHECK(b.size() == 3);
        near(b[0],13.0); near(b[1],15.0); near(b[2],-4.0);
        const Vector back = matrix.transpose()*b;
        near(back[0],30.0); near(back[1],58.0);
    });
    return suite.result();
}
