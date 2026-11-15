#include "scanner.hpp"
#include "config.hpp"
#include "test_support.hpp"

#include <limits>

using namespace NumReconstruction;
using Test::near;
using Test::throws;
constexpr double pi = std::numbers::pi;

int main() {
    Test::Suite suite;
    suite.run("source and central detector at cardinal views", [] {
        const FanBeamGeometry scanner{2.5,5.5,4,3,0.6};
        const Vec2 sources[]{{2.5,0.0},{0.0,2.5},{-2.5,0.0},{0.0,-2.5}};
        const Vec2 detectors[]{{-3.0,0.0},{0.0,-3.0},{3.0,0.0},{0.0,3.0}};
        for (int view = 0; view < 4; ++view) {
            near((scanner.source_position(view)-sources[view]).norm(),0.0);
            near((scanner.detector_position(view,1)-detectors[view]).norm(),0.0);
            near(scanner.ray(view,1).length(),5.5);
        }
    });
    suite.run("detector arc sign and exact oblique positions", [] {
        // Three bin centers are -pi/6,0,pi/6 for this aperture.
        const FanBeamGeometry scanner{2.5,5.5,4,3,pi/4.0};
        near(scanner.fan_angle(0),-pi/6.0);
        near(scanner.fan_angle(2),pi/6.0);
        const auto lower = scanner.detector_position(0,0);
        const auto upper = scanner.detector_position(0,2);
        near(lower.x(),2.5-5.5*std::sqrt(3.0)/2.0);
        near(lower.y(),-2.75);
        near(upper.x(),lower.x()); near(upper.y(),2.75);
        // At pi/2 the positive-fan detector lies to the left of the central one.
        near(scanner.detector_position(1,2).x(),-2.75);
        near(scanner.detector_position(1,2).y(),2.5-5.5*std::sqrt(3.0)/2.0);
    });
    suite.run("full orbit excludes repeated final view", [] {
        const FanBeamGeometry scanner{2.5,5.5,4,2,0.6};
        near(scanner.angle(0),0.0); near(scanner.angle(3),1.5*pi);
        CHECK(scanner.angle(3) < 2.0*pi);
        throws<std::out_of_range>([&] { (void)scanner.angle(4); });
    });
    suite.run("partial orbit nonzero start and single detector", [] {
        const FanBeamGeometry scanner{2.5,5.5,2,1,0.6,pi/2.0,pi};
        near(scanner.angle(0),pi/2.0); near(scanner.angle(1),pi);
        near(scanner.fan_angle(0),0.0);
        near(scanner.detector_position(1,0).x(),3.0);
        const FanBeamGeometry single{2.5,5.5,1,1,0.6,-pi/2.0,pi};
        near(single.angle(0),-pi/2.0);
        near(single.source_position(0).y(),-2.5);
        near(single.detector_position(0,0).y(),3.0);
    });
    suite.run("even and odd detector symmetry and source distance", [] {
        for (int count : {1,2,5,96}) {
            const FanBeamGeometry scanner{2.5,5.5,11,count,0.6};
            for (int k = 0; k < count; ++k) {
                near(scanner.fan_angle(k),-scanner.fan_angle(count-1-k));
                CHECK(std::abs(scanner.fan_angle(k)) < scanner.get_half_fan_angle());
                for (int view = 0; view < 11; ++view) {
                    near(scanner.source_position(view).norm(),2.5);
                    near(scanner.ray(view,k).length(),5.5);
                }
            }
        }
    });
    suite.run("measurement order and inverse mapping", [] {
        const FanBeamGeometry scanner{2.5,5.5,4,3,0.6};
        CHECK(scanner.measurement_count() == 12);
        CHECK(scanner.measurement_index(1,0) == 3);
        CHECK(scanner.measurement_indices(11) == MeasurementIndex{3,2});
        for (int v = 0; v < 4; ++v)
            for (int d = 0; d < 3; ++d)
                CHECK(scanner.measurement_indices(scanner.measurement_index(v,d)) == MeasurementIndex{v,d});
    });
    suite.run("Config adapter and default grid compatibility", [] {
        const auto config = Config::Builder{}.build();
        const auto scanner = FanBeamGeometry::from_config(config);
        near(scanner.get_radius(),2.5); near(scanner.get_detector_distance(),5.5);
        CHECK(scanner.get_num_angles() == 360 && scanner.get_num_detectors() == 256);
        CHECK(scanner.measurement_count() == 92160);
        near(scanner.get_start_angle(),0.0); near(scanner.get_scan_angle(),2.0*pi);
        near(scanner.get_half_fan_angle(),35.0*pi/180.0);
        validate_geometry(Grid2D::from_config(config),scanner);
        const auto changed = Config::Builder{}.set_start_angle(0.3).set_scan_angle(pi)
            .set_num_angles(7).set_num_detectors(9).set_radius(3.0)
            .set_detector_distance(7.0).set_half_fan_angle(0.7).build();
        const auto other = FanBeamGeometry::from_config(changed);
        near(other.angle(0),0.3); near(other.get_scan_angle(),pi);
        near(other.get_radius(),3.0); near(other.get_detector_distance(),7.0);
        near(other.get_half_fan_angle(),0.7);
        CHECK(other.measurement_count() == 63);
    });
    suite.run("invalid scanner parameters", [] {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        const double inf = std::numeric_limits<double>::infinity();
        for (double bad : {0.0,-1.0,nan,inf}) {
            throws<std::invalid_argument>([&] { FanBeamGeometry{bad,5.5,4,3,0.6}; });
            throws<std::invalid_argument>([&] { FanBeamGeometry{2.5,bad,4,3,0.6}; });
            throws<std::invalid_argument>([&] { FanBeamGeometry{2.5,5.5,4,3,bad}; });
            throws<std::invalid_argument>([&] { FanBeamGeometry{2.5,5.5,4,3,0.6,0.0,bad}; });
        }
        throws<std::invalid_argument>([] { FanBeamGeometry{2.5,5.5,0,3,0.6}; });
        throws<std::invalid_argument>([] { FanBeamGeometry{2.5,5.5,4,-1,0.6}; });
        throws<std::invalid_argument>([] { FanBeamGeometry{2.5,5.5,4,3,pi/2.0}; });
        throws<std::invalid_argument>([] { FanBeamGeometry{2.5,5.5,4,3,0.6,0.0,3.0*pi}; });
        throws<std::invalid_argument>([&] { FanBeamGeometry{2.5,5.5,4,3,0.6,nan}; });
        throws<std::invalid_argument>([&] { FanBeamGeometry{2.5,5.5,4,3,0.6,inf}; });
        throws<std::invalid_argument>([] { FanBeamGeometry{1e308,1e308,4,3,0.6}; });
        throws<std::invalid_argument>([] { FanBeamGeometry{2.5,5.5,4,3,0.6,1e100}; });
    });
    suite.run("invalid view detector and flattened indices", [] {
        const FanBeamGeometry scanner{2.5,5.5,4,3,0.6};
        for (int bad : {-1,4}) {
            throws<std::out_of_range>([&] { (void)scanner.angle(bad); });
            throws<std::out_of_range>([&] { (void)scanner.source_position(bad); });
            throws<std::out_of_range>([&] { (void)scanner.detector_position(bad,0); });
            throws<std::out_of_range>([&] { (void)scanner.measurement_index(bad,0); });
        }
        for (int bad : {-1,3}) {
            throws<std::out_of_range>([&] { (void)scanner.fan_angle(bad); });
            throws<std::out_of_range>([&] { (void)scanner.detector_position(0,bad); });
            throws<std::out_of_range>([&] { (void)scanner.ray(0,bad); });
            throws<std::out_of_range>([&] { (void)scanner.measurement_index(0,bad); });
        }
        for (Index bad : {Index{-1},Index{12}})
            throws<std::out_of_range>([&] { (void)scanner.measurement_indices(bad); });
    });
    suite.run("cross-object clearance distance and fan validation", [] {
        const Grid2D grid{64,64,2.0,2.0};
        validate_geometry(grid,FanBeamGeometry{2.5,5.5,120,96,35.0*pi/180.0});
        throws<std::invalid_argument>([&] {
            validate_geometry(grid,FanBeamGeometry{1.0,5.5,4,3,0.8});
        });
        throws<std::invalid_argument>([&] {
            validate_geometry(grid,FanBeamGeometry{std::sqrt(2.0),5.5,4,3,1.0});
        });
        throws<std::invalid_argument>([&] {
            validate_geometry(grid,FanBeamGeometry{2.5,3.0,4,3,0.8});
        });
        throws<std::invalid_argument>([&] {
            validate_geometry(grid,FanBeamGeometry{2.5,5.5,4,3,0.1});
        });
        const Grid2D rectangle{3,2,6.0,2.0};
        validate_geometry(rectangle,FanBeamGeometry{5.0,9.0,4,3,0.8});
    });
    suite.run("scanner rays integrate with grid clipping", [] {
        const Grid2D grid{64,64,2.0,2.0};
        const FanBeamGeometry scanner{2.5,5.5,4,3,0.6};
        for (int v = 0; v < 4; ++v) {
            const auto ray = scanner.ray(v,1);
            const auto interval = clip_ray_to_box(ray,grid.bounds());
            CHECK(interval);
            near((interval->t_exit-interval->t_enter)*ray.length(),2.0);
            CHECK(grid.locate(ray.point_at((interval->t_enter+interval->t_exit)/2.0)));
        }
        // Edge detectors may miss the grid even with a valid scanner geometry.
        const FanBeamGeometry wide{2.5,5.5,1,3,1.4};
        validate_geometry(grid,wide);
        CHECK(!clip_ray_to_box(wide.ray(0,0),grid.bounds()));
    });
    return suite.result();
}
