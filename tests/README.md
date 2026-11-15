# Geometry tests

Requires C++20, CMake >= 3.16 and a locally installed Eigen3 (>= 3.3).
No dependencies are downloaded by CMake. Test executables use a small local
test harness, with checks enabled in both Debug and Release builds.

From the repository root:

```sh
cmake -S . -B output/build -DCMAKE_BUILD_TYPE=Debug
cmake --build output/build --parallel
ctest --test-dir output/build --output-on-failure
```

Use `ctest --test-dir output/build -V` to see individual cases.
`numerical_geometry` is a library target; the placeholder `src/main.cpp` is not
built as an application yet. The projector and reconstruction are not implemented.

Coordinate conventions are documented in
[Coordinate conventions](../docs/coordinate-conventions.md).

Coverage includes analytic coordinates and path lengths, rectangular grid
indexing, boundary ownership, invalid inputs, sparse/vector type integration,
source-detector distances, fan symmetry, partial scans, and grid/scanner adapters.
