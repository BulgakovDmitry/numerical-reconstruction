# Coordinate conventions

- Centered rectangular grid; x right, y up, row zero at the top.
- Row-major image index: `row * pixels_x + col`.
- Pixel ownership: `[left,right) x (bottom,top]`. Right/bottom external edges
  are excluded from `locate`; nonfinite/outside points return `nullopt`.
- `clip_ray_to_box` clips against a **closed** rectangle and a finite segment
  `S+t*(D-S)`, `t in [0,1]`. Point contacts are discarded; paths along rectangle
  edges are retained. Pixel ownership is a separate operation.
- Views exclude the scan interval's final angle; detectors are angular bin
  centers. One detector has fan angle zero.
- The detector arc is centered at the **source**, as in the README formula.
  The README prose currently describes a different center; the implementation
  follows the formula explicitly.
- `validate_geometry` uses sufficient conditions for the circumscribed circle:
  `R > rho`, `d > R+rho`, `half_fan_angle >= asin(rho/R)`. These are conservative
  for partial scans and do not guarantee sufficient discrete measurements.
- `SiddonProjector::trace_ray` follows the pixel ownership rule above: paths
  along internal grid lines are counted on one side; paths along the external
  right/bottom edges are excluded. Misses and point contacts give no entries.

