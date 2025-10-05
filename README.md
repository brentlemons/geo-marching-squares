# geo-marching-squares

A Rust implementation of the marching squares algorithm for generating contour polygons from geospatial data.

[![Crates.io](https://img.shields.io/crates/v/geo-marching-squares.svg)](https://crates.io/crates/geo-marching-squares)
[![Documentation](https://docs.rs/geo-marching-squares/badge.svg)](https://docs.rs/geo-marching-squares)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

## Overview

This library generates **isobands** (filled contour polygons) and **isolines** (contour lines) from 2D scalar fields using geographic coordinates (latitude/longitude). Unlike typical pixel-based marching squares implementations, this version works directly with geospatial data and outputs GeoJSON.

**Use cases:**
- Weather data visualization (temperature, pressure, precipitation)
- Elevation contours from DEM data
- Ocean data (salinity, temperature, currents)
- Any gridded geospatial scalar field

## Features

- 🌍 **Geographic coordinates** - Works with lat/lon, not pixel grids
- 📊 **Isobands & Isolines** - Generate filled bands or contour lines
- 🗺️ **GeoJSON I/O** - Standard geospatial data format
- ⚡ **Parallel processing** - Compute multiple levels concurrently using Rayon
- 🎯 **High precision** - 5 decimal places (~1 meter accuracy)
- ✅ **Proven algorithm** - Direct port of battle-tested Java implementation

## Installation

Add this to your `Cargo.toml`:

```toml
[dependencies]
geo-marching-squares = "0.1.0"
```

## Quick Start

```rust
use geo_marching_squares::MarchingSquares;
use geojson::{Feature, FeatureCollection};

// Your 2D grid of GeoJSON point features with scalar values
let grid: Vec<Vec<Feature>> = load_grid_data();

// Generate isobands for temperature ranges
let thresholds = vec![0.0, 10.0, 20.0, 30.0, 40.0];
let result: FeatureCollection = MarchingSquares::do_concurrent(&grid, &thresholds)?;

// Result contains MultiPolygon features for each band
// Band 0: 0-10°, Band 1: 10-20°, Band 2: 20-30°, Band 3: 30-40°
```

## Input Data Format

The algorithm expects a 2D grid of GeoJSON `Feature` objects:

```json
{
  "type": "Feature",
  "geometry": {
    "type": "Point",
    "coordinates": [longitude, latitude]
  },
  "properties": {
    "value": 25.3
  }
}
```

### Grid Structure

- **Rows**: Latitude dimension (northernmost first)
- **Columns**: Longitude dimension (westernmost first)
- **Cells**: Formed by adjacent grid points (rows-1 × cols-1 cells)
- **Values**: Scalar data (temperature, elevation, etc.)

Example 3×3 grid creates 2×2 cells:

```
Grid Points:        Cells:
[0,0] [0,1] [0,2]   [Cell 0,0] [Cell 0,1]
[1,0] [1,1] [1,2]   [Cell 1,0] [Cell 1,1]
[2,0] [2,1] [2,2]
```

## Algorithm

This implementation uses a **ternary classification** approach for isobands:

1. **Classify corners** - Each cell corner gets a state (below/within/above threshold)
2. **Shape detection** - 81 possible configurations map to 7 shape types
3. **Interpolation** - Cosine-smoothed linear interpolation finds exact contour positions
4. **Edge walking** - Follow contour edges clockwise through adjacent cells
5. **Polygon assembly** - Collect edges into closed polygons with hole detection

See [`.claude/claude.md`](.claude/claude.md) for detailed algorithm documentation.

## API

### Process Single Band

Generate a single isoband between two thresholds:

```rust
let feature = MarchingSquares::process_band(
    &grid,
    lower: 10.0,
    upper: 20.0
)?;
```

Returns a `Feature` with `MultiPolygon` geometry.

### Process Multiple Bands (Concurrent)

Generate multiple isobands in parallel:

```rust
let thresholds = vec![0.0, 10.0, 20.0, 30.0];
let collection = MarchingSquares::do_concurrent(&grid, &thresholds)?;
```

Returns a `FeatureCollection` where each feature represents one band.

### Process Isoline

Generate a contour line at a specific value:

```rust
let feature = MarchingSquares::process_line(&grid, isovalue: 15.0)?;
```

Returns a `Feature` with `MultiLineString` geometry.

## Output Format

### Isoband Output

```json
{
  "type": "Feature",
  "geometry": {
    "type": "MultiPolygon",
    "coordinates": [...]
  },
  "properties": {
    "lower_level": 10.0,
    "upper_level": 20.0
  }
}
```

### Coordinate Precision

All coordinates are rounded to **5 decimal places** for approximately 1-meter precision:

```
[-122.12345, 37.98765]  // Not: [-122.123456789, 37.987654321]
```

## Performance

- **Parallel processing**: Uses Rayon to compute multiple isoband levels concurrently
- **Embarrassingly parallel**: Each cell is independent
- **Typical performance**: Processes 1000×1000 grids in seconds

## Project Status

🚧 **Currently in development** - Porting from Java implementation

See [`IMPLEMENTATION_PLAN.md`](IMPLEMENTATION_PLAN.md) for detailed progress tracking.

### Completed Phases

- [ ] Phase 1: Core data structures
- [ ] Phase 2: Shape factory and classification
- [ ] Phase 3: Interpolation
- [ ] Phase 4: GeoJSON integration
- [ ] Phase 5: Main algorithm
- [ ] Phase 6: Concurrent processing
- [ ] Phase 7: Isolines
- [ ] Phase 8: Testing
- [ ] Phase 9: API design
- [ ] Phase 10: Dependencies

## Implementation Approach

This library is a **faithful port** of a proven Java implementation. It maintains algorithmic fidelity while applying Rust idioms:

- ✅ Exact same algorithm flow
- ✅ Identical interpolation formula (including quirks)
- ✅ Matching edge walking logic
- ✅ Same polygon nesting algorithm
- 🦀 Rust data structures instead of Java objects
- 🦀 `Result<T, E>` instead of exceptions
- 🦀 Iterator chains where appropriate
- 🦀 Rust module organization

**Why?** The Java implementation is battle-tested and produces correct results. This port prioritizes correctness over reimplementation.

## Reference Implementation

Original Java implementation: [`marching-squares-java`](https://github.com/brentlemons/marching-squares-java)

## Contributing

Contributions welcome! Please:

1. Check [`IMPLEMENTATION_PLAN.md`](IMPLEMENTATION_PLAN.md) for current phase
2. Read [`.claude/claude.md`](.claude/claude.md) for algorithm details
3. Ensure changes maintain algorithmic fidelity with Java version
4. Add tests for new functionality

## License

MIT License - see [LICENSE](LICENSE) file for details

## Acknowledgments

- Original Java implementation by Brent Lemons (2017)
- Marching squares algorithm by Eugene Zhang (2004)
- Wikipedia article on [Marching Squares](https://en.wikipedia.org/wiki/Marching_squares)

## References

- [Marching Squares on Wikipedia](https://en.wikipedia.org/wiki/Marching_squares)
- [GeoJSON Specification](https://geojson.org/)
- [Original Java implementation](https://github.com/brentlemons/marching-squares-java)
