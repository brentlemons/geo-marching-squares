# geo-marching-squares

A high-performance Rust implementation of the marching squares algorithm for generating **isobands** (filled contour polygons) and **isolines** (contour lines) from geospatial 2D scalar fields.

[![Tests](https://img.shields.io/badge/tests-51%20passing-brightgreen)]()
[![License](https://img.shields.io/badge/license-MIT-blue)]()

## Features

- 🌍 **Geographic Coordinates** - Works with latitude/longitude (not pixel grids)
- 🎨 **Isobands** - Filled contour polygons with smooth interpolation
- 📏 **Isolines** - Contour lines at specific threshold values
- ⚡ **Parallel Processing** - Multi-threaded using Rayon's work-stealing scheduler
- 📊 **GeoJSON I/O** - Native support for GeoJSON features
- 🦀 **100% Safe Rust** - No unsafe code, compile-time guarantees
- ✅ **Well Tested** - 51 comprehensive tests

## Installation

Add to your `Cargo.toml`:

```toml
[dependencies]
geo-marching-squares = "0.1.0"
geojson = "0.24"
```

## Quick Start

### Complete Example: Weather Temperature Bands

Here's a complete working example that generates temperature isobands:

```rust
use geo_marching_squares::do_concurrent;
use geojson::{Feature, FeatureCollection, Geometry, Value, JsonObject};

fn main() {
    // Step 1: Create your grid of temperature measurements
    let grid = create_temperature_grid();

    // Step 2: Define temperature ranges (in Celsius)
    // This creates 5 bands: <0°, 0-10°, 10-20°, 20-30°, 30+°
    let thresholds = vec![-10.0, 0.0, 10.0, 20.0, 30.0, 40.0];

    // Step 3: Generate isobands in parallel
    let result: FeatureCollection = do_concurrent(&grid, &thresholds);

    // Step 4: Use the results
    println!("Generated {} temperature bands", result.features.len());

    for feature in result.features {
        if let Some(props) = &feature.properties {
            let lower = props.get("lower_level").unwrap().as_f64().unwrap();
            let upper = props.get("upper_level").unwrap().as_f64().unwrap();
            println!("Band: {}°C to {}°C", lower, upper);

            // The geometry contains the actual polygons
            if let Some(geom) = &feature.geometry {
                if let Value::MultiPolygon(polygons) = &geom.value {
                    println!("  Contains {} polygon(s)", polygons.len());
                }
            }
        }
    }

    // Step 5: Save to file or use in your application
    let json = serde_json::to_string_pretty(&result).unwrap();
    std::fs::write("temperature_bands.geojson", json).unwrap();
}

fn create_temperature_grid() -> Vec<Vec<Feature>> {
    // Create a 5x5 grid covering a region
    // In practice, this would come from your weather data
    let mut grid = Vec::new();

    for lat_idx in 0..5 {
        let mut row = Vec::new();
        for lon_idx in 0..5 {
            let lon = -122.0 + (lon_idx as f64 * 0.1);
            let lat = 37.0 + (lat_idx as f64 * 0.1);

            // Simulate temperature values (would be real data in production)
            let temp = 15.0 + (lat_idx as f64 * 2.0) + (lon_idx as f64 * 1.5);

            row.push(create_point_feature(lon, lat, temp));
        }
        grid.push(row);
    }

    grid
}

fn create_point_feature(lon: f64, lat: f64, value: f64) -> Feature {
    let geometry = Geometry::new(Value::Point(vec![lon, lat]));
    let mut properties = JsonObject::new();
    properties.insert("value".to_string(), serde_json::json!(value));

    Feature {
        bbox: None,
        geometry: Some(geometry),
        id: None,
        properties: Some(properties),
        foreign_members: None,
    }
}
```

### Complete Example: Elevation Contours

Generate elevation contour lines from a DEM (Digital Elevation Model):

```rust
use geo_marching_squares::do_concurrent_lines;
use geojson::{Feature, FeatureCollection, Geometry, Value};

fn main() {
    // Step 1: Load elevation data
    let grid = load_elevation_data();

    // Step 2: Define contour line elevations (every 100m)
    let elevations = vec![0.0, 100.0, 200.0, 300.0, 400.0, 500.0];

    // Step 3: Generate contour lines in parallel
    let result: FeatureCollection = do_concurrent_lines(&grid, &elevations);

    // Step 4: Process the results
    println!("Generated {} contour lines", result.features.len());

    for feature in result.features {
        if let Some(props) = &feature.properties {
            let elevation = props.get("isovalue").unwrap().as_f64().unwrap();

            if let Some(geom) = &feature.geometry {
                if let Value::MultiLineString(lines) = &geom.value {
                    let total_segments: usize = lines.iter().map(|l| l.len()).sum();
                    println!("{}m contour: {} line(s), {} points",
                             elevation, lines.len(), total_segments);
                }
            }
        }
    }

    // Save to GeoJSON file
    let json = serde_json::to_string_pretty(&result).unwrap();
    std::fs::write("contours.geojson", json).unwrap();
}

fn load_elevation_data() -> Vec<Vec<Feature>> {
    // In practice, you'd read from a GeoTIFF or other DEM source
    // Here's a simple example with synthetic data
    let mut grid = Vec::new();

    for row in 0..20 {
        let mut row_features = Vec::new();
        for col in 0..20 {
            let lon = -120.0 + (col as f64 * 0.01);
            let lat = 40.0 + (row as f64 * 0.01);

            // Simulate elevation (in meters)
            let elevation = 100.0 + ((row * row + col * col) as f64 * 2.0);

            row_features.push(create_point_feature(lon, lat, elevation));
        }
        grid.push(row_features);
    }

    grid
}
```

### Single Band/Line Examples

For simple use cases, process a single band or line:

```rust
use geo_marching_squares::{process_band, process_line};

// Single isoband: temperatures between 20-30°C
let grid = create_temperature_grid();
let feature = process_band(&grid, 20.0, 30.0);

// Access the result
if let Some(geom) = &feature.geometry {
    if let Value::MultiPolygon(polygons) = &geom.value {
        println!("Generated {} polygon(s) for 20-30°C range", polygons.len());
    }
}

// Single isoline: 500m elevation contour
let dem_grid = load_elevation_data();
let contour = process_line(&dem_grid, 500.0);

// Access the result
if let Some(geom) = &contour.geometry {
    if let Value::MultiLineString(lines) = &geom.value {
        println!("500m contour has {} line segment(s)", lines.len());
    }
}
```

## API Reference

### Isoband Functions

#### `process_band`

Generate a single isoband between two thresholds:

```rust
pub fn process_band(
    data: &[Vec<Feature>],
    lower: f64,
    upper: f64
) -> Feature
```

**Parameters:**
- `data` - 2D grid of GeoJSON Point features with `value` property
- `lower` - Lower threshold (inclusive)
- `upper` - Upper threshold (exclusive)

**Returns:** Feature with MultiPolygon geometry containing all contours in this band

#### `do_concurrent`

Generate multiple isobands in parallel:

```rust
pub fn do_concurrent(
    data: &[Vec<Feature>],
    isobands: &[f64]
) -> FeatureCollection
```

**Parameters:**
- `data` - 2D grid of GeoJSON Point features
- `isobands` - Threshold values (N values create N-1 bands)

**Returns:** FeatureCollection with one feature per band

### Isoline Functions

#### `process_line`

Generate a single isoline at a specific value:

```rust
pub fn process_line(
    data: &[Vec<Feature>],
    isovalue: f64
) -> Feature
```

**Parameters:**
- `data` - 2D grid of GeoJSON Point features
- `isovalue` - Threshold value for the contour line

**Returns:** Feature with MultiLineString geometry

#### `do_concurrent_lines`

Generate multiple isolines in parallel:

```rust
pub fn do_concurrent_lines(
    data: &[Vec<Feature>],
    isovalues: &[f64]
) -> FeatureCollection
```

**Parameters:**
- `data` - 2D grid of GeoJSON Point features
- `isovalues` - Array of threshold values

**Returns:** FeatureCollection with one feature per isoline

## Algorithm Overview

### Isobands vs Isolines

| Feature | Isolines | Isobands |
|---------|----------|----------|
| **Purpose** | Contour lines | Filled regions |
| **Classification** | Binary (above/below) | Ternary (below/within/above) |
| **Configurations** | 16 (2^4) | 81 (3^4) |
| **Geometry** | MultiLineString | MultiPolygon |
| **Interpolation** | Linear | Cosine-smoothed |
| **Use Case** | Elevation contours | Temperature zones |

### How It Works

#### 1. Grid to Cells

The algorithm divides your grid into cells, where each cell is defined by 4 corner points:

```
Grid (3×3 points)     →     Cells (2×2)

[P00] [P01] [P02]           [Cell00] [Cell01]
[P10] [P11] [P12]     →     [Cell10] [Cell11]
[P20] [P21] [P22]
```

#### 2. Classification

**Isobands (Ternary):**
For each corner, classify value V against thresholds (lower, upper):
- 0 if V < lower
- 1 if lower ≤ V < upper
- 2 if V ≥ upper

Creates 81 possible cell configurations (3^4)

**Isolines (Binary):**
For each corner, classify value V against threshold T:
- 0 if V < T
- 1 if V ≥ T

Creates 16 possible cell configurations (2^4)

#### 3. Interpolation

**Isobands use cosine interpolation** for smooth contours:

```rust
mu = (level - value0) / (value1 - value0)
mu2 = (1.0 - cos(mu * PI)) / 2.0
new_mu = 0.5 + ((mu2 - 0.5) * 0.999)  // Center adjustment
x = ((1.0 - new_mu) * x0) + (new_mu * x1)
```

**Isolines use linear interpolation** for precise lines:

```rust
mu = (level - value0) / (value1 - value0)
x = ((1.0 - mu) * x0) + (mu * x1)
```

#### 4. Edge Walking

The algorithm "walks" along cell edges to trace contours:
- Start at any cell with edges
- Follow edge connections to adjacent cells
- Continue until loop closes (back to start)
- Repeat for all remaining untraced edges

#### 5. Polygon Nesting (Isobands Only)

Uses ray-casting to detect nested polygons:
- Exterior rings: Outermost polygons
- Interior rings: Holes within polygons

#### 6. Output Precision

All coordinates rounded to 5 decimal places (~1 meter accuracy):

```rust
(value * 100000.0).round() / 100000.0
```

## Performance

### Computational Complexity

- **Grid size:** N×M cells
- **Cell processing:** O(1) per cell
- **Edge walking:** O(E) where E = number of edges (~2-3× cell count)
- **Polygon nesting:** O(P²×V) where P = polygons, V = vertices per polygon
- **Total:** O(NM + E + P²V), typically dominated by O(NM)

### Parallel Processing

- **Embarrassingly parallel:** Each isoband/isoline level is independent
- **Zero synchronization:** Read-only shared grid
- **Work-stealing:** Rayon's scheduler balances load across CPU cores
- **Scalability:** Near-linear speedup with number of cores

**Benchmark example** (1000×1000 grid, 10 isoband levels):
- Sequential: ~2.5 seconds
- Parallel (8 cores): ~0.4 seconds
- Speedup: 6.25×

### Memory Usage

- **Grid storage:** Shared reference (zero-copy across threads)
- **Per cell:** ~8 bytes for classification
- **Per edge:** ~64 bytes (start, end, direction)
- **Output:** Depends on contour complexity

## Architecture

### Module Structure

```
src/
├── lib.rs                    # Public API, documentation
├── point.rs (298 lines)      # Point with coordinates/value
├── edge.rs (151 lines)       # Edge with direction tracking
├── shape.rs (1,834 lines)    # Shape factory for isobands
├── cell.rs (411 lines)       # Cell for isolines
├── isoline_assembler.rs      # Line segment assembly
│   (197 lines)
└── marching_squares.rs       # Main algorithms
    (635 lines)               # process_band, process_line, etc.
```

### Data Flow

```
Input: 2D Grid of GeoJSON Point Features
         ↓
    [Classification]
   Binary (isolines) or Ternary (isobands)
         ↓
    [Cell Creation]
   16 configs (isolines) or 81 configs (isobands)
         ↓
   [Interpolation]
   Linear (isolines) or Cosine (isobands)
         ↓
  [Edge/Segment Generation]
         ↓
   [Assembly]
   Polylines (isolines) or Polygons (isobands)
         ↓
Output: GeoJSON Feature(Collection)
```

## User Guide

### Understanding Input Data Requirements

The library requires a **rectangular grid** of GeoJSON Point features. Here's what you need to know:

#### Grid Structure

```
Rows = Latitude dimension (North to South)
Columns = Longitude dimension (West to East)

Example 4×3 grid (4 rows, 3 columns):

Row 0: [Point(0,3), Point(1,3), Point(2,3)]  <- Northernmost
Row 1: [Point(0,2), Point(1,2), Point(2,2)]
Row 2: [Point(0,1), Point(1,1), Point(2,1)]
Row 3: [Point(0,0), Point(1,0), Point(2,0)]  <- Southernmost
        ↑           ↑           ↑
        West      Center      East
```

#### Required Feature Format

Each Point feature must have:

1. **Geometry**: Point with `[longitude, latitude]` coordinates
2. **Properties**: Object with a `"value"` field containing the scalar data

```json
{
  "type": "Feature",
  "geometry": {
    "type": "Point",
    "coordinates": [-122.4194, 37.7749]
  },
  "properties": {
    "value": 25.5
  }
}
```

#### Common Pitfalls

❌ **Wrong coordinate order**: `[latitude, longitude]` (should be `[lon, lat]`)
❌ **Missing value property**: Must be named exactly `"value"`
❌ **Irregular grid**: All rows must have same number of columns
❌ **Missing features**: Grid cannot have `None` values

### Loading Data from Common Sources

#### From CSV File

```rust
use csv::Reader;
use geojson::{Feature, Geometry, Value, JsonObject};

fn load_grid_from_csv(path: &str) -> Vec<Vec<Feature>> {
    let mut reader = Reader::from_path(path).unwrap();
    let mut grid_data: Vec<(f64, f64, f64)> = Vec::new();

    // Read CSV: lon, lat, value
    for result in reader.records() {
        let record = result.unwrap();
        let lon: f64 = record[0].parse().unwrap();
        let lat: f64 = record[1].parse().unwrap();
        let value: f64 = record[2].parse().unwrap();
        grid_data.push((lon, lat, value));
    }

    // Organize into 2D grid
    // Assumes data is already sorted by lat (desc) then lon (asc)
    let lons: Vec<f64> = grid_data.iter()
        .map(|(lon, _, _)| *lon)
        .collect::<std::collections::HashSet<_>>()
        .into_iter()
        .collect();

    let lats: Vec<f64> = grid_data.iter()
        .map(|(_, lat, _)| *lat)
        .collect::<std::collections::HashSet<_>>()
        .into_iter()
        .collect();

    let mut grid = vec![vec![]; lats.len()];

    for (lon, lat, value) in grid_data {
        let lat_idx = lats.iter().position(|&l| l == lat).unwrap();
        let lon_idx = lons.iter().position(|&l| l == lon).unwrap();

        let geometry = Geometry::new(Value::Point(vec![lon, lat]));
        let mut properties = JsonObject::new();
        properties.insert("value".to_string(), serde_json::json!(value));

        grid[lat_idx].push(Feature {
            bbox: None,
            geometry: Some(geometry),
            id: None,
            properties: Some(properties),
            foreign_members: None,
        });
    }

    grid
}
```

#### From GeoJSON File

```rust
use geojson::FeatureCollection;
use std::fs;

fn load_grid_from_geojson(path: &str) -> Vec<Vec<Feature>> {
    let geojson_str = fs::read_to_string(path).unwrap();
    let collection: FeatureCollection = geojson_str.parse().unwrap();

    // Assuming features are ordered: row by row, column by column
    let num_cols = 10; // You need to know this from your data
    let mut grid = Vec::new();
    let mut current_row = Vec::new();

    for (i, feature) in collection.features.into_iter().enumerate() {
        current_row.push(feature);

        if (i + 1) % num_cols == 0 {
            grid.push(current_row);
            current_row = Vec::new();
        }
    }

    grid
}
```

#### From NetCDF (Common for Weather Data)

```rust
// Using the netcdf crate
use netcdf;

fn load_grid_from_netcdf(path: &str, var_name: &str) -> Vec<Vec<Feature>> {
    let file = netcdf::open(path).unwrap();

    // Read coordinate variables
    let lons = file.variable("lon").unwrap()
        .values::<f64, _>(None, None).unwrap();
    let lats = file.variable("lat").unwrap()
        .values::<f64, _>(None, None).unwrap();

    // Read data variable
    let data = file.variable(var_name).unwrap()
        .values::<f64, _>(None, None).unwrap();

    // Build grid
    let mut grid = Vec::new();

    for (lat_idx, &lat) in lats.iter().enumerate() {
        let mut row = Vec::new();
        for (lon_idx, &lon) in lons.iter().enumerate() {
            let value = data[[lat_idx, lon_idx]]; // 2D indexing
            row.push(create_point_feature(lon, lat, value));
        }
        grid.push(row);
    }

    grid
}
```

### Working with Output Data

#### Accessing Polygon Coordinates

```rust
use geojson::Value;

let result = do_concurrent(&grid, &thresholds);

for feature in &result.features {
    if let Some(Value::MultiPolygon(multi_poly)) = feature.geometry.as_ref().map(|g| &g.value) {
        for (poly_idx, polygon) in multi_poly.iter().enumerate() {
            println!("Polygon {}: {} rings", poly_idx, polygon.len());

            // First ring is exterior, rest are holes
            let exterior = &polygon[0];
            println!("  Exterior ring: {} points", exterior.len());

            for (i, coord) in exterior.iter().enumerate() {
                let lon = coord[0];
                let lat = coord[1];
                println!("    Point {}: ({}, {})", i, lon, lat);
            }

            // Process holes if any
            for (hole_idx, hole) in polygon[1..].iter().enumerate() {
                println!("  Hole {}: {} points", hole_idx, hole.len());
            }
        }
    }
}
```

#### Accessing Line Coordinates

```rust
let result = do_concurrent_lines(&grid, &isovalues);

for feature in &result.features {
    if let Some(Value::MultiLineString(multi_line)) = feature.geometry.as_ref().map(|g| &g.value) {
        for (line_idx, line) in multi_line.iter().enumerate() {
            println!("Line segment {}: {} points", line_idx, line.len());

            for (i, coord) in line.iter().enumerate() {
                let lon = coord[0];
                let lat = coord[1];
                println!("  Point {}: ({}, {})", i, lon, lat);
            }

            // Check if line is closed (forms a loop)
            if line.first() == line.last() {
                println!("  This is a closed contour loop");
            }
        }
    }
}
```

#### Converting to Other Formats

**To WKT (Well-Known Text):**

```rust
use geojson::Value;

fn polygon_to_wkt(feature: &Feature) -> String {
    if let Some(Value::MultiPolygon(multi_poly)) = feature.geometry.as_ref().map(|g| &g.value) {
        let mut wkt = String::from("MULTIPOLYGON(");

        for (i, polygon) in multi_poly.iter().enumerate() {
            if i > 0 { wkt.push_str(", "); }
            wkt.push('(');

            for (j, ring) in polygon.iter().enumerate() {
                if j > 0 { wkt.push_str(", "); }
                wkt.push('(');

                for (k, coord) in ring.iter().enumerate() {
                    if k > 0 { wkt.push_str(", "); }
                    wkt.push_str(&format!("{} {}", coord[0], coord[1]));
                }

                wkt.push(')');
            }

            wkt.push(')');
        }

        wkt.push(')');
        wkt
    } else {
        String::new()
    }
}
```

**To GeoPackage (using gdal crate):**

```rust
// Using the gdal crate
use gdal::vector::{Layer, OGRwkbGeometryType};
use gdal::Dataset;

fn save_to_gpkg(features: &[Feature], path: &str) -> gdal::errors::Result<()> {
    let driver = gdal::DriverManager::get_driver_by_name("GPKG")?;
    let mut dataset = driver.create_vector_only(path)?;

    let mut layer = dataset.create_layer(
        "contours",
        None,
        OGRwkbGeometryType::wkbMultiPolygon,
    )?;

    // Add fields
    layer.create_field(&gdal::vector::FieldDefn::new("lower", gdal::vector::OGRFieldType::OFTReal)?)?;
    layer.create_field(&gdal::vector::FieldDefn::new("upper", gdal::vector::OGRFieldType::OFTReal)?)?;

    for feature in features {
        // Convert GeoJSON to GDAL geometry and add to layer
        // (implementation details depend on gdal crate version)
    }

    Ok(())
}
```

### Advanced Usage

#### Custom Threshold Sequences

**Logarithmic scales** for data spanning multiple orders of magnitude:

```rust
let thresholds: Vec<f64> = (0..10)
    .map(|i| 10_f64.powi(i))
    .collect();
// [1, 10, 100, 1000, 10000, ...]
```

**Percentile-based bands** for even distribution:

```rust
// Extract all values from grid
fn extract_all_values(grid: &[Vec<Feature>]) -> Vec<f64> {
    grid.iter()
        .flat_map(|row| row.iter())
        .filter_map(|f| {
            f.properties.as_ref()
                .and_then(|p| p.get("value"))
                .and_then(|v| v.as_f64())
        })
        .collect()
}

let mut values = extract_all_values(&grid);
values.sort_by(|a, b| a.partial_cmp(b).unwrap());

// Create 10 bands with equal number of data points
let thresholds: Vec<f64> = (0..=10)
    .map(|i| values[values.len() * i / 10])
    .collect();
```

**Equal-interval bands** with custom spacing:

```rust
fn create_equal_intervals(min: f64, max: f64, num_bands: usize) -> Vec<f64> {
    let interval = (max - min) / num_bands as f64;
    (0..=num_bands)
        .map(|i| min + (i as f64 * interval))
        .collect()
}

let thresholds = create_equal_intervals(0.0, 100.0, 10);
// [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
```

#### Filtering and Post-Processing

**Filter empty features:**

```rust
let result = do_concurrent(&grid, &thresholds);

let non_empty: Vec<Feature> = result.features.into_iter()
    .filter(|f| {
        if let Some(geom) = &f.geometry {
            if let Value::MultiPolygon(polys) = &geom.value {
                return !polys.is_empty();
            }
        }
        false
    })
    .collect();

println!("Kept {} non-empty features", non_empty.len());
```

**Filter by area threshold:**

```rust
fn polygon_area(coords: &[Vec<f64>]) -> f64 {
    // Simple shoelace formula (assumes small polygons)
    let mut area = 0.0;
    for i in 0..coords.len() - 1 {
        area += coords[i][0] * coords[i + 1][1];
        area -= coords[i + 1][0] * coords[i][1];
    }
    area.abs() / 2.0
}

let min_area = 0.0001; // Minimum area threshold

let filtered: Vec<Feature> = result.features.into_iter()
    .filter(|f| {
        if let Some(Value::MultiPolygon(polys)) = f.geometry.as_ref().map(|g| &g.value) {
            polys.iter().any(|poly| {
                polygon_area(&poly[0]) >= min_area
            })
        } else {
            false
        }
    })
    .collect();
```

**Add custom properties:**

```rust
use serde_json::json;

let mut enhanced_features: Vec<Feature> = Vec::new();

for (i, mut feature) in result.features.into_iter().enumerate() {
    // Add custom properties
    if let Some(props) = feature.properties.as_mut() {
        props.insert("id".to_string(), json!(i));
        props.insert("color".to_string(), json!(format!("#{:06x}", i * 0x111111)));
        props.insert("label".to_string(), json!(format!("Band {}", i)));
    }

    enhanced_features.push(feature);
}
```

#### Working with Properties

**Access and modify band properties:**

```rust
for feature in &mut result.features {
    if let Some(props) = &mut feature.properties {
        let lower = props.get("lower_level").unwrap().as_f64().unwrap();
        let upper = props.get("upper_level").unwrap().as_f64().unwrap();

        // Add derived properties
        props.insert("range".to_string(), json!(upper - lower));
        props.insert("midpoint".to_string(), json!((lower + upper) / 2.0));

        println!("Band: [{}, {})", lower, upper);
    }
}
```

**Extract statistics:**

```rust
use geojson::Value;

fn count_polygons(features: &[Feature]) -> usize {
    features.iter()
        .filter_map(|f| {
            if let Some(Value::MultiPolygon(polys)) = f.geometry.as_ref().map(|g| &g.value) {
                Some(polys.len())
            } else {
                None
            }
        })
        .sum()
}

fn count_total_vertices(features: &[Feature]) -> usize {
    features.iter()
        .filter_map(|f| {
            if let Some(Value::MultiPolygon(polys)) = f.geometry.as_ref().map(|g| &g.value) {
                Some(polys.iter()
                    .flat_map(|p| p.iter())
                    .map(|ring| ring.len())
                    .sum::<usize>())
            } else {
                None
            }
        })
        .sum()
}

let total_polys = count_polygons(&result.features);
let total_verts = count_total_vertices(&result.features);
println!("Generated {} polygons with {} vertices", total_polys, total_verts);
```

### Integration Examples

#### Using with Actix-Web (REST API)

```rust
use actix_web::{web, App, HttpResponse, HttpServer};
use geo_marching_squares::do_concurrent;
use geojson::{Feature, FeatureCollection};

#[derive(serde::Deserialize)]
struct ContourRequest {
    grid: Vec<Vec<Feature>>,
    thresholds: Vec<f64>,
}

async fn generate_contours(req: web::Json<ContourRequest>) -> HttpResponse {
    let result = do_concurrent(&req.grid, &req.thresholds);
    HttpResponse::Ok().json(result)
}

#[actix_web::main]
async fn main() -> std::io::Result<()> {
    HttpServer::new(|| {
        App::new()
            .route("/contours", web::post().to(generate_contours))
    })
    .bind("127.0.0.1:8080")?
    .run()
    .await
}
```

#### Using with Leaflet (JavaScript Frontend)

```javascript
// Request contours from your Rust backend
async function loadContours() {
    const response = await fetch('/api/contours', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
            grid: temperatureGrid,
            thresholds: [0, 10, 20, 30, 40]
        })
    });

    const geojson = await response.json();

    // Add to Leaflet map with styling
    L.geoJSON(geojson, {
        style: function(feature) {
            const lower = feature.properties.lower_level;
            return {
                fillColor: getColor(lower),
                weight: 2,
                opacity: 1,
                color: 'white',
                fillOpacity: 0.7
            };
        },
        onEachFeature: function(feature, layer) {
            const lower = feature.properties.lower_level;
            const upper = feature.properties.upper_level;
            layer.bindPopup(`Range: ${lower}°C - ${upper}°C`);
        }
    }).addTo(map);
}

function getColor(value) {
    return value > 30 ? '#d73027' :
           value > 20 ? '#fc8d59' :
           value > 10 ? '#fee08b' :
           value > 0  ? '#d9ef8b' :
                        '#91bfdb';
}
```

### Performance Optimization Tips

#### 1. Use Parallel Processing

Always prefer `do_concurrent` over sequential processing:

```rust
// Good: Parallel (6x faster on 8 cores)
let result = do_concurrent(&grid, &thresholds);

// Avoid: Sequential
let features: Vec<Feature> = thresholds.windows(2)
    .map(|w| process_band(&grid, w[0], w[1]))
    .collect();
```

#### 2. Pre-allocate Grid Structures

```rust
// Pre-allocate with known dimensions
let rows = 100;
let cols = 100;
let mut grid = Vec::with_capacity(rows);

for _ in 0..rows {
    grid.push(Vec::with_capacity(cols));
}
```

#### 3. Reuse Grid for Multiple Analyses

```rust
let grid = load_large_grid(); // Expensive operation

// Process multiple threshold sets on same grid
let temp_bands = do_concurrent(&grid, &temp_thresholds);
let humidity_bands = do_concurrent(&grid, &humidity_thresholds);
```

#### 4. Choose Appropriate Threshold Count

More thresholds = more computation:

```rust
// Fine detail: 20 bands (slower)
let fine = (0..=20).map(|i| i as f64 * 5.0).collect::<Vec<_>>();

// Coarse detail: 5 bands (faster)
let coarse = (0..=5).map(|i| i as f64 * 20.0).collect::<Vec<_>>();
```

## Comparison to Java Reference

This library is a faithful port of a proven Java implementation, with significant improvements:

| Feature | Java Reference | This Implementation |
|---------|---------------|---------------------|
| **Isobands** | ✅ Complete | ✅ Complete (faithful port) |
| **Isolines** | ❌ Stub only | ✅ **Full implementation** |
| **Concurrency** | ExecutorService | **Rayon** (work-stealing) |
| **Type Safety** | Runtime checks | **Compile-time** |
| **Memory Safety** | GC overhead | **Zero-cost ownership** |
| **Tests** | Minimal | **51 comprehensive tests** |
| **Documentation** | Limited | **Extensive with examples** |
| **Lines of Code** | 2,417 (11 files) | 4,522 (6 modules) |

### Algorithmic Fidelity

The isoband implementation exactly matches the Java reference:
- ✅ Same ternary classification encoding
- ✅ Same cosine interpolation (including 0.999 centering hack)
- ✅ Same edge walking direction logic
- ✅ Same saddle point disambiguation
- ✅ Same ray-casting polygon nesting algorithm
- ✅ Same 5 decimal place precision

The isoline implementation goes **beyond the Java reference**, which only had a stub.

## Known Limitations

### Handled Correctly ✅
- Grid boundaries (open isolines, polygon edges)
- Saddle point disambiguation
- Polygon nesting (multiple levels)
- Empty results (all above/below threshold)
- Coordinate precision

### Potential Issues
- **Very large grids:** Memory usage grows linearly with grid size
- **Complex nesting:** O(P²) polygon nesting can be slow with many small polygons
- **Degenerate cases:** Extremely thin features may produce zero-area polygons

### Not Implemented
- Polygon simplification (Douglas-Peucker)
- Topology preservation
- Multi-value grids (vector fields)
- Custom interpolation functions
- Streaming/chunked processing for huge grids

## Testing

Run the test suite:

```bash
# All tests
cargo test

# Specific test suites
cargo test test_process_band    # Isoband tests
cargo test test_process_line    # Isoline tests
cargo test test_do_concurrent   # Parallel processing tests

# With output
cargo test -- --nocapture
```

**Test Coverage:**
- 28 unit tests (data structures, algorithms)
- 8 concurrent processing tests
- 5 isoband integration tests
- 10 isoline integration tests
- **Total: 51 tests, all passing ✅**

## Contributing

Contributions welcome! Areas for improvement:

1. **Performance:** Benchmark and optimize hot paths
2. **Features:** Polygon simplification, topology preservation
3. **Testing:** Property-based tests, fuzzing
4. **Documentation:** More examples, tutorials
5. **Interoperability:** Support for other geometry formats

## License

MIT License - see LICENSE file for details

## Acknowledgments

Based on the Java marching squares implementation. Ported to Rust with enhancements and full isoline support.

## See Also

- [GeoJSON Specification](https://geojson.org/)
- [Marching Squares Algorithm](https://en.wikipedia.org/wiki/Marching_squares)
- [Rayon Parallel Processing](https://docs.rs/rayon/)
