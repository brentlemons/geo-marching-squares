# Implementation Plan: Java to Rust Port

## Status: Phase 2 - Shape Factory (Ready to Start)

This document tracks the progress of porting the Java marching squares implementation to Rust.

---

## Phase 1: Core Data Structures ✅ COMPLETED

### 1.1 Point Structure ✅
- [x] Create `Point` struct with:
  - [x] `x: Option<f64>`, `y: Option<f64>` for coordinates
  - [x] `value: Option<f64>`, `limit: Option<f64>` for threshold tracking
  - [x] `side: Option<Side>` enum (Top, Right, Bottom, Left)
- [x] Implement `PartialEq`, `Eq`, `Hash` for HashMap usage
- [x] Constructors matching Java's two patterns:
  - [x] `new(x: f64, y: f64)` - simple coordinate point
  - [x] `new_with_limit(value: f64, limit: f64, side: Side)` - interpolation marker

### 1.2 Edge Structure ✅
- [x] Create `Edge` struct with:
  - [x] `start: Point`
  - [x] `end: Point`
  - [x] `edge_type: Option<EdgeType>` enum (Outside, Inside)
  - [x] `move_direction: Move` enum (Right, Down, Left, Up, Unknown)
- [x] Simple constructors matching Java versions

### 1.3 Shape Hierarchy ✅
Replace Java inheritance with Rust enum + trait pattern:
- [x] `enum ShapeType` with variants:
  - [x] `Triangle`, `Pentagon`, `Rectangle`, `Trapezoid`, `Hexagon`, `Saddle`, `Square`
- [x] `struct Shape` containing:
  - [x] `shape_type: ShapeType`
  - [x] Common fields: `top_left`, `top_right`, `bottom_right`, `bottom_left` (Point)
  - [x] `value: u8` (ternary cell classification)
  - [x] `points: Vec<Point>` (computed edge points)
  - [x] `edges: HashMap<Point, Edge>`
  - [x] Grid position: `x: usize`, `y: usize`
  - [x] State tracking: `cleared: bool`, `used_edges: usize`
  - [x] Edge flags: `top_edge`, `right_edge`, `bottom_edge`, `left_edge` (bool)
  - [x] Corner values: `tl`, `tr`, `bl`, `br` (f64)
  - [x] Thresholds: `lower`, `upper` (f64)

### 1.4 Cell Structure (for isolines) ✅
- [x] Separate `Cell` struct for isoline algorithm
- [x] Similar structure to Shape but simpler (binary classification)
- [x] May consolidate with Shape later if patterns emerge

**Commit:** `18e930f` - feat: implement core data structures (Point, Edge, Shape, Cell)
**Tests:** 22 passing
**Completed:** 2025-10-05

---

## Phase 2: Shape Factory and Classification ⬜ NOT STARTED

### 2.1 Ternary Classification Function ⬜
```rust
fn classify_corner(value: f64, lower: f64, upper: f64) -> u8 {
    if value < lower { 0 }
    else if value >= upper { 2 }
    else { 1 }
}
```

### 2.2 Shape::create() Factory Function ⬜
- [ ] Static method that:
  - [ ] Takes 4 GeoJSON features + lower/upper thresholds + position + edge flags
  - [ ] Computes ternary value (0-170 range, base-3 encoding)
  - [ ] Matches value to determine ShapeType
  - [ ] Constructs appropriate Shape variant
  - [ ] Returns `Option<Shape>` (None for empty/full cells)

### 2.3 Shape-Specific Edge Construction ⬜
Start with Triangle, then expand:
- [ ] `build_triangle_edges()` - IMPLEMENT FIRST
- [ ] `build_pentagon_edges()`
- [ ] `build_rectangle_edges()`
- [ ] `build_trapezoid_edges()`
- [ ] `build_hexagon_edges()`
- [ ] `build_saddle_edges()` (includes center point disambiguation)
- [ ] `build_square_edges()`

Mirror exact switch/case logic from Java subclasses.

**Commit Message:** `feat: implement shape factory and classification`

---

## Phase 3: Interpolation ⬜ NOT STARTED

### 3.1 Side Blank Detection ⬜
```rust
impl Shape {
    fn is_top_blank(&self) -> bool
    fn is_right_blank(&self) -> bool
    fn is_bottom_blank(&self) -> bool
    fn is_left_blank(&self) -> bool
}
```

### 3.2 Interpolation Methods ⬜
```rust
impl Shape {
    fn interpolate(
        &self,
        level: f64,
        value0: f64,
        value1: f64,
        point0: &Point,
        point1: &Point
    ) -> Point {
        // Cosine interpolation with center adjustment
    }

    fn interpolate_side(&self, level: f64, side: Side) -> Point {
        // Delegates to interpolate() with appropriate corners
    }
}
```

### 3.3 Point Generation ⬜
```rust
impl Shape {
    fn get_points(&self) -> Vec<Point> {
        // Creates 8 potential points (2 per side)
        // Filters nulls, deduplicates
        // Interpolates points with None coordinates
    }
}
```

**Testing:** Test interpolation matches Java (same inputs → same outputs)

**Commit Message:** `feat: implement interpolation with cosine smoothing`

---

## Phase 4: GeoJSON Integration ⬜ NOT STARTED

### 4.1 Input Type Definition ⬜
Define Rust equivalent of `Feature[][]`:
- [ ] Use `geojson` crate (already established in Rust ecosystem)
- [ ] Create helper to extract:
  - [ ] Point coordinates from geometry
  - [ ] Value property from feature properties
- [ ] Input signature: `&[Vec<geojson::Feature>]` or custom wrapper

### 4.2 Coordinate Extraction ⬜
```rust
fn extract_point_data(feature: &geojson::Feature) -> Result<(f64, f64, f64), Error> {
    // Returns (longitude, latitude, value)
}
```

**Commit Message:** `feat: add GeoJSON integration and coordinate extraction`

---

## Phase 5: Main Algorithm - processBand() ⬜ NOT STARTED

### 5.1 Grid to Cells Conversion ⬜
```rust
fn process_band(
    data: &[Vec<geojson::Feature>],
    lower: f64,
    upper: f64
) -> geojson::Feature {
    let rows = data.len();
    let cols = data[0].len();
    let mut cells: Vec<Vec<Option<Shape>>> = vec![vec![None; cols-1]; rows-1];

    // Nested loop to create shapes
    for r in 0..rows-1 {
        for c in 0..cols-1 {
            cells[r][c] = Shape::create(
                &data[r][c],
                &data[r][c+1],
                &data[r+1][c+1],
                &data[r+1][c],
                lower, upper,
                c, r,
                r==0, c+1==cols-1, r+1==rows-1, c==0
            );
        }
    }
}
```

### 5.2 Edge Walking ⬜
Mirror Java's exact algorithm:
```rust
// For each cell, while edges remain:
let mut edges = Vec::new();
let mut current_edge: Option<Edge> = None;
let mut x = c;
let mut y = r;
let mut go_on = true;

while go_on && !cells[y][x].get_edges(...).is_empty() {
    // Get edges from current cell
    // Track current_edge
    // Check if loop closed
    // Update x,y based on move direction
}

// Collect edges into polygon coordinates
```

### 5.3 Polygon Nesting Resolution ⬜
```rust
fn polygon_in_polygon(subject: &Polygon, container: &Polygon) -> bool {
    // Point-in-polygon ray casting algorithm
    // Exact port of Java logic
}
```

Build exterior rings and interior rings (holes) using LinkedList/VecDeque pattern.

**Commit Message:** `feat: implement processBand main algorithm with edge walking`

---

## Phase 6: Concurrent Processing ⬜ NOT STARTED

### 6.1 Task Structure ⬜
```rust
struct BandTask {
    index: usize,
    lower: f64,
    upper: f64,
}

impl BandTask {
    fn execute(&self, data: &[Vec<geojson::Feature>]) -> Result<geojson::Feature, Error> {
        process_band(data, self.lower, self.upper)
    }
}
```

### 6.2 do_concurrent() Implementation ⬜
```rust
fn do_concurrent(
    data: &[Vec<geojson::Feature>],
    isobands: &[f64]
) -> geojson::FeatureCollection {
    use rayon::prelude::*;

    let features: Vec<geojson::Feature> = (0..isobands.len()-1)
        .into_par_iter()
        .map(|i| {
            process_band(data, isobands[i], isobands[i+1])
        })
        .filter(|f| /* has coordinates */)
        .collect();

    geojson::FeatureCollection {
        features,
        bbox: None,
        foreign_members: None,
    }
}
```

**Commit Message:** `feat: add concurrent processing with rayon`

---

## Phase 7: Isolines (process_line) ⬜ NOT STARTED

### 7.1 Binary Classification ⬜
Similar to process_band but simpler:
- [ ] Only 2 states (above/below threshold)
- [ ] 16 possible configurations instead of 81
- [ ] Simpler Cell structure
- [ ] Returns line geometries instead of polygons

Defer this until after isobands work correctly.

**Commit Message:** `feat: implement processLine for isolines`

---

## Phase 8: Testing Strategy

### 8.1 Unit Tests ⬜ ONGOING
- [ ] Test interpolation matches Java (same inputs → same outputs)
- [ ] Test shape classification (value → correct ShapeType)
- [ ] Test edge creation for each shape type
- [ ] Test point-in-polygon

### 8.2 Integration Tests ⬜ NOT STARTED
- [ ] Port/replicate any Java test data
- [ ] Compare GeoJSON output against Java version
- [ ] Coordinate precision (5 decimal places)
- [ ] Polygon winding order

### 8.3 Property-Based Testing ⬜ NOT STARTED
- [ ] Use `proptest` crate
- [ ] Generate random grids, verify:
  - [ ] No panics
  - [ ] All polygons closed
  - [ ] Valid GeoJSON output

**Note:** Tests should be added throughout implementation, not just at the end.

---

## Phase 9: API Design ⬜ NOT STARTED

### 9.1 Public API ⬜
```rust
pub struct MarchingSquares;

impl MarchingSquares {
    pub fn process_band(
        data: &[Vec<geojson::Feature>],
        lower: f64,
        upper: f64
    ) -> Result<geojson::Feature, Error>;

    pub fn process_line(
        data: &[Vec<geojson::Feature>],
        isovalue: f64
    ) -> Result<geojson::Feature, Error>;

    pub fn do_concurrent(
        data: &[Vec<geojson::Feature>],
        isobands: &[f64]
    ) -> Result<geojson::FeatureCollection, Error>;
}
```

### 9.2 Error Handling ⬜
```rust
#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error("Invalid grid dimensions")]
    InvalidGrid,
    #[error("Missing value property in feature")]
    MissingValue,
    #[error("Invalid geometry type")]
    InvalidGeometry,
    // etc.
}
```

**Commit Message:** `feat: finalize public API and error types`

---

## Phase 10: Dependencies ⬜ NOT STARTED

### 10.1 Cargo.toml ⬜
```toml
[dependencies]
geojson = "0.24"          # GeoJSON types
serde = { version = "1.0", features = ["derive"] }
serde_json = "1.0"
rayon = "1.10"             # Parallel processing
thiserror = "1.0"          # Error derive macros

[dev-dependencies]
approx = "0.5"             # Floating point comparison in tests
proptest = "1.0"           # Property-based testing
```

**Commit Message:** `chore: add project dependencies`

---

## Implementation Order

1. ✅ **Phase 10** - Dependencies (Cargo.toml) - DO THIS FIRST
2. **Phase 1** - Data structures (Point, Edge, Shape enum)
3. **Phase 2** - Shape factory (just Triangle first to validate approach)
4. **Phase 3** - Interpolation (test against known Java outputs)
5. **Phase 2 (complete)** - All shape types
6. **Phase 4** - GeoJSON integration
7. **Phase 5** - Main algorithm (process_band)
8. **Phase 8.1** - Unit tests as we go
9. **Phase 8.2** - Integration test with real data
10. **Phase 6** - Concurrent processing
11. **Phase 7** - Isolines (if needed)
12. **Phase 9** - Finalize public API
13. **Phase 8.3** - Property-based tests

---

## Key Rust Idioms to Apply

1. **Ownership**: Pass features by reference, return owned GeoJSON
2. **Error handling**: `Result<T, Error>` instead of exceptions
3. **Options**: Use `Option<T>` for nullable values
4. **Iterators**: Prefer iterator chains over explicit loops where clear
5. **Pattern matching**: Replace Java switch/case with Rust match
6. **Traits**: Implement common traits (Debug, Clone, PartialEq, etc.)
7. **Modules**: Organize as `src/point.rs`, `src/edge.rs`, `src/shape.rs`, etc.
8. **Documentation**: /// doc comments on all public items

---

## Non-Negotiable Constraints

- ✅ Algorithm flow MUST match Java exactly
- ✅ Interpolation formula MUST be identical (including the 0.999 hack)
- ✅ Edge walking direction logic MUST be preserved
- ✅ Saddle point disambiguation MUST use same average calculation
- ✅ Output precision MUST be 5 decimal places
- ✅ Polygon nesting algorithm MUST be same ray-casting approach
- ✅ No algorithmic "improvements" or optimizations that change behavior

---

## Session Handoff Checklist

Before ending a session and starting a new one:

- [ ] Update this document with ✅ for completed items
- [ ] Update phase status (⬜ → 🔄 → ✅)
- [ ] Commit all changes with appropriate message
- [ ] Push to GitHub
- [ ] Update `.claude/claude.md` with current progress
- [ ] Note any blockers or decisions needed in next session
