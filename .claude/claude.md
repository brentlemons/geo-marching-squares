# geo-marching-squares

## Project Status

**Current Phase:** ALL PHASES COMPLETE ✅
**Last Updated:** 2025-10-05
**Last Commit:** `d2b9697` - Phase 7 complete (full isoline implementation)

See `IMPLEMENTATION_PLAN.md` for detailed phase tracking.

### Phase 1 Complete ✅
- ✅ Point, Edge, Shape, Cell structures implemented
- ✅ 22 tests passing
- ✅ Hash/Eq traits on Point for HashMap usage
- ✅ All constructors and getters working

### Phase 2 Complete ✅
- ✅ Shape::create() factory with ternary classification (all 81 configurations)
- ✅ All 7 shape types with edge construction (Triangle, Pentagon, Rectangle, Trapezoid, Hexagon, Saddle, Square)
- ✅ All 14 saddle configurations with center-average disambiguation
- ✅ Cosine interpolation with 0.999 centering hack
- ✅ GeoJSON feature extraction integrated
- ✅ 23 tests passing
- ✅ 1,901 lines (more concise than Java's 2,417)

### Phases 3-4 Complete ✅
- ✅ Interpolation merged into Phase 2
- ✅ GeoJSON integration merged into Phase 2

### Phase 5 Complete ✅
- ✅ process_band() main algorithm (grid-to-cells conversion, edge walking, polygon assembly)
- ✅ polygon_in_polygon() ray-casting for nesting resolution
- ✅ GeoJSON Feature output with MultiPolygon geometry
- ✅ Coordinate rounding (5 decimal places)
- ✅ 30 tests passing (25 unit + 5 integration)
- ✅ ~3,087 total lines

### Phase 6 Complete ✅
- ✅ do_concurrent() function with Rayon parallel processing
- ✅ has_coordinates() helper for filtering empty features
- ✅ Work-stealing thread pool (superior to Java's ExecutorService)
- ✅ Enhanced documentation with concurrent examples
- ✅ 38 tests passing (25 unit + 8 concurrent + 5 process_band)
- ✅ ~3,517 total lines

### Phase 7 Complete ✅
- ✅ Cell::create() with binary classification (16 configurations)
- ✅ Linear interpolation for isolines (no cosine smoothing)
- ✅ Line segment generation for all cases including saddles
- ✅ IsolineAssembler for segment-to-polyline assembly
- ✅ process_line() for single isoline generation
- ✅ do_concurrent_lines() for parallel isoline processing
- ✅ Support for closed loops and open segments
- ✅ 51 tests passing (28 unit + 8 concurrent + 5 process_band + 10 process_line)
- ✅ ~4,522 total lines

**Note:** Phase 7 goes beyond the Java reference (which had incomplete processLine())

---

## Project Overview

This Rust library implements the marching squares algorithm for generating contour polygons (isolines and isobands) from geospatial 2D scalar fields. It is a direct port of a proven Java implementation, maintaining algorithmic fidelity while applying Rust idioms.

### Key Features

- **Geographic coordinates:** Works with latitude/longitude (not pixel coordinates)
- **Isobands & Isolines:** Generate filled contour bands or contour lines
- **GeoJSON I/O:** Input and output GeoJSON features
- **Parallel processing:** Compute multiple isoband levels concurrently
- **Exact port:** Matches Java reference implementation behavior

---

## Quick Start for New Sessions

1. **Check current phase:** See status at top of this document
2. **Review implementation plan:** Open `IMPLEMENTATION_PLAN.md`
3. **Check last commits:** `git log --oneline -5`
4. **Review reference architecture:** See below
5. **Java reference code:** Located at `~/source/marching-squares-java`

---

## Project Structure

```
geo-marching-squares/
├── .claude/
│   └── claude.md           # This file - session continuity
├── src/
│   ├── lib.rs              # Library entry point
│   ├── point.rs            # Point data structure
│   ├── edge.rs             # Edge data structure
│   ├── shape.rs            # Shape enum and implementations
│   ├── cell.rs             # Cell for isolines
│   ├── marching_squares.rs # Main algorithm
│   └── error.rs            # Error types
├── tests/
│   ├── integration_tests.rs
│   └── test_data/          # Sample GeoJSON grids
├── Cargo.toml
├── IMPLEMENTATION_PLAN.md  # Phase-by-phase tracking
└── README.md               # Public-facing documentation
```

---

## Implementation Constraints

**CRITICAL:** This is a faithful port, not a rewrite.

- ✅ Algorithm flow MUST match Java exactly
- ✅ Interpolation formula MUST be identical (including the 0.999 hack)
- ✅ Edge walking direction logic MUST be preserved
- ✅ Saddle point disambiguation MUST use same average calculation
- ✅ Output precision MUST be 5 decimal places
- ✅ Polygon nesting algorithm MUST be same ray-casting approach
- ✅ No algorithmic "improvements" or optimizations that change behavior

---

## Reference Architecture

This document describes the reference architecture from the original Java implementation that this Rust library is porting.

## Overview

The marching squares algorithm generates contour polygons (isolines and isobands) from a 2D scalar field using geographic coordinates (latitude/longitude). Unlike typical pixel-based implementations, this version works directly with geospatial data.

## Data Structure

**Input:** 2D array of `Feature[][]` containing GeoJSON points
- Each feature has:
  - Geographic coordinates (longitude, latitude)
  - Scalar value property
- Grid forms cells from adjacent points
- Each cell defined by 4 corners: topLeft, topRight, bottomRight, bottomLeft

## Algorithm Flow

### 1. Ternary Classification (Isobands)

For each cell corner value, classify into one of three states:
- **0** = value < lower threshold
- **1** = lower ≤ value < upper threshold
- **2** = value ≥ upper threshold

This creates **81 possible cell configurations** (3^4 combinations)

### 2. Shape Creation

Factory pattern via `Shape.create()` returns specific shape types based on cell configuration:
- **Triangle** (8 configurations)
- **Pentagon** (24 configurations)
- **Rectangle** (12 configurations)
- **Trapezoid** (8 configurations)
- **Hexagon** (12 configurations)
- **Saddle** (14 configurations)
- **Square** (1 configuration)

Each shape determines its edge points by walking clockwise around the cell and checking which sides contain contour crossings.

### 3. Edge Point Determination

For each of the 4 cell sides (TOP, RIGHT, BOTTOM, LEFT):
- **isBlank()** check: Does the contour cross this side?
  - Side is blank if both endpoints are on same side of threshold
  - Side has crossing if endpoints straddle the threshold
- Creates 8 potential points (2 per side, for lower and upper thresholds)
- Filters to only points needed for this specific cell configuration

### 4. Interpolation

**Linear interpolation with cosine smoothing** finds precise contour crossing positions:

```java
mu = (level - value0) / (value1 - value0)
mu2 = (1.0 - cos(mu * PI)) / 2.0
newMu = 0.5 + ((mu2 - 0.5) * 0.999)  // Center adjustment
x = ((1.0 - newMu) * point0.x) + (newMu * point1.x)
y = ((1.0 - newMu) * point0.y) + (newMu * point1.y)
```

Applied to each side (TOP, RIGHT, BOTTOM, LEFT) where contour crosses.

### 5. Edge Creation and Clockwise Walking

Each cell creates **edges** that connect interpolated points:
- **Edge** structure:
  - `start`: Point where edge begins
  - `end`: Point where edge ends
  - `move`: Direction to next cell (RIGHT, DOWN, LEFT, UP, UNK)
- Edges stored in HashMap keyed by start point
- Each shape subclass defines edges by walking clockwise around its perimeter
- Edges at cell boundaries indicate which adjacent cell to visit next

### 6. Polygon Assembly

Main algorithm follows edges to form closed polygons:

1. Iterate through all cells in grid
2. For each uncleared cell, start edge walking:
   - Get edges from current cell starting at previous edge's end point
   - Follow the `move` direction to adjacent cell
   - Continue until loop closes (current end = first start)
3. Collect all edge points into polygon coordinate list
4. Handle polygon-in-polygon relationships:
   - Determine if polygons are nested (holes/interior rings)
   - Use point-in-polygon test
   - Exterior rings vs interior rings in GeoJSON Polygon

### 7. Saddle Point Disambiguation

Saddle configurations (where opposite corners have same classification) are ambiguous:
- Calculate cell center average: `avg = (tl + tr + br + bl) / 4.0`
- Compare average to thresholds to determine contour routing:
  - If `avg >= upper`: One routing pattern
  - If `lower <= avg < upper`: Different routing pattern
  - If `avg < lower`: Third routing pattern
- Determines whether contours connect diagonally or form separate loops

### 8. GeoJSON Output

Returns **GeoJSON FeatureCollection**:
- Each Feature contains a **MultiPolygon** geometry
- Coordinates in `[longitude, latitude]` format
- Properties include:
  - `lower_level`: Lower threshold value
  - `upper_level`: Upper threshold value
- Precision: 5 decimal places (approximately 1 meter accuracy)

## Key Characteristics

- **Embarrassingly parallel**: Each cell can be processed independently
- **Edge reuse**: Interpolation results can be cached for shared cell boundaries
- **Clockwise convention**: All edges walk clockwise to ensure consistent polygon winding
- **Grid independence**: Works with any rectangular grid of geographic points
- **Concurrent processing**: Multiple isoband levels can be computed in parallel

## Cell Coordinate System

Cells are indexed by `(row, column)` where:
- Row 0 is at the top (northernmost)
- Column 0 is at the left (westernmost)
- Cells created from grid points: `cell[r][c]` uses points from rows `r` to `r+1`, columns `c` to `c+1`

## Edge Movement Directions

When an edge exits a cell, the `move` direction indicates the next cell:
- **RIGHT**: `column++` (move east)
- **DOWN**: `row++` (move south)
- **LEFT**: `column--` (move west)
- **UP**: `row--` (move north)
- **UNK**: Unknown/error state (should not occur in valid configurations)
