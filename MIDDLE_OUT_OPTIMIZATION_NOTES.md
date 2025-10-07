# Middle-Out Polygon Nesting Optimization

## Inspiration
Based on the "Pied Piper middle-out algorithm" concept - instead of post-processing all polygons with O(N²) comparisons, build the containment hierarchy during edge walking by recursively processing polygons "inside-out".

## The Problem We're Trying to Solve

### Baseline Performance (60761eb - LAST WORKING VERSION)
- **Algorithm**: Post-processing O(N²) polygon nesting
- **Performance**: 39.5s total for 16 bands, 86% spent in polygon nesting
- **Band 11 Example**: 13.6s nesting for 4,556 polygons
- **Approach**:
  1. Walk entire grid, collect all polygons as flat list (exterior rings only)
  2. Post-process: for each polygon, test against all previously processed polygons
  3. Use polygon_in_polygon() to determine containment
  4. Add as hole to container or mark as exterior
  5. Handle reprocessing when polygon is inside subject

**This version WORKS CORRECTLY but is slow due to O(N²) comparisons**

## The Recursive "Inside-Out" Approach

### Core Concept
When edge walking closes a polygon loop:
1. **Step into** the polygon from the boundary
2. **Flood fill** to find all cells inside the polygon
3. **Recursively process** any unprocessed cells with shapes (these are holes!)
4. Build polygon structure [exterior, hole1, hole2, ...] immediately
5. Mark all processed cells so we don't reprocess them

### Expected Benefits
- **Performance**: O(grid size) instead of O(N²) polygon comparisons
- **Natural hierarchy**: Holes discovered during parent polygon processing
- **Handles nesting**: Recursive calls handle donuts with islands automatically
- **No reprocessing**: Each cell processed exactly once

## Implementation History

### Commit 9b00692 - Initial Recursive Implementation
**What we built**:
- `step_into_polygon()`: Determine interior cell based on last edge direction
  - Right edge → step down (interior below)
  - Down edge → step left (interior to left)
  - Left edge → step up (interior above)
  - Up edge → step right (interior to right)
- `flood_fill_interior()`: BFS to find all cells inside boundary
- `walk_polygon_recursive()`: Main recursive function
- Replaced entire nesting algorithm in `process_band_from_cells()`

**Result**: INFINITE LOOP - flood fill expanded across entire grid

### Commit 13babce - Fix Infinite Expansion
**Problem identified**: Flood fill had no bounds - expanded to entire 1.9M cell grid

**Solution attempted**: Add point-in-polygon test
```rust
// Check if cell center is actually inside the polygon
let cell_center = [nr as f64 + 0.5, nc as f64 + 0.5];
if !point_in_ring(&cell_center, polygon_ring) {
    continue;  // Stop flood fill at polygon boundary
}
```

**Result**: Didn't crash, but found almost no holes (6 out of 2,538 polygons)

### Commit 193a327 - Debug Logging
**Added diagnostics**:
- Log interior cells found per polygon
- Log each hole being processed
- Log final hole count per polygon

**Purpose**: Understand why so few holes were detected

### Commit ed6196f - Remove Point-in-Polygon (BROKEN)
**Problem discovered**: Coordinate space mismatch!
- Testing grid indices `[10.5, 20.5]` (r+0.5, c+0.5)
- Against geographic coordinates `[-87.59587, 52.01399]` (lat/lon)
- These will NEVER match!

**Solution attempted**: Remove point-in-polygon test, trust boundary_cells

**Result**: INFINITE LOOP AGAIN - flood fill expands across entire grid

## The Fundamental Challenge

### Why Flood Fill Fails

**The coordinate space problem**:
- Polygon rings contain geographic coordinates (lat/lon after edge transformation)
- Cell indices are in grid space (row, col)
- Cannot use point-in-polygon to test grid cells against geographic polygon

**The boundary problem**:
- `boundary_cells` = cells visited during edge walking
- These cells don't form a complete perimeter in grid space
- A polygon spanning 100 cells might only touch 20 boundary cells during edge walking
- Flood fill leaks through gaps in the boundary

### What We Need

A way to determine "which grid cells are inside a polygon" that:
1. Works in grid coordinate space (row, col indices)
2. Properly bounds the flood fill
3. Doesn't require O(N²) comparisons
4. Handles the fact that edge coordinates are geographic but cells are grid-based

## Current State (193a327)

**Status**: Has recursive infrastructure + debug logging, but flood fill bounded by point-in-polygon

**Behavior**:
- Finds very few holes (6 out of 2,538 polygons)
- Due to coordinate mismatch in point-in-polygon test
- BUT doesn't crash or infinitely expand

**Use case**: Safe to run, produces valid (though incomplete) GeoJSON

## Next Steps to Explore

### Option 1: Fix Coordinate Space for Point-in-Polygon
- Keep polygon coordinates in grid space during nesting
- Only convert to geographic at final output
- Test `[r+0.5, c+0.5]` against grid-space polygon

### Option 2: Use Scanline in Grid Space
- Implement scanline algorithm to determine interior cells
- Works on grid indices, not geographic coordinates
- For each row, find entry/exit points of polygon
- Mark cells between entry/exit as interior

### Option 3: Track Interior Cells During Edge Walking
- As we walk polygon edges, mark which cells are "left" vs "right"
- Use winding direction to determine interior side
- Build interior cell list during edge walk

### Option 4: Revert to O(N²) Post-Processing
- Accept the performance cost
- Commit 60761eb works correctly with swap_remove optimization
- 13.6s for 4,556 polygons is acceptable for many use cases

## GeoJSON Winding Rules

**RFC 7946 Specification**:
- **Exterior rings**: Counter-clockwise (CCW) - right-hand rule
- **Holes**: Clockwise (CW)
- **Nested holes**: Alternate CCW → CW → CCW → CW

**Important**: Our recursive algorithm doesn't currently track depth to reverse winding direction for holes!

## Test Data

**Single Band Test** (test123.geojson):
- Total polygons: 2,538
- Polygons with holes: 6 (using 193a327)
- Max depth: 2 rings (1 exterior + 1 hole)
- Expected: Many more holes based on visual inspection

**Grid Size**: 1799 × 1059 = 1,904,541 cells

## Key Insights

1. **The recursive approach is theoretically sound** - if we can correctly identify interior cells
2. **Coordinate space mismatch** is the core blocker - edges are geographic, cells are grid-based
3. **Flood fill needs proper bounds** - boundary_cells alone insufficient
4. **60761eb is the last known working version** - O(N²) but correct
5. **Performance gain would be massive** - 13.6s → potentially <1s for nesting

## Files Modified

- `src/marching_squares.rs`: Main algorithm implementation
  - Lines 354-376: `step_into_polygon()`
  - Lines 378-431: `flood_fill_interior()`
  - Lines 433-452: `point_in_ring()`
  - Lines 454-608: `walk_polygon_recursive()`
  - Lines 634-663: Main loop using recursive approach

## Author Notes

The middle-out concept is brilliant and should work. The challenge is purely technical - finding a way to bound the flood fill in grid coordinate space. Once solved, this could reduce polygon nesting from O(N²) to O(grid size), a massive performance improvement for meteorological contour generation.
