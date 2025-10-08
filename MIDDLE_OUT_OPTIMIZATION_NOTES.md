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

---

# SOLUTION FOUND! 🎉

## Final Solution (Commits 731a7eb + cef1c4b)

After extensive experimentation, we discovered that **winding direction is the perfect classifier** for polygon types, eliminating the need for flood fill entirely.

### The Breakthrough

**Experimental Discovery**:
- Commit 17bb4bf: Filtered to show only clockwise polygons → showed main bodies + islands
- Commit fffdcc6: Filtered to show only counter-clockwise polygons → showed hundreds of holes
- **Conclusion**: Winding direction perfectly identifies polygon type!

**Winding Rule**:
- **Clockwise (CW)** = Exterior rings (main bodies) OR islands (inside holes)
- **Counter-clockwise (CCW)** = Holes (inside exterior rings)

This matches GeoJSON RFC 7946 winding requirements exactly.

### Algorithm: Winding-Based Hierarchy with Bbox Optimization

**Commit 731a7eb - Core Implementation**:

1. **Simplified Edge Walking**:
   - Walk all polygons, return just exterior rings (no flood fill!)
   - Compute winding direction during walk (shoelace formula)
   - ~2,545 raw polygon rings collected

2. **Post-Processing Containment Detection**:
   ```rust
   fn build_polygon_hierarchy(raw_polygons) -> hierarchical_polygons {
       // Compute bbox + winding for each polygon
       for each polygon:
           bbox = BBox::from_ring(ring)
           is_clockwise = compute_winding(ring) == "clockwise"

       // Find containers for CCW polygons (holes)
       for each CCW polygon:
           find smallest CW polygon that contains it (bbox + point-in-polygon)

       // Find containers for CW polygons (islands in holes)
       for each CW polygon:
           check if inside any CCW polygon (island in hole)

       // Build final structure
       root_level_CW_polygons.attach_holes()
       islands.create_separate_polygons()
   }
   ```

3. **Bbox Pre-Filtering**:
   - Fast bbox containment check before expensive point-in-polygon
   - Eliminates most comparisons (~90%+ in practice)
   - Reduces effective complexity from O(N²) to O(N log N) typical case

### Performance Results

**Before (60761eb - O(N²) post-processing)**:
- Total time: 39.5 seconds for 16 bands
- 86% in polygon nesting: ~34 seconds
- Band 11: 13.6s nesting for 4,556 polygons

**After (731a7eb - Winding + bbox optimization)**:
- Total time: 12.5 seconds for 16 bands
- **68% faster overall!**
- **3x faster** than original O(N²) approach
- Correctly handles hundreds of holes (not just 6!)

### Why It Works

1. **Winding direction is deterministic**: Marching squares generates consistent winding based on cell traversal
2. **No coordinate space issues**: Winding computed in geographic space after edge transformation
3. **Bbox optimization**: Most polygons spatially disjoint, bbox eliminates them quickly
4. **Correct hierarchy**: Handles nested structures (holes in bodies, islands in holes, holes in islands)

### Failed Approaches (Learning Notes)

1. **Flood Fill (Commits 9b00692 - ed6196f)**:
   - Problem: Coordinate space mismatch (grid indices vs geographic coordinates)
   - Result: Infinite expansion or false positives

2. **Side-Detection (Commit c6c1a1c)**:
   - Problem: Complex edge orientation logic, wrong assumptions
   - Result: All polygons classified as EXTERIOR

3. **Middle-Out Recursive**:
   - Concept was sound but implementation path was wrong
   - True "middle-out" = winding-based post-processing, not flood fill!

### Code Cleanup (Commit cef1c4b)

Removed all dead code from failed attempts:
- `step_into_polygon()` - flood fill helper
- `flood_fill_interior()` - broken coordinate space logic
- `point_in_ring()` - old implementation
- `detect_edge_side()` - side-detection logic
- Simplified `walk_polygon_recursive()` to just edge walking

Result: -232 lines, zero warnings, cleaner codebase

### Key Commits

- `60761eb`: Last working O(N²) version (baseline: 39.5s)
- `c6c1a1c`: Side-detection attempt (failed)
- `17bb4bf`: Clockwise-only filter experiment (breakthrough!)
- `fffdcc6`: Counter-clockwise-only filter (confirmation!)
- `731a7eb`: Winding-based hierarchy (final solution: 12.5s)
- `cef1c4b`: Code cleanup

### The Lesson

Sometimes the "middle-out" solution isn't about changing the algorithm structure, but about finding the right classification metric. Winding direction was there all along - we just needed to test it experimentally to discover its power.

---

# R-tree Spatial Indexing Attempt (Commit a3f46e4 - REVERTED)

## Motivation (2025-10-08)

After achieving 68% speedup with winding-based hierarchy, profiling showed remaining bottleneck:
- Hole Containment: 249ms for 1,682 point-in-poly checks
- Point-in-polygon tests are ~1000x slower than bbox checks (~148μs vs ~90ns)
- Early exit optimization (98bc26a) reduced bbox checks by 91% but didn't improve wall time
- Attempted R-tree spatial indexing to reduce candidates checked

## Implementation (Commit a3f46e4)

Added `rstar` crate for spatial indexing:

```rust
// Build R-tree indexes
let cw_rtree: RTree<SpatialEntry> = RTree::bulk_load(cw_polygons);
let ccw_rtree: RTree<SpatialEntry> = RTree::bulk_load(ccw_polygons);

// Query for candidates
for hole in ccw_polygons {
    let query_envelope = AABB::from_corners(
        [hole.bbox.min_x, hole.bbox.min_y],
        [hole.bbox.max_x, hole.bbox.max_y],
    );

    for candidate in rtree.locate_in_envelope(&query_envelope) {
        // Check bbox containment + point-in-polygon
    }
}
```

## Results - FAILED ❌

**Log stream**: `api/grib-inspector-api/368568d03cfe46f3a81d740a2013b978`

### Critical Failure: Zero Point-in-Polygon Tests

```
Metadata (bbox+winding): 1.358ms for 2502 polygons (1020 CW, 1482 CCW)
R-tree Build: 402μs (1020 CW polygons indexed)
Hole Containment: 2.634ms (1482 R-tree queries, 116 candidates checked, 0 point-in-poly checks)
CCW R-tree Build: 356μs (1482 CCW polygons indexed)
Island Containment: 201μs (1020 R-tree queries, 1640 candidates checked, 0 point-in-poly checks)
```

**The Bug**:
- 1,482 holes queried R-tree → 116 candidates returned
- **0 point-in-polygon tests performed** (should be ~1,682 based on previous runs)
- Result: No holes matched to containers (hundreds of "no container" warnings)
- GeoJSON output had no nested polygons (holes rendered as separate polygons)

### Root Cause Analysis

**Problem**: `locate_in_envelope()` returns polygons whose bboxes **overlap** the query envelope, not polygons that **contain** it.

```rust
// Query with hole's bbox
let query_envelope = AABB::from_corners(hole.bbox.min, hole.bbox.max);

// R-tree returns CW polygons that OVERLAP hole's bbox
for candidate in rtree.locate_in_envelope(&query_envelope) {
    // But we check: does hole.bbox fit INSIDE candidate.bbox?
    if !hole.bbox.is_inside(&candidate.bbox) {
        continue;  // ← ALL candidates fail this check!
    }
}
```

**Why it fails**:
- Query envelope = hole's exact bbox
- R-tree returns: CW polygons whose bboxes overlap/intersect
- Overlap ≠ Containment
- Bbox containment check fails for all candidates
- No point-in-polygon tests ever run
- No holes matched

**What we needed**:
- Query method that finds "bboxes that contain this bbox"
- OR query with a point instead of envelope
- OR different R-tree query strategy

### Performance Comparison

| Metric | Early Exit (98bc26a) | R-tree (a3f46e4) |
|--------|---------------------|------------------|
| Hole Containment | 249ms | 2.6ms ⚠️ |
| Bbox checks | 235K | 116 candidates |
| Point-in-poly checks | 1,682 | **0** ❌ |
| Holes matched | Correct | **ZERO** ❌ |
| GeoJSON output | Valid | **Broken** ❌ |

**Conclusion**: R-tree was faster but **completely broken** - matched zero holes due to query semantics mismatch.

## Why Previous Spatial Index Attempt Failed (Historical)

User reported: "we tried spatial indexing before and it was 2.7x slower"

Likely reasons:
1. R-tree construction overhead for small polygon counts
2. Query overhead exceeding linear scan for certain distributions
3. Similar containment query issues
4. Poor cache locality vs sequential iteration

## Decision: Revert to 98bc26a

**Reverted** R-tree implementation (commit a3f46e4) back to early exit optimization (98bc26a):
- Keeps working hole/island detection
- Maintains ~12.5s total performance (68% improvement over baseline)
- Simple O(N²) with early exit is correct and fast enough
- Point-in-polygon tests still expensive but unavoidable

**Files kept**: `RECURSIVE_WALKING_ALGORITHM.md` - comprehensive walking algorithm documentation

## Future Optimization Ideas (If Needed)

1. **Faster point-in-polygon algorithm**:
   - Current: Ray casting (~148μs per test)
   - Consider: Winding number algorithm
   - Or: Pre-compute additional spatial hints

2. **Better R-tree query**:
   - Use point containment query instead of envelope query
   - Query with hole's first coordinate, find all CW bboxes containing it
   - Would need custom R-tree query predicate

3. **Accept current performance**:
   - 12.5s for 16 bands processing 4,556 polygons
   - 68% faster than baseline (39.5s)
   - Good enough for meteorological contour generation

4. **Optimize point-in-polygon itself**:
   - Most expensive operation at ~148μs per test
   - 1,682 tests × 148μs = ~249ms (matches profiling!)
   - 1000x slower than bbox checks
   - Direct optimization target

## Key Takeaway

**Spatial indexing can be fast but wrong**. In this case:
- R-tree reduced candidates from 235K to 116 (99.95% reduction!)
- But query semantics mismatch caused zero matches
- Fast broken code is worse than slow correct code
- Simple early exit O(N²) is good enough and correct
