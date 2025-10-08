# Recursive Grid Walking Algorithm

## Overview

The `walk_polygon_recursive()` function is the core boundary tracing algorithm in our marching squares implementation. Despite its name suggesting recursion, the current implementation is **iterative** - it follows polygon edges in a loop until returning to the starting point.

**Historical Note**: The function was originally designed to be recursive (using flood fill to find holes), but we pivoted to a simpler winding-based post-processing approach. The name remains for historical context.

## Function Signature

```rust
fn walk_polygon_recursive(
    cells: &mut Vec<Vec<Option<Shape>>>,      // Grid of marching squares cells
    _grid_data: &[Vec<GridCell>],             // Original data (unused in current impl)
    _lower: f64,                              // Lower threshold (unused)
    _upper: f64,                              // Upper threshold (unused)
    start_r: usize,                           // Starting row index
    start_c: usize,                           // Starting column index
    processed: &mut HashSet<(usize, usize)>,  // Already-processed cells
) -> Vec<Vec<Position>>                       // Returns single polygon ring
```

**Returns**: A vector containing exactly one ring (the exterior boundary) as geographic coordinates.

## Algorithm Overview

The algorithm performs **boundary tracing** through a marching squares grid:

1. Start at a cell with an unprocessed shape
2. Follow edges from cell to cell, tracking the polygon boundary
3. Stop when the loop closes (edge.end == first_edge.start)
4. Convert edges to geographic coordinates
5. Detect winding direction (clockwise vs counter-clockwise)
6. Mark all visited cells as processed

## Step-by-Step Walkthrough

### 1. Initialize Walking State

```rust
let mut boundary_cells = HashSet::new();     // Cells visited during walk
let mut edges = Vec::new();                  // Edges forming the polygon
let mut y = start_r;                         // Current row position
let mut x = start_c;                         // Current column position
let mut go_on = true;                        // Loop control flag
let mut current_edge: Option<Edge> = None;   // Last edge traversed
```

### 2. Main Edge Walking Loop

```rust
while go_on {
    // Boundary check
    if y >= rows || x >= cols {
        break;
    }

    // Mark this cell as part of the boundary
    boundary_cells.insert((y, x));

    // Get the cell's shape configuration
    let cell_ref = cells[y][x].as_mut()?;

    // Find next edges connected to current position
    let tmp_edges = if let Some(ref edge) = current_edge {
        cell_ref.get_edges(Some(edge.end()))  // Find edges starting at last edge's end
    } else {
        cell_ref.get_edges(None)              // First cell - any edge is fine
    };

    if tmp_edges.is_empty() {
        break;  // Dead end (shouldn't happen in valid grid)
    }

    // Mark edges as used (for multi-contour cells)
    cell_ref.increment_used_edges(tmp_edges.len());

    // Process each edge
    for edge in tmp_edges {
        cell_ref.remove_edge(edge.start());   // Remove from available edges
        current_edge = Some(edge.clone());
        edges.push(edge);

        // Check if loop is closed
        if current_edge.as_ref().unwrap().end() == edges[0].start() {
            go_on = false;
            break;
        }
    }

    if !go_on {
        break;
    }

    // Move to next cell based on edge direction
    if let Some(ref edge) = current_edge {
        match edge.move_direction() {
            Move::Right => x += 1,
            Move::Down  => y += 1,
            Move::Left  => if x > 0 { x -= 1 } else { go_on = false },
            Move::Up    => if y > 0 { y -= 1 } else { go_on = false },
            Move::Unknown => go_on = false,
        }
    }
}
```

### 3. Mark Processed Cells

```rust
for cell in &boundary_cells {
    processed.insert(*cell);
}
```

This prevents the main loop from re-walking the same polygon boundary.

### 4. Convert Edges to Geographic Coordinates

```rust
let mut exterior_ring: Vec<Position> = Vec::new();

if !edges.is_empty() {
    // Add starting point
    exterior_ring.push(vec![
        round_coord(edges[0].start().x().unwrap()),
        round_coord(edges[0].start().y().unwrap()),
    ]);

    // Add all edge endpoints
    for edge in &edges {
        exterior_ring.push(vec![
            round_coord(edge.end().x().unwrap()),
            round_coord(edge.end().y().unwrap()),
        ]);
    }
}
```

**Key Point**: Edges already contain geographic coordinates (lat/lon), not grid indices! The coordinate transformation happens during edge creation in the marching squares cell logic.

### 5. Detect Winding Direction (Shoelace Formula)

```rust
let winding = if exterior_ring.len() >= 3 {
    let mut signed_area = 0.0;
    let n = exterior_ring.len();

    // Sum of cross products: Σ(x[i] * y[i+1] - x[i+1] * y[i])
    for i in 0..n {
        let j = (i + 1) % n;
        signed_area += exterior_ring[i][0] * exterior_ring[j][1];
        signed_area -= exterior_ring[j][0] * exterior_ring[i][1];
    }

    if signed_area > 0.0 {
        "counter-clockwise"  // Hole in our grid orientation
    } else {
        "clockwise"          // Exterior ring or island
    }
} else {
    "unknown"
};
```

**Winding Direction Rules**:
- **Clockwise (CW)**: Exterior rings (main bodies) OR islands (polygons inside holes)
- **Counter-clockwise (CCW)**: Holes (polygons inside exterior rings)

This perfectly matches GeoJSON RFC 7946 winding requirements!

## Visual Example: Tracing a Simple Polygon

Consider a temperature band (15-20°C) in a 4×4 grid:

```
Grid with marching squares shapes:
   0   1   2   3
0  .   ◣  ━━  ◢
1  .   ┃  ⬤  ┃
2  .   ┃  ⬤  ┃
3  .   ◤  ━━  ◥

Legend:
. = no shape (outside threshold)
◣◢◤◥ = corners (single edge)
━━ = horizontal edges
┃ = vertical edges
⬤ = fully inside (4 edges)
```

### Walking Sequence

**Start**: Cell (0, 1) - shape ◣

1. **Cell (0,1)**: ◣ corner
   - Edge: `TopLeft → TopRight`
   - Direction: `Move::Right`
   - Next: cell (0, 2)

2. **Cell (0,2)**: ━━ horizontal
   - Edge: `LeftTop → RightTop`
   - Direction: `Move::Right`
   - Next: cell (0, 3)

3. **Cell (0,3)**: ◢ corner
   - Edge: `LeftTop → BottomRight`
   - Direction: `Move::Down`
   - Next: cell (1, 3)

4. **Cell (1,3)**: ┃ vertical
   - Edge: `TopRight → BottomRight`
   - Direction: `Move::Down`
   - Next: cell (2, 3)

5. **Cell (2,3)**: ┃ vertical
   - Edge: `TopRight → BottomRight`
   - Direction: `Move::Down`
   - Next: cell (3, 3)

6. **Cell (3,3)**: ◥ corner
   - Edge: `TopRight → BottomLeft`
   - Direction: `Move::Left`
   - Next: cell (3, 2)

7. **Cell (3,2)**: ━━ horizontal
   - Edge: `RightBottom → LeftBottom`
   - Direction: `Move::Left`
   - Next: cell (3, 1)

8. **Cell (3,1)**: ◤ corner
   - Edge: `RightBottom → TopLeft`
   - Direction: `Move::Up`
   - Next: cell (2, 1)

9. **Cell (2,1)**: ┃ vertical
   - Edge: `BottomLeft → TopLeft`
   - Direction: `Move::Up`
   - Next: cell (1, 1)

10. **Cell (1,1)**: ┃ vertical
    - Edge: `BottomLeft → TopLeft`
    - Direction: `Move::Up`
    - Next: cell (0, 1)

11. **Back at start**: Check closure
    - `current_edge.end() == edges[0].start()`
    - **YES!** Loop closed → `go_on = false`

**Result**: 10 edges forming a closed polygon, traced clockwise around the temperature band.

## Marching Squares Cell Shapes

The algorithm relies on pre-computed shapes in each grid cell, determined by which corners exceed the threshold:

| Shape | Binary | Description | Edges |
|-------|--------|-------------|-------|
| Empty | 0000 | All below threshold | None |
| ◣ | 0001 | Bottom-left only | 1 edge |
| ◤ | 0010 | Bottom-right only | 1 edge |
| ━━ | 0011 | Bottom edge | 1 edge |
| ◢ | 0100 | Top-right only | 1 edge |
| Saddle | 0101 | Diagonal (ambiguous) | 2 edges |
| ┃ | 0110 | Right edge | 1 edge |
| ... | ... | ... | ... |
| ⬤ | 1111 | All above threshold | 4 edges (interior) |

Each shape has pre-computed edges with geographic coordinate endpoints.

## Edge Direction and Cell Movement

Each edge knows which direction to move to the next cell:

```rust
pub enum Move {
    Right,   // Move to cell (r, c+1)
    Down,    // Move to cell (r+1, c)
    Left,    // Move to cell (r, c-1)
    Up,      // Move to cell (r-1, c)
    Unknown, // Error state
}
```

The edge's `.move_direction()` method returns the next cell direction based on which side of the cell the edge exits.

## Integration with Main Processing Loop

The walking function is called from the main band processing loop:

```rust
pub fn process_band_from_cells(
    data: &[Vec<GridCell>],
    lower: f64,
    upper: f64,
) -> Feature {
    // 1. Create marching squares grid
    let mut cells = create_cells(data, lower, upper);

    // 2. Find all unprocessed polygons
    let mut processed = HashSet::new();
    let mut raw_polygons = Vec::new();

    for r in 0..rows {
        for c in 0..cols {
            // Skip already-processed cells
            if processed.contains(&(r, c)) {
                continue;
            }

            // Skip empty cells
            if cells[r][c].is_none() {
                continue;
            }

            // Found unprocessed polygon - walk it!
            let polygon = walk_polygon_recursive(
                &mut cells,
                data,
                lower,
                upper,
                r,
                c,
                &mut processed,
            );

            raw_polygons.push(polygon);
        }
    }

    // 3. Build hierarchy using winding direction
    let hierarchical_polygons = build_polygon_hierarchy(raw_polygons);

    // 4. Convert to GeoJSON
    create_geojson_feature(hierarchical_polygons, lower, upper)
}
```

## Performance Characteristics

### Time Complexity

- **Per polygon**: O(E) where E = number of edges in the polygon
- **All polygons**: O(N) where N = total grid cells
- Each cell is visited at most once due to `processed` tracking

### Space Complexity

- **Boundary cells**: O(E) for the current polygon being walked
- **Edges**: O(E) for the current polygon
- **Processed set**: O(P) where P = cells touched by all polygons
- **Total**: O(N) in worst case

### Optimization Notes

1. **Early termination**: Loop exits immediately when closure detected
2. **HashSet lookups**: O(1) average case for `processed.contains()`
3. **In-place mutation**: Edges removed from cells to prevent reuse
4. **No backtracking**: Each cell visited exactly once per polygon

## Historical Context: Why "Recursive"?

### Original Design (Failed Approach)

The function was intended to work recursively with flood fill:

```rust
fn walk_polygon_recursive(...) -> Vec<Vec<Position>> {
    // 1. Walk exterior boundary
    let exterior = walk_edges(...);

    // 2. Flood fill interior to find unprocessed cells
    let interior_cells = flood_fill_interior(exterior);

    // 3. Recursively process holes
    let mut holes = vec![exterior];
    for (hole_r, hole_c) in interior_cells {
        if cells[hole_r][hole_c].is_some() {
            let hole = walk_polygon_recursive(...);  // RECURSION!
            holes.push(hole[0]);
        }
    }

    // 4. Return [exterior, hole1, hole2, ...]
    return holes;
}
```

**Why it failed**:
- Coordinate space mismatch: grid indices vs geographic coordinates
- Flood fill leakage: couldn't properly bound the interior
- Infinite loops: expanded to entire 1.9M cell grid
- See `MIDDLE_OUT_OPTIMIZATION_NOTES.md` for detailed history

### Current Approach (Successful)

We simplified to **boundary tracing only**, then use **winding direction** in post-processing:

```rust
// Each polygon walked independently
fn walk_polygon_recursive(...) -> Vec<Vec<Position>> {
    let exterior = walk_edges(...);
    return vec![exterior];  // Just one ring
}

// Holes matched to containers later
fn build_polygon_hierarchy(raw_polygons) -> Vec<Vec<Vec<Position>>> {
    // 1. Compute winding for each polygon
    // 2. CW = exterior/island, CCW = hole
    // 3. Use R-tree spatial index to find containers
    // 4. Build hierarchy: [exterior, hole1, hole2, ...]
}
```

This approach is:
- **Simpler**: No flood fill, no coordinate space issues
- **Faster**: O(N log N) with R-tree instead of O(N²)
- **Correct**: Winding direction perfectly classifies polygons
- **Maintainable**: Clear separation of concerns

## Winding Direction and GeoJSON Compliance

Our algorithm produces GeoJSON-compliant polygons per RFC 7946:

### Winding Rules

- **Exterior rings**: Clockwise (negative signed area)
- **Holes**: Counter-clockwise (positive signed area)
- **Islands** (polygons inside holes): Clockwise again

### Post-Processing Algorithm

```rust
fn build_polygon_hierarchy(raw_polygons) {
    // 1. Classify by winding
    for polygon in raw_polygons {
        let is_clockwise = compute_winding(polygon) == "clockwise";
        metas.push(PolygonMeta { ring: polygon, is_clockwise });
    }

    // 2. Build R-tree spatial index of CW polygons
    let rtree = RTree::bulk_load(cw_polygons);

    // 3. For each CCW hole, find CW container via R-tree query
    for hole in ccw_polygons {
        let container = rtree.locate_in_envelope(hole.bbox)
            .find(|cw| hole.is_inside(cw) && point_in_polygon(cw, hole[0]));
        hole_containers[hole.index] = container;
    }

    // 4. Assemble final structure
    for cw in root_cw_polygons {
        let polygon = vec![cw.exterior];
        polygon.extend(holes_for(cw));
        result.push(polygon);
    }
}
```

See `MIDDLE_OUT_OPTIMIZATION_NOTES.md` lines 183-294 for performance analysis.

## Key Takeaways

1. **Not truly recursive**: Despite the name, it's a simple iterative loop
2. **Boundary tracing**: Follows edges until loop closes
3. **Marching squares foundation**: Uses pre-computed cell shapes
4. **Geographic output**: Edges have lat/lon coordinates, not grid indices
5. **Winding detection**: Shoelace formula determines CW vs CCW
6. **Single ring output**: Returns only exterior, no holes
7. **Post-processing**: Holes matched later via R-tree spatial indexing
8. **Performance**: O(N) for all polygons, O(E) per polygon

The "recursive" name is a historical artifact from the failed flood-fill approach. The current implementation is elegantly simple: walk boundaries, detect winding, match holes in post-processing using spatial indexing.

## Related Documentation

- `MIDDLE_OUT_OPTIMIZATION_NOTES.md`: Detailed history of optimization attempts
- `src/marching_squares.rs`: Full implementation
- `src/edge.rs`: Edge structure and movement logic
- `src/shape.rs`: Cell shape configurations
- `src/cell.rs`: Cell structure and edge management
