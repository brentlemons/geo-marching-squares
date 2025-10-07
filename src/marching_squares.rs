//! Main marching squares algorithm implementation
//!
//! This module implements the core algorithm that converts a 2D grid of
//! GeoJSON features into contour polygons (isobands) and contour lines (isolines).

use crate::cell::Cell;
use crate::edge::{Edge, Move};
use crate::isoline_assembler::IsolineAssembler;
use crate::shape::Shape;
use geojson::{Feature, Geometry, JsonObject, Position, Value as GeoValue};
use std::collections::VecDeque;

/// Round a coordinate value to 5 decimal places (approximately 1 meter accuracy)
///
/// Matches Java's BigDecimal.setScale(5, RoundingMode.HALF_UP)
fn round_coord(value: f64) -> f64 {
    (value * 100000.0).round() / 100000.0
}

/// Bounding box for spatial optimization of polygon nesting
#[derive(Debug, Clone)]
struct BBox {
    min_x: f64,
    max_x: f64,
    min_y: f64,
    max_y: f64,
}

impl BBox {
    /// Compute bounding box from polygon ring (exterior ring)
    fn from_ring(ring: &[Position]) -> Self {
        let mut min_x = f64::INFINITY;
        let mut max_x = f64::NEG_INFINITY;
        let mut min_y = f64::INFINITY;
        let mut max_y = f64::NEG_INFINITY;

        for point in ring {
            min_x = min_x.min(point[0]);
            max_x = max_x.max(point[0]);
            min_y = min_y.min(point[1]);
            max_y = max_y.max(point[1]);
        }

        Self { min_x, max_x, min_y, max_y }
    }

    /// Check if this bbox is completely inside another bbox
    /// Fast check before expensive point-in-polygon test
    fn is_inside(&self, other: &BBox) -> bool {
        self.min_x >= other.min_x
            && self.max_x <= other.max_x
            && self.min_y >= other.min_y
            && self.max_y <= other.max_y
    }

    /// Check if two bboxes don't overlap at all
    /// Fast rejection before expensive polygon tests
    fn disjoint(&self, other: &BBox) -> bool {
        self.max_x < other.min_x
            || self.min_x > other.max_x
            || self.max_y < other.min_y
            || self.min_y > other.max_y
    }
}

/// Test if a polygon is completely contained within another polygon
///
/// Uses point-in-polygon ray-casting algorithm.
/// Reference: http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
///
/// Returns true if ALL points of subject are inside container.
fn polygon_in_polygon(subject: &[Vec<Position>], container: &[Vec<Position>]) -> bool {
    let container_points = &container[0]; // Exterior ring
    let subject_points = &subject[0]; // Exterior ring

    for subject_point in subject_points {
        let mut inside = false;
        let mut j = container_points.len() - 1;

        for i in 0..container_points.len() {
            let one = &container_points[i];
            let two = &container_points[j];

            // Ray casting: count edge crossings
            // Position format: [longitude, latitude] = [x, y]
            if ((one[1] > subject_point[1]) != (two[1] > subject_point[1]))
                && (subject_point[0]
                    < (two[0] - one[0]) * (subject_point[1] - one[1]) / (two[1] - one[1])
                        + one[0])
            {
                inside = !inside;
            }

            j = i;
        }

        if !inside {
            return false;
        }
    }

    true
}

/// Process a single isoband level and return a GeoJSON Feature with MultiPolygon geometry
///
/// This is the main marching squares algorithm that:
/// 1. Converts the 2D grid into cells (shapes)
/// 2. Walks edges to form closed polygon rings
/// 3. Resolves polygon nesting (exterior vs interior rings)
/// 4. Returns a Feature with MultiPolygon geometry and level properties
///
/// # Arguments
///
/// * `data` - 2D array of GeoJSON Point features with "value" property
/// * `lower` - Lower threshold for this isoband
/// * `upper` - Upper threshold for this isoband
///
/// # Returns
///
/// A GeoJSON Feature containing:
/// - MultiPolygon geometry with all contour polygons for this band
/// - Properties: "lower_level" and "upper_level"
pub fn process_band(data: &[Vec<Feature>], lower: f64, upper: f64) -> Feature {
    let rows = data.len();
    let cols = data[0].len();

    // Create cells from grid (rows-1 × cols-1)
    let mut cells: Vec<Vec<Option<Shape>>> = Vec::with_capacity(rows - 1);
    for _ in 0..rows - 1 {
        cells.push(vec![None; cols - 1]);
    }

    for r in 0..rows - 1 {
        for c in 0..cols - 1 {
            cells[r][c] = Shape::create(
                &data[r][c],
                &data[r][c + 1],
                &data[r + 1][c + 1],
                &data[r + 1][c],
                lower,
                upper,
                c,
                r,
                r == 0,
                c + 1 == cols - 1,
                r + 1 == rows - 1,
                c == 0,
            );
        }
    }

    let cell_rows = cells.len();
    let cell_cols = cells[0].len();

    let edge_walking_start = std::time::Instant::now();
    let mut hold_polygons: VecDeque<Vec<Vec<Position>>> = VecDeque::new();
    let mut cells_with_shapes = 0;
    let mut total_edges_found = 0;

    // Walk edges to form polygons
    for r in 0..cell_rows {
        for c in 0..cell_cols {
            if let Some(cell) = &mut cells[r][c] {
                cells_with_shapes += 1;
                if !cell.is_cleared() {
                    let mut y = r;
                    let mut x = c;
                    let mut go_on = true;
                    let mut edges = Vec::new();
                    let mut current_edge: Option<Edge> = None;

                    while go_on {
                        let cell_ref = cells[y][x].as_mut().unwrap();

                        let tmp_edges = if let Some(ref edge) = current_edge {
                            cell_ref.get_edges(Some(edge.end()))
                        } else {
                            cell_ref.get_edges(None)
                        };

                        if tmp_edges.is_empty() {
                            break;
                        }

                        cell_ref.increment_used_edges(tmp_edges.len());

                        for edge in tmp_edges {
                            cell_ref.remove_edge(edge.start());
                            current_edge = Some(edge.clone());
                            edges.push(edge);

                            // Check if loop closed
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
                                Move::Down => y += 1,
                                Move::Left => {
                                    if x > 0 {
                                        x -= 1
                                    } else {
                                        go_on = false;
                                    }
                                }
                                Move::Up => {
                                    if y > 0 {
                                        y -= 1
                                    } else {
                                        go_on = false;
                                    }
                                }
                                Move::Unknown => {
                                    go_on = false;
                                    eprintln!("Unknown edge move at ({}, {})", y, x);
                                }
                            }
                        }
                    }

                    // Convert edges to polygon coordinates
                    if !edges.is_empty() {
                        total_edges_found += edges.len();
                        let mut ring: Vec<Position> = Vec::new();

                        // Add first edge's start point
                        ring.push(vec![
                            round_coord(edges[0].start().x().unwrap()),
                            round_coord(edges[0].start().y().unwrap()),
                        ]);

                        // Add all end points
                        for edge in &edges {
                            ring.push(vec![
                                round_coord(edge.end().x().unwrap()),
                                round_coord(edge.end().y().unwrap()),
                            ]);
                        }

                        hold_polygons.push_back(vec![ring]);
                    }
                }
            }
        }
    }

    let edge_walking_elapsed = edge_walking_start.elapsed();
    eprintln!("[geo-marching-squares] ⏱️  Edge Walking: {:?} ({} shapes, {} total edges, {} polygons)",
        edge_walking_elapsed, cells_with_shapes, total_edges_found, hold_polygons.len());

    // Resolve polygon nesting (exterior vs interior rings) with spatial optimization
    let nesting_start = std::time::Instant::now();
    // Pre-compute bounding boxes for all polygons to avoid expensive polygon-in-polygon tests
    let mut polygons: Vec<Vec<Vec<Position>>> = Vec::new();
    let mut polygon_bboxes: Vec<BBox> = Vec::new();

    while let Some(subject) = hold_polygons.pop_front() {
        let subject_bbox = BBox::from_ring(&subject[0]);
        let mut external = true;

        for i in (0..polygons.len()).rev() {
            let polygon_bbox = &polygon_bboxes[i];

            // Fast rejection: if bboxes don't overlap, skip expensive polygon test
            if subject_bbox.disjoint(polygon_bbox) {
                continue;
            }

            let polygon = &polygons[i];
            let mut push_out = false;

            // Case 1: subject is inside polygon
            // Only test if bbox could be inside
            if subject_bbox.is_inside(polygon_bbox) && polygon_in_polygon(&subject, polygon) {
                // Check if subject is in a hole (should be pushed out)
                for hole in polygon.iter().skip(1) {
                    let hole_bbox = BBox::from_ring(hole);
                    if subject_bbox.is_inside(&hole_bbox) && polygon_in_polygon(&subject, std::slice::from_ref(hole)) {
                        push_out = true;
                        break;
                    }
                }

                if !push_out {
                    // Add subject as interior ring (hole) to polygon
                    // Avoid full clone by using Vec::with_capacity + extend
                    let mut new_polygon = Vec::with_capacity(polygon.len() + 1);
                    new_polygon.extend_from_slice(polygon);
                    new_polygon.push(subject[0].clone());
                    polygons[i] = new_polygon;
                    external = false;
                    break;
                }
            }
            // Case 2: polygon is inside subject
            // Only test if bbox could contain polygon
            else if polygon_bbox.is_inside(&subject_bbox) && polygon_in_polygon(polygon, &subject) {
                // Break apart polygon and reprocess
                if polygon.len() > 1 {
                    // Has interior rings - push them back for reprocessing
                    for hole in &polygon[1..] {
                        hold_polygons.push_back(vec![hole.clone()]);
                    }
                    hold_polygons.push_back(vec![polygon[0].clone()]);
                } else {
                    hold_polygons.push_back(polygon.clone());
                }
                // O(1) swap_remove instead of O(N) remove
                // Safe because we iterate in reverse - swapped element is beyond current index
                polygons.swap_remove(i);
                polygon_bboxes.swap_remove(i);
            }
        }

        if external {
            polygons.push(subject);
            polygon_bboxes.push(subject_bbox);
        }
    }
    let nesting_elapsed = nesting_start.elapsed();
    eprintln!("[geo-marching-squares] ⏱️  Polygon Nesting Resolution: {:?} ({} final polygons)",
        nesting_elapsed, polygons.len());

    // Build MultiPolygon geometry
    let multi_polygon = GeoValue::MultiPolygon(polygons);

    // Create Feature with properties
    let mut feature = Feature {
        bbox: None,
        geometry: Some(Geometry::new(multi_polygon)),
        id: None,
        properties: Some(JsonObject::new()),
        foreign_members: None,
    };

    if let Some(props) = feature.properties.as_mut() {
        props.insert("lower_level".to_string(), serde_json::json!(lower));
        props.insert("upper_level".to_string(), serde_json::json!(upper));
    }

    feature
}

/// Determine which cell to step into from the polygon boundary
/// Based on the last edge's move direction, step perpendicular inward
fn step_into_polygon(last_move: &Move, current_r: usize, current_c: usize, rows: usize, cols: usize) -> Option<(usize, usize)> {
    match last_move {
        Move::Right => {
            // Moving right, interior is below (down)
            if current_r + 1 < rows { Some((current_r + 1, current_c)) } else { None }
        }
        Move::Down => {
            // Moving down, interior is to the left
            if current_c > 0 { Some((current_r, current_c - 1)) } else { None }
        }
        Move::Left => {
            // Moving left, interior is above (up)
            if current_r > 0 { Some((current_r - 1, current_c)) } else { None }
        }
        Move::Up => {
            // Moving up, interior is to the right
            if current_c + 1 < cols { Some((current_r, current_c + 1)) } else { None }
        }
        Move::Unknown => None,
    }
}

/// Flood fill to find all cells inside a closed polygon boundary
/// Uses BFS starting from an interior cell with polygon containment check
fn flood_fill_interior(
    cells: &[Vec<Option<Shape>>],
    start: Option<(usize, usize)>,
    boundary_cells: &std::collections::HashSet<(usize, usize)>,
    polygon_ring: &[Position],
) -> std::collections::HashSet<(usize, usize)> {
    use std::collections::{HashSet, VecDeque};

    let mut interior = HashSet::new();
    let mut queue = VecDeque::new();

    if let Some(start_pos) = start {
        queue.push_back(start_pos);
        interior.insert(start_pos);
    }

    let rows = cells.len();
    let cols = if rows > 0 { cells[0].len() } else { 0 };

    while let Some((r, c)) = queue.pop_front() {
        // Add neighbors to queue (4-connected)
        let neighbors = [
            (r.wrapping_sub(1), c), // up
            (r + 1, c),              // down
            (r, c.wrapping_sub(1)),  // left
            (r, c + 1),              // right
        ];

        for (nr, nc) in neighbors {
            // Bounds check
            if nr >= rows || nc >= cols {
                continue;
            }

            // Already visited or on boundary
            if interior.contains(&(nr, nc)) || boundary_cells.contains(&(nr, nc)) {
                continue;
            }

            // Check if cell center is actually inside the polygon
            // Cell center is at (r + 0.5, c + 0.5)
            let cell_center = [nr as f64 + 0.5, nc as f64 + 0.5];
            if !point_in_ring(&cell_center, polygon_ring) {
                continue;
            }

            interior.insert((nr, nc));
            queue.push_back((nr, nc));
        }
    }

    interior
}

/// Check if a point is inside a polygon ring using ray casting
fn point_in_ring(point: &[f64], ring: &[Position]) -> bool {
    let mut inside = false;
    let mut j = ring.len() - 1;

    for i in 0..ring.len() {
        let one = &ring[i];
        let two = &ring[j];

        if ((one[1] > point[1]) != (two[1] > point[1]))
            && (point[0] < (two[0] - one[0]) * (point[1] - one[1]) / (two[1] - one[1]) + one[0])
        {
            inside = !inside;
        }

        j = i;
    }

    inside
}

/// Polygon orientation classification based on side-detection
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PolygonOrientation {
    Exterior,  // Inside band is outside the polygon (typical exterior ring)
    Hole,      // Inside band is inside the polygon (typical hole)
    Uncertain, // Cannot determine from side-detection alone
}

/// Check if a grid cell value is "inside band" (between lower and upper thresholds)
fn is_inside_band(value: f64, lower: f64, upper: f64) -> bool {
    value >= lower && value < upper
}

/// Perform side-detection for a single edge during polygon walking
/// Returns (high_side_score, low_side_score) based on which side has "inside band" values
///
/// Logic per edge direction:
/// - Right edge: high side is below (r+1), low side is above (r)
/// - Down edge: high side is left (c-1), low side is right (c)
/// - Left edge: high side is above (r), low side is below (r+1)
/// - Up edge: high side is right (c), low side is left (c-1)
fn detect_edge_side(
    grid_data: &[Vec<crate::GridCell>],
    lower: f64,
    upper: f64,
    cell_r: usize,
    cell_c: usize,
    move_dir: &Move,
) -> (i32, i32) {
    let grid_rows = grid_data.len();
    let grid_cols = if grid_rows > 0 { grid_data[0].len() } else { 0 };

    match move_dir {
        Move::Right => {
            // Moving right: check cell below (high side) and current cell row (low side)
            let high = if cell_r + 2 < grid_rows && cell_c + 1 < grid_cols {
                is_inside_band(grid_data[cell_r + 2][cell_c + 1].value, lower, upper)
            } else {
                false
            };
            let low = if cell_c + 1 < grid_cols {
                is_inside_band(grid_data[cell_r][cell_c + 1].value, lower, upper)
            } else {
                false
            };
            (if high { 1 } else { 0 }, if low { 1 } else { 0 })
        }
        Move::Down => {
            // Moving down: check cell to left (high side) and current cell column (low side)
            let high = if cell_r + 1 < grid_rows && cell_c > 0 {
                is_inside_band(grid_data[cell_r + 1][cell_c].value, lower, upper)
            } else {
                false
            };
            let low = if cell_r + 1 < grid_rows && cell_c + 1 < grid_cols {
                is_inside_band(grid_data[cell_r + 1][cell_c + 1].value, lower, upper)
            } else {
                false
            };
            (if high { 1 } else { 0 }, if low { 1 } else { 0 })
        }
        Move::Left => {
            // Moving left: check cell above (high side) and current cell row (low side)
            let high = if cell_r > 0 {
                is_inside_band(grid_data[cell_r][cell_c].value, lower, upper)
            } else {
                false
            };
            let low = if cell_r + 2 < grid_rows {
                is_inside_band(grid_data[cell_r + 2][cell_c].value, lower, upper)
            } else {
                false
            };
            (if high { 1 } else { 0 }, if low { 1 } else { 0 })
        }
        Move::Up => {
            // Moving up: check cell to right (high side) and current cell column (low side)
            let high = if cell_r > 0 && cell_c + 1 < grid_cols {
                is_inside_band(grid_data[cell_r][cell_c + 1].value, lower, upper)
            } else {
                false
            };
            let low = if cell_r > 0 && cell_c > 0 {
                is_inside_band(grid_data[cell_r][cell_c].value, lower, upper)
            } else {
                false
            };
            (if high { 1 } else { 0 }, if low { 1 } else { 0 })
        }
        Move::Unknown => (0, 0),
    }
}

/// Recursively walk a polygon and detect its orientation using side-detection
/// Returns a polygon structure: [exterior_ring, hole1, hole2, ...]
fn walk_polygon_recursive(
    cells: &mut Vec<Vec<Option<Shape>>>,
    grid_data: &[Vec<crate::GridCell>],
    lower: f64,
    upper: f64,
    start_r: usize,
    start_c: usize,
    processed: &mut std::collections::HashSet<(usize, usize)>,
) -> Vec<Vec<Position>> {
    use std::collections::HashSet;

    let rows = cells.len();
    let cols = if rows > 0 { cells[0].len() } else { 0 };

    // Track cells we visit while walking this polygon's boundary
    let mut boundary_cells = HashSet::new();
    let mut edges = Vec::new();
    let mut y = start_r;
    let mut x = start_c;
    let mut go_on = true;
    let mut current_edge: Option<Edge> = None;
    let mut last_move = Move::Unknown;

    // Side-detection: track which side of edges has "inside band" values
    let mut high_side_count = 0i32;  // Inside band on "high" side (right/below)
    let mut low_side_count = 0i32;   // Inside band on "low" side (left/above)

    // Walk edges to close the loop
    while go_on {
        if y >= rows || x >= cols {
            break;
        }

        boundary_cells.insert((y, x));

        let cell_ref = match cells[y][x].as_mut() {
            Some(cell) => cell,
            None => break,
        };

        let tmp_edges = if let Some(ref edge) = current_edge {
            cell_ref.get_edges(Some(edge.end()))
        } else {
            cell_ref.get_edges(None)
        };

        if tmp_edges.is_empty() {
            break;
        }

        cell_ref.increment_used_edges(tmp_edges.len());

        for edge in tmp_edges {
            cell_ref.remove_edge(edge.start());
            last_move = *edge.move_direction();

            // Side-detection: accumulate scores for each edge (before edge is moved)
            let (high, low) = detect_edge_side(grid_data, lower, upper, y, x, edge.move_direction());
            high_side_count += high;
            low_side_count += low;

            current_edge = Some(edge.clone());
            edges.push(edge);

            // Check if loop closed
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
                Move::Down => y += 1,
                Move::Left => {
                    if x > 0 {
                        x -= 1
                    } else {
                        go_on = false;
                    }
                }
                Move::Up => {
                    if y > 0 {
                        y -= 1
                    } else {
                        go_on = false;
                    }
                }
                Move::Unknown => {
                    go_on = false;
                }
            }
        }
    }

    // Mark boundary cells as processed
    for cell in &boundary_cells {
        processed.insert(*cell);
    }

    // Convert edges to exterior ring
    let mut exterior_ring: Vec<Position> = Vec::new();
    if !edges.is_empty() {
        exterior_ring.push(vec![
            round_coord(edges[0].start().x().unwrap()),
            round_coord(edges[0].start().y().unwrap()),
        ]);

        for edge in &edges {
            exterior_ring.push(vec![
                round_coord(edge.end().x().unwrap()),
                round_coord(edge.end().y().unwrap()),
            ]);
        }
    }

    // Detect winding direction using shoelace formula (signed area)
    let winding = if exterior_ring.len() >= 3 {
        let mut signed_area = 0.0;
        let n = exterior_ring.len();
        for i in 0..n {
            let j = (i + 1) % n;
            signed_area += exterior_ring[i][0] * exterior_ring[j][1];
            signed_area -= exterior_ring[j][0] * exterior_ring[i][1];
        }
        if signed_area > 0.0 {
            "counter-clockwise"
        } else {
            "clockwise"
        }
    } else {
        "unknown"
    };

    // Classify polygon orientation using side-detection
    let orientation = if high_side_count > low_side_count * 3 / 2 {
        PolygonOrientation::Exterior
    } else if low_side_count > high_side_count * 3 / 2 {
        PolygonOrientation::Hole
    } else {
        PolygonOrientation::Uncertain
    };

    let orientation_str = match orientation {
        PolygonOrientation::Exterior => "EXTERIOR",
        PolygonOrientation::Hole => "HOLE",
        PolygonOrientation::Uncertain => "UNCERTAIN",
    };

    eprintln!("[geo-marching-squares] Polygon at ({}, {}) winding: {}, orientation: {} (high={}, low={}, edges={})",
        start_r, start_c, winding, orientation_str, high_side_count, low_side_count, exterior_ring.len());

    // Step into polygon to find interior starting cell
    let interior_start = step_into_polygon(&last_move, y, x, rows, cols);

    // Flood fill to find all interior cells (with polygon containment check)
    let interior_cells = flood_fill_interior(cells, interior_start, &boundary_cells, &exterior_ring);

    // Debug: log interior cells found
    if !interior_cells.is_empty() {
        eprintln!("[geo-marching-squares] Found {} interior cells for polygon at ({}, {})",
            interior_cells.len(), start_r, start_c);
    }

    // Recursively process any holes (unprocessed cells with shapes inside)
    let mut holes = Vec::new();
    let mut holes_found = 0;
    for (r, c) in interior_cells {
        if !processed.contains(&(r, c)) {
            if let Some(cell) = &cells[r][c] {
                if !cell.is_cleared() {
                    // Found a hole! Recursively walk it
                    holes_found += 1;
                    eprintln!("[geo-marching-squares] Processing hole {} at ({}, {}) inside polygon ({}, {})",
                        holes_found, r, c, start_r, start_c);
                    let hole_polygon = walk_polygon_recursive(cells, grid_data, lower, upper, r, c, processed);
                    // Take only the exterior ring of the hole (index 0)
                    if !hole_polygon.is_empty() {
                        holes.push(hole_polygon[0].clone());
                    }
                }
            }
        }
    }

    if holes_found > 0 {
        eprintln!("[geo-marching-squares] Polygon at ({}, {}) has {} holes",
            start_r, start_c, holes.len());
    }

    // Build polygon: [exterior, hole1, hole2, ...]
    let mut polygon = vec![exterior_ring];
    polygon.extend(holes);
    polygon
}

/// Process a single isoband level from GridCell data (memory-efficient)
///
/// This is the high-performance alternative to `process_band()` for large grids.
/// Uses lightweight GridCell structs instead of full GeoJSON Features,
/// reducing memory usage by ~12x.
///
/// # Arguments
///
/// * `data` - 2D array of GridCell structs
/// * `lower` - Lower threshold for this isoband
/// * `upper` - Upper threshold for this isoband
///
/// # Returns
///
/// A GeoJSON Feature containing:
/// - MultiPolygon geometry with all contour polygons for this band
/// - Properties: "lower_level" and "upper_level"
///
/// # Example
///
/// ```rust,ignore
/// use geo_marching_squares::{process_band_from_cells, GridCell};
///
/// let grid: Vec<Vec<GridCell>> = load_grid_data();
/// let feature = process_band_from_cells(&grid, 10.0, 20.0);
/// ```
pub fn process_band_from_cells(data: &[Vec<crate::GridCell>], lower: f64, upper: f64) -> Feature {
    let band_start = std::time::Instant::now();
    let rows = data.len();
    let cols = data[0].len();

    eprintln!("[geo-marching-squares] process_band_from_cells: Starting band [{}, {}) on {}x{} grid",
        lower, upper, rows, cols);

    // Create cells from grid (rows-1 × cols-1)
    let mut cells: Vec<Vec<Option<Shape>>> = Vec::with_capacity(rows - 1);
    for _ in 0..rows - 1 {
        cells.push(vec![None; cols - 1]);
    }

    let classification_start = std::time::Instant::now();
    eprintln!("[geo-marching-squares] process_band_from_cells: Creating {} cells", (rows - 1) * (cols - 1));
    for r in 0..rows - 1 {
        for c in 0..cols - 1 {
            cells[r][c] = Shape::create_from_cells(
                &data[r][c],
                &data[r][c + 1],
                &data[r + 1][c + 1],
                &data[r + 1][c],
                lower,
                upper,
                c,
                r,
                r == 0,
                c + 1 == cols - 1,
                r + 1 == rows - 1,
                c == 0,
            );
        }
    }
    let classification_elapsed = classification_start.elapsed();

    let cell_rows = cells.len();
    let cell_cols = cells[0].len();

    eprintln!("[geo-marching-squares] ⏱️  Shape Classification: {:?}", classification_elapsed);
    eprintln!("[geo-marching-squares] process_band_from_cells: Cells created, now walking edges");

    // Recursive polygon walking with automatic hole detection
    let recursive_start = std::time::Instant::now();
    let mut polygons: Vec<Vec<Vec<Position>>> = Vec::new();
    let mut processed = std::collections::HashSet::new();

    eprintln!("[geo-marching-squares] process_band_from_cells: Walking edges with recursive hole detection");

    // Walk grid in order, recursively processing polygons and their holes
    for r in 0..cell_rows {
        for c in 0..cell_cols {
            if !processed.contains(&(r, c)) {
                if let Some(cell) = &cells[r][c] {
                    if !cell.is_cleared() {
                        // Found an unprocessed polygon - recursively walk it and its holes
                        let polygon = walk_polygon_recursive(&mut cells, data, lower, upper, r, c, &mut processed);
                        if !polygon.is_empty() && !polygon[0].is_empty() {
                            polygons.push(polygon);
                        }
                    }
                }
            }
        }
    }

    let recursive_elapsed = recursive_start.elapsed();
    eprintln!("[geo-marching-squares] ⏱️  Recursive Edge Walking + Hole Detection: {:?} ({} final polygons)",
        recursive_elapsed, polygons.len());

    eprintln!("[geo-marching-squares] process_band_from_cells: Polygon processing complete, building feature with {} polygons",
        polygons.len());

    // Build MultiPolygon geometry
    let geometry_start = std::time::Instant::now();
    let multi_polygon = GeoValue::MultiPolygon(polygons);

    // Create Feature with properties
    let mut feature = Feature {
        bbox: None,
        geometry: Some(Geometry::new(multi_polygon)),
        id: None,
        properties: Some(JsonObject::new()),
        foreign_members: None,
    };

    if let Some(props) = feature.properties.as_mut() {
        props.insert("lower_level".to_string(), serde_json::json!(lower));
        props.insert("upper_level".to_string(), serde_json::json!(upper));
    }
    let geometry_elapsed = geometry_start.elapsed();
    let total_elapsed = band_start.elapsed();

    eprintln!("[geo-marching-squares] ⏱️  Geometry Building: {:?}", geometry_elapsed);
    eprintln!("[geo-marching-squares] ⏱️  TOTAL BAND TIME: {:?}", total_elapsed);
    eprintln!("[geo-marching-squares] process_band_from_cells: Feature complete for band [{}, {})", lower, upper);

    feature
}

/// Check if a Feature has non-empty MultiPolygon geometry
///
/// Used to filter out empty results from isoband processing.
fn has_coordinates(feature: &Feature) -> bool {
    match &feature.geometry {
        Some(geometry) => match &geometry.value {
            GeoValue::MultiPolygon(polygons) => !polygons.is_empty(),
            _ => false,
        },
        None => false,
    }
}

/// Process multiple isoband levels concurrently using parallel processing
///
/// Takes a grid and a list of threshold values, and computes all isobands
/// in parallel using Rayon's work-stealing thread pool. This is the primary
/// API for generating multiple contour levels efficiently.
///
/// # Arguments
///
/// * `data` - 2D array of GeoJSON Point features with "value" property
/// * `isobands` - Sorted list of threshold values (e.g., [0.0, 10.0, 20.0, 30.0])
///               Creates N-1 isobands from N threshold values
///
/// # Returns
///
/// A GeoJSON FeatureCollection containing all isoband features.
/// Only features with non-empty geometry are included.
/// Features may not be in order (parallel processing is non-deterministic).
///
/// # Example
///
/// ```rust,ignore
/// use geo_marching_squares::do_concurrent;
///
/// let grid = load_grid_data(); // Your 2D grid
/// let thresholds = vec![0.0, 10.0, 20.0, 30.0];
///
/// // Creates 3 isobands: 0-10, 10-20, 20-30
/// let result = do_concurrent(&grid, &thresholds);
/// assert_eq!(result.features.len(), 3); // assuming all have data
/// ```
pub fn do_concurrent(
    data: &[Vec<Feature>],
    isobands: &[f64],
) -> geojson::FeatureCollection {
    use rayon::prelude::*;

    // Process each isoband pair in parallel
    let features: Vec<Feature> = (0..isobands.len() - 1)
        .into_par_iter()
        .map(|i| {
            // Each thread processes one isoband level
            process_band(data, isobands[i], isobands[i + 1])
        })
        .filter(|feature| {
            // Filter out empty features (no geometry)
            has_coordinates(feature)
        })
        .collect();

    // Build and return FeatureCollection
    geojson::FeatureCollection {
        bbox: None,
        foreign_members: None,
        features,
    }
}

/// Process multiple isoband levels from GridCell data concurrently (memory-efficient)
///
/// This is the high-performance alternative to `do_concurrent()` for large grids.
/// Uses lightweight GridCell structs instead of full GeoJSON Features,
/// reducing memory usage by ~12x.
///
/// # Arguments
///
/// * `data` - 2D array of GridCell structs
/// * `isobands` - Sorted list of threshold values (e.g., [0.0, 10.0, 20.0, 30.0])
///               Creates N-1 isobands from N threshold values
///
/// # Returns
///
/// A GeoJSON FeatureCollection containing all isoband features.
/// Only features with non-empty geometry are included.
///
/// # Memory Savings
///
/// For a 1799×1059 HRRR grid (1.9M points):
/// - GridCell array: ~46 MB
/// - Feature array: ~380-570 MB
/// - **Savings: ~12x reduction**
///
/// # Example
///
/// ```rust,ignore
/// use geo_marching_squares::{do_concurrent_from_cells, GridCell};
///
/// let grid: Vec<Vec<GridCell>> = load_grid_data();
/// let thresholds = vec![0.0, 10.0, 20.0, 30.0];
///
/// // Creates 3 isobands: 0-10, 10-20, 20-30
/// let result = do_concurrent_from_cells(&grid, &thresholds);
/// ```
pub fn do_concurrent_from_cells(
    data: &[Vec<crate::GridCell>],
    isobands: &[f64],
) -> geojson::FeatureCollection {
    use rayon::prelude::*;

    eprintln!("[geo-marching-squares] do_concurrent_from_cells: Starting with {} thresholds, {} bands to generate",
        isobands.len(), isobands.len().saturating_sub(1));
    eprintln!("[geo-marching-squares] Grid size: {}x{} = {} cells",
        data.len(),
        if data.is_empty() { 0 } else { data[0].len() },
        data.len() * if data.is_empty() { 0 } else { data[0].len() });

    // Process each isoband pair in parallel
    let features: Vec<Feature> = (0..isobands.len() - 1)
        .into_par_iter()
        .map(|i| {
            eprintln!("[geo-marching-squares] Band {}/{}: Processing [{}, {})",
                i + 1, isobands.len() - 1, isobands[i], isobands[i + 1]);
            // Each thread processes one isoband level
            let result = process_band_from_cells(data, isobands[i], isobands[i + 1]);
            eprintln!("[geo-marching-squares] Band {}/{}: Completed", i + 1, isobands.len() - 1);
            result
        })
        .filter(|feature| {
            // Filter out empty features (no geometry)
            has_coordinates(feature)
        })
        .collect();

    eprintln!("[geo-marching-squares] do_concurrent_from_cells: Completed, generated {} features", features.len());

    // Build and return FeatureCollection
    geojson::FeatureCollection {
        bbox: None,
        foreign_members: None,
        features,
    }
}

/// Check if a Feature has non-empty MultiLineString geometry
///
/// Used to filter out empty results from isoline processing.
fn has_line_coordinates(feature: &Feature) -> bool {
    match &feature.geometry {
        Some(geometry) => match &geometry.value {
            GeoValue::MultiLineString(lines) => !lines.is_empty(),
            _ => false,
        },
        None => false,
    }
}

/// Process a single isoline (contour line) level
///
/// Generates contour lines at the specified threshold value using binary
/// classification. Returns a Feature with MultiLineString geometry containing
/// all contour line segments for this level.
///
/// Unlike isobands which use ternary classification and cosine interpolation,
/// isolines use binary classification and linear interpolation.
///
/// # Arguments
///
/// * `data` - 2D array of GeoJSON Point features with "value" property
/// * `isovalue` - Threshold value for the contour line
///
/// # Returns
///
/// A GeoJSON Feature containing:
/// - MultiLineString geometry with all contour lines
/// - Property: "isovalue" with the threshold value
///
/// # Example
///
/// ```rust,ignore
/// use geo_marching_squares::process_line;
///
/// let grid = load_grid_data();
/// let feature = process_line(&grid, 15.0);
/// // Returns Feature with MultiLineString geometry
/// ```
pub fn process_line(data: &[Vec<Feature>], isovalue: f64) -> Feature {
    let rows = data.len();
    let cols = data[0].len();

    // Create cells from grid (binary classification)
    let mut cells: Vec<Vec<Option<Cell>>> = Vec::with_capacity(rows - 1);
    for _ in 0..rows - 1 {
        cells.push(vec![None; cols - 1]);
    }

    for r in 0..rows - 1 {
        for c in 0..cols - 1 {
            cells[r][c] = Cell::create(
                &data[r][c],
                &data[r][c + 1],
                &data[r + 1][c + 1],
                &data[r + 1][c],
                isovalue,
            );
        }
    }

    // Assemble line segments into polylines
    let mut assembler = IsolineAssembler::new();

    for r in 0..cells.len() {
        for c in 0..cells[0].len() {
            if let Some(cell) = &cells[r][c] {
                let segments = cell.get_line_segments();
                assembler.add_cell_segments(r, c, segments);
            }
        }
    }

    let polylines = assembler.assemble();

    // Build MultiLineString geometry
    let multi_linestring = GeoValue::MultiLineString(polylines);

    // Create Feature with properties
    let mut feature = Feature {
        bbox: None,
        geometry: Some(Geometry::new(multi_linestring)),
        id: None,
        properties: Some(JsonObject::new()),
        foreign_members: None,
    };

    if let Some(props) = feature.properties.as_mut() {
        props.insert("isovalue".to_string(), serde_json::json!(isovalue));
    }

    feature
}

/// Process a single isoline level from GridCell data (memory-efficient)
///
/// This is the high-performance alternative to `process_line()` for large grids.
/// Uses lightweight GridCell structs instead of full GeoJSON Features,
/// reducing memory usage by ~12x.
///
/// # Arguments
///
/// * `data` - 2D array of GridCell structs
/// * `isovalue` - Threshold value for the contour line
///
/// # Returns
///
/// A GeoJSON Feature containing:
/// - MultiLineString geometry with all contour lines
/// - Property: "isovalue" with the threshold value
///
/// # Example
///
/// ```rust,ignore
/// use geo_marching_squares::{process_line_from_cells, GridCell};
///
/// let grid: Vec<Vec<GridCell>> = load_grid_data();
/// let feature = process_line_from_cells(&grid, 15.0);
/// ```
pub fn process_line_from_cells(data: &[Vec<crate::GridCell>], isovalue: f64) -> Feature {
    let rows = data.len();
    let cols = data[0].len();

    // Create cells from grid (binary classification)
    let mut cells: Vec<Vec<Option<Cell>>> = Vec::with_capacity(rows - 1);
    for _ in 0..rows - 1 {
        cells.push(vec![None; cols - 1]);
    }

    for r in 0..rows - 1 {
        for c in 0..cols - 1 {
            cells[r][c] = Cell::create_from_cells(
                &data[r][c],
                &data[r][c + 1],
                &data[r + 1][c + 1],
                &data[r + 1][c],
                isovalue,
            );
        }
    }

    // Assemble line segments into polylines
    let mut assembler = IsolineAssembler::new();

    for r in 0..cells.len() {
        for c in 0..cells[0].len() {
            if let Some(cell) = &cells[r][c] {
                let segments = cell.get_line_segments();
                assembler.add_cell_segments(r, c, segments);
            }
        }
    }

    let polylines = assembler.assemble();

    // Build MultiLineString geometry
    let multi_linestring = GeoValue::MultiLineString(polylines);

    // Create Feature with properties
    let mut feature = Feature {
        bbox: None,
        geometry: Some(Geometry::new(multi_linestring)),
        id: None,
        properties: Some(JsonObject::new()),
        foreign_members: None,
    };

    if let Some(props) = feature.properties.as_mut() {
        props.insert("isovalue".to_string(), serde_json::json!(isovalue));
    }

    feature
}

/// Process multiple isoline levels concurrently using parallel processing
///
/// Takes a grid and a list of threshold values, and computes all isolines
/// in parallel using Rayon's work-stealing thread pool.
///
/// # Arguments
///
/// * `data` - 2D array of GeoJSON Point features with "value" property
/// * `isovalues` - List of threshold values for contour lines
///
/// # Returns
///
/// A GeoJSON FeatureCollection containing all isoline features.
/// Only features with non-empty geometry are included.
///
/// # Example
///
/// ```rust,ignore
/// use geo_marching_squares::do_concurrent_lines;
///
/// let grid = load_grid_data();
/// let isovalues = vec![0.0, 10.0, 20.0, 30.0];
///
/// // Creates 4 isolines at values 0, 10, 20, and 30
/// let result = do_concurrent_lines(&grid, &isovalues);
/// ```
pub fn do_concurrent_lines(
    data: &[Vec<Feature>],
    isovalues: &[f64],
) -> geojson::FeatureCollection {
    use rayon::prelude::*;

    let features: Vec<Feature> = isovalues
        .par_iter()
        .map(|&isovalue| process_line(data, isovalue))
        .filter(|feature| has_line_coordinates(feature))
        .collect();

    geojson::FeatureCollection {
        bbox: None,
        foreign_members: None,
        features,
    }
}

/// Process multiple isoline levels from GridCell data concurrently (memory-efficient)
///
/// This is the high-performance alternative to `do_concurrent_lines()` for large grids.
/// Uses lightweight GridCell structs instead of full GeoJSON Features,
/// reducing memory usage by ~12x.
///
/// # Arguments
///
/// * `data` - 2D array of GridCell structs
/// * `isovalues` - List of threshold values for contour lines
///
/// # Returns
///
/// A GeoJSON FeatureCollection containing all isoline features.
/// Only features with non-empty geometry are included.
///
/// # Memory Savings
///
/// For a 1799×1059 HRRR grid (1.9M points):
/// - GridCell array: ~46 MB
/// - Feature array: ~380-570 MB
/// - **Savings: ~12x reduction**
///
/// # Example
///
/// ```rust,ignore
/// use geo_marching_squares::{do_concurrent_lines_from_cells, GridCell};
///
/// let grid: Vec<Vec<GridCell>> = load_grid_data();
/// let isovalues = vec![0.0, 10.0, 20.0, 30.0];
///
/// // Creates 4 isolines at values 0, 10, 20, and 30
/// let result = do_concurrent_lines_from_cells(&grid, &isovalues);
/// ```
pub fn do_concurrent_lines_from_cells(
    data: &[Vec<crate::GridCell>],
    isovalues: &[f64],
) -> geojson::FeatureCollection {
    use rayon::prelude::*;

    let features: Vec<Feature> = isovalues
        .par_iter()
        .map(|&isovalue| process_line_from_cells(data, isovalue))
        .filter(|feature| has_line_coordinates(feature))
        .collect();

    geojson::FeatureCollection {
        bbox: None,
        foreign_members: None,
        features,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_round_coord() {
        assert_eq!(round_coord(1.234567), 1.23457);
        assert_eq!(round_coord(1.234564), 1.23456);
        assert_eq!(round_coord(-122.123456789), -122.12346);
        assert_eq!(round_coord(45.678901234), 45.6789);
    }

    #[test]
    fn test_polygon_in_polygon_simple() {
        // Square container: (0,0) -> (10,0) -> (10,10) -> (0,10)
        let container = vec![vec![
            vec![0.0, 0.0],
            vec![10.0, 0.0],
            vec![10.0, 10.0],
            vec![0.0, 10.0],
            vec![0.0, 0.0],
        ]];

        // Small square inside: (2,2) -> (8,2) -> (8,8) -> (2,8)
        let inside = vec![vec![
            vec![2.0, 2.0],
            vec![8.0, 2.0],
            vec![8.0, 8.0],
            vec![2.0, 8.0],
            vec![2.0, 2.0],
        ]];

        // Small square outside: (12,12) -> (18,12) -> (18,18) -> (12,18)
        let outside = vec![vec![
            vec![12.0, 12.0],
            vec![18.0, 12.0],
            vec![18.0, 18.0],
            vec![12.0, 18.0],
            vec![12.0, 12.0],
        ]];

        assert!(polygon_in_polygon(&inside, &container));
        assert!(!polygon_in_polygon(&outside, &container));
        assert!(!polygon_in_polygon(&container, &inside));
    }
}
