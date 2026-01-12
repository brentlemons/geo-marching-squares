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

/// Default coordinate precision (5 decimal places, ~1.1m at equator)
pub const DEFAULT_PRECISION: u32 = 5;

/// Round a coordinate value to specified decimal places
///
/// # Arguments
///
/// * `value` - The coordinate value to round
/// * `precision` - Number of decimal places (1-6 recommended)
///   - 1 = ~11 km, 2 = ~1.1 km, 3 = ~110 m, 4 = ~11 m, 5 = ~1.1 m, 6 = ~0.11 m
///
/// Matches Java's BigDecimal.setScale(n, RoundingMode.HALF_UP)
fn round_coord_with_precision(value: f64, precision: u32) -> f64 {
    let factor = 10_f64.powi(precision as i32);
    (value * factor).round() / factor
}

/// Round a coordinate value to 5 decimal places (approximately 1 meter accuracy)
///
/// Matches Java's BigDecimal.setScale(5, RoundingMode.HALF_UP)
fn round_coord(value: f64) -> f64 {
    round_coord_with_precision(value, DEFAULT_PRECISION)
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

        Self {
            min_x,
            max_x,
            min_y,
            max_y,
        }
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
            // Position format: [x, y]
            if ((one[1] > subject_point[1]) != (two[1] > subject_point[1]))
                && (subject_point[0]
                    < (two[0] - one[0]) * (subject_point[1] - one[1]) / (two[1] - one[1]) + one[0])
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

    let mut hold_polygons: VecDeque<Vec<Vec<Position>>> = VecDeque::new();

    // Walk edges to form polygons
    for r in 0..cell_rows {
        for c in 0..cell_cols {
            if let Some(cell) = &mut cells[r][c] {
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
                                }
                            }
                        }
                    }

                    // Convert edges to polygon coordinates
                    if !edges.is_empty() {
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

    // Resolve polygon nesting (exterior vs interior rings) with spatial optimization
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
                    if subject_bbox.is_inside(&hole_bbox)
                        && polygon_in_polygon(&subject, std::slice::from_ref(hole))
                    {
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
            else if polygon_bbox.is_inside(&subject_bbox) && polygon_in_polygon(polygon, &subject)
            {
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

    // Orient polygons to follow RFC 7946 right-hand rule
    let oriented_polygons = orient_polygons(polygons);

    // Build MultiPolygon geometry
    let multi_polygon = GeoValue::MultiPolygon(oriented_polygons);

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

/// Walk a polygon starting from a cell and return its exterior ring
/// Containment hierarchy will be determined in post-processing based on winding direction
fn walk_polygon_recursive_with_precision(
    cells: &mut Vec<Vec<Option<Shape>>>,
    _grid_data: &[Vec<crate::GridCell>],
    _lower: f64,
    _upper: f64,
    start_r: usize,
    start_c: usize,
    processed: &mut std::collections::HashSet<(usize, usize)>,
    precision: u32,
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
            round_coord_with_precision(edges[0].start().x().unwrap(), precision),
            round_coord_with_precision(edges[0].start().y().unwrap(), precision),
        ]);

        for edge in &edges {
            exterior_ring.push(vec![
                round_coord_with_precision(edge.end().x().unwrap(), precision),
                round_coord_with_precision(edge.end().y().unwrap(), precision),
            ]);
        }
    }

    vec![exterior_ring]
}

/// Metadata for a polygon during hierarchy building
#[derive(Debug, Clone)]
struct PolygonMeta {
    ring: Vec<Position>,
    bbox: BBox,
    is_clockwise: bool,
}

/// Build polygon hierarchy using winding direction and bbox-optimized containment detection
///
/// Takes a flat list of polygon rings and organizes them into proper GeoJSON MultiPolygon structure:
/// - Clockwise rings (exterior/islands) become separate polygons or nested within holes
/// - Counter-clockwise rings (holes) become holes in their containing clockwise polygon
///
/// Algorithm:
/// 1. Compute bbox and winding for each polygon
/// 2. For each counter-clockwise polygon (hole), find its clockwise container
/// 3. For each clockwise polygon at root level, attach its holes
/// 4. Nested clockwise polygons (islands in holes) become separate polygons
pub fn build_polygon_hierarchy(raw_polygons: Vec<Vec<Vec<Position>>>) -> Vec<Vec<Vec<Position>>> {
    if raw_polygons.is_empty() {
        return vec![];
    }

    // Step 1: Build metadata for all polygons (just exterior rings)
    let metas: Vec<PolygonMeta> = raw_polygons
        .into_iter()
        .filter_map(|mut poly| {
            if poly.is_empty() || poly[0].is_empty() {
                return None;
            }
            let ring = poly.remove(0);
            let bbox = BBox::from_ring(&ring);
            let is_clockwise = compute_winding(&ring) == "clockwise";
            Some(PolygonMeta {
                ring,
                bbox,
                is_clockwise,
            })
        })
        .collect();

    // Step 2: Find container for each counter-clockwise polygon (holes)
    let mut hole_containers: Vec<Option<usize>> = vec![None; metas.len()];

    for i in 0..metas.len() {
        if !metas[i].is_clockwise {
            for j in 0..metas.len() {
                if i == j || !metas[j].is_clockwise {
                    continue;
                }
                if !metas[i].bbox.is_inside(&metas[j].bbox) {
                    continue;
                }
                if !metas[i].ring.is_empty()
                    && polygon_contains_point(&metas[j].ring, &metas[i].ring[0])
                {
                    hole_containers[i] = Some(j);
                    break;
                }
            }
        }
    }

    // Step 3: For each clockwise polygon, find if it's nested inside a hole
    let mut cw_containers: Vec<Option<usize>> = vec![None; metas.len()];

    for i in 0..metas.len() {
        if metas[i].is_clockwise {
            for j in 0..metas.len() {
                if i == j || metas[j].is_clockwise {
                    continue;
                }
                if !metas[i].bbox.is_inside(&metas[j].bbox) {
                    continue;
                }
                if !metas[i].ring.is_empty()
                    && polygon_contains_point(&metas[j].ring, &metas[i].ring[0])
                {
                    cw_containers[i] = Some(j);
                    break;
                }
            }
        }
    }

    // Step 4: Build final polygon structure
    let mut result: Vec<Vec<Vec<Position>>> = Vec::new();

    // Root-level clockwise polygons
    for i in 0..metas.len() {
        if metas[i].is_clockwise && cw_containers[i].is_none() {
            let mut polygon = vec![metas[i].ring.clone()];
            for j in 0..metas.len() {
                if !metas[j].is_clockwise && hole_containers[j] == Some(i) {
                    polygon.push(metas[j].ring.clone());
                }
            }
            result.push(polygon);
        }
    }

    // Islands (clockwise polygons inside holes)
    for i in 0..metas.len() {
        if metas[i].is_clockwise && cw_containers[i].is_some() {
            let mut polygon = vec![metas[i].ring.clone()];
            for j in 0..metas.len() {
                if !metas[j].is_clockwise && hole_containers[j] == Some(i) {
                    polygon.push(metas[j].ring.clone());
                }
            }
            result.push(polygon);
        }
    }

    // Orphan counter-clockwise polygons (holes with no container)
    for i in 0..metas.len() {
        if !metas[i].is_clockwise && hole_containers[i].is_none() {
            result.push(vec![metas[i].ring.clone()]);
        }
    }

    result
}

/// Compute winding direction for a ring (clockwise or counter-clockwise)
fn compute_winding(ring: &[Position]) -> &'static str {
    if ring.len() < 3 {
        return "unknown";
    }

    // Shoelace formula (signed area)
    let mut signed_area = 0.0;
    let n = ring.len();
    for i in 0..n {
        let j = (i + 1) % n;
        signed_area += ring[i][0] * ring[j][1];
        signed_area -= ring[j][0] * ring[i][1];
    }

    if signed_area > 0.0 {
        "counter-clockwise"
    } else {
        "clockwise"
    }
}

/// Check if a polygon contains a point using ray casting algorithm
fn polygon_contains_point(ring: &[Position], point: &Position) -> bool {
    if ring.len() < 3 {
        return false;
    }

    let mut inside = false;
    let px = point[0];
    let py = point[1];

    let n = ring.len();
    let mut j = n - 1;

    for i in 0..n {
        let xi = ring[i][0];
        let yi = ring[i][1];
        let xj = ring[j][0];
        let yj = ring[j][1];

        if ((yi > py) != (yj > py)) && (px < (xj - xi) * (py - yi) / (yj - yi) + xi) {
            inside = !inside;
        }

        j = i;
    }

    inside
}

/// Compute signed area of a ring using the shoelace formula
/// Positive = counter-clockwise, Negative = clockwise
fn signed_area(ring: &[Position]) -> f64 {
    if ring.len() < 3 {
        return 0.0;
    }

    let mut area = 0.0;
    let n = ring.len();
    for i in 0..n {
        let j = (i + 1) % n;
        area += ring[i][0] * ring[j][1];
        area -= ring[j][0] * ring[i][1];
    }
    area / 2.0
}

/// Check if a ring is counter-clockwise (positive signed area)
fn is_ccw(ring: &[Position]) -> bool {
    signed_area(ring) > 0.0
}

/// Reverse the winding order of a ring
fn reverse_ring(ring: &mut [Position]) {
    ring.reverse();
}

/// Orient polygons to follow RFC 7946 right-hand rule:
/// - Exterior rings: counter-clockwise (CCW)
/// - Interior rings (holes): clockwise (CW)
///
/// This ensures GeoJSON compliance and correct rendering in mapping libraries.
pub fn orient_polygons(polygons: Vec<Vec<Vec<Position>>>) -> Vec<Vec<Vec<Position>>> {
    polygons
        .into_iter()
        .map(|mut polygon| {
            for (i, ring) in polygon.iter_mut().enumerate() {
                if ring.len() < 3 {
                    continue;
                }

                let ccw = is_ccw(ring);

                if i == 0 {
                    // Exterior ring should be CCW
                    if !ccw {
                        reverse_ring(ring);
                    }
                } else {
                    // Interior rings (holes) should be CW
                    if ccw {
                        reverse_ring(ring);
                    }
                }
            }
            polygon
        })
        .collect()
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
    process_band_from_cells_with_precision(data, lower, upper, DEFAULT_PRECISION)
}

/// Process a single isoband level with configurable coordinate precision
///
/// # Arguments
///
/// * `data` - 2D array of GridCell structs
/// * `lower` - Lower threshold for this isoband
/// * `upper` - Upper threshold for this isoband
/// * `precision` - Number of decimal places for coordinates (1-6 recommended)
pub fn process_band_from_cells_with_precision(
    data: &[Vec<crate::GridCell>],
    lower: f64,
    upper: f64,
    precision: u32,
) -> Feature {
    let rows = data.len();
    let cols = data[0].len();

    // Create cells from grid (rows-1 × cols-1)
    let mut cells: Vec<Vec<Option<Shape>>> = Vec::with_capacity(rows - 1);
    for _ in 0..rows - 1 {
        cells.push(vec![None; cols - 1]);
    }

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

    let cell_rows = cells.len();
    let cell_cols = cells[0].len();

    // Walk grid in order, recursively processing polygons and their holes
    let mut polygons: Vec<Vec<Vec<Position>>> = Vec::new();
    let mut processed = std::collections::HashSet::new();

    for r in 0..cell_rows {
        for c in 0..cell_cols {
            if !processed.contains(&(r, c)) {
                if let Some(cell) = &cells[r][c] {
                    if !cell.is_cleared() {
                        let polygon = walk_polygon_recursive_with_precision(
                            &mut cells,
                            data,
                            lower,
                            upper,
                            r,
                            c,
                            &mut processed,
                            precision,
                        );
                        if !polygon.is_empty() && !polygon[0].is_empty() {
                            polygons.push(polygon);
                        }
                    }
                }
            }
        }
    }

    // Build containment hierarchy and orient polygons
    let hierarchical_polygons = build_polygon_hierarchy(polygons);
    let oriented_polygons = orient_polygons(hierarchical_polygons);

    // Build MultiPolygon geometry
    let multi_polygon = GeoValue::MultiPolygon(oriented_polygons);

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
pub fn do_concurrent(data: &[Vec<Feature>], isobands: &[f64]) -> geojson::FeatureCollection {
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
    do_concurrent_from_cells_with_precision(data, isobands, DEFAULT_PRECISION)
}

/// Process multiple isoband levels with configurable coordinate precision
///
/// # Arguments
///
/// * `data` - 2D array of GridCell structs
/// * `isobands` - Sorted list of threshold values
/// * `precision` - Number of decimal places for coordinates (1-6 recommended)
///   - 1 = ~11 km, 2 = ~1.1 km, 3 = ~110 m, 4 = ~11 m, 5 = ~1.1 m, 6 = ~0.11 m
pub fn do_concurrent_from_cells_with_precision(
    data: &[Vec<crate::GridCell>],
    isobands: &[f64],
    precision: u32,
) -> geojson::FeatureCollection {
    let features: Vec<Feature> = (0..isobands.len() - 1)
        .map(|i| {
            process_band_from_cells_with_precision(data, isobands[i], isobands[i + 1], precision)
        })
        .filter(|feature| has_coordinates(feature))
        .collect();

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
    process_line_from_cells_with_precision(data, isovalue, DEFAULT_PRECISION)
}

/// Process a single isoline with configurable coordinate precision
pub fn process_line_from_cells_with_precision(
    data: &[Vec<crate::GridCell>],
    isovalue: f64,
    precision: u32,
) -> Feature {
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

    let polylines = assembler.assemble_with_precision(precision);

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
pub fn do_concurrent_lines(data: &[Vec<Feature>], isovalues: &[f64]) -> geojson::FeatureCollection {
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
    do_concurrent_lines_from_cells_with_precision(data, isovalues, DEFAULT_PRECISION)
}

/// Generate multiple isolines with configurable coordinate precision
///
/// # Arguments
///
/// * `data` - 2D array of GridCell structs
/// * `isovalues` - Array of contour level values
/// * `precision` - Number of decimal places for coordinates (1-6 recommended)
///   - 1 = ~11 km, 2 = ~1.1 km, 3 = ~110 m, 4 = ~11 m, 5 = ~1.1 m, 6 = ~0.11 m
pub fn do_concurrent_lines_from_cells_with_precision(
    data: &[Vec<crate::GridCell>],
    isovalues: &[f64],
    precision: u32,
) -> geojson::FeatureCollection {
    use rayon::prelude::*;

    let features: Vec<Feature> = isovalues
        .par_iter()
        .map(|&isovalue| process_line_from_cells_with_precision(data, isovalue, precision))
        .filter(|feature| has_line_coordinates(feature))
        .collect();

    geojson::FeatureCollection {
        bbox: None,
        foreign_members: None,
        features,
    }
}

// =============================================================================
// Flat Array Processing (for Lambda/high-performance use cases)
// =============================================================================

/// A polygon with exterior ring and optional holes, using f32 grid-space coordinates.
///
/// This is the output format for `process_band_flat()`, designed for efficient
/// serialization and transmission (e.g., in Lambda responses).
///
/// Coordinates are in grid-space:
/// - x = column index (0 to width-1, with sub-pixel interpolation)
/// - y = row index (0 to height-1, with sub-pixel interpolation)
///
/// Winding order follows RFC 7946:
/// - Exterior ring: counter-clockwise (CCW)
/// - Holes: clockwise (CW)
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct ContourPolygon {
    /// Exterior ring (CCW winding in grid space)
    pub exterior: Vec<[f32; 2]>,

    /// Interior rings / holes (CW winding in grid space)
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub holes: Vec<Vec<[f32; 2]>>,
}

/// Timing metrics for contour generation phases.
///
/// All times are in milliseconds (f64).
#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct ContourMetrics {
    /// Time to build GridCell array from flat values
    pub grid_build_ms: f64,
    /// Time to classify cells into marching squares shapes
    pub classification_ms: f64,
    /// Time to walk edges and trace polygon boundaries
    pub edge_walking_ms: f64,
    /// Time to determine polygon hierarchy (holes vs exteriors)
    pub hierarchy_ms: f64,
    /// Time to orient polygons per RFC 7946
    pub orientation_ms: f64,
    /// Time to convert to final ContourPolygon format
    pub conversion_ms: f64,
    /// Total time for all phases
    pub total_ms: f64,
    /// Number of raw polygons before hierarchy resolution
    pub raw_polygon_count: usize,
    /// Number of final polygons after hierarchy resolution
    pub final_polygon_count: usize,
}

/// Process a single isoband from a flat f32 array
///
/// Returns polygons in grid-space coordinates (f32).
/// Coordinates represent positions in the grid where:
/// - x = column index (0 to width-1, with sub-pixel interpolation)
/// - y = row index (0 to height-1, with sub-pixel interpolation)
///
/// # Arguments
///
/// * `values` - Flat row-major array of grid values (f32, e.g., from Zarr)
/// * `width` - Number of columns in the grid
/// * `height` - Number of rows in the grid
/// * `lower` - Lower threshold (inclusive)
/// * `upper` - Upper threshold (exclusive)
/// * `precision` - Decimal places for coordinate rounding (default 5)
///
/// # Returns
///
/// Vector of polygons, each with exterior ring and optional holes.
/// Winding order follows RFC 7946 (exterior CCW, holes CW).
///
/// # Panics
///
/// Panics if `values.len() != width * height`
///
/// # Example
///
/// ```rust
/// use geo_marching_squares::process_band_flat;
///
/// // 3x3 grid with a high-value region in the bottom-right
/// let values: Vec<f32> = vec![
///     5.0, 5.0, 5.0,
///     5.0, 15.0, 15.0,
///     5.0, 15.0, 15.0,
/// ];
///
/// let polygons = process_band_flat(&values, 3, 3, 10.0, 20.0, 5);
/// // Returns contour polygon(s) for values in [10, 20)
/// ```
pub fn process_band_flat(
    values: &[f32],
    width: usize,
    height: usize,
    lower: f32,
    upper: f32,
    precision: u32,
) -> Vec<ContourPolygon> {
    assert_eq!(
        values.len(),
        width * height,
        "Grid size mismatch: expected {}x{}={}, got {}",
        width,
        height,
        width * height,
        values.len()
    );

    let rows = height;
    let cols = width;

    // Build a 2D GridCell array with grid indices as coordinates
    // Row 0 = top of grid, y increases downward
    let grid: Vec<Vec<crate::GridCell>> = (0..rows)
        .map(|r| {
            (0..cols)
                .map(|c| {
                    let idx = r * cols + c;
                    crate::GridCell {
                        x: c as f64,
                        y: r as f64,
                        value: values[idx] as f64,
                    }
                })
                .collect()
        })
        .collect();

    // Use existing process_band_from_cells_with_precision to generate GeoJSON
    let feature =
        process_band_from_cells_with_precision(&grid, lower as f64, upper as f64, precision);

    // Convert GeoJSON MultiPolygon to ContourPolygon format with f32 coords
    convert_feature_to_contour_polygons(feature, precision)
}

/// Convert a GeoJSON Feature with MultiPolygon geometry to ContourPolygon format
fn convert_feature_to_contour_polygons(feature: Feature, precision: u32) -> Vec<ContourPolygon> {
    let geometry = match feature.geometry {
        Some(g) => g,
        None => return vec![],
    };

    let multi_polygon = match geometry.value {
        GeoValue::MultiPolygon(mp) => mp,
        _ => return vec![],
    };

    let factor = 10_f32.powi(precision as i32);

    multi_polygon
        .into_iter()
        .map(|polygon_rings| {
            let exterior = polygon_rings
                .get(0)
                .map(|ring| {
                    ring.iter()
                        .map(|pos| {
                            [
                                (pos[0] as f32 * factor).round() / factor,
                                (pos[1] as f32 * factor).round() / factor,
                            ]
                        })
                        .collect()
                })
                .unwrap_or_default();

            let holes: Vec<Vec<[f32; 2]>> = polygon_rings
                .iter()
                .skip(1)
                .map(|ring| {
                    ring.iter()
                        .map(|pos| {
                            [
                                (pos[0] as f32 * factor).round() / factor,
                                (pos[1] as f32 * factor).round() / factor,
                            ]
                        })
                        .collect()
                })
                .collect();

            ContourPolygon { exterior, holes }
        })
        .collect()
}

/// Process a single isoband from a flat f32 array, returning detailed timing metrics.
///
/// This is the metrics-enabled version of `process_band_flat()`. Use this when you
/// need detailed performance data for each phase of contour generation.
///
/// # Arguments
///
/// * `values` - Flat row-major array of grid values (f32, e.g., from Zarr)
/// * `width` - Number of columns in the grid
/// * `height` - Number of rows in the grid
/// * `lower` - Lower threshold (inclusive)
/// * `upper` - Upper threshold (exclusive)
/// * `precision` - Decimal places for coordinate rounding (default 5)
///
/// # Returns
///
/// A tuple of (polygons, metrics) where:
/// - `polygons`: Vector of polygons with exterior ring and optional holes
/// - `metrics`: Timing data for each processing phase
///
/// # Example
///
/// ```rust,ignore
/// use geo_marching_squares::process_band_flat_with_metrics;
///
/// let values: Vec<f32> = fetch_zarr_data();
/// let (polygons, metrics) = process_band_flat_with_metrics(&values, 1799, 1059, 10.0, 20.0, 5);
///
/// println!("Generated {} polygons", polygons.len());
/// println!("Classification: {:.2}ms", metrics.classification_ms);
/// println!("Edge walking: {:.2}ms", metrics.edge_walking_ms);
/// println!("Hierarchy: {:.2}ms", metrics.hierarchy_ms);
/// ```
pub fn process_band_flat_with_metrics(
    values: &[f32],
    width: usize,
    height: usize,
    lower: f32,
    upper: f32,
    precision: u32,
) -> (Vec<ContourPolygon>, ContourMetrics) {
    use std::time::Instant;

    let total_start = Instant::now();

    assert_eq!(
        values.len(),
        width * height,
        "Grid size mismatch: expected {}x{}={}, got {}",
        width,
        height,
        width * height,
        values.len()
    );

    let rows = height;
    let cols = width;

    // Phase 1: Build GridCell array
    let grid_start = Instant::now();
    let grid: Vec<Vec<crate::GridCell>> = (0..rows)
        .map(|r| {
            (0..cols)
                .map(|c| {
                    let idx = r * cols + c;
                    crate::GridCell {
                        x: c as f64,
                        y: r as f64,
                        value: values[idx] as f64,
                    }
                })
                .collect()
        })
        .collect();
    let grid_build_ms = grid_start.elapsed().as_secs_f64() * 1000.0;

    // Phase 2: Shape classification (create cells)
    let classification_start = Instant::now();
    let mut cells: Vec<Vec<Option<Shape>>> = Vec::with_capacity(rows - 1);
    for _ in 0..rows - 1 {
        cells.push(vec![None; cols - 1]);
    }

    for r in 0..rows - 1 {
        for c in 0..cols - 1 {
            cells[r][c] = Shape::create_from_cells(
                &grid[r][c],
                &grid[r][c + 1],
                &grid[r + 1][c + 1],
                &grid[r + 1][c],
                lower as f64,
                upper as f64,
                c,
                r,
                r == 0,
                c + 1 == cols - 1,
                r + 1 == rows - 1,
                c == 0,
            );
        }
    }
    let classification_ms = classification_start.elapsed().as_secs_f64() * 1000.0;

    let cell_rows = cells.len();
    let cell_cols = cells[0].len();

    // Phase 3: Edge walking
    let edge_walking_start = Instant::now();
    let mut polygons: Vec<Vec<Vec<Position>>> = Vec::new();
    let mut processed = std::collections::HashSet::new();

    for r in 0..cell_rows {
        for c in 0..cell_cols {
            if !processed.contains(&(r, c)) {
                if let Some(cell) = &cells[r][c] {
                    if !cell.is_cleared() {
                        let polygon = walk_polygon_recursive_with_precision(
                            &mut cells,
                            &grid,
                            lower as f64,
                            upper as f64,
                            r,
                            c,
                            &mut processed,
                            precision,
                        );
                        if !polygon.is_empty() && !polygon[0].is_empty() {
                            polygons.push(polygon);
                        }
                    }
                }
            }
        }
    }
    let edge_walking_ms = edge_walking_start.elapsed().as_secs_f64() * 1000.0;
    let raw_polygon_count = polygons.len();

    // Phase 4: Build polygon hierarchy
    let hierarchy_start = Instant::now();
    let hierarchical_polygons = build_polygon_hierarchy(polygons);
    let hierarchy_ms = hierarchy_start.elapsed().as_secs_f64() * 1000.0;

    // Phase 5: Orient polygons
    let orientation_start = Instant::now();
    let oriented_polygons = orient_polygons(hierarchical_polygons);
    let orientation_ms = orientation_start.elapsed().as_secs_f64() * 1000.0;
    let final_polygon_count = oriented_polygons.len();

    // Phase 6: Convert to ContourPolygon format
    let conversion_start = Instant::now();
    let multi_polygon = GeoValue::MultiPolygon(oriented_polygons);
    let feature = Feature {
        bbox: None,
        geometry: Some(Geometry::new(multi_polygon)),
        id: None,
        properties: Some(JsonObject::new()),
        foreign_members: None,
    };
    let contour_polygons = convert_feature_to_contour_polygons(feature, precision);
    let conversion_ms = conversion_start.elapsed().as_secs_f64() * 1000.0;

    let total_ms = total_start.elapsed().as_secs_f64() * 1000.0;

    let metrics = ContourMetrics {
        grid_build_ms,
        classification_ms,
        edge_walking_ms,
        hierarchy_ms,
        orientation_ms,
        conversion_ms,
        total_ms,
        raw_polygon_count,
        final_polygon_count,
    };

    (contour_polygons, metrics)
}

/// Process multiple isoband levels from a flat f32 array in parallel
///
/// This is the parallel version of `process_band_flat()`, processing all bands
/// concurrently using Rayon's work-stealing thread pool.
///
/// # Arguments
///
/// * `values` - Flat row-major array of grid values (f32)
/// * `width` - Number of columns in the grid
/// * `height` - Number of rows in the grid
/// * `thresholds` - Array of threshold values (N thresholds = N-1 bands)
/// * `precision` - Decimal places for coordinate rounding
///
/// # Returns
///
/// Vector of (lower, upper, polygons) tuples, one per band.
/// Empty bands (no polygons) are filtered out.
///
/// # Example
///
/// ```rust,ignore
/// use geo_marching_squares::do_concurrent_flat;
///
/// let values: Vec<f32> = load_zarr_data();
/// let thresholds = vec![0.0, 10.0, 20.0, 30.0];
///
/// // Process 3 bands in parallel: 0-10, 10-20, 20-30
/// let results = do_concurrent_flat(&values, 1799, 1059, &thresholds, 5);
/// ```
pub fn do_concurrent_flat(
    values: &[f32],
    width: usize,
    height: usize,
    thresholds: &[f32],
    precision: u32,
) -> Vec<(f32, f32, Vec<ContourPolygon>)> {
    use rayon::prelude::*;

    if thresholds.len() < 2 {
        return vec![];
    }

    // Generate band intervals from thresholds
    let bands: Vec<(f32, f32)> = thresholds.windows(2).map(|w| (w[0], w[1])).collect();

    // Process all bands in parallel
    bands
        .into_par_iter()
        .map(|(lower, upper)| {
            let polygons = process_band_flat(values, width, height, lower, upper, precision);
            (lower, upper, polygons)
        })
        .filter(|(_, _, polygons)| !polygons.is_empty())
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_round_coord() {
        // Test with default precision (5 decimals)
        assert_eq!(round_coord_with_precision(1.234567, 5), 1.23457);
        assert_eq!(round_coord_with_precision(1.234564, 5), 1.23456);
        assert_eq!(round_coord_with_precision(-122.123456789, 5), -122.12346);
        assert_eq!(round_coord_with_precision(45.678901234, 5), 45.6789);
        // Test with lower precision
        assert_eq!(round_coord_with_precision(45.678901234, 3), 45.679);
        assert_eq!(round_coord_with_precision(45.678901234, 2), 45.68);
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

    #[test]
    fn test_signed_area_ccw() {
        // Counter-clockwise square (positive area)
        let ccw_ring = vec![
            vec![0.0, 0.0],
            vec![10.0, 0.0],
            vec![10.0, 10.0],
            vec![0.0, 10.0],
            vec![0.0, 0.0],
        ];
        assert!(signed_area(&ccw_ring) > 0.0);
        assert!(is_ccw(&ccw_ring));
    }

    #[test]
    fn test_signed_area_cw() {
        // Clockwise square (negative area)
        let cw_ring = vec![
            vec![0.0, 0.0],
            vec![0.0, 10.0],
            vec![10.0, 10.0],
            vec![10.0, 0.0],
            vec![0.0, 0.0],
        ];
        assert!(signed_area(&cw_ring) < 0.0);
        assert!(!is_ccw(&cw_ring));
    }

    #[test]
    fn test_orient_polygons_fixes_exterior() {
        // Clockwise exterior (wrong) - should be fixed to CCW
        let cw_exterior = vec![
            vec![0.0, 0.0],
            vec![0.0, 10.0],
            vec![10.0, 10.0],
            vec![10.0, 0.0],
            vec![0.0, 0.0],
        ];
        let polygons = vec![vec![cw_exterior]];
        let oriented = orient_polygons(polygons);

        // After orientation, exterior should be CCW
        assert!(is_ccw(&oriented[0][0]));
    }

    #[test]
    fn test_orient_polygons_fixes_hole() {
        // CCW exterior (correct) with CCW hole (wrong - should be CW)
        let ccw_exterior = vec![
            vec![0.0, 0.0],
            vec![10.0, 0.0],
            vec![10.0, 10.0],
            vec![0.0, 10.0],
            vec![0.0, 0.0],
        ];
        let ccw_hole = vec![
            vec![2.0, 2.0],
            vec![8.0, 2.0],
            vec![8.0, 8.0],
            vec![2.0, 8.0],
            vec![2.0, 2.0],
        ];
        let polygons = vec![vec![ccw_exterior, ccw_hole]];
        let oriented = orient_polygons(polygons);

        // After orientation: exterior CCW, hole CW
        assert!(is_ccw(&oriented[0][0])); // exterior is CCW
        assert!(!is_ccw(&oriented[0][1])); // hole is CW
    }

    #[test]
    fn test_orient_polygons_preserves_correct() {
        // Already correct: CCW exterior, CW hole
        let ccw_exterior = vec![
            vec![0.0, 0.0],
            vec![10.0, 0.0],
            vec![10.0, 10.0],
            vec![0.0, 10.0],
            vec![0.0, 0.0],
        ];
        let cw_hole = vec![
            vec![2.0, 2.0],
            vec![2.0, 8.0],
            vec![8.0, 8.0],
            vec![8.0, 2.0],
            vec![2.0, 2.0],
        ];
        let original_exterior = ccw_exterior.clone();
        let original_hole = cw_hole.clone();

        let polygons = vec![vec![ccw_exterior, cw_hole]];
        let oriented = orient_polygons(polygons);

        // Should be unchanged
        assert_eq!(oriented[0][0], original_exterior);
        assert_eq!(oriented[0][1], original_hole);
    }

    // =========================================================================
    // Tests for flat array processing (process_band_flat, do_concurrent_flat)
    // =========================================================================

    #[test]
    fn test_process_band_flat_simple() {
        // 3x3 grid with a high-value region in the bottom-right
        // Grid layout (values):
        //   5   5   5
        //   5  15  15
        //   5  15  15
        let values: Vec<f32> = vec![5.0, 5.0, 5.0, 5.0, 15.0, 15.0, 5.0, 15.0, 15.0];

        let polygons = process_band_flat(&values, 3, 3, 10.0, 20.0, 5);

        // Should produce at least one polygon for the 10-20 band
        assert!(!polygons.is_empty(), "Should have at least one polygon");

        // Check that exterior ring has points
        assert!(
            !polygons[0].exterior.is_empty(),
            "Exterior ring should have points"
        );
    }

    #[test]
    fn test_process_band_flat_empty_result() {
        // All values below threshold
        let values: Vec<f32> = vec![5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0];

        let polygons = process_band_flat(&values, 3, 3, 10.0, 20.0, 5);

        // Should produce no polygons when all values are below threshold
        assert!(
            polygons.is_empty(),
            "Should have no polygons for uniform grid below threshold"
        );
    }

    #[test]
    fn test_process_band_flat_all_above() {
        // All values above threshold
        let values: Vec<f32> = vec![25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0];

        let polygons = process_band_flat(&values, 3, 3, 10.0, 20.0, 5);

        // Should produce no polygons when all values are above threshold
        assert!(
            polygons.is_empty(),
            "Should have no polygons for uniform grid above threshold"
        );
    }

    #[test]
    fn test_process_band_flat_grid_coords() {
        // 3x3 grid with values that create a contour
        let values: Vec<f32> = vec![5.0, 5.0, 5.0, 5.0, 15.0, 15.0, 5.0, 15.0, 15.0];

        let polygons = process_band_flat(&values, 3, 3, 10.0, 20.0, 5);

        assert!(!polygons.is_empty());

        // Verify coordinates are in grid space (0-2 range for a 3x3 grid)
        for poly in &polygons {
            for point in &poly.exterior {
                assert!(
                    point[0] >= 0.0 && point[0] <= 2.0,
                    "X coordinate {} should be in grid range [0, 2]",
                    point[0]
                );
                assert!(
                    point[1] >= 0.0 && point[1] <= 2.0,
                    "Y coordinate {} should be in grid range [0, 2]",
                    point[1]
                );
            }
        }
    }

    #[test]
    fn test_process_band_flat_f32_precision() {
        // Verify f32 coordinates maintain precision
        let values: Vec<f32> = vec![5.0, 5.0, 5.0, 5.0, 15.0, 15.0, 5.0, 15.0, 15.0];

        let polygons = process_band_flat(&values, 3, 3, 10.0, 20.0, 3);

        assert!(!polygons.is_empty());

        // Check that coordinates are rounded to 3 decimal places
        for poly in &polygons {
            for point in &poly.exterior {
                // Multiply by 1000, round, check if it's close to an integer
                let x_scaled = point[0] * 1000.0;
                let y_scaled = point[1] * 1000.0;
                assert!(
                    (x_scaled - x_scaled.round()).abs() < 0.001,
                    "X coordinate should be rounded to 3 decimal places"
                );
                assert!(
                    (y_scaled - y_scaled.round()).abs() < 0.001,
                    "Y coordinate should be rounded to 3 decimal places"
                );
            }
        }
    }

    #[test]
    fn test_do_concurrent_flat_multiple_bands() {
        // 5x5 grid with gradient values
        let mut values: Vec<f32> = Vec::with_capacity(25);
        for r in 0..5 {
            for c in 0..5 {
                values.push((r * 5 + c) as f32 * 2.0); // Values from 0 to 48
            }
        }

        let thresholds = vec![0.0, 10.0, 20.0, 30.0, 40.0];
        let results = do_concurrent_flat(&values, 5, 5, &thresholds, 5);

        // Should have at least some bands with polygons
        assert!(
            !results.is_empty(),
            "Should have at least one band with polygons"
        );

        // Each result should have valid lower/upper thresholds
        for (lower, upper, polygons) in &results {
            assert!(lower < upper, "Lower should be less than upper");
            assert!(
                !polygons.is_empty(),
                "Each returned band should have polygons"
            );
        }
    }

    #[test]
    fn test_do_concurrent_flat_filters_empty() {
        // Uniform grid - no bands will have polygons
        let values: Vec<f32> = vec![5.0; 25]; // 5x5 grid of 5.0

        let thresholds = vec![10.0, 20.0, 30.0];
        let results = do_concurrent_flat(&values, 5, 5, &thresholds, 5);

        // All bands should be filtered out (values are all below thresholds)
        assert!(results.is_empty(), "Should filter out empty bands");
    }

    #[test]
    fn test_contour_polygon_serde() {
        let polygon = ContourPolygon {
            exterior: vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.0, 0.0]],
            holes: vec![vec![
                [0.2, 0.2],
                [0.2, 0.8],
                [0.8, 0.8],
                [0.8, 0.2],
                [0.2, 0.2],
            ]],
        };

        // Serialize to JSON
        let json = serde_json::to_string(&polygon).unwrap();

        // Deserialize back
        let deserialized: ContourPolygon = serde_json::from_str(&json).unwrap();

        assert_eq!(polygon.exterior.len(), deserialized.exterior.len());
        assert_eq!(polygon.holes.len(), deserialized.holes.len());
    }

    #[test]
    fn test_contour_polygon_empty_holes_skipped() {
        let polygon = ContourPolygon {
            exterior: vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.0, 0.0]],
            holes: vec![], // Empty holes
        };

        // Serialize to JSON
        let json = serde_json::to_string(&polygon).unwrap();

        // Empty holes should be skipped in serialization
        assert!(
            !json.contains("holes"),
            "Empty holes should be skipped: {}",
            json
        );
    }

    #[test]
    #[should_panic(expected = "Grid size mismatch")]
    fn test_process_band_flat_size_mismatch() {
        let values: Vec<f32> = vec![5.0; 10]; // 10 values
        process_band_flat(&values, 3, 3, 10.0, 20.0, 5); // Expects 9 values
    }
}
