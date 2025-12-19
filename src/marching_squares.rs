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

    eprintln!("[geo-marching-squares] Polygon at ({}, {}) winding: {}, edges: {}",
        start_r, start_c, winding, exterior_ring.len());

    // Return just the exterior ring - containment will be determined in post-processing
    // No more flood fill or recursive hole processing
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
fn build_polygon_hierarchy(raw_polygons: Vec<Vec<Vec<Position>>>) -> Vec<Vec<Vec<Position>>> {
    let hierarchy_start = std::time::Instant::now();

    if raw_polygons.is_empty() {
        return vec![];
    }

    // Step 1: Build metadata for all polygons (just exterior rings)
    let metadata_start = std::time::Instant::now();
    let metas: Vec<PolygonMeta> = raw_polygons
        .into_iter()
        .enumerate()
        .filter_map(|(_idx, mut poly)| {
            if poly.is_empty() || poly[0].is_empty() {
                return None;
            }
            let ring = poly.remove(0); // Take just the exterior ring
            let bbox = BBox::from_ring(&ring);
            let is_clockwise = compute_winding(&ring) == "clockwise";
            Some(PolygonMeta {
                ring,
                bbox,
                is_clockwise,
            })
        })
        .collect();
    let metadata_elapsed = metadata_start.elapsed();

    let cw_count = metas.iter().filter(|m| m.is_clockwise).count();
    let ccw_count = metas.iter().filter(|m| !m.is_clockwise).count();
    eprintln!("[geo-marching-squares]   ⏱️  Metadata (bbox+winding): {:?} for {} polygons ({} CW, {} CCW)",
        metadata_elapsed, metas.len(), cw_count, ccw_count);

    // Step 2: Find container for each counter-clockwise polygon (holes)
    let hole_containment_start = std::time::Instant::now();
    let mut hole_containers: Vec<Option<usize>> = vec![None; metas.len()];
    let mut bbox_checks = 0usize;
    let mut pip_checks = 0usize;

    for i in 0..metas.len() {
        if !metas[i].is_clockwise {
            // This is a hole - find its clockwise container (early exit on first match)
            let mut container: Option<usize> = None;

            for j in 0..metas.len() {
                if i == j || !metas[j].is_clockwise {
                    continue; // Skip self and other holes
                }

                bbox_checks += 1;
                // Quick bbox check first
                if !metas[i].bbox.is_inside(&metas[j].bbox) {
                    continue;
                }

                pip_checks += 1;
                // Point-in-polygon check - test first point of hole against container
                if metas[i].ring.len() > 0 && polygon_contains_point(&metas[j].ring, &metas[i].ring[0]) {
                    // Found a container! Early exit instead of finding smallest
                    container = Some(j);
                    break;
                }
            }

            hole_containers[i] = container;
            if container.is_none() {
                eprintln!("[geo-marching-squares]   WARNING: Counter-clockwise polygon {} has no container", i);
            }
        }
    }
    let hole_containment_elapsed = hole_containment_start.elapsed();
    eprintln!("[geo-marching-squares]   ⏱️  Hole Containment: {:?} ({} bbox checks, {} point-in-poly checks)",
        hole_containment_elapsed, bbox_checks, pip_checks);

    // Step 3: For each clockwise polygon, find if it's nested inside a hole
    let island_containment_start = std::time::Instant::now();
    let mut cw_containers: Vec<Option<usize>> = vec![None; metas.len()];
    let mut island_bbox_checks = 0usize;
    let mut island_pip_checks = 0usize;

    for i in 0..metas.len() {
        if metas[i].is_clockwise {
            // Check if this clockwise polygon is inside a counter-clockwise polygon (island in hole)
            for j in 0..metas.len() {
                if i == j || metas[j].is_clockwise {
                    continue;
                }

                island_bbox_checks += 1;
                if !metas[i].bbox.is_inside(&metas[j].bbox) {
                    continue;
                }

                island_pip_checks += 1;
                if metas[i].ring.len() > 0 && polygon_contains_point(&metas[j].ring, &metas[i].ring[0]) {
                    cw_containers[i] = Some(j);
                    break; // Found container
                }
            }
        }
    }
    let island_containment_elapsed = island_containment_start.elapsed();
    eprintln!("[geo-marching-squares]   ⏱️  Island Containment: {:?} ({} bbox checks, {} point-in-poly checks)",
        island_containment_elapsed, island_bbox_checks, island_pip_checks);

    // Step 4: Build final polygon structure
    let assembly_start = std::time::Instant::now();
    // Root-level clockwise polygons (not inside any hole) become top-level multipolygon entries
    let mut result: Vec<Vec<Vec<Position>>> = Vec::new();

    for i in 0..metas.len() {
        if metas[i].is_clockwise && cw_containers[i].is_none() {
            // This is a root-level clockwise polygon
            let mut polygon = vec![metas[i].ring.clone()];

            // Attach its holes (counter-clockwise polygons that have this as container)
            for j in 0..metas.len() {
                if !metas[j].is_clockwise && hole_containers[j] == Some(i) {
                    polygon.push(metas[j].ring.clone());
                }
            }

            result.push(polygon);
        }
    }

    // Also add islands (clockwise polygons inside holes) as separate polygons
    for i in 0..metas.len() {
        if metas[i].is_clockwise && cw_containers[i].is_some() {
            // This is an island - add as separate polygon with its own holes
            let mut polygon = vec![metas[i].ring.clone()];

            // Find holes inside this island
            for j in 0..metas.len() {
                if !metas[j].is_clockwise && hole_containers[j] == Some(i) {
                    polygon.push(metas[j].ring.clone());
                }
            }

            result.push(polygon);
        }
    }

    // Handle orphan counter-clockwise polygons (holes with no container)
    // These can occur at grid boundaries or in test cases with small grids
    // Treat them as standalone polygons (they'll appear as "inverted" shapes)
    for i in 0..metas.len() {
        if !metas[i].is_clockwise && hole_containers[i].is_none() {
            // Lone CCW polygon - add as standalone
            result.push(vec![metas[i].ring.clone()]);
        }
    }
    let assembly_elapsed = assembly_start.elapsed();
    let hierarchy_total = hierarchy_start.elapsed();

    eprintln!("[geo-marching-squares]   ⏱️  Polygon Assembly: {:?} ({} output polygons)",
        assembly_elapsed, result.len());
    eprintln!("[geo-marching-squares]   ⏱️  TOTAL HIERARCHY: {:?}", hierarchy_total);

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
                        let polygon = walk_polygon_recursive_with_precision(&mut cells, data, lower, upper, r, c, &mut processed, precision);
                        if !polygon.is_empty() && !polygon[0].is_empty() {
                            polygons.push(polygon);
                        }
                    }
                }
            }
        }
    }

    let recursive_elapsed = recursive_start.elapsed();
    eprintln!("[geo-marching-squares] ⏱️  Edge Walking: {:?} ({} raw polygons)",
        recursive_elapsed, polygons.len());

    // Build containment hierarchy using winding direction + bbox optimization
    let nesting_start = std::time::Instant::now();
    let hierarchical_polygons = build_polygon_hierarchy(polygons);
    let nesting_elapsed = nesting_start.elapsed();
    eprintln!("[geo-marching-squares] ⏱️  Polygon Nesting: {:?} ({} final polygons)",
        nesting_elapsed, hierarchical_polygons.len());

    eprintln!("[geo-marching-squares] process_band_from_cells: Polygon processing complete, building feature with {} polygons",
        hierarchical_polygons.len());

    // Build MultiPolygon geometry
    let geometry_start = std::time::Instant::now();
    let multi_polygon = GeoValue::MultiPolygon(hierarchical_polygons);

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
    // Note: rayon removed for sequential processing test

    eprintln!("[geo-marching-squares] do_concurrent_from_cells: Starting with {} thresholds, {} bands to generate (precision={})",
        isobands.len(), isobands.len().saturating_sub(1), precision);
    eprintln!("[geo-marching-squares] Grid size: {}x{} = {} cells",
        data.len(),
        if data.is_empty() { 0 } else { data[0].len() },
        data.len() * if data.is_empty() { 0 } else { data[0].len() });

    // Process each isoband pair sequentially (testing performance)
    let features: Vec<Feature> = (0..isobands.len() - 1)
        .into_iter()  // Changed from into_par_iter() for sequential processing
        .map(|i| {
            eprintln!("[geo-marching-squares] Band {}/{}: Processing [{}, {})",
                i + 1, isobands.len() - 1, isobands[i], isobands[i + 1]);
            // Each iteration processes one isoband level
            let result = process_band_from_cells_with_precision(data, isobands[i], isobands[i + 1], precision);
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
}
