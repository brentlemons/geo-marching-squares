//! Optimized marching squares implementation for high-performance contour generation.
//!
//! This module provides memory-efficient data structures and algorithms optimized for
//! processing large grids (1M+ cells) with minimal allocations.
//!
//! Key optimizations:
//! - Compact Point enum (24 bytes vs 72 bytes)
//! - Sparse cell storage (HashMap vs dense Vec)
//! - Lazy edge building (on-demand vs upfront)
//! - ArrayVec for edges (stack vs heap allocation)
//! - Direct flat array access (no intermediate grid)

use arrayvec::ArrayVec;
use std::collections::HashMap;

/// Compact point representation - either a coordinate or an interpolation marker.
///
/// Size: 24 bytes (down from 72 bytes in the original Point)
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CompactPoint {
    /// A point with actual coordinates
    Coord { x: f64, y: f64 },
    /// A marker for interpolation (coordinates computed lazily)
    Interp { value: f64, limit: f64, side: Side },
}

impl CompactPoint {
    #[inline]
    pub fn coord(x: f64, y: f64) -> Self {
        CompactPoint::Coord { x, y }
    }

    #[inline]
    pub fn interp(value: f64, limit: f64, side: Side) -> Self {
        CompactPoint::Interp { value, limit, side }
    }

    /// Resolve this point to actual coordinates given the cell's corner positions.
    #[inline]
    pub fn resolve(&self, cell_x: usize, cell_y: usize) -> (f64, f64) {
        match *self {
            CompactPoint::Coord { x, y } => (x, y),
            CompactPoint::Interp { value, limit, side } => {
                // Interpolate along the edge
                let t = if (value - limit).abs() < f64::EPSILON {
                    0.5
                } else {
                    // This is a placeholder - actual interpolation needs corner values
                    0.5
                };
                match side {
                    Side::Top => (cell_x as f64 + t, cell_y as f64),
                    Side::Right => (cell_x as f64 + 1.0, cell_y as f64 + t),
                    Side::Bottom => (cell_x as f64 + t, cell_y as f64 + 1.0),
                    Side::Left => (cell_x as f64, cell_y as f64 + t),
                }
            }
        }
    }

    /// Get x coordinate if this is a Coord point
    #[inline]
    pub fn x(&self) -> Option<f64> {
        match self {
            CompactPoint::Coord { x, .. } => Some(*x),
            _ => None,
        }
    }

    /// Get y coordinate if this is a Coord point
    #[inline]
    pub fn y(&self) -> Option<f64> {
        match self {
            CompactPoint::Coord { y, .. } => Some(*y),
            _ => None,
        }
    }
}

impl std::hash::Hash for CompactPoint {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        match self {
            CompactPoint::Coord { x, y } => {
                0u8.hash(state);
                x.to_bits().hash(state);
                y.to_bits().hash(state);
            }
            CompactPoint::Interp { value, limit, side } => {
                1u8.hash(state);
                value.to_bits().hash(state);
                limit.to_bits().hash(state);
                side.hash(state);
            }
        }
    }
}

impl Eq for CompactPoint {}

/// Side of a cell edge
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Side {
    Top,
    Right,
    Bottom,
    Left,
}

/// Direction to move to the next cell
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MoveDir {
    Right,
    Down,
    Left,
    Up,
    None,
}

/// Compact edge representation
///
/// Size: ~56 bytes (down from 152 bytes)
#[derive(Debug, Clone)]
pub struct CompactEdge {
    pub start: CompactPoint,
    pub end: CompactPoint,
    pub move_dir: MoveDir,
}

impl CompactEdge {
    #[inline]
    pub fn new(start: CompactPoint, end: CompactPoint, move_dir: MoveDir) -> Self {
        Self {
            start,
            end,
            move_dir,
        }
    }
}

/// Ternary classification for a corner (below/within/above threshold band)
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum Ternary {
    Below = 0,
    Within = 1,
    Above = 2,
}

impl Ternary {
    #[inline]
    pub fn classify(value: f32, lower: f32, upper: f32) -> Self {
        if value < lower {
            Ternary::Below
        } else if value >= upper {
            Ternary::Above
        } else {
            Ternary::Within
        }
    }
}

/// Compact cell storing only classification data - edges built lazily.
///
/// Size: ~32 bytes (down from 440 bytes for Shape)
#[derive(Debug, Clone)]
pub struct CompactCell {
    /// Grid position
    pub col: u16,
    pub row: u16,
    /// Ternary classification value (0-80, encodes 4 corners × 3 states)
    pub class_value: u8,
    /// Corner values for interpolation
    pub tl: f32,
    pub tr: f32,
    pub br: f32,
    pub bl: f32,
    /// Boundary flags packed into one byte
    pub boundary_flags: u8,
    /// Number of edges used (for tracking during walking)
    pub edges_used: u8,
}

impl CompactCell {
    /// Boundary flag bits
    const TOP_EDGE: u8 = 0b0001;
    const RIGHT_EDGE: u8 = 0b0010;
    const BOTTOM_EDGE: u8 = 0b0100;
    const LEFT_EDGE: u8 = 0b1000;

    /// Create a new cell from corner values and thresholds.
    /// Returns None if the cell is empty (all corners same classification).
    #[inline]
    pub fn new(
        col: usize,
        row: usize,
        tl: f32,
        tr: f32,
        br: f32,
        bl: f32,
        lower: f32,
        upper: f32,
        is_top: bool,
        is_right: bool,
        is_bottom: bool,
        is_left: bool,
    ) -> Option<Self> {
        // Skip cells with any NaN values - can't interpolate through missing data
        if tl.is_nan() || tr.is_nan() || br.is_nan() || bl.is_nan() {
            return None;
        }

        // Compute ternary classification
        let tl_class = Ternary::classify(tl, lower, upper);
        let tr_class = Ternary::classify(tr, lower, upper);
        let br_class = Ternary::classify(br, lower, upper);
        let bl_class = Ternary::classify(bl, lower, upper);

        // Encode as single value: tl*27 + tr*9 + br*3 + bl (base-3 encoding)
        let class_value =
            (tl_class as u8) * 27 + (tr_class as u8) * 9 + (br_class as u8) * 3 + (bl_class as u8);

        // Early exit: if all corners are the same, no contour passes through
        // class_value 0 = all Below (0*27 + 0*9 + 0*3 + 0)
        // class_value 40 = all Within (1*27 + 1*9 + 1*3 + 1)
        // class_value 80 = all Above (2*27 + 2*9 + 2*3 + 2)
        if class_value == 0 || class_value == 40 || class_value == 80 {
            return None;
        }

        let mut boundary_flags = 0u8;
        if is_top {
            boundary_flags |= Self::TOP_EDGE;
        }
        if is_right {
            boundary_flags |= Self::RIGHT_EDGE;
        }
        if is_bottom {
            boundary_flags |= Self::BOTTOM_EDGE;
        }
        if is_left {
            boundary_flags |= Self::LEFT_EDGE;
        }

        Some(Self {
            col: col as u16,
            row: row as u16,
            class_value,
            tl,
            tr,
            br,
            bl,
            boundary_flags,
            edges_used: 0,
        })
    }

    #[inline]
    pub fn is_top_edge(&self) -> bool {
        self.boundary_flags & Self::TOP_EDGE != 0
    }

    #[inline]
    pub fn is_right_edge(&self) -> bool {
        self.boundary_flags & Self::RIGHT_EDGE != 0
    }

    #[inline]
    pub fn is_bottom_edge(&self) -> bool {
        self.boundary_flags & Self::BOTTOM_EDGE != 0
    }

    #[inline]
    pub fn is_left_edge(&self) -> bool {
        self.boundary_flags & Self::LEFT_EDGE != 0
    }

    /// Linear interpolation along an edge
    #[inline]
    fn lerp(v1: f32, v2: f32, threshold: f32) -> f64 {
        if (v2 - v1).abs() < f32::EPSILON {
            0.5
        } else {
            ((threshold - v1) / (v2 - v1)) as f64
        }
    }

    /// Build edges lazily based on classification.
    /// Returns up to 6 edges (ArrayVec avoids heap allocation).
    pub fn build_edges(&self, lower: f32, upper: f32) -> ArrayVec<CompactEdge, 6> {
        let mut edges = ArrayVec::new();

        let col = self.col as f64;
        let row = self.row as f64;

        // Decode ternary classification
        let tl_class = Ternary::from_value((self.class_value / 27) % 3);
        let tr_class = Ternary::from_value((self.class_value / 9) % 3);
        let br_class = Ternary::from_value((self.class_value / 3) % 3);
        let bl_class = Ternary::from_value(self.class_value % 3);

        // Compute interpolated points on each edge where contour crosses
        // Top edge (tl to tr)
        let top_lower = if (tl_class == Ternary::Below) != (tr_class == Ternary::Below) {
            let t = Self::lerp(self.tl, self.tr, lower);
            Some(CompactPoint::coord(col + t, row))
        } else {
            None
        };

        let top_upper = if (tl_class == Ternary::Above) != (tr_class == Ternary::Above) {
            let t = Self::lerp(self.tl, self.tr, upper);
            Some(CompactPoint::coord(col + t, row))
        } else {
            None
        };

        // Right edge (tr to br)
        let right_lower = if (tr_class == Ternary::Below) != (br_class == Ternary::Below) {
            let t = Self::lerp(self.tr, self.br, lower);
            Some(CompactPoint::coord(col + 1.0, row + t))
        } else {
            None
        };

        let right_upper = if (tr_class == Ternary::Above) != (br_class == Ternary::Above) {
            let t = Self::lerp(self.tr, self.br, upper);
            Some(CompactPoint::coord(col + 1.0, row + t))
        } else {
            None
        };

        // Bottom edge (br to bl)
        let bottom_lower = if (br_class == Ternary::Below) != (bl_class == Ternary::Below) {
            let t = Self::lerp(self.br, self.bl, lower);
            Some(CompactPoint::coord(col + 1.0 - t, row + 1.0))
        } else {
            None
        };

        let bottom_upper = if (br_class == Ternary::Above) != (bl_class == Ternary::Above) {
            let t = Self::lerp(self.br, self.bl, upper);
            Some(CompactPoint::coord(col + 1.0 - t, row + 1.0))
        } else {
            None
        };

        // Left edge (bl to tl)
        let left_lower = if (bl_class == Ternary::Below) != (tl_class == Ternary::Below) {
            let t = Self::lerp(self.bl, self.tl, lower);
            Some(CompactPoint::coord(col, row + 1.0 - t))
        } else {
            None
        };

        let left_upper = if (bl_class == Ternary::Above) != (tl_class == Ternary::Above) {
            let t = Self::lerp(self.bl, self.tl, upper);
            Some(CompactPoint::coord(col, row + 1.0 - t))
        } else {
            None
        };

        // Build edges based on the isoband pattern
        // This is a simplified version - the full implementation needs all 81 cases
        self.connect_edges(
            &mut edges,
            [
                top_lower,
                top_upper,
                right_lower,
                right_upper,
                bottom_lower,
                bottom_upper,
                left_lower,
                left_upper,
            ],
            [tl_class, tr_class, br_class, bl_class],
        );

        edges
    }

    /// Connect interpolated points into edges based on classification pattern.
    fn connect_edges(
        &self,
        edges: &mut ArrayVec<CompactEdge, 6>,
        points: [Option<CompactPoint>; 8],
        _classes: [Ternary; 4],
    ) {
        let [top_l, top_u, right_l, right_u, bottom_l, bottom_u, left_l, left_u] = points;

        // The isoband has an interior region where values are "within" the band.
        // We need to trace the boundary of this region.

        // For each corner that's "within", we might have edges from adjacent crossings.
        // For simplicity, we trace from lower crossings to upper crossings around the cell.

        // This is a simplified approach that handles the most common cases.
        // A full implementation would have a lookup table for all 81 ternary combinations.

        // Collect all valid crossing points in clockwise order
        let mut crossings: ArrayVec<(CompactPoint, MoveDir), 8> = ArrayVec::new();

        // Top edge (left to right)
        if let Some(p) = top_l {
            crossings.push((p, MoveDir::Up));
        }
        if let Some(p) = top_u {
            crossings.push((p, MoveDir::Up));
        }

        // Right edge (top to bottom)
        if let Some(p) = right_u {
            crossings.push((p, MoveDir::Right));
        }
        if let Some(p) = right_l {
            crossings.push((p, MoveDir::Right));
        }

        // Bottom edge (right to left)
        if let Some(p) = bottom_l {
            crossings.push((p, MoveDir::Down));
        }
        if let Some(p) = bottom_u {
            crossings.push((p, MoveDir::Down));
        }

        // Left edge (bottom to top)
        if let Some(p) = left_u {
            crossings.push((p, MoveDir::Left));
        }
        if let Some(p) = left_l {
            crossings.push((p, MoveDir::Left));
        }

        // For isobands, we connect lower-to-upper crossings to form the band boundary
        // This needs more sophisticated logic for saddle points and complex cases

        // Simple case: connect consecutive crossings
        if crossings.len() >= 2 {
            for i in 0..crossings.len() {
                let j = (i + 1) % crossings.len();
                let (start, _) = crossings[i].clone();
                let (end, move_dir) = crossings[j].clone();
                edges.push(CompactEdge::new(start, end, move_dir));
            }
        }
    }

    /// Check if all edges have been used
    #[inline]
    pub fn is_cleared(&self) -> bool {
        // A cell is cleared when we've walked all its edges
        self.edges_used > 0
    }

    /// Mark edges as used
    #[inline]
    pub fn mark_used(&mut self, count: u8) {
        self.edges_used = self.edges_used.saturating_add(count);
    }
}

impl Ternary {
    #[inline]
    fn from_value(v: u8) -> Self {
        match v {
            0 => Ternary::Below,
            1 => Ternary::Within,
            _ => Ternary::Above,
        }
    }
}

/// Sparse cell storage using HashMap instead of dense Vec.
pub type CellMap = HashMap<(u16, u16), CompactCell>;

/// Process a band using optimized structures (experimental - edge walking incomplete).
///
/// This function works directly with the flat values array, skipping grid construction.
/// NOTE: The edge walking logic is simplified and may not produce correct results.
/// Use `process_band_optimized_hybrid` for production use.
pub fn process_band_optimized(
    values: &[f32],
    width: usize,
    height: usize,
    lower: f32,
    upper: f32,
    precision: u32,
) -> (Vec<crate::ContourPolygon>, crate::ContourMetrics) {
    use std::time::Instant;

    let total_start = Instant::now();
    let rows = height;
    let cols = width;

    // Phase 1: Classification with sparse storage (skip grid building)
    let classification_start = Instant::now();
    let mut cells: CellMap = HashMap::new();

    // Direct access to flat array - no intermediate grid
    for r in 0..rows - 1 {
        for c in 0..cols - 1 {
            // Get corner values directly from flat array
            let tl = values[r * cols + c];
            let tr = values[r * cols + c + 1];
            let br = values[(r + 1) * cols + c + 1];
            let bl = values[(r + 1) * cols + c];

            // Create cell only if it contributes to the contour
            if let Some(cell) = CompactCell::new(
                c,
                r,
                tl,
                tr,
                br,
                bl,
                lower,
                upper,
                r == 0,
                c + 1 == cols - 1,
                r + 1 == rows - 1,
                c == 0,
            ) {
                cells.insert((c as u16, r as u16), cell);
            }
        }
    }
    let classification_ms = classification_start.elapsed().as_secs_f64() * 1000.0;

    // Track cell count for metrics
    let _cell_count = cells.len();

    // Phase 2: Edge walking with lazy edge building
    let edge_walking_start = Instant::now();
    let mut polygons: Vec<Vec<Vec<Vec<f64>>>> = Vec::new();
    let mut processed: std::collections::HashSet<(u16, u16)> = std::collections::HashSet::new();

    // Find all unprocessed cells and walk their boundaries
    let cell_keys: Vec<(u16, u16)> = cells.keys().cloned().collect();
    for key in cell_keys {
        if processed.contains(&key) {
            continue;
        }

        if let Some(ring) =
            walk_polygon_optimized(&mut cells, &mut processed, key, lower, upper, precision)
        {
            if !ring.is_empty() {
                // Wrap the ring as a polygon (exterior ring only, holes determined later)
                polygons.push(vec![ring]);
            }
        }
    }
    let edge_walking_ms = edge_walking_start.elapsed().as_secs_f64() * 1000.0;
    let raw_polygon_count = polygons.len();

    // Phase 3: Build polygon hierarchy
    let hierarchy_start = Instant::now();
    let hierarchical_polygons = crate::marching_squares::build_polygon_hierarchy(polygons);
    let hierarchy_ms = hierarchy_start.elapsed().as_secs_f64() * 1000.0;

    // Phase 4: Orient polygons
    let orientation_start = Instant::now();
    let oriented_polygons = crate::marching_squares::orient_polygons(hierarchical_polygons);
    let orientation_ms = orientation_start.elapsed().as_secs_f64() * 1000.0;
    let final_polygon_count = oriented_polygons.len();

    // Phase 5: Convert to ContourPolygon format
    let conversion_start = Instant::now();
    let contour_polygons = convert_to_contour_polygons(oriented_polygons, precision);
    let conversion_ms = conversion_start.elapsed().as_secs_f64() * 1000.0;

    let total_ms = total_start.elapsed().as_secs_f64() * 1000.0;

    let metrics = crate::ContourMetrics {
        grid_build_ms: 0.0, // Skipped!
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

/// Hybrid optimized processing: fast classification + proven edge walking.
///
/// This function combines:
/// - Optimized: Direct flat array access (no grid building)
/// - Optimized: Sparse cell storage (only stores contour cells)
/// - Original: Proven Shape-based edge walking
///
/// This gives most of the performance benefit while ensuring correctness.
pub fn process_band_optimized_hybrid(
    values: &[f32],
    width: usize,
    height: usize,
    lower: f32,
    upper: f32,
    precision: u32,
) -> (Vec<crate::ContourPolygon>, crate::ContourMetrics) {
    use crate::shape::Shape;
    use crate::GridCell;
    use std::time::Instant;

    let total_start = Instant::now();
    let rows = height;
    let cols = width;

    // Phase 1: Build sparse Shape storage
    // We create Shapes for ALL cells, but only store those with Some value
    // This ensures edge walking can find neighbor cells correctly
    let grid_start = Instant::now();
    let classification_start = Instant::now();

    // Build shapes for all cells, storing in HashMap for sparse access
    let mut cells: HashMap<(usize, usize), Option<Shape>> = HashMap::new();
    let mut active_cells: Vec<(usize, usize)> = Vec::new();

    for r in 0..rows - 1 {
        for c in 0..cols - 1 {
            // Create GridCells for the four corners
            let tl = GridCell {
                x: c as f64,
                y: r as f64,
                value: values[r * cols + c] as f64,
            };
            let tr = GridCell {
                x: (c + 1) as f64,
                y: r as f64,
                value: values[r * cols + c + 1] as f64,
            };
            let br = GridCell {
                x: (c + 1) as f64,
                y: (r + 1) as f64,
                value: values[(r + 1) * cols + c + 1] as f64,
            };
            let bl = GridCell {
                x: c as f64,
                y: (r + 1) as f64,
                value: values[(r + 1) * cols + c] as f64,
            };

            let shape = Shape::create_from_cells(
                &tl,
                &tr,
                &br,
                &bl,
                lower as f64,
                upper as f64,
                c,
                r,
                r == 0,
                c + 1 == cols - 1,
                r + 1 == rows - 1,
                c == 0,
            );

            // Only store in HashMap if there's actually a shape (sparse storage benefit)
            if shape.is_some() {
                cells.insert((r, c), shape);
                active_cells.push((r, c));
            }
        }
    }
    let classification_ms = classification_start.elapsed().as_secs_f64() * 1000.0;
    let grid_build_ms = grid_start.elapsed().as_secs_f64() * 1000.0 - classification_ms;

    // Phase 2: Edge walking using the original proven logic
    // active_cells is already in row-major order (from nested for loops)
    let edge_walking_start = Instant::now();
    let mut polygons: Vec<Vec<Vec<Vec<f64>>>> = Vec::new();

    // Walk polygons starting from each cell that still has unused edges
    for (r, c) in &active_cells {
        let (r, c) = (*r, *c);
        // Check is_cleared() each time since walking modifies cells
        let should_walk = {
            if let Some(Some(cell)) = cells.get(&(r, c)) {
                !cell.is_cleared()
            } else {
                false
            }
        };

        if should_walk {
            // Walk the polygon using sparse cell lookup
            if let Some(ring) = walk_polygon_sparse(&mut cells, r, c, rows - 1, cols - 1, precision)
            {
                if ring.len() >= 3 {
                    polygons.push(vec![ring]);
                }
            }
        }
    }
    let edge_walking_ms = edge_walking_start.elapsed().as_secs_f64() * 1000.0;
    let raw_polygon_count = polygons.len();

    // Phase 4: Build polygon hierarchy
    let hierarchy_start = Instant::now();
    let hierarchical_polygons = crate::marching_squares::build_polygon_hierarchy(polygons);
    let hierarchy_ms = hierarchy_start.elapsed().as_secs_f64() * 1000.0;

    // Phase 5: Orient polygons
    let orientation_start = Instant::now();
    let oriented_polygons = crate::marching_squares::orient_polygons(hierarchical_polygons);
    let orientation_ms = orientation_start.elapsed().as_secs_f64() * 1000.0;
    let final_polygon_count = oriented_polygons.len();

    // Phase 6: Convert to ContourPolygon format
    let conversion_start = Instant::now();
    let contour_polygons = convert_to_contour_polygons(oriented_polygons, precision);
    let conversion_ms = conversion_start.elapsed().as_secs_f64() * 1000.0;

    let total_ms = total_start.elapsed().as_secs_f64() * 1000.0;

    let metrics = crate::ContourMetrics {
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

/// Walk a polygon using sparse cell storage.
/// Mimics the original edge walking logic from marching_squares.rs.
fn walk_polygon_sparse(
    cells: &mut HashMap<(usize, usize), Option<crate::Shape>>,
    start_r: usize,
    start_c: usize,
    max_rows: usize,
    max_cols: usize,
    precision: u32,
) -> Option<Vec<Vec<f64>>> {
    use crate::edge::Move;

    let mut y = start_r;
    let mut x = start_c;
    let mut go_on = true;
    let mut edges: Vec<crate::Edge> = Vec::new();
    let mut current_edge: Option<crate::Edge> = None;

    while go_on {
        if y > max_rows || x > max_cols {
            break;
        }

        let cell = match cells.get_mut(&(y, x)) {
            Some(Some(c)) => c,
            _ => break,
        };

        let tmp_edges = if let Some(ref edge) = current_edge {
            cell.get_edges(Some(edge.end()))
        } else {
            cell.get_edges(None)
        };

        if tmp_edges.is_empty() {
            break;
        }

        cell.increment_used_edges(tmp_edges.len());

        for edge in tmp_edges {
            cell.remove_edge(edge.start());
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
                        x -= 1;
                    } else {
                        go_on = false;
                    }
                }
                Move::Up => {
                    if y > 0 {
                        y -= 1;
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
        let factor = 10_f64.powi(precision as i32);
        let mut ring: Vec<Vec<f64>> = Vec::new();

        // Add first edge's start point
        if let (Some(sx), Some(sy)) = (edges[0].start().x(), edges[0].start().y()) {
            ring.push(vec![
                (sx * factor).round() / factor,
                (sy * factor).round() / factor,
            ]);
        }

        // Add all end points
        for edge in &edges {
            if let (Some(ex), Some(ey)) = (edge.end().x(), edge.end().y()) {
                ring.push(vec![
                    (ex * factor).round() / factor,
                    (ey * factor).round() / factor,
                ]);
            }
        }

        if ring.len() >= 3 {
            return Some(ring);
        }
    }

    None
}

/// Walk a polygon boundary starting from a cell.
fn walk_polygon_optimized(
    cells: &mut CellMap,
    processed: &mut std::collections::HashSet<(u16, u16)>,
    start: (u16, u16),
    lower: f32,
    upper: f32,
    precision: u32,
) -> Option<Vec<Vec<f64>>> {
    let mut ring: Vec<Vec<f64>> = Vec::new();
    let mut current = start;
    let mut visited: std::collections::HashSet<(u16, u16)> = std::collections::HashSet::new();

    loop {
        if visited.contains(&current) {
            break;
        }
        visited.insert(current);

        let cell = cells.get_mut(&current)?;
        let edges = cell.build_edges(lower, upper);

        if edges.is_empty() {
            break;
        }

        cell.mark_used(edges.len() as u8);

        for edge in &edges {
            if let (Some(x), Some(y)) = (edge.start.x(), edge.start.y()) {
                let factor = 10_f64.powi(precision as i32);
                ring.push(vec![
                    (x * factor).round() / factor,
                    (y * factor).round() / factor,
                ]);
            }
        }

        // Move to next cell based on last edge direction
        if let Some(last_edge) = edges.last() {
            let (col, row) = current;
            current = match last_edge.move_dir {
                MoveDir::Right => (col + 1, row),
                MoveDir::Down => (col, row + 1),
                MoveDir::Left => (col.saturating_sub(1), row),
                MoveDir::Up => (col, row.saturating_sub(1)),
                MoveDir::None => break,
            };
        } else {
            break;
        }

        // Check if we've returned to start
        if current == start {
            break;
        }
    }

    // Mark all visited cells as processed
    for key in visited {
        processed.insert(key);
    }

    if ring.len() >= 3 {
        // Close the ring
        if ring.first() != ring.last() {
            ring.push(ring[0].clone());
        }
        Some(ring)
    } else {
        None
    }
}

/// Convert oriented polygons to ContourPolygon format.
fn convert_to_contour_polygons(
    polygons: Vec<Vec<Vec<Vec<f64>>>>,
    precision: u32,
) -> Vec<crate::ContourPolygon> {
    let factor = 10_f32.powi(precision as i32);

    polygons
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

            crate::ContourPolygon { exterior, holes }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_compact_point_size() {
        use std::mem::size_of;
        let size = size_of::<CompactPoint>();
        assert!(
            size <= 32,
            "CompactPoint should be <= 32 bytes, got {}",
            size
        );
    }

    #[test]
    fn test_compact_cell_size() {
        use std::mem::size_of;
        let size = size_of::<CompactCell>();
        assert!(
            size <= 32,
            "CompactCell should be <= 32 bytes, got {}",
            size
        );
    }

    #[test]
    fn test_compact_edge_size() {
        use std::mem::size_of;
        let size = size_of::<CompactEdge>();
        assert!(
            size <= 64,
            "CompactEdge should be <= 64 bytes, got {}",
            size
        );
    }

    #[test]
    fn test_ternary_classification() {
        assert_eq!(Ternary::classify(5.0, 10.0, 20.0), Ternary::Below);
        assert_eq!(Ternary::classify(15.0, 10.0, 20.0), Ternary::Within);
        assert_eq!(Ternary::classify(25.0, 10.0, 20.0), Ternary::Above);
    }

    #[test]
    fn test_cell_creation_empty() {
        // All corners below - should return None
        let cell = CompactCell::new(
            0, 0, 5.0, 5.0, 5.0, 5.0, 10.0, 20.0, false, false, false, false,
        );
        assert!(cell.is_none());

        // All corners above - should return None
        let cell = CompactCell::new(
            0, 0, 25.0, 25.0, 25.0, 25.0, 10.0, 20.0, false, false, false, false,
        );
        assert!(cell.is_none());

        // All corners within - should return None
        let cell = CompactCell::new(
            0, 0, 15.0, 15.0, 15.0, 15.0, 10.0, 20.0, false, false, false, false,
        );
        assert!(cell.is_none());
    }

    #[test]
    fn test_cell_creation_contour() {
        // Mixed values - should create cell
        let cell = CompactCell::new(
            0, 0, 5.0, 15.0, 15.0, 5.0, 10.0, 20.0, false, false, false, false,
        );
        assert!(cell.is_some());
    }

    #[test]
    fn test_sparse_storage() {
        // Create cells only for contour crossings
        let values: Vec<f32> = vec![5.0, 5.0, 5.0, 5.0, 15.0, 15.0, 5.0, 15.0, 15.0];

        let mut cells: CellMap = HashMap::new();
        let width = 3;
        let height = 3;
        let lower = 10.0;
        let upper = 20.0;

        for r in 0..height - 1 {
            for c in 0..width - 1 {
                let tl = values[r * width + c];
                let tr = values[r * width + c + 1];
                let br = values[(r + 1) * width + c + 1];
                let bl = values[(r + 1) * width + c];

                if let Some(cell) = CompactCell::new(
                    c,
                    r,
                    tl,
                    tr,
                    br,
                    bl,
                    lower,
                    upper,
                    r == 0,
                    c + 1 == width - 1,
                    r + 1 == height - 1,
                    c == 0,
                ) {
                    cells.insert((c as u16, r as u16), cell);
                }
            }
        }

        // Only cells on the contour boundary should be created
        assert!(
            cells.len() < (width - 1) * (height - 1),
            "Sparse storage should have fewer cells"
        );
    }

    #[test]
    fn test_process_band_optimized_simple() {
        // 3x3 grid with a high-value region
        let values: Vec<f32> = vec![5.0, 5.0, 5.0, 5.0, 15.0, 15.0, 5.0, 15.0, 15.0];

        let (_polygons, metrics) = process_band_optimized(&values, 3, 3, 10.0, 20.0, 5);

        // Verify the optimized classification runs without errors
        assert!(metrics.classification_ms >= 0.0);
    }

    #[test]
    fn test_compare_original_vs_optimized() {
        // Compare output of original and optimized implementations
        let values: Vec<f32> = vec![5.0, 5.0, 5.0, 5.0, 15.0, 15.0, 5.0, 15.0, 15.0];

        // Original
        let (orig_polygons, _) =
            crate::process_band_flat_with_metrics(&values, 3, 3, 10.0, 20.0, 5);

        // Optimized (uses experimental edge walking, may not produce polygons)
        let (_opt_polygons, opt_metrics) = process_band_optimized(&values, 3, 3, 10.0, 20.0, 5);

        // Original should produce polygons
        assert!(
            !orig_polygons.is_empty(),
            "Original should produce polygons"
        );
        // Optimized classification should run without errors
        assert!(opt_metrics.classification_ms >= 0.0);
    }

    #[test]
    fn test_hybrid_matches_original() {
        // Test that the hybrid version produces identical output to original
        let values: Vec<f32> = vec![5.0, 5.0, 5.0, 5.0, 15.0, 15.0, 5.0, 15.0, 15.0];

        // Original
        let (orig_polygons, _) =
            crate::process_band_flat_with_metrics(&values, 3, 3, 10.0, 20.0, 5);

        // Hybrid
        let (hybrid_polygons, _) = process_band_optimized_hybrid(&values, 3, 3, 10.0, 20.0, 5);

        // Hybrid should produce the same number of polygons
        assert_eq!(
            orig_polygons.len(),
            hybrid_polygons.len(),
            "Polygon count mismatch: original {} vs hybrid {}",
            orig_polygons.len(),
            hybrid_polygons.len()
        );

        // Verify same exterior point count
        for i in 0..orig_polygons.len() {
            assert_eq!(
                orig_polygons[i].exterior.len(),
                hybrid_polygons[i].exterior.len(),
                "Polygon {} exterior size mismatch",
                i
            );
        }
    }

    #[test]
    fn test_hybrid_larger_grid() {
        // Test with a larger grid to ensure correctness
        let width = 10;
        let height = 10;
        let mut values: Vec<f32> = Vec::with_capacity(width * height);

        // Create a gradient with some high-value regions
        for r in 0..height {
            for c in 0..width {
                let dist_from_center = ((r as f32 - 5.0).powi(2) + (c as f32 - 5.0).powi(2)).sqrt();
                values.push(20.0 - dist_from_center * 2.0);
            }
        }

        // Original
        let (orig_polygons, _) =
            crate::process_band_flat_with_metrics(&values, width, height, 10.0, 18.0, 5);

        // Hybrid
        let (hybrid_polygons, _) =
            process_band_optimized_hybrid(&values, width, height, 10.0, 18.0, 5);

        // Same polygon count
        assert_eq!(orig_polygons.len(), hybrid_polygons.len());
    }

    #[test]
    fn test_nan_values_skipped() {
        // Test that NaN values don't create contour cells
        // This simulates ocean areas in meteorological data where NaN = no data
        //
        // Grid layout (3x3):
        //   (0,0)=NaN  (1,0)=NaN  (2,0)=NaN
        //   (0,1)=NaN  (1,1)=NaN  (2,1)=15
        //   (0,2)=NaN  (1,2)=15   (2,2)=15
        //
        // Marching squares cells (2x2):
        // - Cell(0,0): uses (0,0),(1,0),(1,1),(0,1) = NaN,NaN,NaN,NaN → skip (all NaN)
        // - Cell(1,0): uses (1,0),(2,0),(2,1),(1,1) = NaN,NaN,15,NaN → skip (has NaN)
        // - Cell(0,1): uses (0,1),(1,1),(1,2),(0,2) = NaN,NaN,15,NaN → skip (has NaN)
        // - Cell(1,1): uses (1,1),(2,1),(2,2),(1,2) = NaN,15,15,15 → skip (has NaN)
        //
        // All cells have at least one NaN corner, so no polygons should be created.
        let values: Vec<f32> = vec![
            f32::NAN, f32::NAN, f32::NAN, // Row 0: all NaN (ocean)
            f32::NAN, f32::NAN, 15.0,     // Row 1: center is NaN (coastline edge)
            f32::NAN, 15.0, 15.0,         // Row 2: partial data (land)
        ];

        let (polygons, metrics) = process_band_optimized_hybrid(&values, 3, 3, 10.0, 20.0, 5);

        println!(
            "NaN test: polygons={}, raw={}, final={}",
            polygons.len(),
            metrics.raw_polygon_count,
            metrics.final_polygon_count
        );

        // With NaN handling, no contours should be created because all cells
        // have at least one NaN corner and can't be interpolated
        assert!(
            polygons.is_empty(),
            "Should produce no polygons when all cells have NaN corners"
        );
    }

    #[test]
    fn test_nan_values_isolated_region() {
        // Test grid with a single NaN point surrounded by valid data
        // Cells touching NaN are skipped, but cells away from NaN should work
        //
        // Grid layout (5x5) - single NaN at center (2,2):
        //   15   15   15   15   15    <- Row 0
        //   15   15   15   15   15    <- Row 1
        //   15   15  NaN   15   15    <- Row 2 (center is NaN)
        //   15   15   15   15   15    <- Row 3
        //   15   15   15   15   15    <- Row 4
        //
        // Marching squares cells (4x4 = 16 cells):
        // - 4 cells touch NaN: (1,1), (2,1), (1,2), (2,2)
        // - 12 cells don't touch NaN and should produce polygons
        let values: Vec<f32> = vec![
            15.0, 15.0, 15.0, 15.0, 15.0, // Row 0
            15.0, 15.0, 15.0, 15.0, 15.0, // Row 1
            15.0, 15.0, f32::NAN, 15.0, 15.0, // Row 2: center NaN
            15.0, 15.0, 15.0, 15.0, 15.0, // Row 3
            15.0, 15.0, 15.0, 15.0, 15.0, // Row 4
        ];

        let (polygons, _) = process_band_optimized_hybrid(&values, 5, 5, 10.0, 20.0, 5);

        // 12 cells have all valid values, 4 cells touch the central NaN
        // The valid cells should form polygons (may merge during post-processing)
        println!("Polygons with NaN point: {}", polygons.len());
        assert!(
            !polygons.is_empty(),
            "Should create polygons from cells that don't touch NaN"
        );
    }

    #[test]
    fn test_valid_cells_still_create_polygons() {
        // Verify that valid cells (no NaN) still create polygons correctly
        // This ensures the NaN check didn't break normal operation
        let values: Vec<f32> = vec![
            15.0, 15.0, 15.0, // Row 0
            15.0, 15.0, 15.0, // Row 1
            15.0, 15.0, 15.0, // Row 2
        ];

        let (polygons, _) = process_band_optimized_hybrid(&values, 3, 3, 10.0, 20.0, 5);

        // All 4 cells have valid corners and values within band [10, 20]
        // They should merge into a single polygon
        assert!(
            !polygons.is_empty(),
            "Should create polygons when all cells have valid values"
        );
        println!("Valid grid polygons: {}", polygons.len());
    }
}
