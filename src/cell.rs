use crate::edge::Edge;
use crate::point::{Point, Side};
use std::collections::HashMap;

/// Represents a line segment within a cell for isoline generation
#[derive(Debug, Clone)]
pub struct LineSegment {
    pub start: Point,
    pub end: Point,
    pub start_side: Side,
    pub end_side: Side,
}

/// A cell for isoline (contour line) generation
///
/// Unlike Shape which uses ternary classification for isobands,
/// Cell uses binary classification (above/below threshold) for isolines.
/// This results in 16 possible configurations instead of 81.
#[derive(Debug, Clone)]
pub struct Cell {
    /// Corner points
    top_left: Point,
    top_right: Point,
    bottom_right: Point,
    bottom_left: Point,

    /// Binary classification value (0-15)
    /// Bit pattern: tl(8) | tr(4) | br(2) | bl(1)
    value: u8,

    /// Computed edge points (after interpolation)
    points: Vec<Point>,

    /// Edges for the contour line
    edges: HashMap<Point, Edge>,

    /// Threshold value for the isoline
    isovalue: f64,

    /// Corner scalar values
    tl: f64,
    tr: f64,
    bl: f64,
    br: f64,
}

impl Cell {
    /// Create a new cell with binary classification
    ///
    /// Classifies each corner as 0 (below isovalue) or 1 (above isovalue)
    /// and computes the 4-bit cell configuration value.
    pub fn new(
        top_left: Point,
        top_right: Point,
        bottom_right: Point,
        bottom_left: Point,
        isovalue: f64,
        tl: f64,
        tr: f64,
        bl: f64,
        br: f64,
    ) -> Self {
        // Binary classification: 0 if below, 1 if above/equal
        // Bit positions: TL=8, TR=4, BR=2, BL=1
        let mut value = 0u8;
        value |= if tl < isovalue { 0 } else { 8 };
        value |= if tr < isovalue { 0 } else { 4 };
        value |= if br < isovalue { 0 } else { 2 };
        value |= if bl < isovalue { 0 } else { 1 };

        Self {
            top_left,
            top_right,
            bottom_right,
            bottom_left,
            value,
            points: Vec::new(),
            edges: HashMap::new(),
            isovalue,
            tl,
            tr,
            bl,
            br,
        }
    }

    /// Get the binary classification value (0-15)
    pub fn value(&self) -> u8 {
        self.value
    }

    /// Get the computed edge points
    pub fn points(&self) -> &[Point] {
        &self.points
    }

    /// Get all edges
    pub fn edges(&self) -> &HashMap<Point, Edge> {
        &self.edges
    }

    /// Get the isovalue threshold
    pub fn isovalue(&self) -> f64 {
        self.isovalue
    }

    /// Get corner values
    pub fn corner_values(&self) -> (f64, f64, f64, f64) {
        (self.tl, self.tr, self.bl, self.br)
    }

    /// Get corner points
    pub fn corner_points(&self) -> (&Point, &Point, &Point, &Point) {
        (
            &self.top_left,
            &self.top_right,
            &self.bottom_right,
            &self.bottom_left,
        )
    }

    /// Set the computed points (used after interpolation)
    pub fn set_points(&mut self, points: Vec<Point>) {
        self.points = points;
    }

    /// Insert an edge into the edge map
    pub fn insert_edge(&mut self, start: Point, edge: Edge) {
        self.edges.insert(start, edge);
    }

    /// Create a cell from GeoJSON features using binary classification
    ///
    /// This is the primary way to create cells. It:
    /// 1. Extracts coordinates and values from GeoJSON features
    /// 2. Computes binary classification
    /// 3. Returns None for empty/full cells (0, 15)
    /// 4. Generates interpolated points and line segments
    pub fn create(
        top_left_feature: &geojson::Feature,
        top_right_feature: &geojson::Feature,
        bottom_right_feature: &geojson::Feature,
        bottom_left_feature: &geojson::Feature,
        isovalue: f64,
    ) -> Option<Self> {
        // Extract coordinates
        let tl_coords = Self::extract_coords(top_left_feature)?;
        let tr_coords = Self::extract_coords(top_right_feature)?;
        let br_coords = Self::extract_coords(bottom_right_feature)?;
        let bl_coords = Self::extract_coords(bottom_left_feature)?;

        // Extract values
        let tl_val = Self::extract_value(top_left_feature)?;
        let tr_val = Self::extract_value(top_right_feature)?;
        let br_val = Self::extract_value(bottom_right_feature)?;
        let bl_val = Self::extract_value(bottom_left_feature)?;

        // Binary classification
        let mut value = 0u8;
        value |= if tl_val < isovalue { 0 } else { 8 };
        value |= if tr_val < isovalue { 0 } else { 4 };
        value |= if br_val < isovalue { 0 } else { 2 };
        value |= if bl_val < isovalue { 0 } else { 1 };

        // Skip empty cases (0 = all below, 15 = all above)
        if value == 0 || value == 15 {
            return None;
        }

        // Create cell
        let mut cell = Self::new(
            Point::new(tl_coords[0], tl_coords[1]),
            Point::new(tr_coords[0], tr_coords[1]),
            Point::new(br_coords[0], br_coords[1]),
            Point::new(bl_coords[0], bl_coords[1]),
            isovalue,
            tl_val,
            tr_val,
            bl_val,
            br_val,
        );

        // Generate interpolated points
        cell.generate_points();

        Some(cell)
    }

    /// Create a cell from GridCell data (memory-efficient alternative)
    ///
    /// This is the high-performance alternative to `create()` for large grids.
    /// Uses lightweight GridCell structs instead of full GeoJSON Features,
    /// reducing memory usage by ~12x.
    ///
    /// # Arguments
    ///
    /// * `top_left` - Top-left grid cell
    /// * `top_right` - Top-right grid cell
    /// * `bottom_right` - Bottom-right grid cell
    /// * `bottom_left` - Bottom-left grid cell
    /// * `isovalue` - Threshold value for the isoline
    pub fn create_from_cells(
        top_left: &crate::GridCell,
        top_right: &crate::GridCell,
        bottom_right: &crate::GridCell,
        bottom_left: &crate::GridCell,
        isovalue: f64,
    ) -> Option<Self> {
        let tl_val = top_left.value;
        let tr_val = top_right.value;
        let br_val = bottom_right.value;
        let bl_val = bottom_left.value;

        // Binary classification
        let mut value = 0u8;
        value |= if tl_val < isovalue { 0 } else { 8 };
        value |= if tr_val < isovalue { 0 } else { 4 };
        value |= if br_val < isovalue { 0 } else { 2 };
        value |= if bl_val < isovalue { 0 } else { 1 };

        // Skip empty cases (0 = all below, 15 = all above)
        if value == 0 || value == 15 {
            return None;
        }

        // Create cell
        let mut cell = Self::new(
            Point::new(top_left.x, top_left.y),
            Point::new(top_right.x, top_right.y),
            Point::new(bottom_right.x, bottom_right.y),
            Point::new(bottom_left.x, bottom_left.y),
            isovalue,
            tl_val,
            tr_val,
            bl_val,
            br_val,
        );

        // Generate interpolated points
        cell.generate_points();

        Some(cell)
    }

    /// Extract coordinates from a GeoJSON Point feature
    fn extract_coords(feature: &geojson::Feature) -> Option<Vec<f64>> {
        if let Some(geojson::Geometry {
            value: geojson::Value::Point(coords),
            ..
        }) = &feature.geometry
        {
            Some(coords.clone())
        } else {
            None
        }
    }

    /// Extract the "value" property from a GeoJSON feature
    fn extract_value(feature: &geojson::Feature) -> Option<f64> {
        feature.properties.as_ref()?.get("value")?.as_f64()
    }

    /// Check if top side has no contour crossing
    fn is_top_blank(&self) -> bool {
        ((self.tl >= self.isovalue) && (self.tr >= self.isovalue))
            || ((self.tl < self.isovalue) && (self.tr < self.isovalue))
    }

    /// Check if right side has no contour crossing
    fn is_right_blank(&self) -> bool {
        ((self.tr >= self.isovalue) && (self.br >= self.isovalue))
            || ((self.tr < self.isovalue) && (self.br < self.isovalue))
    }

    /// Check if bottom side has no contour crossing
    fn is_bottom_blank(&self) -> bool {
        ((self.bl >= self.isovalue) && (self.br >= self.isovalue))
            || ((self.bl < self.isovalue) && (self.br < self.isovalue))
    }

    /// Check if left side has no contour crossing
    fn is_left_blank(&self) -> bool {
        ((self.tl >= self.isovalue) && (self.bl >= self.isovalue))
            || ((self.tl < self.isovalue) && (self.bl < self.isovalue))
    }

    /// Linear interpolation for isolines (no cosine smoothing)
    fn interpolate(level: f64, value0: f64, value1: f64, point0: &Point, point1: &Point) -> Point {
        let mu = (level - value0) / (value1 - value0);
        let x = ((1.0 - mu) * point0.x().unwrap()) + (mu * point1.x().unwrap());
        let y = ((1.0 - mu) * point0.y().unwrap()) + (mu * point1.y().unwrap());
        Point::new(x, y)
    }

    /// Interpolate a point on a specific side
    fn interpolate_side(&self, side: Side) -> Point {
        match side {
            Side::Top => Self::interpolate(
                self.isovalue,
                self.tl,
                self.tr,
                &self.top_left,
                &self.top_right,
            ),
            Side::Right => Self::interpolate(
                self.isovalue,
                self.tr,
                self.br,
                &self.top_right,
                &self.bottom_right,
            ),
            Side::Bottom => Self::interpolate(
                self.isovalue,
                self.bl,
                self.br,
                &self.bottom_left,
                &self.bottom_right,
            ),
            Side::Left => Self::interpolate(
                self.isovalue,
                self.tl,
                self.bl,
                &self.top_left,
                &self.bottom_left,
            ),
        }
    }

    /// Generate interpolated points for sides that cross the threshold
    fn generate_points(&mut self) {
        let mut points = Vec::new();

        if !self.is_top_blank() {
            points.push((Side::Top, self.interpolate_side(Side::Top)));
        }
        if !self.is_right_blank() {
            points.push((Side::Right, self.interpolate_side(Side::Right)));
        }
        if !self.is_bottom_blank() {
            points.push((Side::Bottom, self.interpolate_side(Side::Bottom)));
        }
        if !self.is_left_blank() {
            points.push((Side::Left, self.interpolate_side(Side::Left)));
        }

        self.points = points.into_iter().map(|(_, p)| p).collect();
    }

    /// Get line segments for this cell based on its configuration
    pub fn get_line_segments(&self) -> Vec<LineSegment> {
        let mut segments = Vec::new();

        match self.value {
            0 | 15 => {} // No segments

            1 => segments.push(self.create_segment(Side::Left, Side::Bottom)),
            2 => segments.push(self.create_segment(Side::Bottom, Side::Right)),
            3 => segments.push(self.create_segment(Side::Left, Side::Right)),
            4 => segments.push(self.create_segment(Side::Right, Side::Top)),
            6 => segments.push(self.create_segment(Side::Bottom, Side::Top)),
            7 => segments.push(self.create_segment(Side::Left, Side::Top)),
            8 => segments.push(self.create_segment(Side::Top, Side::Left)),
            9 => segments.push(self.create_segment(Side::Top, Side::Bottom)),
            11 => segments.push(self.create_segment(Side::Top, Side::Right)),
            12 => segments.push(self.create_segment(Side::Right, Side::Left)),
            13 => segments.push(self.create_segment(Side::Bottom, Side::Left)),
            14 => segments.push(self.create_segment(Side::Right, Side::Bottom)),

            // Saddle cases - two separate segments
            5 => {
                segments.push(self.create_segment(Side::Left, Side::Bottom));
                segments.push(self.create_segment(Side::Top, Side::Right));
            }
            10 => {
                segments.push(self.create_segment(Side::Left, Side::Top));
                segments.push(self.create_segment(Side::Bottom, Side::Right));
            }

            _ => unreachable!(),
        }

        segments
    }

    /// Create a line segment between two sides
    fn create_segment(&self, start_side: Side, end_side: Side) -> LineSegment {
        let start = self.get_point_for_side(start_side).clone();
        let end = self.get_point_for_side(end_side).clone();

        LineSegment {
            start,
            end,
            start_side,
            end_side,
        }
    }

    /// Get the interpolated point for a specific side
    fn get_point_for_side(&self, side: Side) -> &Point {
        // Points are stored in order: Top, Right, Bottom, Left (if they exist)
        // We need to find which index corresponds to this side
        let mut point_index = 0;

        if !self.is_top_blank() && side == Side::Top {
            return &self.points[point_index];
        }
        if !self.is_top_blank() {
            point_index += 1;
        }

        if !self.is_right_blank() && side == Side::Right {
            return &self.points[point_index];
        }
        if !self.is_right_blank() {
            point_index += 1;
        }

        if !self.is_bottom_blank() && side == Side::Bottom {
            return &self.points[point_index];
        }
        if !self.is_bottom_blank() {
            point_index += 1;
        }

        if !self.is_left_blank() && side == Side::Left {
            return &self.points[point_index];
        }

        panic!("Point not found for side {:?}", side);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cell_binary_classification_all_below() {
        let tl = Point::new(0.0, 2.0);
        let tr = Point::new(1.0, 2.0);
        let br = Point::new(1.0, 1.0);
        let bl = Point::new(0.0, 1.0);

        // All values below threshold
        let cell = Cell::new(tl, tr, br, bl, 10.0, 5.0, 5.0, 5.0, 5.0);

        // Binary: 0000 = 0
        assert_eq!(cell.value(), 0);
    }

    #[test]
    fn test_cell_binary_classification_all_above() {
        let tl = Point::new(0.0, 2.0);
        let tr = Point::new(1.0, 2.0);
        let br = Point::new(1.0, 1.0);
        let bl = Point::new(0.0, 1.0);

        // All values above threshold
        let cell = Cell::new(tl, tr, br, bl, 10.0, 15.0, 15.0, 15.0, 15.0);

        // Binary: 1111 = 15
        assert_eq!(cell.value(), 15);
    }

    #[test]
    fn test_cell_binary_classification_mixed() {
        let tl = Point::new(0.0, 2.0);
        let tr = Point::new(1.0, 2.0);
        let br = Point::new(1.0, 1.0);
        let bl = Point::new(0.0, 1.0);

        // Values: tl=5, tr=15, br=15, bl=5, isovalue=10
        let cell = Cell::new(tl, tr, br, bl, 10.0, 5.0, 15.0, 5.0, 15.0);

        // Binary: tl(0), tr(1), br(1), bl(0) = 0110 = 6
        assert_eq!(cell.value(), 6);
    }

    #[test]
    fn test_cell_getters() {
        let tl = Point::new(0.0, 2.0);
        let tr = Point::new(1.0, 2.0);
        let br = Point::new(1.0, 1.0);
        let bl = Point::new(0.0, 1.0);

        let cell = Cell::new(tl, tr, br, bl, 10.0, 5.0, 15.0, 5.0, 15.0);

        assert_eq!(cell.isovalue(), 10.0);
        assert_eq!(cell.corner_values(), (5.0, 15.0, 5.0, 15.0));
        assert_eq!(cell.points().len(), 0); // No points computed yet
    }

    #[test]
    fn test_cell_edge_cases() {
        let tl = Point::new(0.0, 2.0);
        let tr = Point::new(1.0, 2.0);
        let br = Point::new(1.0, 1.0);
        let bl = Point::new(0.0, 1.0);

        // Corner exactly on threshold - should be classified as "above"
        let cell = Cell::new(tl, tr, br, bl, 10.0, 10.0, 15.0, 5.0, 15.0);

        // Binary: tl(1), tr(1), br(1), bl(0) = 1110 = 14
        assert_eq!(cell.value(), 14);
    }
}
