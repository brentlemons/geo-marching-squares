use crate::edge::Edge;
use crate::point::Point;
use std::collections::HashMap;

/// A cell for isoline (contour line) generation
///
/// Unlike Shape which uses ternary classification for isobands,
/// Cell uses binary classification (above/below threshold) for isolines.
/// This results in 16 possible configurations instead of 81.
#[derive(Debug)]
pub struct Cell {
    /// Corner points (geographic coordinates)
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
        (&self.top_left, &self.top_right, &self.bottom_right, &self.bottom_left)
    }

    /// Set the computed points (used after interpolation)
    pub fn set_points(&mut self, points: Vec<Point>) {
        self.points = points;
    }

    /// Insert an edge into the edge map
    pub fn insert_edge(&mut self, start: Point, edge: Edge) {
        self.edges.insert(start, edge);
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
