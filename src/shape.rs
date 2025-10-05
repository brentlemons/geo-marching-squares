use crate::edge::Edge;
use crate::point::Point;
use std::collections::HashMap;

/// The geometric shape type formed by a cell configuration
///
/// Based on ternary classification (below/within/above threshold),
/// there are 81 possible cell configurations that map to these 7 shape types.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ShapeType {
    Triangle,
    Pentagon,
    Rectangle,
    Trapezoid,
    Hexagon,
    Saddle,
    Square,
}

/// A shape represents a single cell's contribution to an isoband
///
/// Each shape knows its type, corner points, computed edge points,
/// and the edges that form its contribution to the contour polygon.
#[derive(Debug)]
pub struct Shape {
    /// The geometric type of this shape
    shape_type: ShapeType,

    /// Corner points (geographic coordinates)
    top_left: Point,
    top_right: Point,
    bottom_right: Point,
    bottom_left: Point,

    /// Ternary classification value (0-170)
    /// Computed from corner values vs thresholds
    value: u8,

    /// Computed edge points (after interpolation)
    points: Vec<Point>,

    /// Edges keyed by start point for efficient lookup
    edges: HashMap<Point, Edge>,

    /// Grid position (column, row)
    x: usize,
    y: usize,

    /// Has this cell been fully processed?
    cleared: bool,

    /// How many edges have been used?
    used_edges: usize,

    /// Is this cell on a grid boundary?
    top_edge: bool,
    right_edge: bool,
    bottom_edge: bool,
    left_edge: bool,

    /// Corner scalar values
    tl: f64,
    tr: f64,
    bl: f64,
    br: f64,

    /// Isoband thresholds
    lower: f64,
    upper: f64,
}

impl Shape {
    /// Create a new shape
    ///
    /// This is a low-level constructor. Typically shapes are created via
    /// the factory method in Phase 2.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        shape_type: ShapeType,
        top_left: Point,
        top_right: Point,
        bottom_right: Point,
        bottom_left: Point,
        value: u8,
        lower: f64,
        upper: f64,
        x: usize,
        y: usize,
        top_edge: bool,
        right_edge: bool,
        bottom_edge: bool,
        left_edge: bool,
        tl: f64,
        tr: f64,
        bl: f64,
        br: f64,
    ) -> Self {
        Self {
            shape_type,
            top_left,
            top_right,
            bottom_right,
            bottom_left,
            value,
            points: Vec::new(),
            edges: HashMap::new(),
            x,
            y,
            cleared: false,
            used_edges: 0,
            top_edge,
            right_edge,
            bottom_edge,
            left_edge,
            tl,
            tr,
            bl,
            br,
            lower,
            upper,
        }
    }

    /// Get the shape type
    pub fn shape_type(&self) -> ShapeType {
        self.shape_type
    }

    /// Get the ternary classification value
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

    /// Get edges starting from a specific point
    ///
    /// Returns edges in order, following the contour around the cell.
    /// If start is None, begins at the first available edge.
    pub fn get_edges(&self, start: Option<&Point>) -> Vec<Edge> {
        if self.edges.len() <= 1 {
            return self.edges.values().cloned().collect();
        }

        let mut result = Vec::new();
        let mut current = match start {
            Some(s) => s.clone(),
            None => {
                // Find first point that has an edge
                self.points
                    .iter()
                    .find(|p| self.edges.contains_key(p))
                    .cloned()
                    .unwrap_or_else(|| self.points[0].clone())
            }
        };

        while self.edges.contains_key(&current) && result.len() < self.edges.len() {
            let edge = self.edges.get(&current).unwrap().clone();
            current = edge.end().clone();
            result.push(edge);
        }

        result
    }

    /// Check if this shape has been fully processed
    pub fn is_cleared(&self) -> bool {
        self.cleared
    }

    /// Get grid x position (column)
    pub fn x(&self) -> usize {
        self.x
    }

    /// Get grid y position (row)
    pub fn y(&self) -> usize {
        self.y
    }

    /// Get top edge boundary flag
    pub fn is_top_edge(&self) -> bool {
        self.top_edge
    }

    /// Get right edge boundary flag
    pub fn is_right_edge(&self) -> bool {
        self.right_edge
    }

    /// Get bottom edge boundary flag
    pub fn is_bottom_edge(&self) -> bool {
        self.bottom_edge
    }

    /// Get left edge boundary flag
    pub fn is_left_edge(&self) -> bool {
        self.left_edge
    }

    /// Increment the count of used edges
    ///
    /// When all edges have been used, mark the shape as cleared.
    pub fn increment_used_edges(&mut self, count: usize) {
        self.used_edges += count;
        if self.used_edges >= self.edges.len() {
            self.cleared = true;
        }
    }

    /// Remove an edge from the edge map
    pub fn remove_edge(&mut self, key: &Point) {
        self.edges.remove(key);
    }

    /// Set the computed points (used after interpolation)
    pub fn set_points(&mut self, points: Vec<Point>) {
        self.points = points;
    }

    /// Insert an edge into the edge map
    pub fn insert_edge(&mut self, start: Point, edge: Edge) {
        self.edges.insert(start, edge);
    }

    /// Get corner values
    pub fn corner_values(&self) -> (f64, f64, f64, f64) {
        (self.tl, self.tr, self.bl, self.br)
    }

    /// Get thresholds
    pub fn thresholds(&self) -> (f64, f64) {
        (self.lower, self.upper)
    }

    /// Get corner points
    pub fn corner_points(&self) -> (&Point, &Point, &Point, &Point) {
        (&self.top_left, &self.top_right, &self.bottom_right, &self.bottom_left)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_shape_creation() {
        let tl = Point::new(0.0, 2.0);
        let tr = Point::new(1.0, 2.0);
        let br = Point::new(1.0, 1.0);
        let bl = Point::new(0.0, 1.0);

        let shape = Shape::new(
            ShapeType::Triangle,
            tl,
            tr,
            br,
            bl,
            1, // value
            10.0, // lower
            20.0, // upper
            0, // x
            0, // y
            true, false, false, false, // edges
            5.0, 15.0, 5.0, 15.0, // corner values
        );

        assert_eq!(shape.shape_type(), ShapeType::Triangle);
        assert_eq!(shape.value(), 1);
        assert_eq!(shape.x(), 0);
        assert_eq!(shape.y(), 0);
        assert!(!shape.is_cleared());
    }

    #[test]
    fn test_shape_edge_tracking() {
        let tl = Point::new(0.0, 2.0);
        let tr = Point::new(1.0, 2.0);
        let br = Point::new(1.0, 1.0);
        let bl = Point::new(0.0, 1.0);

        let mut shape = Shape::new(
            ShapeType::Rectangle,
            tl,
            tr,
            br,
            bl,
            5,
            10.0,
            20.0,
            0,
            0,
            false, false, false, false,
            5.0, 15.0, 15.0, 5.0,
        );

        assert_eq!(shape.used_edges, 0);
        assert!(!shape.is_cleared());

        // Simulate adding edges
        shape.increment_used_edges(2);
        assert_eq!(shape.used_edges, 2);
    }
}
