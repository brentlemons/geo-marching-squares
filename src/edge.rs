use crate::point::Point;

/// Type of edge (currently unused but mirroring Java structure)
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EdgeType {
    Outside,
    Inside,
}

/// Direction to move to the next cell when following this edge
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Move {
    Right, // column++
    Down,  // row++
    Left,  // column--
    Up,    // row--
    Unknown,
}

/// An edge connecting two points in a contour
///
/// Edges form the boundaries of isoband polygons. Each edge knows its start and
/// end points, and which direction to move to find the next adjacent cell.
#[derive(Debug, Clone, PartialEq)]
pub struct Edge {
    start: Point,
    end: Point,
    edge_type: Option<EdgeType>,
    move_direction: Move,
}

impl Edge {
    /// Create a basic edge with unknown move direction
    pub fn new(start: Point, end: Point) -> Self {
        Self {
            start,
            end,
            edge_type: None,
            move_direction: Move::Unknown,
        }
    }

    /// Create an edge with a specified move direction
    ///
    /// This is the most common constructor. The move direction indicates
    /// which adjacent cell to visit when following the contour.
    pub fn new_with_move(start: Point, end: Point, move_direction: Move) -> Self {
        Self {
            start,
            end,
            edge_type: None,
            move_direction,
        }
    }

    /// Create an edge with a specified edge type
    pub fn new_with_type(start: Point, end: Point, edge_type: EdgeType) -> Self {
        Self {
            start,
            end,
            edge_type: Some(edge_type),
            move_direction: Move::Unknown,
        }
    }

    /// Get the start point of this edge
    pub fn start(&self) -> &Point {
        &self.start
    }

    /// Get the end point of this edge
    pub fn end(&self) -> &Point {
        &self.end
    }

    /// Get the move direction for this edge
    pub fn move_direction(&self) -> &Move {
        &self.move_direction
    }

    /// Get the edge type if set
    pub fn edge_type(&self) -> Option<&EdgeType> {
        self.edge_type.as_ref()
    }

    /// Set the move direction (builder pattern)
    pub fn set_move(mut self, move_direction: Move) -> Self {
        self.move_direction = move_direction;
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_edge_basic_constructor() {
        let start = Point::new(0.0, 0.0);
        let end = Point::new(1.0, 1.0);
        let edge = Edge::new(start.clone(), end.clone());

        assert_eq!(edge.start(), &start);
        assert_eq!(edge.end(), &end);
        assert_eq!(edge.move_direction(), &Move::Unknown);
        assert_eq!(edge.edge_type(), None);
    }

    #[test]
    fn test_edge_with_move() {
        let start = Point::new(0.0, 0.0);
        let end = Point::new(1.0, 0.0);
        let edge = Edge::new_with_move(start, end, Move::Right);

        assert_eq!(edge.move_direction(), &Move::Right);
    }

    #[test]
    fn test_edge_with_type() {
        let start = Point::new(0.0, 0.0);
        let end = Point::new(1.0, 0.0);
        let edge = Edge::new_with_type(start, end, EdgeType::Outside);

        assert_eq!(edge.edge_type(), Some(&EdgeType::Outside));
        assert_eq!(edge.move_direction(), &Move::Unknown);
    }

    #[test]
    fn test_edge_set_move() {
        let start = Point::new(0.0, 0.0);
        let end = Point::new(1.0, 0.0);
        let edge = Edge::new(start, end).set_move(Move::Down);

        assert_eq!(edge.move_direction(), &Move::Down);
    }
}
