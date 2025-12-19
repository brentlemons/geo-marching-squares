//! Isoline segment assembly module
//!
//! This module handles assembling individual cell line segments into complete
//! polylines (either open or closed contour lines).

use crate::cell::LineSegment;
use crate::point::{Point, Side};
use geojson::Position;

/// Default coordinate precision
const DEFAULT_PRECISION: u32 = 5;

/// Round a coordinate value to specified decimal places
fn round_coord_with_precision(value: f64, precision: u32) -> f64 {
    let factor = 10_f64.powi(precision as i32);
    (value * factor).round() / factor
}


/// Check if two points are equal within a small epsilon
fn points_equal(p1: &Point, p2: &Point) -> bool {
    const EPSILON: f64 = 1e-9;
    (p1.x().unwrap() - p2.x().unwrap()).abs() < EPSILON
        && (p1.y().unwrap() - p2.y().unwrap()).abs() < EPSILON
}

/// Assembles cell line segments into complete polylines
pub struct IsolineAssembler {
    segments: Vec<(usize, usize, LineSegment)>, // (row, col, segment)
}

impl IsolineAssembler {
    /// Create a new assembler
    pub fn new() -> Self {
        Self {
            segments: Vec::new(),
        }
    }

    /// Add segments from a cell at position (r, c)
    pub fn add_cell_segments(&mut self, r: usize, c: usize, segments: Vec<LineSegment>) {
        for segment in segments {
            self.segments.push((r, c, segment));
        }
    }

    /// Assemble all segments into polylines
    pub fn assemble(&mut self) -> Vec<Vec<Position>> {
        self.assemble_with_precision(DEFAULT_PRECISION)
    }

    /// Assemble all segments into polylines with configurable precision
    pub fn assemble_with_precision(&mut self, precision: u32) -> Vec<Vec<Position>> {
        let mut polylines = Vec::new();
        let mut used = vec![false; self.segments.len()];

        for start_idx in 0..self.segments.len() {
            if used[start_idx] {
                continue;
            }

            // Start a new polyline
            let polyline = self.trace_polyline_with_precision(start_idx, &mut used, precision);
            if !polyline.is_empty() {
                polylines.push(polyline);
            }
        }

        polylines
    }

    /// Trace a single polyline with configurable precision
    fn trace_polyline_with_precision(
        &self,
        start_idx: usize,
        used: &mut [bool],
        precision: u32,
    ) -> Vec<Position> {
        let mut polyline = Vec::new();
        let mut current_idx = start_idx;

        used[current_idx] = true;
        let (_, _, segment) = &self.segments[current_idx];

        // Add start point
        polyline.push(vec![
            round_coord_with_precision(segment.start.x().unwrap(), precision),
            round_coord_with_precision(segment.start.y().unwrap(), precision),
        ]);

        // Add end point
        polyline.push(vec![
            round_coord_with_precision(segment.end.x().unwrap(), precision),
            round_coord_with_precision(segment.end.y().unwrap(), precision),
        ]);

        // Try to extend polyline forward
        loop {
            let (current_r, current_c, current_segment) = &self.segments[current_idx];
            let last_point = &current_segment.end;
            let last_side = &current_segment.end_side;

            // Find next connecting segment
            if let Some(next_idx) =
                self.find_connecting_segment(*current_r, *current_c, last_point, last_side, used)
            {
                used[next_idx] = true;
                let next_segment = &self.segments[next_idx].2;

                polyline.push(vec![
                    round_coord_with_precision(next_segment.end.x().unwrap(), precision),
                    round_coord_with_precision(next_segment.end.y().unwrap(), precision),
                ]);

                // Check if loop closed
                if points_equal(&next_segment.end, &self.segments[start_idx].2.start) {
                    break; // Closed loop
                }

                current_idx = next_idx;
            } else {
                break; // No more connections - open polyline
            }
        }

        polyline
    }

    /// Find a segment that connects to the given point/side
    fn find_connecting_segment(
        &self,
        current_r: usize,
        current_c: usize,
        point: &Point,
        side: &Side,
        used: &[bool],
    ) -> Option<usize> {
        // Calculate which neighbor cell to check based on exit side
        let (neighbor_r, neighbor_c) = match side {
            Side::Top => {
                if current_r == 0 {
                    return None; // At grid boundary
                }
                (current_r - 1, current_c)
            }
            Side::Right => (current_r, current_c + 1),
            Side::Bottom => (current_r + 1, current_c),
            Side::Left => {
                if current_c == 0 {
                    return None; // At grid boundary
                }
                (current_r, current_c - 1)
            }
        };

        // Find unused segment in neighbor cell that starts at this point
        for (idx, (r, c, segment)) in self.segments.iter().enumerate() {
            if used[idx] {
                continue;
            }

            if *r == neighbor_r && *c == neighbor_c {
                if points_equal(&segment.start, point) {
                    return Some(idx);
                }
            }
        }

        None
    }
}

impl Default for IsolineAssembler {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::point::Point;

    #[test]
    fn test_round_coord() {
        // Test with default precision (5 decimals)
        assert_eq!(round_coord_with_precision(1.234567, 5), 1.23457);
        assert_eq!(round_coord_with_precision(1.234564, 5), 1.23456);
        assert_eq!(round_coord_with_precision(-122.123456789, 5), -122.12346);
        // Test with lower precision
        assert_eq!(round_coord_with_precision(1.234567, 3), 1.235);
        assert_eq!(round_coord_with_precision(-122.123456, 2), -122.12);
    }

    #[test]
    fn test_points_equal() {
        let p1 = Point::new(1.0, 2.0);
        let p2 = Point::new(1.0, 2.0);
        let p3 = Point::new(1.0001, 2.0);

        assert!(points_equal(&p1, &p2));
        assert!(!points_equal(&p1, &p3));
    }

    #[test]
    fn test_assembler_single_segment() {
        let mut assembler = IsolineAssembler::new();

        let segment = LineSegment {
            start: Point::new(0.0, 0.0),
            end: Point::new(1.0, 1.0),
            start_side: Side::Left,
            end_side: Side::Right,
        };

        assembler.add_cell_segments(0, 0, vec![segment]);

        let polylines = assembler.assemble();

        assert_eq!(polylines.len(), 1);
        assert_eq!(polylines[0].len(), 2);
        assert_eq!(polylines[0][0], vec![0.0, 0.0]);
        assert_eq!(polylines[0][1], vec![1.0, 1.0]);
    }
}
