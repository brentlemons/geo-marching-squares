//! Main marching squares algorithm implementation
//!
//! This module implements the core algorithm that converts a 2D grid of
//! GeoJSON features into contour polygons (isobands).

use crate::edge::{Edge, Move};
use crate::shape::Shape;
use geojson::{Feature, Geometry, JsonObject, Position, Value as GeoValue};
use std::collections::VecDeque;

/// Round a coordinate value to 5 decimal places (approximately 1 meter accuracy)
///
/// Matches Java's BigDecimal.setScale(5, RoundingMode.HALF_UP)
fn round_coord(value: f64) -> f64 {
    (value * 100000.0).round() / 100000.0
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
                                    eprintln!("Unknown edge move at ({}, {})", y, x);
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

    // Resolve polygon nesting (exterior vs interior rings)
    let mut polygons: Vec<Vec<Vec<Position>>> = Vec::new();

    while let Some(subject) = hold_polygons.pop_front() {
        let mut external = true;

        for i in (0..polygons.len()).rev() {
            let polygon = &polygons[i];
            let mut push_out = false;

            // Case 1: subject is inside polygon
            if polygon_in_polygon(&subject, polygon) {
                // Check if subject is in a hole (should be pushed out)
                for j in 1..polygon.len() {
                    let hole = &polygon[j..j + 1];
                    if polygon_in_polygon(&subject, &[hole[0].clone()]) {
                        push_out = true;
                        break;
                    }
                }

                if !push_out {
                    // Add subject as interior ring (hole) to polygon
                    let mut new_polygon = polygon.clone();
                    new_polygon.push(subject[0].clone());
                    polygons[i] = new_polygon;
                    external = false;
                    break;
                }
            }
            // Case 2: polygon is inside subject
            else if polygon_in_polygon(polygon, &subject) {
                // Break apart polygon and reprocess
                if polygon.len() > 1 {
                    // Has interior rings
                    for j in 1..polygon.len() {
                        hold_polygons.push_back(vec![polygon[j].clone()]);
                    }
                    hold_polygons.push_back(vec![polygon[0].clone()]);
                } else {
                    hold_polygons.push_back(polygon.clone());
                }
                polygons.remove(i);
            }
        }

        if external {
            polygons.push(subject);
        }
    }

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
