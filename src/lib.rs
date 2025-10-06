//! # geo-marching-squares
//!
//! A Rust implementation of the marching squares algorithm for generating contour polygons
//! from geospatial data.
//!
//! This library generates **isobands** (filled contour polygons) and **isolines** (contour lines)
//! from 2D scalar fields using geographic coordinates (latitude/longitude).
//!
//! ## Examples
//!
//! ### Concurrent Processing (Recommended)
//!
//! Process multiple isoband levels in parallel for maximum performance:
//!
//! ```rust,ignore
//! use geo_marching_squares::do_concurrent;
//! use geojson::{Feature, FeatureCollection};
//!
//! // Your 2D grid of GeoJSON point features with scalar values
//! let grid: Vec<Vec<Feature>> = load_grid_data();
//!
//! // Define isoband thresholds (creates N-1 bands from N thresholds)
//! let thresholds = vec![0.0, 10.0, 20.0, 30.0, 40.0];
//!
//! // Generate all isobands in parallel using Rayon
//! let result: FeatureCollection = do_concurrent(&grid, &thresholds);
//!
//! // Result contains 4 isobands: 0-10, 10-20, 20-30, 30-40
//! println!("Generated {} isoband features", result.features.len());
//! ```
//!
//! ### Single Isoband Processing
//!
//! Process a single isoband level:
//!
//! ```rust,ignore
//! use geo_marching_squares::process_band;
//! use geojson::Feature;
//!
//! let grid: Vec<Vec<Feature>> = load_grid_data();
//!
//! // Generate contour for values between 10 and 20
//! let feature: Feature = process_band(&grid, 10.0, 20.0);
//!
//! // Feature contains MultiPolygon geometry with all contours for this band
//! ```
//!
//! ## Isolines (Contour Lines)
//!
//! Generate contour lines at specific threshold values:
//!
//! ```rust,ignore
//! use geo_marching_squares::{process_line, do_concurrent_lines};
//!
//! let grid = load_grid_data();
//!
//! // Single isoline at value 15.0
//! let feature = process_line(&grid, 15.0);
//! // Returns Feature with MultiLineString geometry
//!
//! // Multiple isolines in parallel
//! let isovalues = vec![0.0, 10.0, 20.0, 30.0];
//! let collection = do_concurrent_lines(&grid, &isovalues);
//! // Returns FeatureCollection with 4 isoline features
//! ```
//!
//! ### Isolines vs Isobands
//!
//! | Feature | Isolines | Isobands |
//! |---------|----------|----------|
//! | **Classification** | Binary (above/below) | Ternary (below/within/above) |
//! | **Configurations** | 16 possible (2^4) | 81 possible (3^4) |
//! | **Geometry** | MultiLineString | MultiPolygon |
//! | **Interpolation** | Linear | Cosine-smoothed |
//! | **Output** | Contour lines | Filled regions |
//!
//! ## Performance
//!
//! - **Parallel processing**: Uses Rayon's work-stealing thread pool
//! - **Embarrassingly parallel**: Each level computed independently
//! - **Zero-copy**: Input grid shared across all threads (read-only)
//! - **Efficient**: Rust's zero-cost abstractions ensure optimal performance

mod point;
mod edge;
mod shape;
mod cell;
mod isoline_assembler;
mod marching_squares;

pub use point::{Point, Side};
pub use edge::{Edge, EdgeType, Move};
pub use shape::{Shape, ShapeType};
pub use cell::{Cell, LineSegment};
pub use marching_squares::{process_band, do_concurrent, process_line, do_concurrent_lines};

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    #[test]
    fn test_point_creation() {
        let p1 = Point::new(10.0, 20.0);
        assert_eq!(p1.x(), Some(10.0));
        assert_eq!(p1.y(), Some(20.0));

        let p2 = Point::new_with_limit(15.0, 10.0, Side::Top);
        assert_eq!(p2.x(), None);
        assert_eq!(p2.value(), Some(15.0));
    }

    #[test]
    fn test_edge_creation() {
        let start = Point::new(0.0, 0.0);
        let end = Point::new(1.0, 1.0);
        let edge = Edge::new_with_move(start, end, Move::Right);

        assert_eq!(edge.move_direction(), &Move::Right);
    }

    #[test]
    fn test_cell_binary_classification() {
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
            1,    // value
            10.0, // lower
            20.0, // upper
            0,    // x
            0,    // y
            true, false, false, false, // edges
            5.0, 15.0, 5.0, 15.0, // corner values
        );

        assert_eq!(shape.shape_type(), ShapeType::Triangle);
        assert_eq!(shape.value(), 1);
    }

    #[test]
    fn test_point_hashmap_usage() {
        let mut map = HashMap::new();
        let p1 = Point::new(5.0, 10.0);
        let p2 = Point::new(5.0, 10.0);

        map.insert(p1.clone(), "value1");
        map.insert(p2.clone(), "value2");

        // Should overwrite since p1 == p2
        assert_eq!(map.len(), 1);
        assert_eq!(map.get(&p1), Some(&"value2"));
    }

    #[test]
    fn test_shape_types_all_variants() {
        // Ensure all shape type variants exist
        let types = vec![
            ShapeType::Triangle,
            ShapeType::Pentagon,
            ShapeType::Rectangle,
            ShapeType::Trapezoid,
            ShapeType::Hexagon,
            ShapeType::Saddle,
            ShapeType::Square,
        ];

        assert_eq!(types.len(), 7);
    }

    #[test]
    fn test_move_directions() {
        // Ensure all move directions exist
        let moves = vec![
            Move::Right,
            Move::Down,
            Move::Left,
            Move::Up,
            Move::Unknown,
        ];

        assert_eq!(moves.len(), 5);
    }
}
