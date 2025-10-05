use crate::edge::{Edge, Move};
use crate::point::{Point, Side};
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
    /// Create a shape from GeoJSON features using the factory pattern
    ///
    /// This is the primary way to create shapes. It:
    /// 1. Extracts coordinates and values from GeoJSON features
    /// 2. Computes ternary classification
    /// 3. Determines appropriate ShapeType
    /// 4. Builds edges for that shape type
    /// 5. Returns None for empty/full cells
    #[allow(clippy::too_many_arguments)]
    pub fn create(
        top_left_feature: &geojson::Feature,
        top_right_feature: &geojson::Feature,
        bottom_right_feature: &geojson::Feature,
        bottom_left_feature: &geojson::Feature,
        lower: f64,
        upper: f64,
        x: usize,
        y: usize,
        top_edge: bool,
        right_edge: bool,
        bottom_edge: bool,
        left_edge: bool,
    ) -> Option<Self> {
        // Extract coordinates
        let tl_coords = Self::extract_coords(top_left_feature)?;
        let tr_coords = Self::extract_coords(top_right_feature)?;
        let br_coords = Self::extract_coords(bottom_right_feature)?;
        let bl_coords = Self::extract_coords(bottom_left_feature)?;

        // Extract values
        let tl = Self::extract_value(top_left_feature)?;
        let tr = Self::extract_value(top_right_feature)?;
        let br = Self::extract_value(bottom_right_feature)?;
        let bl = Self::extract_value(bottom_left_feature)?;

        // Create corner points
        let top_left = Point::new(tl_coords.0, tl_coords.1);
        let top_right = Point::new(tr_coords.0, tr_coords.1);
        let bottom_right = Point::new(br_coords.0, br_coords.1);
        let bottom_left = Point::new(bl_coords.0, bl_coords.1);

        // Compute ternary classification value using bitwise encoding
        // This matches Java's exact logic in Shape.java lines 52-55
        let mut value = 0u8;
        value |= if tl < lower { 0 } else if tl >= upper { 128 } else { 64 };
        value |= if tr < lower { 0 } else if tr >= upper { 32 } else { 16 };
        value |= if br < lower { 0 } else if br >= upper { 8 } else { 4 };
        value |= if bl < lower { 0 } else if bl >= upper { 2 } else { 1 };

        // Determine shape type from value (matches Java's switch statement)
        let shape_type = match value {
            0 | 170 => return None, // Empty or full cell

            169 | 166 | 154 | 106 | 1 | 4 | 16 | 64 => ShapeType::Triangle,

            101 | 149 | 86 | 89 | 69 | 21 | 84 | 81 | 96 | 24 | 6 | 129 | 74 | 146 | 164 | 41 |
            66 | 144 | 36 | 9 | 104 | 26 | 134 | 161 => ShapeType::Pentagon,

            5 | 20 | 80 | 65 | 165 | 150 | 90 | 105 | 160 | 130 | 10 | 40 => ShapeType::Rectangle,

            168 | 2 | 162 | 8 | 138 | 32 | 42 | 128 => ShapeType::Trapezoid,

            37 | 133 | 148 | 22 | 82 | 88 | 73 | 97 | 145 | 25 | 70 | 100 => ShapeType::Hexagon,

            153 | 102 | 68 | 17 | 136 | 34 | 152 | 18 | 137 | 33 | 98 | 72 | 38 | 132 => ShapeType::Saddle,

            85 => ShapeType::Square,

            _ => {
                eprintln!("Unknown shape value: {}", value);
                return None;
            }
        };

        let mut shape = Self {
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
        };

        // Generate points and build edges
        shape.points = shape.get_points();
        shape.build_edges();

        Some(shape)
    }

    /// Extract coordinates from a GeoJSON feature
    fn extract_coords(feature: &geojson::Feature) -> Option<(f64, f64)> {
        if let Some(geojson::Geometry { value: geojson::Value::Point(coords), .. }) = &feature.geometry {
            Some((coords[0], coords[1]))
        } else {
            None
        }
    }

    /// Extract value property from a GeoJSON feature
    fn extract_value(feature: &geojson::Feature) -> Option<f64> {
        feature.property("value")
            .and_then(|v| v.as_f64())
    }

    /// Build edges for this shape based on its type and value
    ///
    /// This implements the exact edge construction logic from the Java subclasses
    fn build_edges(&mut self) {
        match self.shape_type {
            ShapeType::Triangle => self.build_triangle_edges(),
            ShapeType::Pentagon => self.build_pentagon_edges(),
            ShapeType::Rectangle => self.build_rectangle_edges(),
            ShapeType::Trapezoid => self.build_trapezoid_edges(),
            ShapeType::Hexagon => self.build_hexagon_edges(),
            ShapeType::Saddle => self.build_saddle_edges(),
            ShapeType::Square => self.build_square_edges(),
        }
    }

    /// Build edges for Triangle shapes (8 configurations)
    fn build_triangle_edges(&mut self) {
        let pts = &self.points;
        match self.value {
            169 | 1 => {
                // 2221 | 0001
                if self.bottom_edge {
                    self.edges.insert(pts[0].clone(), Edge::new_with_move(pts[0].clone(), pts[1].clone(), Move::Left));
                }
                if self.left_edge {
                    self.edges.insert(pts[1].clone(), Edge::new(pts[1].clone(), pts[2].clone()));
                }
                self.edges.insert(pts[2].clone(), Edge::new_with_move(pts[2].clone(), pts[0].clone(), Move::Down));
            }
            166 | 4 => {
                // 2212 | 0010
                if self.right_edge {
                    self.edges.insert(pts[0].clone(), Edge::new_with_move(pts[0].clone(), pts[1].clone(), Move::Down));
                }
                if self.bottom_edge {
                    self.edges.insert(pts[1].clone(), Edge::new(pts[1].clone(), pts[2].clone()));
                }
                self.edges.insert(pts[2].clone(), Edge::new_with_move(pts[2].clone(), pts[0].clone(), Move::Right));
            }
            154 | 16 => {
                // 2122 | 0100
                if self.right_edge {
                    self.edges.insert(pts[0].clone(), Edge::new(pts[0].clone(), pts[1].clone()));
                }
                self.edges.insert(pts[1].clone(), Edge::new_with_move(pts[1].clone(), pts[2].clone(), Move::Up));
                if self.top_edge {
                    self.edges.insert(pts[2].clone(), Edge::new_with_move(pts[2].clone(), pts[0].clone(), Move::Right));
                }
            }
            106 | 64 => {
                // 1222 | 1000
                self.edges.insert(pts[0].clone(), Edge::new_with_move(pts[0].clone(), pts[1].clone(), Move::Left));
                if self.left_edge {
                    self.edges.insert(pts[1].clone(), Edge::new_with_move(pts[1].clone(), pts[2].clone(), Move::Up));
                }
                if self.top_edge {
                    self.edges.insert(pts[2].clone(), Edge::new(pts[2].clone(), pts[0].clone()));
                }
            }
            _ => {}
        }
    }

    /// Build edges for Pentagon shapes (24 configurations)
    fn build_pentagon_edges(&mut self) {
        let pts = &self.points;
        match self.value {
            101 | 69 => {
                self.edges.insert(pts[0].clone(), Edge::new_with_move(pts[0].clone(), pts[1].clone(), Move::Right));
                if self.right_edge {
                    self.edges.insert(pts[1].clone(), Edge::new_with_move(pts[1].clone(), pts[2].clone(), Move::Down));
                }
                if self.bottom_edge {
                    self.edges.insert(pts[2].clone(), Edge::new_with_move(pts[2].clone(), pts[3].clone(), Move::Left));
                }
                if self.left_edge {
                    self.edges.insert(pts[3].clone(), Edge::new_with_move(pts[3].clone(), pts[4].clone(), Move::Up));
                }
                if self.top_edge {
                    self.edges.insert(pts[4].clone(), Edge::new(pts[4].clone(), pts[0].clone()));
                }
            }
            149 | 21 => {
                if self.right_edge {
                    self.edges.insert(pts[0].clone(), Edge::new_with_move(pts[0].clone(), pts[1].clone(), Move::Down));
                }
                if self.bottom_edge {
                    self.edges.insert(pts[1].clone(), Edge::new_with_move(pts[1].clone(), pts[2].clone(), Move::Left));
                }
                if self.left_edge {
                    self.edges.insert(pts[2].clone(), Edge::new(pts[2].clone(), pts[3].clone()));
                }
                self.edges.insert(pts[3].clone(), Edge::new_with_move(pts[3].clone(), pts[4].clone(), Move::Up));
                if self.top_edge {
                    self.edges.insert(pts[4].clone(), Edge::new_with_move(pts[4].clone(), pts[0].clone(), Move::Right));
                }
            }
            86 | 84 => {
                if self.right_edge {
                    self.edges.insert(pts[0].clone(), Edge::new_with_move(pts[0].clone(), pts[1].clone(), Move::Down));
                }
                if self.bottom_edge {
                    self.edges.insert(pts[1].clone(), Edge::new(pts[1].clone(), pts[2].clone()));
                }
                self.edges.insert(pts[2].clone(), Edge::new_with_move(pts[2].clone(), pts[3].clone(), Move::Left));
                if self.left_edge {
                    self.edges.insert(pts[3].clone(), Edge::new_with_move(pts[3].clone(), pts[4].clone(), Move::Up));
                }
                if self.top_edge {
                    self.edges.insert(pts[4].clone(), Edge::new_with_move(pts[4].clone(), pts[0].clone(), Move::Right));
                }
            }
            89 | 81 => {
                if self.right_edge {
                    self.edges.insert(pts[0].clone(), Edge::new(pts[0].clone(), pts[1].clone()));
                }
                self.edges.insert(pts[1].clone(), Edge::new_with_move(pts[1].clone(), pts[2].clone(), Move::Down));
                if self.bottom_edge {
                    self.edges.insert(pts[2].clone(), Edge::new_with_move(pts[2].clone(), pts[3].clone(), Move::Left));
                }
                if self.left_edge {
                    self.edges.insert(pts[3].clone(), Edge::new_with_move(pts[3].clone(), pts[4].clone(), Move::Up));
                }
                if self.top_edge {
                    self.edges.insert(pts[4].clone(), Edge::new_with_move(pts[4].clone(), pts[0].clone(), Move::Right));
                }
            }
            96 | 74 => {
                self.edges.insert(pts[0].clone(), Edge::new_with_move(pts[0].clone(), pts[1].clone(), Move::Right));
                if self.right_edge {
                    self.edges.insert(pts[1].clone(), Edge::new(pts[1].clone(), pts[2].clone()));
                }
                self.edges.insert(pts[2].clone(), Edge::new_with_move(pts[2].clone(), pts[3].clone(), Move::Left));
                if self.left_edge {
                    self.edges.insert(pts[3].clone(), Edge::new_with_move(pts[3].clone(), pts[4].clone(), Move::Up));
                }
                if self.top_edge {
                    self.edges.insert(pts[4].clone(), Edge::new(pts[4].clone(), pts[0].clone()));
                }
            }
            24 | 146 => {
                if self.right_edge {
                    self.edges.insert(pts[0].clone(), Edge::new(pts[0].clone(), pts[1].clone()));
                }
                self.edges.insert(pts[1].clone(), Edge::new_with_move(pts[1].clone(), pts[2].clone(), Move::Down));
                if self.bottom_edge {
                    self.edges.insert(pts[2].clone(), Edge::new(pts[2].clone(), pts[3].clone()));
                }
                self.edges.insert(pts[3].clone(), Edge::new_with_move(pts[3].clone(), pts[4].clone(), Move::Up));
                if self.top_edge {
                    self.edges.insert(pts[4].clone(), Edge::new_with_move(pts[4].clone(), pts[0].clone(), Move::Right));
                }
            }
            6 | 164 => {
                if self.right_edge {
                    self.edges.insert(pts[0].clone(), Edge::new_with_move(pts[0].clone(), pts[1].clone(), Move::Down));
                }
                if self.bottom_edge {
                    self.edges.insert(pts[1].clone(), Edge::new(pts[1].clone(), pts[2].clone()));
                }
                self.edges.insert(pts[2].clone(), Edge::new_with_move(pts[2].clone(), pts[3].clone(), Move::Left));
                if self.left_edge {
                    self.edges.insert(pts[3].clone(), Edge::new(pts[3].clone(), pts[4].clone()));
                }
                self.edges.insert(pts[4].clone(), Edge::new_with_move(pts[4].clone(), pts[0].clone(), Move::Right));
            }
            129 | 41 => {
                self.edges.insert(pts[0].clone(), Edge::new_with_move(pts[0].clone(), pts[1].clone(), Move::Down));
                if self.bottom_edge {
                    self.edges.insert(pts[1].clone(), Edge::new_with_move(pts[1].clone(), pts[2].clone(), Move::Left));
                }
                if self.left_edge {
                    self.edges.insert(pts[2].clone(), Edge::new(pts[2].clone(), pts[3].clone()));
                }
                self.edges.insert(pts[3].clone(), Edge::new_with_move(pts[3].clone(), pts[4].clone(), Move::Up));
                if self.top_edge {
                    self.edges.insert(pts[4].clone(), Edge::new(pts[4].clone(), pts[0].clone()));
                }
            }
            66 | 104 => {
                self.edges.insert(pts[0].clone(), Edge::new_with_move(pts[0].clone(), pts[1].clone(), Move::Down));
                if self.bottom_edge {
                    self.edges.insert(pts[1].clone(), Edge::new(pts[1].clone(), pts[2].clone()));
                }
                self.edges.insert(pts[2].clone(), Edge::new_with_move(pts[2].clone(), pts[3].clone(), Move::Left));
                if self.left_edge {
                    self.edges.insert(pts[3].clone(), Edge::new_with_move(pts[3].clone(), pts[4].clone(), Move::Up));
                }
                if self.top_edge {
                    self.edges.insert(pts[4].clone(), Edge::new(pts[4].clone(), pts[0].clone()));
                }
            }
            144 | 26 => {
                if self.right_edge {
                    self.edges.insert(pts[0].clone(), Edge::new(pts[0].clone(), pts[1].clone()));
                }
                self.edges.insert(pts[1].clone(), Edge::new_with_move(pts[1].clone(), pts[2].clone(), Move::Left));
                if self.left_edge {
                    self.edges.insert(pts[2].clone(), Edge::new(pts[2].clone(), pts[3].clone()));
                }
                self.edges.insert(pts[3].clone(), Edge::new_with_move(pts[3].clone(), pts[4].clone(), Move::Up));
                if self.top_edge {
                    self.edges.insert(pts[4].clone(), Edge::new_with_move(pts[4].clone(), pts[0].clone(), Move::Right));
                }
            }
            36 | 134 => {
                self.edges.insert(pts[0].clone(), Edge::new_with_move(pts[0].clone(), pts[1].clone(), Move::Right));
                if self.right_edge {
                    self.edges.insert(pts[1].clone(), Edge::new_with_move(pts[1].clone(), pts[2].clone(), Move::Down));
                }
                if self.bottom_edge {
                    self.edges.insert(pts[2].clone(), Edge::new(pts[2].clone(), pts[3].clone()));
                }
                self.edges.insert(pts[3].clone(), Edge::new_with_move(pts[3].clone(), pts[4].clone(), Move::Up));
                if self.top_edge {
                    self.edges.insert(pts[4].clone(), Edge::new(pts[4].clone(), pts[0].clone()));
                }
            }
            9 | 161 => {
                if self.right_edge {
                    self.edges.insert(pts[0].clone(), Edge::new(pts[0].clone(), pts[1].clone()));
                }
                self.edges.insert(pts[1].clone(), Edge::new_with_move(pts[1].clone(), pts[2].clone(), Move::Down));
                if self.bottom_edge {
                    self.edges.insert(pts[2].clone(), Edge::new_with_move(pts[2].clone(), pts[3].clone(), Move::Left));
                }
                if self.left_edge {
                    self.edges.insert(pts[3].clone(), Edge::new(pts[3].clone(), pts[4].clone()));
                }
                self.edges.insert(pts[4].clone(), Edge::new_with_move(pts[4].clone(), pts[0].clone(), Move::Right));
            }
            _ => {}
        }
    }

    /// Build edges for Rectangle shapes (12 configurations)
    fn build_rectangle_edges(&mut self) {
        let pts = &self.points;
        match self.value {
            5 | 165 => {
                if self.right_edge {
                    self.edges.insert(pts[0].clone(), Edge::new_with_move(pts[0].clone(), pts[1].clone(), Move::Down));
                }
                if self.bottom_edge {
                    self.edges.insert(pts[1].clone(), Edge::new_with_move(pts[1].clone(), pts[2].clone(), Move::Left));
                }
                if self.left_edge {
                    self.edges.insert(pts[2].clone(), Edge::new(pts[2].clone(), pts[3].clone()));
                }
                self.edges.insert(pts[3].clone(), Edge::new_with_move(pts[3].clone(), pts[0].clone(), Move::Right));
            }
            20 | 150 => {
                if self.right_edge {
                    self.edges.insert(pts[0].clone(), Edge::new_with_move(pts[0].clone(), pts[1].clone(), Move::Down));
                }
                if self.bottom_edge {
                    self.edges.insert(pts[1].clone(), Edge::new(pts[1].clone(), pts[2].clone()));
                }
                self.edges.insert(pts[2].clone(), Edge::new_with_move(pts[2].clone(), pts[3].clone(), Move::Up));
                if self.top_edge {
                    self.edges.insert(pts[3].clone(), Edge::new_with_move(pts[3].clone(), pts[0].clone(), Move::Right));
                }
            }
            80 | 90 => {
                if self.right_edge {
                    self.edges.insert(pts[0].clone(), Edge::new(pts[0].clone(), pts[1].clone()));
                }
                self.edges.insert(pts[1].clone(), Edge::new_with_move(pts[1].clone(), pts[2].clone(), Move::Left));
                if self.left_edge {
                    self.edges.insert(pts[2].clone(), Edge::new_with_move(pts[2].clone(), pts[3].clone(), Move::Up));
                }
                if self.top_edge {
                    self.edges.insert(pts[3].clone(), Edge::new_with_move(pts[3].clone(), pts[0].clone(), Move::Right));
                }
            }
            65 | 105 => {
                self.edges.insert(pts[0].clone(), Edge::new_with_move(pts[0].clone(), pts[1].clone(), Move::Down));
                if self.bottom_edge {
                    self.edges.insert(pts[1].clone(), Edge::new_with_move(pts[1].clone(), pts[2].clone(), Move::Left));
                }
                if self.left_edge {
                    self.edges.insert(pts[2].clone(), Edge::new_with_move(pts[2].clone(), pts[3].clone(), Move::Up));
                }
                if self.top_edge {
                    self.edges.insert(pts[3].clone(), Edge::new(pts[3].clone(), pts[0].clone()));
                }
            }
            160 | 10 => {
                if self.right_edge {
                    self.edges.insert(pts[0].clone(), Edge::new(pts[0].clone(), pts[1].clone()));
                }
                self.edges.insert(pts[1].clone(), Edge::new_with_move(pts[1].clone(), pts[2].clone(), Move::Left));
                if self.left_edge {
                    self.edges.insert(pts[2].clone(), Edge::new(pts[2].clone(), pts[3].clone()));
                }
                self.edges.insert(pts[3].clone(), Edge::new_with_move(pts[3].clone(), pts[0].clone(), Move::Right));
            }
            130 | 40 => {
                self.edges.insert(pts[0].clone(), Edge::new_with_move(pts[0].clone(), pts[1].clone(), Move::Down));
                if self.bottom_edge {
                    self.edges.insert(pts[1].clone(), Edge::new(pts[1].clone(), pts[2].clone()));
                }
                self.edges.insert(pts[2].clone(), Edge::new_with_move(pts[2].clone(), pts[3].clone(), Move::Up));
                if self.top_edge {
                    self.edges.insert(pts[3].clone(), Edge::new(pts[3].clone(), pts[0].clone()));
                }
            }
            _ => {}
        }
    }

    /// Build edges for Trapezoid shapes (8 configurations)
    fn build_trapezoid_edges(&mut self) {
        let pts = &self.points;
        match self.value {
            168 | 2 => {
                if self.bottom_edge {
                    self.edges.insert(pts[0].clone(), Edge::new(pts[0].clone(), pts[1].clone()));
                }
                self.edges.insert(pts[1].clone(), Edge::new_with_move(pts[1].clone(), pts[2].clone(), Move::Left));
                if self.left_edge {
                    self.edges.insert(pts[2].clone(), Edge::new(pts[2].clone(), pts[3].clone()));
                }
                self.edges.insert(pts[3].clone(), Edge::new_with_move(pts[3].clone(), pts[0].clone(), Move::Down));
            }
            162 | 8 => {
                if self.right_edge {
                    self.edges.insert(pts[0].clone(), Edge::new(pts[0].clone(), pts[1].clone()));
                }
                self.edges.insert(pts[1].clone(), Edge::new_with_move(pts[1].clone(), pts[2].clone(), Move::Down));
                if self.bottom_edge {
                    self.edges.insert(pts[2].clone(), Edge::new(pts[2].clone(), pts[3].clone()));
                }
                self.edges.insert(pts[3].clone(), Edge::new_with_move(pts[3].clone(), pts[0].clone(), Move::Right));
            }
            138 | 32 => {
                self.edges.insert(pts[0].clone(), Edge::new_with_move(pts[0].clone(), pts[1].clone(), Move::Right));
                if self.right_edge {
                    self.edges.insert(pts[1].clone(), Edge::new(pts[1].clone(), pts[2].clone()));
                }
                self.edges.insert(pts[2].clone(), Edge::new_with_move(pts[2].clone(), pts[3].clone(), Move::Up));
                if self.top_edge {
                    self.edges.insert(pts[3].clone(), Edge::new(pts[3].clone(), pts[0].clone()));
                }
            }
            42 | 128 => {
                self.edges.insert(pts[0].clone(), Edge::new_with_move(pts[0].clone(), pts[1].clone(), Move::Left));
                if self.left_edge {
                    self.edges.insert(pts[1].clone(), Edge::new(pts[1].clone(), pts[2].clone()));
                }
                self.edges.insert(pts[2].clone(), Edge::new_with_move(pts[2].clone(), pts[3].clone(), Move::Up));
                if self.top_edge {
                    self.edges.insert(pts[3].clone(), Edge::new(pts[3].clone(), pts[0].clone()));
                }
            }
            _ => {}
        }
    }

    /// Build edges for Hexagon shapes (12 configurations)
    fn build_hexagon_edges(&mut self) {
        let pts = &self.points;
        match self.value {
            37 | 133 => {
                self.edges.insert(pts[0].clone(), Edge::new_with_move(pts[0].clone(), pts[1].clone(), Move::Right));
                if self.right_edge {
                    self.edges.insert(pts[1].clone(), Edge::new_with_move(pts[1].clone(), pts[2].clone(), Move::Down));
                }
                if self.bottom_edge {
                    self.edges.insert(pts[2].clone(), Edge::new_with_move(pts[2].clone(), pts[3].clone(), Move::Left));
                }
                if self.left_edge {
                    self.edges.insert(pts[3].clone(), Edge::new(pts[3].clone(), pts[4].clone()));
                }
                self.edges.insert(pts[4].clone(), Edge::new_with_move(pts[4].clone(), pts[5].clone(), Move::Up));
                if self.top_edge {
                    self.edges.insert(pts[5].clone(), Edge::new(pts[5].clone(), pts[0].clone()));
                }
            }
            148 | 22 => {
                if self.right_edge {
                    self.edges.insert(pts[0].clone(), Edge::new_with_move(pts[0].clone(), pts[1].clone(), Move::Down));
                }
                if self.bottom_edge {
                    self.edges.insert(pts[1].clone(), Edge::new(pts[1].clone(), pts[2].clone()));
                }
                self.edges.insert(pts[2].clone(), Edge::new_with_move(pts[2].clone(), pts[3].clone(), Move::Left));
                if self.left_edge {
                    self.edges.insert(pts[3].clone(), Edge::new(pts[3].clone(), pts[4].clone()));
                }
                self.edges.insert(pts[4].clone(), Edge::new_with_move(pts[4].clone(), pts[5].clone(), Move::Up));
                if self.top_edge {
                    self.edges.insert(pts[5].clone(), Edge::new_with_move(pts[5].clone(), pts[0].clone(), Move::Right));
                }
            }
            82 | 88 => {
                if self.right_edge {
                    self.edges.insert(pts[0].clone(), Edge::new(pts[0].clone(), pts[1].clone()));
                }
                self.edges.insert(pts[1].clone(), Edge::new_with_move(pts[1].clone(), pts[2].clone(), Move::Down));
                if self.bottom_edge {
                    self.edges.insert(pts[2].clone(), Edge::new(pts[2].clone(), pts[3].clone()));
                }
                self.edges.insert(pts[3].clone(), Edge::new_with_move(pts[3].clone(), pts[4].clone(), Move::Left));
                if self.left_edge {
                    self.edges.insert(pts[4].clone(), Edge::new_with_move(pts[4].clone(), pts[5].clone(), Move::Up));
                }
                if self.top_edge {
                    self.edges.insert(pts[5].clone(), Edge::new_with_move(pts[5].clone(), pts[0].clone(), Move::Right));
                }
            }
            73 | 97 => {
                self.edges.insert(pts[0].clone(), Edge::new_with_move(pts[0].clone(), pts[1].clone(), Move::Right));
                if self.right_edge {
                    self.edges.insert(pts[1].clone(), Edge::new(pts[1].clone(), pts[2].clone()));
                }
                self.edges.insert(pts[2].clone(), Edge::new_with_move(pts[2].clone(), pts[3].clone(), Move::Down));
                if self.bottom_edge {
                    self.edges.insert(pts[3].clone(), Edge::new_with_move(pts[3].clone(), pts[4].clone(), Move::Left));
                }
                if self.left_edge {
                    self.edges.insert(pts[4].clone(), Edge::new_with_move(pts[4].clone(), pts[5].clone(), Move::Up));
                }
                if self.top_edge {
                    self.edges.insert(pts[5].clone(), Edge::new(pts[5].clone(), pts[0].clone()));
                }
            }
            145 | 25 => {
                if self.right_edge {
                    self.edges.insert(pts[0].clone(), Edge::new(pts[0].clone(), pts[1].clone()));
                }
                self.edges.insert(pts[1].clone(), Edge::new_with_move(pts[1].clone(), pts[2].clone(), Move::Down));
                if self.bottom_edge {
                    self.edges.insert(pts[2].clone(), Edge::new_with_move(pts[2].clone(), pts[3].clone(), Move::Left));
                }
                if self.left_edge {
                    self.edges.insert(pts[3].clone(), Edge::new(pts[3].clone(), pts[4].clone()));
                }
                self.edges.insert(pts[4].clone(), Edge::new_with_move(pts[4].clone(), pts[5].clone(), Move::Up));
                if self.top_edge {
                    self.edges.insert(pts[5].clone(), Edge::new_with_move(pts[5].clone(), pts[0].clone(), Move::Right));
                }
            }
            70 | 100 => {
                self.edges.insert(pts[0].clone(), Edge::new_with_move(pts[0].clone(), pts[1].clone(), Move::Right));
                if self.right_edge {
                    self.edges.insert(pts[1].clone(), Edge::new_with_move(pts[1].clone(), pts[2].clone(), Move::Down));
                }
                if self.bottom_edge {
                    self.edges.insert(pts[2].clone(), Edge::new(pts[2].clone(), pts[3].clone()));
                }
                self.edges.insert(pts[3].clone(), Edge::new_with_move(pts[3].clone(), pts[4].clone(), Move::Left));
                if self.left_edge {
                    self.edges.insert(pts[4].clone(), Edge::new_with_move(pts[4].clone(), pts[5].clone(), Move::Up));
                }
                if self.top_edge {
                    self.edges.insert(pts[5].clone(), Edge::new(pts[5].clone(), pts[0].clone()));
                }
            }
            _ => {}
        }
    }

    /// Build edges for Saddle shapes (14 configurations)
    ///
    /// Saddle points require disambiguation using cell center average
    fn build_saddle_edges(&mut self) {
        let average = (self.tl + self.tr + self.br + self.bl) / 4.0;

        // Saddle shapes are complex and have different point configurations
        // We need to manually create points for saddles instead of using get_points()
        match self.value {
            153 => self.build_saddle_153(average),
            102 => self.build_saddle_102(average),
            68 => self.build_saddle_68(average),
            17 => self.build_saddle_17(average),
            136 => self.build_saddle_136(average),
            34 => self.build_saddle_34(average),
            152 => self.build_saddle_152(average),
            18 => self.build_saddle_18(average),
            137 => self.build_saddle_137(average),
            33 => self.build_saddle_33(average),
            98 => self.build_saddle_98(average),
            72 => self.build_saddle_72(average),
            38 => self.build_saddle_38(average),
            132 => self.build_saddle_132(average),
            _ => {}
        }
    }

    /// Build edges for Square shape (1 configuration)
    fn build_square_edges(&mut self) {
        let pts = &self.points;
        if self.right_edge {
            self.edges.insert(pts[0].clone(), Edge::new_with_move(pts[0].clone(), pts[1].clone(), Move::Down));
        }
        if self.bottom_edge {
            self.edges.insert(pts[1].clone(), Edge::new_with_move(pts[1].clone(), pts[2].clone(), Move::Left));
        }
        if self.left_edge {
            self.edges.insert(pts[2].clone(), Edge::new_with_move(pts[2].clone(), pts[3].clone(), Move::Up));
        }
        if self.top_edge {
            self.edges.insert(pts[3].clone(), Edge::new_with_move(pts[3].clone(), pts[0].clone(), Move::Right));
        }
    }

    // Saddle-specific implementations - these are very complex
    // For now, implementing the simplest case (153)
    // The rest will follow the exact same Java pattern

    fn build_saddle_153(&mut self, average: f64) {
        if average >= self.upper {
            let pt0 = self.interpolate(self.upper, Side::Right);
            let pt1 = self.interpolate(self.upper, Side::Top);
            let pt2 = self.top_right.clone();

            self.edges.insert(pt0.clone(), Edge::new_with_move(pt0.clone(), pt1.clone(), Move::Up));
            if self.top_edge {
                self.edges.insert(pt1.clone(), Edge::new_with_move(pt1.clone(), pt2.clone(), Move::Right));
            }
            if self.right_edge {
                self.edges.insert(pt2.clone(), Edge::new(pt2, pt0.clone()));
            }

            let pt3 = self.interpolate(self.upper, Side::Left);
            let pt4 = self.interpolate(self.upper, Side::Bottom);
            let pt5 = self.bottom_left.clone();

            self.edges.insert(pt3.clone(), Edge::new_with_move(pt3.clone(), pt4.clone(), Move::Down));
            if self.bottom_edge {
                self.edges.insert(pt4.clone(), Edge::new_with_move(pt4.clone(), pt5.clone(), Move::Left));
            }
            if self.left_edge {
                self.edges.insert(pt5.clone(), Edge::new(pt5, pt3));
            }
        } else if average >= self.lower {
            let pt0 = self.interpolate(self.upper, Side::Right);
            let pt1 = self.interpolate(self.upper, Side::Bottom);
            let pt2 = self.bottom_left.clone();
            let pt3 = self.interpolate(self.upper, Side::Left);
            let pt4 = self.interpolate(self.upper, Side::Top);
            let pt5 = self.top_right.clone();

            self.edges.insert(pt0.clone(), Edge::new_with_move(pt0.clone(), pt1.clone(), Move::Down));
            if self.bottom_edge {
                self.edges.insert(pt1.clone(), Edge::new_with_move(pt1.clone(), pt2.clone(), Move::Left));
            }
            if self.left_edge {
                self.edges.insert(pt2.clone(), Edge::new(pt2, pt3.clone()));
            }
            self.edges.insert(pt3.clone(), Edge::new_with_move(pt3, pt4.clone(), Move::Up));
            if self.top_edge {
                self.edges.insert(pt4.clone(), Edge::new_with_move(pt4, pt5.clone(), Move::Right));
            }
            if self.right_edge {
                self.edges.insert(pt5.clone(), Edge::new(pt5, pt0));
            }
        }
    }

    fn build_saddle_102(&mut self, average: f64) {
        // Case 1212
        if average >= self.upper {
            let pt0 = self.interpolate(self.upper, Side::Top);
            let pt1 = self.interpolate(self.upper, Side::Left);
            let pt2 = self.top_left.clone();

            self.edges.insert(pt0.clone(), Edge::new_with_move(pt0.clone(), pt1.clone(), Move::Left));
            if self.left_edge {
                self.edges.insert(pt1.clone(), Edge::new_with_move(pt1.clone(), pt2.clone(), Move::Up));
            }
            if self.top_edge {
                self.edges.insert(pt2.clone(), Edge::new(pt2, pt0));
            }

            let pt3 = self.interpolate(self.upper, Side::Bottom);
            let pt4 = self.interpolate(self.upper, Side::Right);
            let pt5 = self.bottom_right.clone();

            self.edges.insert(pt3.clone(), Edge::new_with_move(pt3.clone(), pt4.clone(), Move::Right));
            if self.right_edge {
                self.edges.insert(pt4.clone(), Edge::new_with_move(pt4.clone(), pt5.clone(), Move::Down));
            }
            if self.bottom_edge {
                self.edges.insert(pt5.clone(), Edge::new(pt5, pt3));
            }
        } else if average >= self.lower {
            let pt0 = self.interpolate(self.upper, Side::Top);
            let pt1 = self.interpolate(self.upper, Side::Right);
            let pt2 = self.bottom_right.clone();
            let pt3 = self.interpolate(self.upper, Side::Bottom);
            let pt4 = self.interpolate(self.upper, Side::Left);
            let pt5 = self.top_left.clone();

            self.edges.insert(pt0.clone(), Edge::new_with_move(pt0.clone(), pt1.clone(), Move::Right));
            if self.right_edge {
                self.edges.insert(pt1.clone(), Edge::new_with_move(pt1.clone(), pt2.clone(), Move::Down));
            }
            if self.bottom_edge {
                self.edges.insert(pt2.clone(), Edge::new(pt2, pt3.clone()));
            }
            self.edges.insert(pt3.clone(), Edge::new_with_move(pt3, pt4.clone(), Move::Left));
            if self.left_edge {
                self.edges.insert(pt4.clone(), Edge::new_with_move(pt4.clone(), pt5.clone(), Move::Up));
            }
            if self.top_edge {
                self.edges.insert(pt5.clone(), Edge::new(pt5, pt0));
            }
        }
    }

    fn build_saddle_68(&mut self, average: f64) {
        // Case 1010
        if average < self.lower {
            let pt0 = self.interpolate(self.lower, Side::Top);
            let pt1 = self.interpolate(self.lower, Side::Left);
            let pt2 = self.top_left.clone();
            let pt3 = self.interpolate(self.lower, Side::Bottom);
            let pt4 = self.interpolate(self.lower, Side::Right);
            let pt5 = self.bottom_right.clone();

            self.edges.insert(pt0.clone(), Edge::new_with_move(pt0.clone(), pt1.clone(), Move::Left));
            if self.left_edge {
                self.edges.insert(pt1.clone(), Edge::new_with_move(pt1.clone(), pt2.clone(), Move::Up));
            }
            if self.top_edge {
                self.edges.insert(pt2.clone(), Edge::new(pt2, pt0));
            }

            self.edges.insert(pt3.clone(), Edge::new_with_move(pt3.clone(), pt4.clone(), Move::Right));
            if self.right_edge {
                self.edges.insert(pt4.clone(), Edge::new_with_move(pt4.clone(), pt5.clone(), Move::Down));
            }
            if self.bottom_edge {
                self.edges.insert(pt5.clone(), Edge::new(pt5, pt3));
            }
        } else if average >= self.lower && average < self.upper {
            let pt0 = self.interpolate(self.lower, Side::Top);
            let pt1 = self.interpolate(self.lower, Side::Right);
            let pt2 = self.bottom_right.clone();
            let pt3 = self.interpolate(self.lower, Side::Bottom);
            let pt4 = self.interpolate(self.lower, Side::Left);
            let pt5 = self.top_left.clone();

            self.edges.insert(pt0.clone(), Edge::new_with_move(pt0.clone(), pt1.clone(), Move::Right));
            if self.right_edge {
                self.edges.insert(pt1.clone(), Edge::new_with_move(pt1.clone(), pt2.clone(), Move::Down));
            }
            if self.bottom_edge {
                self.edges.insert(pt2.clone(), Edge::new(pt2, pt3.clone()));
            }
            self.edges.insert(pt3.clone(), Edge::new_with_move(pt3, pt4.clone(), Move::Left));
            if self.left_edge {
                self.edges.insert(pt4.clone(), Edge::new_with_move(pt4.clone(), pt5.clone(), Move::Up));
            }
            if self.top_edge {
                self.edges.insert(pt5.clone(), Edge::new(pt5, pt0));
            }
        }
    }

    fn build_saddle_17(&mut self, average: f64) {
        // Case 0101
        if average < self.lower {
            let pt0 = self.interpolate(self.lower, Side::Right);
            let pt1 = self.interpolate(self.lower, Side::Top);
            let pt2 = self.top_right.clone();

            self.edges.insert(pt0.clone(), Edge::new_with_move(pt0.clone(), pt1.clone(), Move::Up));
            if self.top_edge {
                self.edges.insert(pt1.clone(), Edge::new_with_move(pt1.clone(), pt2.clone(), Move::Right));
            }
            if self.right_edge {
                self.edges.insert(pt2.clone(), Edge::new(pt2, pt0));
            }

            let pt3 = self.interpolate(self.lower, Side::Left);
            let pt4 = self.interpolate(self.lower, Side::Bottom);
            let pt5 = self.bottom_left.clone();

            self.edges.insert(pt3.clone(), Edge::new_with_move(pt3.clone(), pt4.clone(), Move::Down));
            if self.bottom_edge {
                self.edges.insert(pt4.clone(), Edge::new_with_move(pt4.clone(), pt5.clone(), Move::Left));
            }
            if self.left_edge {
                self.edges.insert(pt5.clone(), Edge::new(pt5, pt3));
            }
        } else if average >= self.lower && average < self.upper {
            let pt0 = self.interpolate(self.lower, Side::Right);
            let pt1 = self.interpolate(self.lower, Side::Bottom);
            let pt2 = self.bottom_left.clone();
            let pt3 = self.interpolate(self.lower, Side::Left);
            let pt4 = self.interpolate(self.lower, Side::Top);
            let pt5 = self.top_right.clone();

            self.edges.insert(pt0.clone(), Edge::new_with_move(pt0.clone(), pt1.clone(), Move::Down));
            if self.bottom_edge {
                self.edges.insert(pt1.clone(), Edge::new_with_move(pt1.clone(), pt2.clone(), Move::Left));
            }
            if self.left_edge {
                self.edges.insert(pt2.clone(), Edge::new(pt2, pt3.clone()));
            }

            self.edges.insert(pt3.clone(), Edge::new_with_move(pt3, pt4.clone(), Move::Up));
            if self.top_edge {
                self.edges.insert(pt4.clone(), Edge::new_with_move(pt4.clone(), pt5.clone(), Move::Right));
            }
            if self.right_edge {
                self.edges.insert(pt5.clone(), Edge::new(pt5, pt0));
            }
        }
    }

    fn build_saddle_136(&mut self, average: f64) {
        // Case 2020 - 8 points
        if average < self.lower {
            let pt0 = self.interpolate(self.lower, Side::Top);
            let pt1 = self.interpolate(self.lower, Side::Left);
            let pt2 = self.interpolate(self.upper, Side::Left);
            let pt3 = self.interpolate(self.upper, Side::Top);

            self.edges.insert(pt0.clone(), Edge::new_with_move(pt0.clone(), pt1.clone(), Move::Left));
            if self.left_edge {
                self.edges.insert(pt1.clone(), Edge::new(pt1, pt2.clone()));
            }
            self.edges.insert(pt2.clone(), Edge::new_with_move(pt2, pt3.clone(), Move::Up));
            if self.top_edge {
                self.edges.insert(pt3.clone(), Edge::new(pt3, pt0));
            }

            let pt4 = self.interpolate(self.upper, Side::Right);
            let pt5 = self.interpolate(self.upper, Side::Bottom);
            let pt6 = self.interpolate(self.lower, Side::Bottom);
            let pt7 = self.interpolate(self.lower, Side::Right);

            self.edges.insert(pt4.clone(), Edge::new_with_move(pt4.clone(), pt5.clone(), Move::Down));
            if self.bottom_edge {
                self.edges.insert(pt5.clone(), Edge::new(pt5, pt6.clone()));
            }
            self.edges.insert(pt6.clone(), Edge::new_with_move(pt6, pt7.clone(), Move::Right));
            if self.right_edge {
                self.edges.insert(pt7.clone(), Edge::new(pt7, pt4));
            }
        } else if average >= self.upper {
            let pt0 = self.interpolate(self.lower, Side::Top);
            let pt1 = self.interpolate(self.lower, Side::Right);
            let pt2 = self.interpolate(self.upper, Side::Right);
            let pt3 = self.interpolate(self.upper, Side::Top);

            self.edges.insert(pt0.clone(), Edge::new_with_move(pt0.clone(), pt1.clone(), Move::Right));
            if self.right_edge {
                self.edges.insert(pt1.clone(), Edge::new(pt1, pt2.clone()));
            }
            self.edges.insert(pt2.clone(), Edge::new_with_move(pt2, pt3.clone(), Move::Up));
            if self.top_edge {
                self.edges.insert(pt3.clone(), Edge::new(pt3, pt0));
            }

            let pt4 = self.interpolate(self.lower, Side::Bottom);
            let pt5 = self.interpolate(self.lower, Side::Left);
            let pt6 = self.interpolate(self.upper, Side::Left);
            let pt7 = self.interpolate(self.upper, Side::Bottom);

            self.edges.insert(pt4.clone(), Edge::new_with_move(pt4.clone(), pt5.clone(), Move::Left));
            if self.left_edge {
                self.edges.insert(pt5.clone(), Edge::new(pt5, pt6.clone()));
            }
            self.edges.insert(pt6.clone(), Edge::new_with_move(pt6, pt7.clone(), Move::Down));
            if self.bottom_edge {
                self.edges.insert(pt7.clone(), Edge::new(pt7, pt4));
            }
        } else {
            let pt0 = self.interpolate(self.lower, Side::Top);
            let pt1 = self.interpolate(self.lower, Side::Right);
            let pt2 = self.interpolate(self.upper, Side::Right);
            let pt3 = self.interpolate(self.upper, Side::Bottom);
            let pt4 = self.interpolate(self.lower, Side::Bottom);
            let pt5 = self.interpolate(self.lower, Side::Left);
            let pt6 = self.interpolate(self.upper, Side::Left);
            let pt7 = self.interpolate(self.upper, Side::Top);

            self.edges.insert(pt0.clone(), Edge::new_with_move(pt0.clone(), pt1.clone(), Move::Right));
            if self.right_edge {
                self.edges.insert(pt1.clone(), Edge::new(pt1, pt2.clone()));
            }
            self.edges.insert(pt2.clone(), Edge::new_with_move(pt2, pt3.clone(), Move::Down));
            if self.bottom_edge {
                self.edges.insert(pt3.clone(), Edge::new(pt3, pt4.clone()));
            }
            self.edges.insert(pt4.clone(), Edge::new_with_move(pt4, pt5.clone(), Move::Left));
            if self.left_edge {
                self.edges.insert(pt5.clone(), Edge::new(pt5, pt6.clone()));
            }
            self.edges.insert(pt6.clone(), Edge::new_with_move(pt6, pt7.clone(), Move::Up));
            if self.top_edge {
                self.edges.insert(pt7.clone(), Edge::new(pt7, pt0));
            }
        }
    }

    fn build_saddle_34(&mut self, average: f64) {
        // Case 0202 - 8 points
        if average >= self.upper {
            let pt0 = self.interpolate(self.upper, Side::Top);
            let pt1 = self.interpolate(self.upper, Side::Left);
            let pt2 = self.interpolate(self.lower, Side::Left);
            let pt3 = self.interpolate(self.lower, Side::Top);

            self.edges.insert(pt0.clone(), Edge::new_with_move(pt0.clone(), pt1.clone(), Move::Left));
            if self.left_edge {
                self.edges.insert(pt1.clone(), Edge::new(pt1, pt2.clone()));
            }
            self.edges.insert(pt2.clone(), Edge::new_with_move(pt2, pt3.clone(), Move::Up));
            if self.top_edge {
                self.edges.insert(pt3.clone(), Edge::new(pt3, pt0));
            }

            let pt4 = self.interpolate(self.lower, Side::Right);
            let pt5 = self.interpolate(self.lower, Side::Bottom);
            let pt6 = self.interpolate(self.upper, Side::Bottom);
            let pt7 = self.interpolate(self.upper, Side::Right);

            self.edges.insert(pt4.clone(), Edge::new_with_move(pt4.clone(), pt5.clone(), Move::Down));
            if self.bottom_edge {
                self.edges.insert(pt5.clone(), Edge::new(pt5, pt6.clone()));
            }
            self.edges.insert(pt6.clone(), Edge::new_with_move(pt6, pt7.clone(), Move::Right));
            if self.right_edge {
                self.edges.insert(pt7.clone(), Edge::new(pt7, pt4));
            }
        } else if average < self.lower {
            let pt0 = self.interpolate(self.upper, Side::Top);
            let pt1 = self.interpolate(self.upper, Side::Right);
            let pt2 = self.interpolate(self.lower, Side::Right);
            let pt3 = self.interpolate(self.lower, Side::Top);

            self.edges.insert(pt0.clone(), Edge::new_with_move(pt0.clone(), pt1.clone(), Move::Right));
            if self.right_edge {
                self.edges.insert(pt1.clone(), Edge::new(pt1, pt2.clone()));
            }
            self.edges.insert(pt2.clone(), Edge::new_with_move(pt2, pt3.clone(), Move::Up));
            if self.top_edge {
                self.edges.insert(pt3.clone(), Edge::new(pt3, pt0));
            }

            let pt4 = self.interpolate(self.upper, Side::Bottom);
            let pt5 = self.interpolate(self.upper, Side::Left);
            let pt6 = self.interpolate(self.lower, Side::Left);
            let pt7 = self.interpolate(self.lower, Side::Bottom);

            self.edges.insert(pt4.clone(), Edge::new_with_move(pt4.clone(), pt5.clone(), Move::Left));
            if self.left_edge {
                self.edges.insert(pt5.clone(), Edge::new(pt5, pt6.clone()));
            }
            self.edges.insert(pt6.clone(), Edge::new_with_move(pt6, pt7.clone(), Move::Down));
            if self.bottom_edge {
                self.edges.insert(pt7.clone(), Edge::new(pt7, pt4));
            }
        } else {
            let pt0 = self.interpolate(self.upper, Side::Top);
            let pt1 = self.interpolate(self.upper, Side::Right);
            let pt2 = self.interpolate(self.lower, Side::Right);
            let pt3 = self.interpolate(self.lower, Side::Bottom);
            let pt4 = self.interpolate(self.upper, Side::Bottom);
            let pt5 = self.interpolate(self.upper, Side::Left);
            let pt6 = self.interpolate(self.lower, Side::Left);
            let pt7 = self.interpolate(self.lower, Side::Top);

            self.edges.insert(pt0.clone(), Edge::new_with_move(pt0.clone(), pt1.clone(), Move::Right));
            if self.right_edge {
                self.edges.insert(pt1.clone(), Edge::new(pt1, pt2.clone()));
            }
            self.edges.insert(pt2.clone(), Edge::new_with_move(pt2, pt3.clone(), Move::Down));
            if self.bottom_edge {
                self.edges.insert(pt3.clone(), Edge::new(pt3, pt4.clone()));
            }
            self.edges.insert(pt4.clone(), Edge::new_with_move(pt4, pt5.clone(), Move::Left));
            if self.left_edge {
                self.edges.insert(pt5.clone(), Edge::new(pt5, pt6.clone()));
            }
            self.edges.insert(pt6.clone(), Edge::new_with_move(pt6, pt7.clone(), Move::Up));
            if self.top_edge {
                self.edges.insert(pt7.clone(), Edge::new(pt7, pt0));
            }
        }
    }

    fn build_saddle_152(&mut self, average: f64) {
        // Case 2120 - 7 points
        if average < self.lower || average >= self.upper {
            let pt0 = self.interpolate(self.upper, Side::Right);
            let pt1 = self.interpolate(self.upper, Side::Top);
            let pt2 = self.top_right.clone();

            self.edges.insert(pt0.clone(), Edge::new_with_move(pt0.clone(), pt1.clone(), Move::Up));
            if self.top_edge {
                self.edges.insert(pt1.clone(), Edge::new_with_move(pt1.clone(), pt2.clone(), Move::Right));
            }
            if self.right_edge {
                self.edges.insert(pt2.clone(), Edge::new(pt2, pt0));
            }

            let pt3 = self.interpolate(self.lower, Side::Bottom);
            let pt4 = self.interpolate(self.lower, Side::Left);
            let pt5 = self.interpolate(self.upper, Side::Left);
            let pt6 = self.interpolate(self.upper, Side::Bottom);

            self.edges.insert(pt3.clone(), Edge::new_with_move(pt3.clone(), pt4.clone(), Move::Left));
            if self.left_edge {
                self.edges.insert(pt4.clone(), Edge::new(pt4, pt5.clone()));
            }
            self.edges.insert(pt5.clone(), Edge::new_with_move(pt5, pt6.clone(), Move::Down));
            if self.bottom_edge {
                self.edges.insert(pt6.clone(), Edge::new(pt6, pt3));
            }
        } else {
            let pt0 = self.interpolate(self.upper, Side::Right);
            let pt1 = self.interpolate(self.upper, Side::Bottom);
            let pt2 = self.interpolate(self.lower, Side::Bottom);
            let pt3 = self.interpolate(self.lower, Side::Left);
            let pt4 = self.interpolate(self.upper, Side::Left);
            let pt5 = self.interpolate(self.upper, Side::Top);
            let pt6 = self.top_right.clone();

            self.edges.insert(pt0.clone(), Edge::new_with_move(pt0.clone(), pt1.clone(), Move::Down));
            if self.bottom_edge {
                self.edges.insert(pt1.clone(), Edge::new(pt1, pt2.clone()));
            }
            self.edges.insert(pt2.clone(), Edge::new_with_move(pt2, pt3.clone(), Move::Left));
            if self.left_edge {
                self.edges.insert(pt3.clone(), Edge::new(pt3, pt4.clone()));
            }
            self.edges.insert(pt4.clone(), Edge::new_with_move(pt4, pt5.clone(), Move::Up));
            if self.top_edge {
                self.edges.insert(pt5.clone(), Edge::new_with_move(pt5.clone(), pt6.clone(), Move::Right));
            }
            if self.right_edge {
                self.edges.insert(pt6.clone(), Edge::new(pt6, pt0));
            }
        }
    }

    fn build_saddle_18(&mut self, average: f64) {
        // Case 0102 - 7 points
        if average < self.lower || average >= self.upper {
            let pt0 = self.interpolate(self.lower, Side::Right);
            let pt1 = self.interpolate(self.lower, Side::Top);
            let pt2 = self.top_right.clone();

            self.edges.insert(pt0.clone(), Edge::new_with_move(pt0.clone(), pt1.clone(), Move::Up));
            if self.top_edge {
                self.edges.insert(pt1.clone(), Edge::new_with_move(pt1.clone(), pt2.clone(), Move::Right));
            }
            if self.right_edge {
                self.edges.insert(pt2.clone(), Edge::new(pt2, pt0));
            }

            let pt3 = self.interpolate(self.upper, Side::Bottom);
            let pt4 = self.interpolate(self.upper, Side::Left);
            let pt5 = self.interpolate(self.lower, Side::Left);
            let pt6 = self.interpolate(self.lower, Side::Bottom);

            self.edges.insert(pt3.clone(), Edge::new_with_move(pt3.clone(), pt4.clone(), Move::Left));
            if self.left_edge {
                self.edges.insert(pt4.clone(), Edge::new(pt4, pt5.clone()));
            }
            self.edges.insert(pt5.clone(), Edge::new_with_move(pt5, pt6.clone(), Move::Down));
            if self.bottom_edge {
                self.edges.insert(pt6.clone(), Edge::new(pt6, pt3));
            }
        } else {
            let pt0 = self.interpolate(self.lower, Side::Right);
            let pt1 = self.interpolate(self.lower, Side::Bottom);
            let pt2 = self.interpolate(self.upper, Side::Bottom);
            let pt3 = self.interpolate(self.upper, Side::Left);
            let pt4 = self.interpolate(self.lower, Side::Left);
            let pt5 = self.interpolate(self.lower, Side::Top);
            let pt6 = self.top_right.clone();

            self.edges.insert(pt0.clone(), Edge::new_with_move(pt0.clone(), pt1.clone(), Move::Down));
            if self.bottom_edge {
                self.edges.insert(pt1.clone(), Edge::new(pt1, pt2.clone()));
            }
            self.edges.insert(pt2.clone(), Edge::new_with_move(pt2, pt3.clone(), Move::Left));
            if self.left_edge {
                self.edges.insert(pt3.clone(), Edge::new(pt3, pt4.clone()));
            }
            self.edges.insert(pt4.clone(), Edge::new_with_move(pt4, pt5.clone(), Move::Up));
            if self.top_edge {
                self.edges.insert(pt5.clone(), Edge::new_with_move(pt5.clone(), pt6.clone(), Move::Right));
            }
            if self.right_edge {
                self.edges.insert(pt6.clone(), Edge::new(pt6, pt0));
            }
        }
    }

    fn build_saddle_137(&mut self, average: f64) {
        // Case 2021 - 7 points
        if average < self.lower || average >= self.upper {
            let pt0 = self.interpolate(self.lower, Side::Top);
            let pt1 = self.interpolate(self.lower, Side::Right);
            let pt2 = self.interpolate(self.upper, Side::Right);
            let pt3 = self.interpolate(self.upper, Side::Top);

            self.edges.insert(pt0.clone(), Edge::new_with_move(pt0.clone(), pt1.clone(), Move::Right));
            if self.right_edge {
                self.edges.insert(pt1.clone(), Edge::new(pt1, pt2.clone()));
            }
            self.edges.insert(pt2.clone(), Edge::new_with_move(pt2, pt3.clone(), Move::Up));
            if self.top_edge {
                self.edges.insert(pt3.clone(), Edge::new(pt3, pt0));
            }

            let pt4 = self.interpolate(self.upper, Side::Left);
            let pt5 = self.interpolate(self.upper, Side::Bottom);
            let pt6 = self.bottom_left.clone();

            self.edges.insert(pt4.clone(), Edge::new_with_move(pt4.clone(), pt5.clone(), Move::Down));
            if self.bottom_edge {
                self.edges.insert(pt5.clone(), Edge::new_with_move(pt5.clone(), pt6.clone(), Move::Left));
            }
            if self.left_edge {
                self.edges.insert(pt6.clone(), Edge::new(pt6, pt4));
            }
        } else {
            let pt0 = self.interpolate(self.lower, Side::Top);
            let pt1 = self.interpolate(self.lower, Side::Right);
            let pt2 = self.interpolate(self.upper, Side::Right);
            let pt3 = self.interpolate(self.upper, Side::Bottom);
            let pt4 = self.bottom_left.clone();
            let pt5 = self.interpolate(self.upper, Side::Left);
            let pt6 = self.interpolate(self.upper, Side::Top);

            self.edges.insert(pt0.clone(), Edge::new_with_move(pt0.clone(), pt1.clone(), Move::Right));
            if self.right_edge {
                self.edges.insert(pt1.clone(), Edge::new(pt1, pt2.clone()));
            }
            self.edges.insert(pt2.clone(), Edge::new_with_move(pt2, pt3.clone(), Move::Down));
            if self.bottom_edge {
                self.edges.insert(pt3.clone(), Edge::new_with_move(pt3.clone(), pt4.clone(), Move::Left));
            }
            if self.left_edge {
                self.edges.insert(pt4.clone(), Edge::new(pt4, pt5.clone()));
            }
            self.edges.insert(pt5.clone(), Edge::new_with_move(pt5, pt6.clone(), Move::Up));
            if self.top_edge {
                self.edges.insert(pt6.clone(), Edge::new(pt6, pt0));
            }
        }
    }

    fn build_saddle_33(&mut self, average: f64) {
        // Case 0201 - 7 points
        if average < self.lower || average >= self.upper {
            let pt0 = self.interpolate(self.upper, Side::Top);
            let pt1 = self.interpolate(self.upper, Side::Right);
            let pt2 = self.interpolate(self.lower, Side::Right);
            let pt3 = self.interpolate(self.lower, Side::Top);

            self.edges.insert(pt0.clone(), Edge::new_with_move(pt0.clone(), pt1.clone(), Move::Right));
            if self.right_edge {
                self.edges.insert(pt1.clone(), Edge::new(pt1, pt2.clone()));
            }
            self.edges.insert(pt2.clone(), Edge::new_with_move(pt2, pt3.clone(), Move::Up));
            if self.top_edge {
                self.edges.insert(pt3.clone(), Edge::new(pt3, pt0));
            }

            let pt4 = self.interpolate(self.lower, Side::Left);
            let pt5 = self.interpolate(self.lower, Side::Bottom);
            let pt6 = self.bottom_left.clone();

            self.edges.insert(pt4.clone(), Edge::new_with_move(pt4.clone(), pt5.clone(), Move::Down));
            if self.bottom_edge {
                self.edges.insert(pt5.clone(), Edge::new_with_move(pt5.clone(), pt6.clone(), Move::Left));
            }
            if self.left_edge {
                self.edges.insert(pt6.clone(), Edge::new(pt6, pt4));
            }
        } else {
            let pt0 = self.interpolate(self.upper, Side::Top);
            let pt1 = self.interpolate(self.upper, Side::Right);
            let pt2 = self.interpolate(self.lower, Side::Right);
            let pt3 = self.interpolate(self.lower, Side::Bottom);
            let pt4 = self.bottom_left.clone();
            let pt5 = self.interpolate(self.lower, Side::Left);
            let pt6 = self.interpolate(self.lower, Side::Top);

            self.edges.insert(pt0.clone(), Edge::new_with_move(pt0.clone(), pt1.clone(), Move::Right));
            if self.right_edge {
                self.edges.insert(pt1.clone(), Edge::new(pt1, pt2.clone()));
            }
            self.edges.insert(pt2.clone(), Edge::new_with_move(pt2, pt3.clone(), Move::Down));
            if self.bottom_edge {
                self.edges.insert(pt3.clone(), Edge::new_with_move(pt3.clone(), pt4.clone(), Move::Left));
            }
            if self.left_edge {
                self.edges.insert(pt4.clone(), Edge::new(pt4, pt5.clone()));
            }
            self.edges.insert(pt5.clone(), Edge::new_with_move(pt5, pt6.clone(), Move::Up));
            if self.top_edge {
                self.edges.insert(pt6.clone(), Edge::new(pt6, pt0));
            }
        }
    }

    fn build_saddle_98(&mut self, average: f64) {
        // Case 1202 - 7 points
        if average < self.lower || average >= self.upper {
            let pt0 = self.interpolate(self.upper, Side::Top);
            let pt1 = self.interpolate(self.upper, Side::Left);
            let pt2 = self.top_left.clone();

            self.edges.insert(pt0.clone(), Edge::new_with_move(pt0.clone(), pt1.clone(), Move::Left));
            if self.left_edge {
                self.edges.insert(pt1.clone(), Edge::new_with_move(pt1.clone(), pt2.clone(), Move::Up));
            }
            if self.top_edge {
                self.edges.insert(pt2.clone(), Edge::new(pt2, pt0));
            }

            let pt3 = self.interpolate(self.lower, Side::Right);
            let pt4 = self.interpolate(self.lower, Side::Bottom);
            let pt5 = self.interpolate(self.upper, Side::Bottom);
            let pt6 = self.interpolate(self.upper, Side::Right);

            self.edges.insert(pt3.clone(), Edge::new_with_move(pt3.clone(), pt4.clone(), Move::Down));
            if self.bottom_edge {
                self.edges.insert(pt4.clone(), Edge::new(pt4, pt5.clone()));
            }
            self.edges.insert(pt5.clone(), Edge::new_with_move(pt5, pt6.clone(), Move::Right));
            if self.right_edge {
                self.edges.insert(pt6.clone(), Edge::new(pt6, pt3));
            }
        } else {
            let pt0 = self.interpolate(self.upper, Side::Top);
            let pt1 = self.interpolate(self.upper, Side::Right);
            let pt2 = self.interpolate(self.lower, Side::Right);
            let pt3 = self.interpolate(self.lower, Side::Bottom);
            let pt4 = self.interpolate(self.upper, Side::Bottom);
            let pt5 = self.interpolate(self.upper, Side::Left);
            let pt6 = self.top_left.clone();

            self.edges.insert(pt0.clone(), Edge::new_with_move(pt0.clone(), pt1.clone(), Move::Right));
            if self.right_edge {
                self.edges.insert(pt1.clone(), Edge::new(pt1, pt2.clone()));
            }
            self.edges.insert(pt2.clone(), Edge::new_with_move(pt2, pt3.clone(), Move::Down));
            if self.bottom_edge {
                self.edges.insert(pt3.clone(), Edge::new(pt3, pt4.clone()));
            }
            self.edges.insert(pt4.clone(), Edge::new_with_move(pt4, pt5.clone(), Move::Left));
            if self.left_edge {
                self.edges.insert(pt5.clone(), Edge::new_with_move(pt5.clone(), pt6.clone(), Move::Up));
            }
            if self.top_edge {
                self.edges.insert(pt6.clone(), Edge::new(pt6, pt0));
            }
        }
    }

    fn build_saddle_72(&mut self, average: f64) {
        // Case 1020 - 7 points
        if average < self.lower || average >= self.upper {
            let pt0 = self.interpolate(self.lower, Side::Top);
            let pt1 = self.interpolate(self.lower, Side::Left);
            let pt2 = self.top_left.clone();

            self.edges.insert(pt0.clone(), Edge::new_with_move(pt0.clone(), pt1.clone(), Move::Left));
            if self.left_edge {
                self.edges.insert(pt1.clone(), Edge::new_with_move(pt1.clone(), pt2.clone(), Move::Up));
            }
            if self.top_edge {
                self.edges.insert(pt2.clone(), Edge::new(pt2, pt0));
            }

            let pt3 = self.interpolate(self.upper, Side::Right);
            let pt4 = self.interpolate(self.upper, Side::Bottom);
            let pt5 = self.interpolate(self.lower, Side::Bottom);
            let pt6 = self.interpolate(self.lower, Side::Right);

            self.edges.insert(pt3.clone(), Edge::new_with_move(pt3.clone(), pt4.clone(), Move::Down));
            if self.bottom_edge {
                self.edges.insert(pt4.clone(), Edge::new(pt4, pt5.clone()));
            }
            self.edges.insert(pt5.clone(), Edge::new_with_move(pt5, pt6.clone(), Move::Right));
            if self.right_edge {
                self.edges.insert(pt6.clone(), Edge::new(pt6, pt3));
            }
        } else {
            let pt0 = self.interpolate(self.lower, Side::Top);
            let pt1 = self.interpolate(self.lower, Side::Right);
            let pt2 = self.interpolate(self.upper, Side::Right);
            let pt3 = self.interpolate(self.upper, Side::Bottom);
            let pt4 = self.interpolate(self.lower, Side::Bottom);
            let pt5 = self.interpolate(self.lower, Side::Left);
            let pt6 = self.top_left.clone();

            self.edges.insert(pt0.clone(), Edge::new_with_move(pt0.clone(), pt1.clone(), Move::Right));
            if self.right_edge {
                self.edges.insert(pt1.clone(), Edge::new(pt1, pt2.clone()));
            }
            self.edges.insert(pt2.clone(), Edge::new_with_move(pt2, pt3.clone(), Move::Down));
            if self.bottom_edge {
                self.edges.insert(pt3.clone(), Edge::new(pt3, pt4.clone()));
            }
            self.edges.insert(pt4.clone(), Edge::new_with_move(pt4, pt5.clone(), Move::Left));
            if self.left_edge {
                self.edges.insert(pt5.clone(), Edge::new_with_move(pt5.clone(), pt6.clone(), Move::Up));
            }
            if self.top_edge {
                self.edges.insert(pt6.clone(), Edge::new(pt6, pt0));
            }
        }
    }

    fn build_saddle_38(&mut self, average: f64) {
        // Case 0212 - 7 points
        if average < self.lower || average >= self.upper {
            let pt0 = self.interpolate(self.upper, Side::Top);
            let pt1 = self.interpolate(self.upper, Side::Left);
            let pt2 = self.interpolate(self.lower, Side::Left);
            let pt3 = self.interpolate(self.lower, Side::Top);

            self.edges.insert(pt0.clone(), Edge::new_with_move(pt0.clone(), pt1.clone(), Move::Left));
            if self.left_edge {
                self.edges.insert(pt1.clone(), Edge::new(pt1, pt2.clone()));
            }
            self.edges.insert(pt2.clone(), Edge::new_with_move(pt2, pt3.clone(), Move::Up));
            if self.top_edge {
                self.edges.insert(pt3.clone(), Edge::new(pt3, pt0));
            }

            let pt4 = self.interpolate(self.upper, Side::Bottom);
            let pt5 = self.interpolate(self.upper, Side::Right);
            let pt6 = self.bottom_right.clone();

            self.edges.insert(pt4.clone(), Edge::new_with_move(pt4.clone(), pt5.clone(), Move::Right));
            if self.right_edge {
                self.edges.insert(pt5.clone(), Edge::new_with_move(pt5.clone(), pt6.clone(), Move::Down));
            }
            if self.bottom_edge {
                self.edges.insert(pt6.clone(), Edge::new(pt6, pt4));
            }
        } else {
            let pt0 = self.interpolate(self.upper, Side::Top);
            let pt1 = self.interpolate(self.upper, Side::Right);
            let pt2 = self.bottom_right.clone();
            let pt3 = self.interpolate(self.upper, Side::Bottom);
            let pt4 = self.interpolate(self.upper, Side::Left);
            let pt5 = self.interpolate(self.lower, Side::Left);
            let pt6 = self.interpolate(self.lower, Side::Top);

            self.edges.insert(pt0.clone(), Edge::new_with_move(pt0.clone(), pt1.clone(), Move::Right));
            if self.right_edge {
                self.edges.insert(pt1.clone(), Edge::new_with_move(pt1.clone(), pt2.clone(), Move::Down));
            }
            if self.bottom_edge {
                self.edges.insert(pt2.clone(), Edge::new(pt2, pt3.clone()));
            }
            self.edges.insert(pt3.clone(), Edge::new_with_move(pt3, pt4.clone(), Move::Left));
            if self.left_edge {
                self.edges.insert(pt4.clone(), Edge::new(pt4, pt5.clone()));
            }
            self.edges.insert(pt5.clone(), Edge::new_with_move(pt5, pt6.clone(), Move::Up));
            if self.top_edge {
                self.edges.insert(pt6.clone(), Edge::new(pt6, pt0));
            }
        }
    }

    fn build_saddle_132(&mut self, average: f64) {
        // Case 2010 - 7 points
        if average < self.lower || average >= self.upper {
            let pt0 = self.interpolate(self.lower, Side::Top);
            let pt1 = self.interpolate(self.lower, Side::Left);
            let pt2 = self.interpolate(self.upper, Side::Left);
            let pt3 = self.interpolate(self.upper, Side::Top);

            self.edges.insert(pt0.clone(), Edge::new_with_move(pt0.clone(), pt1.clone(), Move::Left));
            if self.left_edge {
                self.edges.insert(pt1.clone(), Edge::new(pt1, pt2.clone()));
            }
            self.edges.insert(pt2.clone(), Edge::new_with_move(pt2, pt3.clone(), Move::Up));
            if self.top_edge {
                self.edges.insert(pt3.clone(), Edge::new(pt3, pt0));
            }

            let pt4 = self.interpolate(self.lower, Side::Bottom);
            let pt5 = self.interpolate(self.lower, Side::Right);
            let pt6 = self.bottom_right.clone();

            self.edges.insert(pt4.clone(), Edge::new_with_move(pt4.clone(), pt5.clone(), Move::Right));
            if self.right_edge {
                self.edges.insert(pt5.clone(), Edge::new_with_move(pt5.clone(), pt6.clone(), Move::Down));
            }
            if self.bottom_edge {
                self.edges.insert(pt6.clone(), Edge::new(pt6, pt4));
            }
        } else {
            let pt0 = self.interpolate(self.lower, Side::Top);
            let pt1 = self.interpolate(self.lower, Side::Right);
            let pt2 = self.bottom_right.clone();
            let pt3 = self.interpolate(self.lower, Side::Bottom);
            let pt4 = self.interpolate(self.lower, Side::Left);
            let pt5 = self.interpolate(self.upper, Side::Left);
            let pt6 = self.interpolate(self.upper, Side::Top);

            self.edges.insert(pt0.clone(), Edge::new_with_move(pt0.clone(), pt1.clone(), Move::Right));
            if self.right_edge {
                self.edges.insert(pt1.clone(), Edge::new_with_move(pt1.clone(), pt2.clone(), Move::Down));
            }
            if self.bottom_edge {
                self.edges.insert(pt2.clone(), Edge::new(pt2, pt3.clone()));
            }
            self.edges.insert(pt3.clone(), Edge::new_with_move(pt3, pt4.clone(), Move::Left));
            if self.left_edge {
                self.edges.insert(pt4.clone(), Edge::new(pt4, pt5.clone()));
            }
            self.edges.insert(pt5.clone(), Edge::new_with_move(pt5, pt6.clone(), Move::Up));
            if self.top_edge {
                self.edges.insert(pt6.clone(), Edge::new(pt6, pt0));
            }
        }
    }

    /// Check if top side is blank (no contour crossing)
    fn is_top_blank(&self) -> bool {
        (self.tl >= self.upper && self.tr >= self.upper) || (self.tl < self.lower && self.tr < self.lower)
    }

    /// Check if right side is blank (no contour crossing)
    fn is_right_blank(&self) -> bool {
        (self.tr >= self.upper && self.br >= self.upper) || (self.tr < self.lower && self.br < self.lower)
    }

    /// Check if bottom side is blank (no contour crossing)
    fn is_bottom_blank(&self) -> bool {
        (self.bl >= self.upper && self.br >= self.upper) || (self.bl < self.lower && self.br < self.lower)
    }

    /// Check if left side is blank (no contour crossing)
    fn is_left_blank(&self) -> bool {
        (self.tl >= self.upper && self.bl >= self.upper) || (self.tl < self.lower && self.bl < self.lower)
    }

    /// Generate the 8 potential edge points for this cell
    ///
    /// This matches Java's getPoints() logic exactly (Shape.java lines 226-246)
    fn get_points(&self) -> Vec<Point> {
        let mut eight_points = Vec::new();

        // Top-right corner, top side
        eight_points.push(if self.is_top_blank() {
            None
        } else if self.tr >= self.upper {
            Some(Point::new_with_limit(self.tr, self.upper, Side::Top))
        } else if self.tr < self.lower {
            Some(Point::new_with_limit(self.tr, self.lower, Side::Top))
        } else {
            Some(self.top_right.clone())
        });

        // Top-right corner, right side
        eight_points.push(if self.is_right_blank() {
            None
        } else if self.tr >= self.upper {
            Some(Point::new_with_limit(self.tr, self.upper, Side::Right))
        } else if self.tr < self.lower {
            Some(Point::new_with_limit(self.tr, self.lower, Side::Right))
        } else {
            Some(self.top_right.clone())
        });

        // Bottom-right corner, right side
        eight_points.push(if self.is_right_blank() {
            None
        } else if self.br >= self.upper {
            Some(Point::new_with_limit(self.br, self.upper, Side::Right))
        } else if self.br < self.lower {
            Some(Point::new_with_limit(self.br, self.lower, Side::Right))
        } else {
            Some(self.bottom_right.clone())
        });

        // Bottom-right corner, bottom side
        eight_points.push(if self.is_bottom_blank() {
            None
        } else if self.br >= self.upper {
            Some(Point::new_with_limit(self.br, self.upper, Side::Bottom))
        } else if self.br < self.lower {
            Some(Point::new_with_limit(self.br, self.lower, Side::Bottom))
        } else {
            Some(self.bottom_right.clone())
        });

        // Bottom-left corner, bottom side
        eight_points.push(if self.is_bottom_blank() {
            None
        } else if self.bl >= self.upper {
            Some(Point::new_with_limit(self.bl, self.upper, Side::Bottom))
        } else if self.bl < self.lower {
            Some(Point::new_with_limit(self.bl, self.lower, Side::Bottom))
        } else {
            Some(self.bottom_left.clone())
        });

        // Bottom-left corner, left side
        eight_points.push(if self.is_left_blank() {
            None
        } else if self.bl >= self.upper {
            Some(Point::new_with_limit(self.bl, self.upper, Side::Left))
        } else if self.bl < self.lower {
            Some(Point::new_with_limit(self.bl, self.lower, Side::Left))
        } else {
            Some(self.bottom_left.clone())
        });

        // Top-left corner, left side
        eight_points.push(if self.is_left_blank() {
            None
        } else if self.tl >= self.upper {
            Some(Point::new_with_limit(self.tl, self.upper, Side::Left))
        } else if self.tl < self.lower {
            Some(Point::new_with_limit(self.tl, self.lower, Side::Left))
        } else {
            Some(self.top_left.clone())
        });

        // Top-left corner, top side
        eight_points.push(if self.is_top_blank() {
            None
        } else if self.tl >= self.upper {
            Some(Point::new_with_limit(self.tl, self.upper, Side::Top))
        } else if self.tl < self.lower {
            Some(Point::new_with_limit(self.tl, self.lower, Side::Top))
        } else {
            Some(self.top_left.clone())
        });

        // Filter out nulls and deduplicate
        let mut slim: Vec<Point> = eight_points.into_iter().flatten().collect();

        // Remove duplicates (Java uses distinct())
        let mut seen = HashMap::new();
        slim.retain(|p| seen.insert(p.clone(), ()).is_none());

        // Interpolate points that need it
        for point in &mut slim {
            if point.x().is_none() && point.y().is_none() {
                if let (Some(limit), Some(side)) = (point.limit(), point.side()) {
                    *point = self.interpolate(limit, side);
                }
            }
        }

        slim
    }

    /// Interpolate a point on a cell side
    ///
    /// Matches Java's interpolate(double level, Side side) method (Shape.java lines 513-525)
    fn interpolate(&self, level: f64, side: Side) -> Point {
        match side {
            Side::Top => self.interpolate_between(level, self.tl, self.tr, &self.top_left, &self.top_right),
            Side::Right => self.interpolate_between(level, self.tr, self.br, &self.top_right, &self.bottom_right),
            Side::Bottom => self.interpolate_between(level, self.bl, self.br, &self.bottom_left, &self.bottom_right),
            Side::Left => self.interpolate_between(level, self.tl, self.bl, &self.top_left, &self.bottom_left),
        }
    }

    /// Interpolate between two points using cosine smoothing
    ///
    /// This is the EXACT formula from Java (Shape.java lines 492-511)
    /// including the 0.999 "hack" for centering
    fn interpolate_between(&self, level: f64, value0: f64, value1: f64, point0: &Point, point1: &Point) -> Point {
        let mu = (level - value0) / (value1 - value0);
        let mu2 = (1.0 - (mu * std::f64::consts::PI).cos()) / 2.0;

        let center_diff = (mu2 - 0.5) * 0.999;
        let new_mu = 0.5 + center_diff;

        let x = ((1.0 - new_mu) * point0.x().unwrap()) + (new_mu * point1.x().unwrap());
        let y = ((1.0 - new_mu) * point0.y().unwrap()) + (new_mu * point1.y().unwrap());

        Point::new(x, y)
    }

    /// Create a new shape (low-level constructor for testing)
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

    // Getters
    pub fn shape_type(&self) -> ShapeType {
        self.shape_type
    }

    pub fn value(&self) -> u8 {
        self.value
    }

    pub fn points(&self) -> &[Point] {
        &self.points
    }

    pub fn edges(&self) -> &HashMap<Point, Edge> {
        &self.edges
    }

    pub fn get_edges(&self, start: Option<&Point>) -> Vec<Edge> {
        if self.edges.len() <= 1 {
            return self.edges.values().cloned().collect();
        }

        let mut result = Vec::new();
        let mut current = match start {
            Some(s) => s.clone(),
            None => {
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

    pub fn is_cleared(&self) -> bool {
        self.cleared
    }

    pub fn x(&self) -> usize {
        self.x
    }

    pub fn y(&self) -> usize {
        self.y
    }

    pub fn is_top_edge(&self) -> bool {
        self.top_edge
    }

    pub fn is_right_edge(&self) -> bool {
        self.right_edge
    }

    pub fn is_bottom_edge(&self) -> bool {
        self.bottom_edge
    }

    pub fn is_left_edge(&self) -> bool {
        self.left_edge
    }

    pub fn increment_used_edges(&mut self, count: usize) {
        self.used_edges += count;
        if self.used_edges >= self.edges.len() {
            self.cleared = true;
        }
    }

    pub fn remove_edge(&mut self, key: &Point) {
        self.edges.remove(key);
    }

    pub fn set_points(&mut self, points: Vec<Point>) {
        self.points = points;
    }

    pub fn insert_edge(&mut self, start: Point, edge: Edge) {
        self.edges.insert(start, edge);
    }

    pub fn corner_values(&self) -> (f64, f64, f64, f64) {
        (self.tl, self.tr, self.bl, self.br)
    }

    pub fn thresholds(&self) -> (f64, f64) {
        (self.lower, self.upper)
    }

    pub fn corner_points(&self) -> (&Point, &Point, &Point, &Point) {
        (&self.top_left, &self.top_right, &self.bottom_right, &self.bottom_left)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ternary_classification() {
        // Test value 1 (0001): bl=1, others=0
        let tl = 5.0;
        let tr = 5.0;
        let br = 5.0;
        let bl = 15.0;
        let lower = 10.0;
        let upper = 20.0;

        let mut value = 0u8;
        value |= if tl < lower { 0 } else if tl >= upper { 128 } else { 64 };
        value |= if tr < lower { 0 } else if tr >= upper { 32 } else { 16 };
        value |= if br < lower { 0 } else if br >= upper { 8 } else { 4 };
        value |= if bl < lower { 0 } else if bl >= upper { 2 } else { 1 };

        assert_eq!(value, 1);
    }

    #[test]
    fn test_ternary_classification_85() {
        // Test value 85 (1111): all corners in middle band
        let tl = 15.0;
        let tr = 15.0;
        let br = 15.0;
        let bl = 15.0;
        let lower = 10.0;
        let upper = 20.0;

        let mut value = 0u8;
        value |= if tl < lower { 0 } else if tl >= upper { 128 } else { 64 };
        value |= if tr < lower { 0 } else if tr >= upper { 32 } else { 16 };
        value |= if br < lower { 0 } else if br >= upper { 8 } else { 4 };
        value |= if bl < lower { 0 } else if bl >= upper { 2 } else { 1 };

        assert_eq!(value, 85); // 64 + 16 + 4 + 1
    }

    #[test]
    fn test_interpolation() {
        let shape = Shape::new(
            ShapeType::Triangle,
            Point::new(0.0, 2.0),
            Point::new(1.0, 2.0),
            Point::new(1.0, 1.0),
            Point::new(0.0, 1.0),
            1,
            10.0,
            20.0,
            0,
            0,
            false, false, false, false,
            5.0, 15.0, 15.0, 5.0,  // tl, tr, br, bl - fixed: br should be 15, not 5
        );

        // Interpolate on the bottom side between bl=5 and br=15 at level=10
        let point = shape.interpolate(10.0, Side::Bottom);

        // Should be between (0.0, 1.0) and (1.0, 1.0)
        assert!(point.x().unwrap() > 0.0 && point.x().unwrap() < 1.0);
        assert_eq!(point.y().unwrap(), 1.0); // Y is constant on bottom edge
    }
}
