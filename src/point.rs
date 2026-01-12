use std::hash::{Hash, Hasher};

/// Side of a cell where an interpolated point lies
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Side {
    Top,
    Right,
    Bottom,
    Left,
}

/// A point in 2D space, either with actual coordinates or as a placeholder for interpolation
///
/// Points can represent two different concepts:
/// 1. Actual coordinates (x, y populated)
/// 2. Markers for interpolation (value, limit, side populated)
#[derive(Debug, Clone)]
pub struct Point {
    x: Option<f64>,
    y: Option<f64>,
    value: Option<f64>,
    limit: Option<f64>,
    side: Option<Side>,
}

impl Point {
    /// Create a point with actual coordinates
    pub fn new(x: f64, y: f64) -> Self {
        Self {
            x: Some(x),
            y: Some(y),
            value: None,
            limit: None,
            side: None,
        }
    }

    /// Create a point as an interpolation marker
    ///
    /// This point needs interpolation to determine actual coordinates.
    /// Used when a cell corner value equals a threshold boundary.
    pub fn new_with_limit(value: f64, limit: f64, side: Side) -> Self {
        Self {
            x: None,
            y: None,
            value: Some(value),
            limit: Some(limit),
            side: Some(side),
        }
    }

    /// Get x coordinate
    pub fn x(&self) -> Option<f64> {
        self.x
    }

    /// Get y coordinate
    pub fn y(&self) -> Option<f64> {
        self.y
    }

    /// Get the value used for interpolation
    pub fn value(&self) -> Option<f64> {
        self.value
    }

    /// Get the threshold limit
    pub fn limit(&self) -> Option<f64> {
        self.limit
    }

    /// Get the side of the cell
    pub fn side(&self) -> Option<Side> {
        self.side
    }
}

impl PartialEq for Point {
    fn eq(&self, other: &Self) -> bool {
        // Compare floats with option wrapping
        let x_eq = match (self.x, other.x) {
            (Some(a), Some(b)) => a == b,
            (None, None) => true,
            _ => false,
        };
        let y_eq = match (self.y, other.y) {
            (Some(a), Some(b)) => a == b,
            (None, None) => true,
            _ => false,
        };
        let value_eq = match (self.value, other.value) {
            (Some(a), Some(b)) => a == b,
            (None, None) => true,
            _ => false,
        };
        let limit_eq = match (self.limit, other.limit) {
            (Some(a), Some(b)) => a == b,
            (None, None) => true,
            _ => false,
        };

        x_eq && y_eq && value_eq && limit_eq && self.side == other.side
    }
}

impl Eq for Point {}

impl Hash for Point {
    fn hash<H: Hasher>(&self, state: &mut H) {
        // Hash option-wrapped f64s by converting to bits
        if let Some(x) = self.x {
            x.to_bits().hash(state);
        } else {
            0u64.hash(state);
        }
        if let Some(y) = self.y {
            y.to_bits().hash(state);
        } else {
            0u64.hash(state);
        }
        if let Some(value) = self.value {
            value.to_bits().hash(state);
        } else {
            0u64.hash(state);
        }
        if let Some(limit) = self.limit {
            limit.to_bits().hash(state);
        } else {
            0u64.hash(state);
        }
        self.side.hash(state);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    #[test]
    fn test_point_with_coordinates() {
        let p = Point::new(10.5, 20.3);
        assert_eq!(p.x(), Some(10.5));
        assert_eq!(p.y(), Some(20.3));
        assert_eq!(p.value(), None);
        assert_eq!(p.limit(), None);
        assert_eq!(p.side(), None);
    }

    #[test]
    fn test_point_with_limit() {
        let p = Point::new_with_limit(15.0, 10.0, Side::Top);
        assert_eq!(p.x(), None);
        assert_eq!(p.y(), None);
        assert_eq!(p.value(), Some(15.0));
        assert_eq!(p.limit(), Some(10.0));
        assert_eq!(p.side(), Some(Side::Top));
    }

    #[test]
    fn test_point_equality() {
        let p1 = Point::new(10.0, 20.0);
        let p2 = Point::new(10.0, 20.0);
        let p3 = Point::new(10.0, 21.0);

        assert_eq!(p1, p2);
        assert_ne!(p1, p3);
    }

    #[test]
    fn test_point_in_hashmap() {
        let mut map = HashMap::new();
        let p1 = Point::new(5.0, 10.0);
        let p2 = Point::new(5.0, 10.0);

        map.insert(p1.clone(), "value1");
        map.insert(p2.clone(), "value2");

        // Should overwrite since p1 == p2
        assert_eq!(map.len(), 1);
        assert_eq!(map.get(&p1), Some(&"value2"));
    }
}
