//! Lightweight grid cell representation for efficient memory usage
//!
//! This module provides a minimal representation of grid points that is
//! much more memory-efficient than full GeoJSON Features.

/// A lightweight grid cell containing geographic coordinates and a scalar value
///
/// This struct is designed for efficient memory usage with large grids.
/// At 24 bytes per cell (3 × f64), it uses ~12x less memory than a full
/// GeoJSON Feature with all its metadata.
///
/// # Memory Comparison
///
/// - `GridCell`: 24 bytes (3 × f64)
/// - `geojson::Feature`: ~200-300 bytes (with geometry, properties, options)
///
/// For a 1799×1059 HRRR grid (1.9M points):
/// - GridCell array: ~46 MB
/// - Feature array: ~380-570 MB
///
/// # Example
///
/// ```
/// use geo_marching_squares::GridCell;
///
/// let cell = GridCell {
///     lon: -122.4194,
///     lat: 37.7749,
///     value: 15.5,
/// };
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GridCell {
    /// Longitude (x-coordinate) in decimal degrees
    pub lon: f64,

    /// Latitude (y-coordinate) in decimal degrees
    pub lat: f64,

    /// Scalar value at this grid point
    pub value: f64,
}

impl GridCell {
    /// Create a new GridCell
    ///
    /// # Arguments
    ///
    /// * `lon` - Longitude in decimal degrees
    /// * `lat` - Latitude in decimal degrees
    /// * `value` - Scalar value at this point
    ///
    /// # Example
    ///
    /// ```
    /// use geo_marching_squares::GridCell;
    ///
    /// let cell = GridCell::new(-122.4194, 37.7749, 15.5);
    /// assert_eq!(cell.lon, -122.4194);
    /// assert_eq!(cell.lat, 37.7749);
    /// assert_eq!(cell.value, 15.5);
    /// ```
    pub fn new(lon: f64, lat: f64, value: f64) -> Self {
        Self { lon, lat, value }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_grid_cell_creation() {
        let cell = GridCell::new(-122.0, 37.0, 290.5);
        assert_eq!(cell.lon, -122.0);
        assert_eq!(cell.lat, 37.0);
        assert_eq!(cell.value, 290.5);
    }

    #[test]
    fn test_grid_cell_struct_literal() {
        let cell = GridCell {
            lon: -74.0,
            lat: 40.7,
            value: 285.3,
        };
        assert_eq!(cell.lon, -74.0);
        assert_eq!(cell.lat, 40.7);
        assert_eq!(cell.value, 285.3);
    }

    #[test]
    fn test_grid_cell_copy() {
        let cell1 = GridCell::new(-122.0, 37.0, 290.5);
        let cell2 = cell1; // Copy
        assert_eq!(cell1, cell2);
    }

    #[test]
    fn test_grid_cell_size() {
        use std::mem::size_of;
        // Should be exactly 24 bytes (3 × f64)
        assert_eq!(size_of::<GridCell>(), 24);
    }
}
