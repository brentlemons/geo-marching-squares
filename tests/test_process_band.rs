use geo_marching_squares::process_band;
use geojson::{Feature, Geometry, JsonObject, Value};

/// Helper function to create a GeoJSON Point feature with a value
fn create_point_feature(lon: f64, lat: f64, value: f64) -> Feature {
    let geometry = Geometry::new(Value::Point(vec![lon, lat]));
    let mut properties = JsonObject::new();
    properties.insert("value".to_string(), serde_json::json!(value));

    Feature {
        bbox: None,
        geometry: Some(geometry),
        id: None,
        properties: Some(properties),
        foreign_members: None,
    }
}

#[test]
fn test_process_band_simple_3x3_grid() {
    // Create a simple 3x3 grid with values that form a clear band
    // Grid layout (values):
    //   5   5   5
    //   5  15  15
    //   5  15  15
    //
    // With lower=10, upper=20, we should get a contour around the 15 values

    let grid = vec![
        vec![
            create_point_feature(0.0, 2.0, 5.0),
            create_point_feature(1.0, 2.0, 5.0),
            create_point_feature(2.0, 2.0, 5.0),
        ],
        vec![
            create_point_feature(0.0, 1.0, 5.0),
            create_point_feature(1.0, 1.0, 15.0),
            create_point_feature(2.0, 1.0, 15.0),
        ],
        vec![
            create_point_feature(0.0, 0.0, 5.0),
            create_point_feature(1.0, 0.0, 15.0),
            create_point_feature(2.0, 0.0, 15.0),
        ],
    ];

    let result = process_band(&grid, 10.0, 20.0);

    // Verify the result has properties
    assert!(result.properties.is_some());
    let props = result.properties.unwrap();
    assert_eq!(props.get("lower_level").unwrap(), &serde_json::json!(10.0));
    assert_eq!(props.get("upper_level").unwrap(), &serde_json::json!(20.0));

    // Verify geometry exists and is MultiPolygon
    assert!(result.geometry.is_some());
    let geometry = result.geometry.unwrap();
    match geometry.value {
        Value::MultiPolygon(polygons) => {
            // Should have at least one polygon
            assert!(!polygons.is_empty(), "Expected at least one polygon");

            // Each polygon should have at least an exterior ring
            for polygon in &polygons {
                assert!(!polygon.is_empty(), "Polygon should have exterior ring");
                assert!(
                    polygon[0].len() >= 3,
                    "Exterior ring should have at least 3 points"
                );

                // Verify polygon is closed (first point == last point)
                let first = &polygon[0][0];
                let last = &polygon[0][polygon[0].len() - 1];
                assert_eq!(first, last, "Polygon should be closed");

                // Verify coordinates are rounded to 5 decimal places
                for point in &polygon[0] {
                    let rounded_x = (point[0] * 100000.0).round() / 100000.0;
                    let rounded_y = (point[1] * 100000.0).round() / 100000.0;
                    assert_eq!(point[0], rounded_x, "X coordinate should be rounded");
                    assert_eq!(point[1], rounded_y, "Y coordinate should be rounded");
                }
            }
        }
        _ => panic!("Expected MultiPolygon geometry"),
    }
}

#[test]
fn test_process_band_all_below_threshold() {
    // All values below lower threshold - should produce empty MultiPolygon
    let grid = vec![
        vec![
            create_point_feature(0.0, 1.0, 5.0),
            create_point_feature(1.0, 1.0, 5.0),
        ],
        vec![
            create_point_feature(0.0, 0.0, 5.0),
            create_point_feature(1.0, 0.0, 5.0),
        ],
    ];

    let result = process_band(&grid, 10.0, 20.0);

    assert!(result.geometry.is_some());
    let geometry = result.geometry.unwrap();
    match geometry.value {
        Value::MultiPolygon(polygons) => {
            assert!(
                polygons.is_empty(),
                "Should have no polygons when all values below threshold"
            );
        }
        _ => panic!("Expected MultiPolygon geometry"),
    }
}

#[test]
fn test_process_band_all_above_threshold() {
    // All values above upper threshold - should produce empty MultiPolygon
    let grid = vec![
        vec![
            create_point_feature(0.0, 1.0, 25.0),
            create_point_feature(1.0, 1.0, 25.0),
        ],
        vec![
            create_point_feature(0.0, 0.0, 25.0),
            create_point_feature(1.0, 0.0, 25.0),
        ],
    ];

    let result = process_band(&grid, 10.0, 20.0);

    assert!(result.geometry.is_some());
    let geometry = result.geometry.unwrap();
    match geometry.value {
        Value::MultiPolygon(polygons) => {
            assert!(
                polygons.is_empty(),
                "Should have no polygons when all values above threshold"
            );
        }
        _ => panic!("Expected MultiPolygon geometry"),
    }
}

#[test]
fn test_process_band_all_within_band() {
    // All values within the band - should produce a square polygon (case 85: 1111)
    let grid = vec![
        vec![
            create_point_feature(0.0, 1.0, 15.0),
            create_point_feature(1.0, 1.0, 15.0),
        ],
        vec![
            create_point_feature(0.0, 0.0, 15.0),
            create_point_feature(1.0, 0.0, 15.0),
        ],
    ];

    let result = process_band(&grid, 10.0, 20.0);

    assert!(result.geometry.is_some());
    let geometry = result.geometry.unwrap();
    match geometry.value {
        Value::MultiPolygon(polygons) => {
            assert!(
                !polygons.is_empty(),
                "Should have a square polygon when all values within band"
            );
            // Should have exactly one polygon (the square)
            assert_eq!(polygons.len(), 1);
            // The exterior ring should have 5 points (4 corners + closing point)
            assert_eq!(
                polygons[0][0].len(),
                5,
                "Square should have 5 points (4 corners + closing)"
            );
        }
        _ => panic!("Expected MultiPolygon geometry"),
    }
}

#[test]
fn test_process_band_single_cell_triangle() {
    // 2x2 grid forming a single cell with triangle shape
    // Bottom-left corner is within band, others below
    let grid = vec![
        vec![
            create_point_feature(0.0, 1.0, 5.0),
            create_point_feature(1.0, 1.0, 5.0),
        ],
        vec![
            create_point_feature(0.0, 0.0, 15.0),
            create_point_feature(1.0, 0.0, 5.0),
        ],
    ];

    let result = process_band(&grid, 10.0, 20.0);

    assert!(result.geometry.is_some());
    let geometry = result.geometry.unwrap();
    match geometry.value {
        Value::MultiPolygon(polygons) => {
            assert!(
                !polygons.is_empty(),
                "Should have at least one polygon for triangle"
            );
        }
        _ => panic!("Expected MultiPolygon geometry"),
    }
}
