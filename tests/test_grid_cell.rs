use geo_marching_squares::{do_concurrent_from_cells, do_concurrent_lines_from_cells, GridCell};

/// Create a simple test grid with GridCell structs
fn create_test_grid() -> Vec<Vec<GridCell>> {
    // Create a 4x4 grid with a gradient pattern
    // Values increase from bottom-left (0.0) to top-right (15.0)
    let mut grid = Vec::new();

    for row in 0..4 {
        let mut row_vec = Vec::new();
        for col in 0..4 {
            let lon = -122.0 + (col as f64 * 0.1);
            let lat = 37.0 + (row as f64 * 0.1);
            let value = (row * 4 + col) as f64;

            row_vec.push(GridCell { lon, lat, value });
        }
        grid.push(row_vec);
    }

    grid
}

#[test]
fn test_grid_cell_isobands() {
    let grid = create_test_grid();
    let thresholds = vec![0.0, 5.0, 10.0, 15.0];

    let result = do_concurrent_from_cells(&grid, &thresholds);

    // Should create 3 isobands: 0-5, 5-10, 10-15
    assert!(!result.features.is_empty(), "Should generate some features");

    // Verify feature properties
    for feature in &result.features {
        let props = feature.properties.as_ref().unwrap();
        assert!(props.contains_key("lower_level"));
        assert!(props.contains_key("upper_level"));
    }
}

#[test]
fn test_grid_cell_isolines() {
    let grid = create_test_grid();
    let isovalues = vec![5.0, 10.0];

    let result = do_concurrent_lines_from_cells(&grid, &isovalues);

    // Should create isolines at values 5.0 and 10.0
    assert!(!result.features.is_empty(), "Should generate some features");

    // Verify feature properties
    for feature in &result.features {
        let props = feature.properties.as_ref().unwrap();
        assert!(props.contains_key("isovalue"));
    }
}

#[test]
fn test_grid_cell_memory_size() {
    use std::mem::size_of;

    // Verify GridCell is exactly 24 bytes (3 × f64)
    assert_eq!(size_of::<GridCell>(), 24);
}

#[test]
fn test_grid_cell_large_grid_simulation() {
    // Simulate HRRR-scale grid dimensions (smaller for test speed)
    let rows = 100;
    let cols = 100;

    let mut grid = Vec::with_capacity(rows);
    for row in 0..rows {
        let mut row_vec = Vec::with_capacity(cols);
        for col in 0..cols {
            let lon = -122.0 + (col as f64 * 0.01);
            let lat = 37.0 + (row as f64 * 0.01);
            let value = ((row + col) as f64).sin() * 10.0 + 15.0;

            row_vec.push(GridCell { lon, lat, value });
        }
        grid.push(row_vec);
    }

    let thresholds = vec![10.0, 15.0, 20.0];
    let result = do_concurrent_from_cells(&grid, &thresholds);

    assert!(!result.features.is_empty(), "Should generate features from large grid");
}

#[test]
fn test_grid_cell_coordinates_preserved() {
    let grid = create_test_grid();
    let thresholds = vec![5.0, 10.0];

    let result = do_concurrent_from_cells(&grid, &thresholds);

    // Verify that generated coordinates are within expected lon/lat bounds
    for feature in &result.features {
        if let Some(geojson::Geometry {
            value: geojson::Value::MultiPolygon(polygons),
            ..
        }) = &feature.geometry
        {
            for polygon in polygons {
                for ring in polygon {
                    for coord in ring {
                        // Longitude should be between -122.0 and -121.7
                        assert!(coord[0] >= -122.1 && coord[0] <= -121.6);
                        // Latitude should be between 37.0 and 37.3
                        assert!(coord[1] >= 36.9 && coord[1] <= 37.4);
                    }
                }
            }
        }
    }
}

#[test]
fn test_grid_cell_vs_feature_compatibility() {
    // This test verifies that GridCell and Feature APIs produce similar results
    use geojson::{Feature, Geometry, JsonObject, Value};
    use geo_marching_squares::do_concurrent;

    // Create small test grid
    let grid_size = 5;
    let mut grid_cells = Vec::new();
    let mut grid_features = Vec::new();

    for row in 0..grid_size {
        let mut row_cells = Vec::new();
        let mut row_features = Vec::new();

        for col in 0..grid_size {
            let lon = -122.0 + (col as f64 * 0.1);
            let lat = 37.0 + (row as f64 * 0.1);
            let value = (row + col) as f64;

            // GridCell version
            row_cells.push(GridCell { lon, lat, value });

            // Feature version
            let geometry = Geometry::new(Value::Point(vec![lon, lat]));
            let mut properties = JsonObject::new();
            properties.insert("value".to_string(), serde_json::json!(value));

            row_features.push(Feature {
                bbox: None,
                geometry: Some(geometry),
                id: None,
                properties: Some(properties),
                foreign_members: None,
            });
        }

        grid_cells.push(row_cells);
        grid_features.push(row_features);
    }

    let thresholds = vec![0.0, 3.0, 6.0];

    let result_cells = do_concurrent_from_cells(&grid_cells, &thresholds);
    let result_features = do_concurrent(&grid_features, &thresholds);

    // Both should generate the same number of features
    assert_eq!(
        result_cells.features.len(),
        result_features.features.len(),
        "GridCell and Feature APIs should produce same number of features"
    );
}
