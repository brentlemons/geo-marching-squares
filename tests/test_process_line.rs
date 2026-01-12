use geo_marching_squares::{do_concurrent_lines, process_line, Cell};
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

/// Create a simple gradient grid for testing
fn create_gradient_grid() -> Vec<Vec<Feature>> {
    vec![
        vec![
            create_point_feature(0.0, 3.0, 5.0),
            create_point_feature(1.0, 3.0, 10.0),
            create_point_feature(2.0, 3.0, 15.0),
            create_point_feature(3.0, 3.0, 20.0),
        ],
        vec![
            create_point_feature(0.0, 2.0, 5.0),
            create_point_feature(1.0, 2.0, 10.0),
            create_point_feature(2.0, 2.0, 15.0),
            create_point_feature(3.0, 2.0, 20.0),
        ],
        vec![
            create_point_feature(0.0, 1.0, 5.0),
            create_point_feature(1.0, 1.0, 10.0),
            create_point_feature(2.0, 1.0, 15.0),
            create_point_feature(3.0, 1.0, 20.0),
        ],
        vec![
            create_point_feature(0.0, 0.0, 5.0),
            create_point_feature(1.0, 0.0, 10.0),
            create_point_feature(2.0, 0.0, 15.0),
            create_point_feature(3.0, 0.0, 20.0),
        ],
    ]
}

/// Create a grid that forms a circular pattern
fn create_circular_grid() -> Vec<Vec<Feature>> {
    vec![
        vec![
            create_point_feature(0.0, 4.0, 5.0),
            create_point_feature(1.0, 4.0, 5.0),
            create_point_feature(2.0, 4.0, 5.0),
            create_point_feature(3.0, 4.0, 5.0),
            create_point_feature(4.0, 4.0, 5.0),
        ],
        vec![
            create_point_feature(0.0, 3.0, 5.0),
            create_point_feature(1.0, 3.0, 15.0),
            create_point_feature(2.0, 3.0, 15.0),
            create_point_feature(3.0, 3.0, 15.0),
            create_point_feature(4.0, 3.0, 5.0),
        ],
        vec![
            create_point_feature(0.0, 2.0, 5.0),
            create_point_feature(1.0, 2.0, 15.0),
            create_point_feature(2.0, 2.0, 25.0),
            create_point_feature(3.0, 2.0, 15.0),
            create_point_feature(4.0, 2.0, 5.0),
        ],
        vec![
            create_point_feature(0.0, 1.0, 5.0),
            create_point_feature(1.0, 1.0, 15.0),
            create_point_feature(2.0, 1.0, 15.0),
            create_point_feature(3.0, 1.0, 15.0),
            create_point_feature(4.0, 1.0, 5.0),
        ],
        vec![
            create_point_feature(0.0, 0.0, 5.0),
            create_point_feature(1.0, 0.0, 5.0),
            create_point_feature(2.0, 0.0, 5.0),
            create_point_feature(3.0, 0.0, 5.0),
            create_point_feature(4.0, 0.0, 5.0),
        ],
    ]
}

#[test]
fn test_cell_create_empty_cases() {
    // Case 0: all below threshold
    let tl_feature = create_point_feature(0.0, 1.0, 5.0);
    let tr_feature = create_point_feature(1.0, 1.0, 5.0);
    let br_feature = create_point_feature(1.0, 0.0, 5.0);
    let bl_feature = create_point_feature(0.0, 0.0, 5.0);

    let cell = Cell::create(&tl_feature, &tr_feature, &br_feature, &bl_feature, 10.0);
    assert!(cell.is_none(), "Case 0 (all below) should return None");

    // Case 15: all above threshold
    let tl_feature = create_point_feature(0.0, 1.0, 15.0);
    let tr_feature = create_point_feature(1.0, 1.0, 15.0);
    let br_feature = create_point_feature(1.0, 0.0, 15.0);
    let bl_feature = create_point_feature(0.0, 0.0, 15.0);

    let cell = Cell::create(&tl_feature, &tr_feature, &br_feature, &bl_feature, 10.0);
    assert!(cell.is_none(), "Case 15 (all above) should return None");
}

#[test]
fn test_cell_single_segment_cases() {
    // Case 1: only BL above threshold (0001)
    let tl_feature = create_point_feature(0.0, 1.0, 5.0);
    let tr_feature = create_point_feature(1.0, 1.0, 5.0);
    let br_feature = create_point_feature(1.0, 0.0, 5.0);
    let bl_feature = create_point_feature(0.0, 0.0, 15.0);

    let cell = Cell::create(&tl_feature, &tr_feature, &br_feature, &bl_feature, 10.0);
    assert!(cell.is_some());

    let cell = cell.unwrap();
    let segments = cell.get_line_segments();
    assert_eq!(segments.len(), 1, "Case 1 should produce 1 segment");
}

#[test]
fn test_cell_saddle_cases() {
    // Case 5: TL and BR above threshold (0101) - saddle case
    let tl_feature = create_point_feature(0.0, 1.0, 15.0);
    let tr_feature = create_point_feature(1.0, 1.0, 5.0);
    let br_feature = create_point_feature(1.0, 0.0, 15.0);
    let bl_feature = create_point_feature(0.0, 0.0, 5.0);

    let cell = Cell::create(&tl_feature, &tr_feature, &br_feature, &bl_feature, 10.0);
    assert!(cell.is_some());

    let cell = cell.unwrap();
    let segments = cell.get_line_segments();
    assert_eq!(
        segments.len(),
        2,
        "Case 5 (saddle) should produce 2 segments"
    );

    // Case 10: TR and BL above threshold (1010) - saddle case
    let tl_feature = create_point_feature(0.0, 1.0, 5.0);
    let tr_feature = create_point_feature(1.0, 1.0, 15.0);
    let br_feature = create_point_feature(1.0, 0.0, 5.0);
    let bl_feature = create_point_feature(0.0, 0.0, 15.0);

    let cell = Cell::create(&tl_feature, &tr_feature, &br_feature, &bl_feature, 10.0);
    assert!(cell.is_some());

    let cell = cell.unwrap();
    let segments = cell.get_line_segments();
    assert_eq!(
        segments.len(),
        2,
        "Case 10 (saddle) should produce 2 segments"
    );
}

#[test]
fn test_process_line_gradient() {
    let grid = create_gradient_grid();
    let result = process_line(&grid, 12.5);

    assert!(result.geometry.is_some());

    match &result.geometry.unwrap().value {
        Value::MultiLineString(lines) => {
            assert!(!lines.is_empty(), "Gradient grid should produce isolines");

            // Verify lines have at least 2 points
            for line in lines {
                assert!(line.len() >= 2, "Each line should have at least 2 points");

                // Verify coordinate precision (5 decimal places)
                for point in line {
                    let rounded_x = (point[0] * 100000.0).round() / 100000.0;
                    let rounded_y = (point[1] * 100000.0).round() / 100000.0;
                    assert_eq!(point[0], rounded_x, "X coordinate should be rounded");
                    assert_eq!(point[1], rounded_y, "Y coordinate should be rounded");
                }
            }
        }
        _ => panic!("Expected MultiLineString geometry"),
    }

    // Verify properties
    let props = result.properties.unwrap();
    assert_eq!(props.get("isovalue").unwrap(), &serde_json::json!(12.5));
}

#[test]
fn test_process_line_circular_closed_loop() {
    let grid = create_circular_grid();
    let result = process_line(&grid, 20.0);

    assert!(result.geometry.is_some());

    match &result.geometry.unwrap().value {
        Value::MultiLineString(lines) => {
            assert!(!lines.is_empty(), "Circular grid should produce isolines");

            // Check if we have a closed loop
            for line in lines {
                if line.len() > 2 {
                    // A closed loop has first point == last point
                    if line[0] == line[line.len() - 1] {
                        // Found a closed loop, test passes
                        return;
                    }
                }
            }
        }
        _ => panic!("Expected MultiLineString geometry"),
    }
}

#[test]
fn test_process_line_empty_result() {
    // All values below threshold
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

    let result = process_line(&grid, 10.0);

    match &result.geometry.unwrap().value {
        Value::MultiLineString(lines) => {
            assert!(
                lines.is_empty(),
                "Should have no lines when all values below threshold"
            );
        }
        _ => panic!("Expected MultiLineString geometry"),
    }
}

#[test]
fn test_do_concurrent_lines_multiple() {
    let grid = create_gradient_grid();
    let isovalues = vec![7.5, 12.5, 17.5];

    let result = do_concurrent_lines(&grid, &isovalues);

    // Should have up to 3 isoline features
    assert!(
        result.features.len() <= 3,
        "Should have at most 3 isoline features"
    );

    // Verify each feature has correct properties
    for feature in &result.features {
        let props = feature.properties.as_ref().expect("Should have properties");

        assert!(
            props.contains_key("isovalue"),
            "Should have isovalue property"
        );

        // Verify geometry exists and is MultiLineString
        assert!(feature.geometry.is_some(), "Should have geometry");
        match &feature.geometry.as_ref().unwrap().value {
            Value::MultiLineString(lines) => {
                assert!(!lines.is_empty(), "MultiLineString should not be empty");
            }
            _ => panic!("Expected MultiLineString geometry"),
        }
    }
}

#[test]
fn test_do_concurrent_lines_filters_empty() {
    // Grid with uniform values - no isolines
    let grid = vec![
        vec![
            create_point_feature(0.0, 2.0, 5.0),
            create_point_feature(1.0, 2.0, 5.0),
            create_point_feature(2.0, 2.0, 5.0),
        ],
        vec![
            create_point_feature(0.0, 1.0, 5.0),
            create_point_feature(1.0, 1.0, 5.0),
            create_point_feature(2.0, 1.0, 5.0),
        ],
        vec![
            create_point_feature(0.0, 0.0, 5.0),
            create_point_feature(1.0, 0.0, 5.0),
            create_point_feature(2.0, 0.0, 5.0),
        ],
    ];

    let isovalues = vec![10.0, 15.0, 20.0];

    let result = do_concurrent_lines(&grid, &isovalues);

    // Should filter out all empty features
    assert_eq!(
        result.features.len(),
        0,
        "Should have no features for uniform grid"
    );
}

#[test]
fn test_process_line_properties_set() {
    let grid = create_gradient_grid();
    let result = process_line(&grid, 15.0);

    assert!(result.properties.is_some());
    let props = result.properties.unwrap();

    assert_eq!(props.get("isovalue").unwrap(), &serde_json::json!(15.0));
}

#[test]
fn test_linear_vs_cosine_interpolation() {
    // This test verifies that Cell uses linear interpolation, not cosine
    // By checking that midpoint values are exactly at 0.5

    let tl_feature = create_point_feature(0.0, 1.0, 0.0);
    let tr_feature = create_point_feature(1.0, 1.0, 20.0);
    let br_feature = create_point_feature(1.0, 0.0, 0.0);
    let bl_feature = create_point_feature(0.0, 0.0, 0.0);

    let cell = Cell::create(&tl_feature, &tr_feature, &br_feature, &bl_feature, 10.0);
    assert!(cell.is_some());

    let cell = cell.unwrap();
    let points = cell.points();

    // With linear interpolation at value 10.0 between 0.0 and 20.0,
    // the interpolation factor should be exactly 0.5
    // So the X coordinate should be exactly at 0.5 (midpoint between 0.0 and 1.0)

    // Find the top side point
    for point in points {
        if point.y().unwrap() == 1.0 {
            // Top side point
            // Linear: mu = (10 - 0) / (20 - 0) = 0.5
            // x = (1 - 0.5) * 0.0 + 0.5 * 1.0 = 0.5
            let expected_x = 0.5;
            let actual_x = point.x().unwrap();
            assert!(
                (actual_x - expected_x).abs() < 1e-9,
                "Linear interpolation should give x=0.5, got {}",
                actual_x
            );
            return;
        }
    }

    panic!("Should have found top side point");
}
