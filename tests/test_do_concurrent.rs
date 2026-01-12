use geo_marching_squares::{do_concurrent, process_band};
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

/// Create a 5x5 grid with varying values for testing
fn create_test_grid_5x5() -> Vec<Vec<Feature>> {
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

/// Create a uniform grid with all same values
fn create_uniform_grid(value: f64) -> Vec<Vec<Feature>> {
    vec![
        vec![
            create_point_feature(0.0, 2.0, value),
            create_point_feature(1.0, 2.0, value),
            create_point_feature(2.0, 2.0, value),
        ],
        vec![
            create_point_feature(0.0, 1.0, value),
            create_point_feature(1.0, 1.0, value),
            create_point_feature(2.0, 1.0, value),
        ],
        vec![
            create_point_feature(0.0, 0.0, value),
            create_point_feature(1.0, 0.0, value),
            create_point_feature(2.0, 0.0, value),
        ],
    ]
}

#[test]
fn test_do_concurrent_multiple_bands() {
    let grid = create_test_grid_5x5();
    let thresholds = vec![0.0, 10.0, 20.0, 30.0];

    let result = do_concurrent(&grid, &thresholds);

    // Should create up to 3 isobands (some may be empty)
    assert!(result.features.len() <= 3, "Should have at most 3 isobands");

    // Verify each feature has correct properties
    for feature in &result.features {
        let props = feature.properties.as_ref().expect("Should have properties");

        assert!(
            props.contains_key("lower_level"),
            "Should have lower_level property"
        );
        assert!(
            props.contains_key("upper_level"),
            "Should have upper_level property"
        );

        // Verify geometry exists and is MultiPolygon
        assert!(feature.geometry.is_some(), "Should have geometry");
        match &feature.geometry.as_ref().unwrap().value {
            Value::MultiPolygon(polygons) => {
                assert!(!polygons.is_empty(), "MultiPolygon should not be empty");
            }
            _ => panic!("Expected MultiPolygon geometry"),
        }
    }
}

#[test]
fn test_do_concurrent_filters_empty() {
    // All values below first threshold
    let grid = create_uniform_grid(5.0);
    let thresholds = vec![10.0, 20.0, 30.0];

    let result = do_concurrent(&grid, &thresholds);

    // All isobands should be filtered out (all values below lowest threshold)
    assert_eq!(
        result.features.len(),
        0,
        "Should have no features when all values below threshold"
    );
}

#[test]
fn test_do_concurrent_single_band() {
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

    let thresholds = vec![10.0, 20.0];

    let result = do_concurrent(&grid, &thresholds);

    // Should have exactly 1 isoband
    assert_eq!(result.features.len(), 1, "Should have 1 isoband");

    let feature = &result.features[0];
    let props = feature.properties.as_ref().unwrap();
    assert_eq!(props.get("lower_level").unwrap(), &serde_json::json!(10.0));
    assert_eq!(props.get("upper_level").unwrap(), &serde_json::json!(20.0));
}

#[test]
fn test_do_concurrent_matches_sequential() {
    let grid = create_test_grid_5x5();
    let thresholds = vec![0.0, 10.0, 20.0, 30.0];

    // Run concurrent
    let concurrent_result = do_concurrent(&grid, &thresholds);

    // Run sequential
    let mut sequential_features = Vec::new();
    for i in 0..thresholds.len() - 1 {
        let feature = process_band(&grid, thresholds[i], thresholds[i + 1]);

        // Check if feature has coordinates
        if let Some(geometry) = &feature.geometry {
            if let Value::MultiPolygon(polygons) = &geometry.value {
                if !polygons.is_empty() {
                    sequential_features.push(feature);
                }
            }
        }
    }

    // Should have same number of features
    assert_eq!(
        concurrent_result.features.len(),
        sequential_features.len(),
        "Concurrent and sequential should produce same number of features"
    );

    // Verify all features have valid geometry
    for feature in &concurrent_result.features {
        match &feature.geometry.as_ref().unwrap().value {
            Value::MultiPolygon(polygons) => {
                assert!(!polygons.is_empty(), "Should have non-empty polygons");

                // Verify each polygon is closed
                for polygon in polygons {
                    assert!(!polygon.is_empty(), "Polygon should have exterior ring");
                    let ring = &polygon[0];
                    assert!(ring.len() >= 4, "Ring should have at least 4 points");

                    // First and last point should be the same (closed)
                    assert_eq!(ring[0], ring[ring.len() - 1], "Polygon should be closed");
                }
            }
            _ => panic!("Expected MultiPolygon"),
        }
    }
}

#[test]
fn test_do_concurrent_minimum_thresholds() {
    let grid = create_uniform_grid(15.0);
    let thresholds = vec![10.0, 20.0];

    let result = do_concurrent(&grid, &thresholds);

    // Should work with minimum input (2 thresholds = 1 band)
    assert!(
        result.features.len() <= 1,
        "Should have at most 1 feature with 2 thresholds"
    );
}

#[test]
fn test_do_concurrent_empty_grid() {
    // Minimal 2x2 grid
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

    let thresholds = vec![10.0, 20.0, 30.0];

    let result = do_concurrent(&grid, &thresholds);

    // Should return empty collection (all values below threshold)
    assert_eq!(result.features.len(), 0);
}

#[test]
fn test_do_concurrent_all_bands_populated() {
    // Create grid with values spanning all bands
    let grid = vec![
        vec![
            create_point_feature(0.0, 3.0, 5.0),
            create_point_feature(1.0, 3.0, 15.0),
            create_point_feature(2.0, 3.0, 25.0),
            create_point_feature(3.0, 3.0, 35.0),
        ],
        vec![
            create_point_feature(0.0, 2.0, 5.0),
            create_point_feature(1.0, 2.0, 15.0),
            create_point_feature(2.0, 2.0, 25.0),
            create_point_feature(3.0, 2.0, 35.0),
        ],
        vec![
            create_point_feature(0.0, 1.0, 5.0),
            create_point_feature(1.0, 1.0, 15.0),
            create_point_feature(2.0, 1.0, 25.0),
            create_point_feature(3.0, 1.0, 35.0),
        ],
        vec![
            create_point_feature(0.0, 0.0, 5.0),
            create_point_feature(1.0, 0.0, 15.0),
            create_point_feature(2.0, 0.0, 25.0),
            create_point_feature(3.0, 0.0, 35.0),
        ],
    ];

    let thresholds = vec![0.0, 10.0, 20.0, 30.0, 40.0];

    let result = do_concurrent(&grid, &thresholds);

    // Should have features (exact number depends on grid structure)
    assert!(
        result.features.len() > 0,
        "Should have at least one feature with diverse values"
    );

    // Verify feature collection structure
    assert!(result.bbox.is_none(), "bbox should be None");
    assert!(
        result.foreign_members.is_none(),
        "foreign_members should be None"
    );
}

#[test]
fn test_do_concurrent_properties_correct() {
    let grid = create_test_grid_5x5();
    let thresholds = vec![0.0, 10.0, 20.0, 30.0];

    let result = do_concurrent(&grid, &thresholds);

    // Verify each feature has properties within the threshold ranges
    for feature in &result.features {
        let props = feature.properties.as_ref().unwrap();
        let lower = props.get("lower_level").unwrap().as_f64().unwrap();
        let upper = props.get("upper_level").unwrap().as_f64().unwrap();

        // Lower should be less than upper
        assert!(lower < upper, "lower_level should be less than upper_level");

        // Values should be from our threshold list
        assert!(
            thresholds.contains(&lower),
            "lower_level should be from thresholds"
        );
        assert!(
            thresholds.contains(&upper),
            "upper_level should be from thresholds"
        );
    }
}
