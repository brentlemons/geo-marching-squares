# zarr-contour Lambda Implementation Plan

**Status**: Ready for Implementation
**Created**: 2025-01-11
**Lambda Name**: `zarr-contour`

---

## Executive Summary

Convert geo-marching-squares into an AWS Lambda that generates contour polygons from Zarr weather grids. The Lambda operates in two modes:

- **Coordinator**: Receives request with thresholds, fans out to workers, assembles results
- **Worker**: Fetches Zarr data, generates single band contour, returns polygons

### Expected Performance

| Metric | Current (grib-inspector) | With zarr-contour |
|--------|--------------------------|-------------------|
| 13 bands total | ~8,000ms | ~1,100ms |
| Contour generation | 6,907ms (sequential) | ~530ms (parallel) |
| Speedup | 1x | **~7x** |

---

## Part 1: Architecture Overview

### Request Flow

```
┌─────────────────────────────────────────────────────────────────────┐
│                     GRIB-INSPECTOR API                               │
│  - Receives contour request from client                             │
│  - Looks up grid metadata from OpenSearch                           │
│  - Constructs coordinator request                                   │
│  - Invokes zarr-contour Lambda (coordinator mode)                   │
│  - Receives grid-space polygons                                     │
│  - Applies PROJ transformation to lat/lon                           │
│  - Applies clipping, simplification                                 │
│  - Returns GeoJSON FeatureCollection                                │
└─────────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────────┐
│                    ZARR-CONTOUR LAMBDA                               │
│                    (Coordinator Mode)                                │
├─────────────────────────────────────────────────────────────────────┤
│  1. Parse coordinator request                                        │
│  2. Generate band intervals from thresholds                         │
│  3. Invoke N worker Lambdas (synchronous, parallel)                 │
│  4. Collect all worker responses                                    │
│  5. Assemble and gzip response                                      │
│  6. Return to caller                                                │
└─────────────────────────────────────────────────────────────────────┘
                              │
              ┌───────────────┼───────────────┐
              ▼               ▼               ▼
┌───────────────────┐ ┌───────────────────┐ ┌───────────────────┐
│   ZARR-CONTOUR    │ │   ZARR-CONTOUR    │ │   ZARR-CONTOUR    │
│   (Worker Mode)   │ │   (Worker Mode)   │ │   (Worker Mode)   │
├───────────────────┤ ├───────────────────┤ ├───────────────────┤
│ 1. Fetch Zarr     │ │ 1. Fetch Zarr     │ │ 1. Fetch Zarr     │
│ 2. Convert to f32 │ │ 2. Convert to f32 │ │ 2. Convert to f32 │
│ 3. Run contour    │ │ 3. Run contour    │ │ 3. Run contour    │
│ 4. Orient polygons│ │ 4. Orient polygons│ │ 4. Orient polygons│
│ 5. Gzip & return  │ │ 5. Gzip & return  │ │ 5. Gzip & return  │
└───────────────────┘ └───────────────────┘ └───────────────────┘
        │                     │                     │
        └─────────────────────┴─────────────────────┘
                              │
                      ~530ms per worker
                      (all run in parallel)
```

---

## Part 2: Data Structures

### 2.1 Request Types

```rust
use serde::{Deserialize, Serialize};

/// Top-level request - mode determines which variant is used
#[derive(Debug, Deserialize)]
#[serde(tag = "mode")]
pub enum ContourRequest {
    #[serde(rename = "coordinator")]
    Coordinator(CoordinatorRequest),

    #[serde(rename = "worker")]
    Worker(WorkerRequest),
}

/// Coordinator request - receives thresholds, fans out to workers
#[derive(Debug, Deserialize)]
pub struct CoordinatorRequest {
    /// S3 bucket containing Zarr data
    pub zarr_bucket: String,

    /// S3 key/prefix for the Zarr array
    pub zarr_key: String,

    /// Grid dimensions
    pub width: usize,
    pub height: usize,

    /// Grid coordinate system (for output metadata)
    pub x_origin: f64,
    pub y_origin: f64,
    pub x_step: f64,
    pub y_step: f64,

    /// Threshold values (N thresholds = N-1 bands)
    /// Example: [0, 10, 20, 30] creates bands 0-10, 10-20, 20-30
    pub thresholds: Vec<f32>,

    /// Coordinate precision (decimal places, default 5)
    #[serde(default = "default_precision")]
    pub precision: u32,

    /// Name of this Lambda function (for invoking workers)
    pub function_name: String,
}

fn default_precision() -> u32 { 5 }

/// Worker request - processes a single band
#[derive(Debug, Serialize, Deserialize)]
pub struct WorkerRequest {
    /// S3 bucket containing Zarr data
    pub zarr_bucket: String,

    /// S3 key/prefix for the Zarr array
    pub zarr_key: String,

    /// Grid dimensions
    pub width: usize,
    pub height: usize,

    /// Grid coordinate system
    pub x_origin: f64,
    pub y_origin: f64,
    pub x_step: f64,
    pub y_step: f64,

    /// Band thresholds
    pub lower: f32,
    pub upper: f32,

    /// Coordinate precision (decimal places)
    pub precision: u32,
}
```

### 2.2 Response Types

```rust
/// Coordinator response - collection of all bands
#[derive(Debug, Serialize)]
pub struct CoordinatorResponse {
    /// Processing statistics
    pub stats: CoordinatorStats,

    /// All band results
    pub bands: Vec<BandResult>,
}

#[derive(Debug, Serialize)]
pub struct CoordinatorStats {
    pub total_ms: f64,
    pub num_bands: usize,
    pub num_workers: usize,
}

/// Worker response - single band result
#[derive(Debug, Serialize, Deserialize)]
pub struct WorkerResponse {
    /// Processing statistics
    pub stats: WorkerStats,

    /// Band result
    pub band: BandResult,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct WorkerStats {
    pub zarr_fetch_ms: f64,
    pub contour_ms: f64,
    pub total_ms: f64,
}

/// Single band's polygon data (grid-space coordinates)
#[derive(Debug, Serialize, Deserialize)]
pub struct BandResult {
    /// Lower threshold (inclusive)
    pub lower: f32,

    /// Upper threshold (exclusive)
    pub upper: f32,

    /// Polygons in grid-space coordinates
    /// Structure: [polygon][ring][point][x,y]
    /// - polygon[0] is exterior ring (CCW in grid space)
    /// - polygon[1..] are holes (CW in grid space)
    pub polygons: Vec<Polygon>,
}

/// A polygon with exterior ring and optional holes
/// Using f32 for grid-space coordinates to reduce payload size
#[derive(Debug, Serialize, Deserialize)]
pub struct Polygon {
    /// Exterior ring (CCW winding in grid space)
    pub exterior: Vec<[f32; 2]>,

    /// Interior rings / holes (CW winding in grid space)
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub holes: Vec<Vec<[f32; 2]>>,
}
```

---

## Part 3: Lambda Handler

### 3.1 Main Entry Point

```rust
// src/bin/lambda.rs

use lambda_runtime::{service_fn, Error, LambdaEvent};
use zarr_contour::{ContourRequest, handle_coordinator, handle_worker};

#[tokio::main]
async fn main() -> Result<(), Error> {
    // Initialize tracing
    tracing_subscriber::fmt()
        .with_env_filter(tracing_subscriber::EnvFilter::from_default_env())
        .json()
        .init();

    lambda_runtime::run(service_fn(handler)).await
}

async fn handler(event: LambdaEvent<ContourRequest>) -> Result<Vec<u8>, Error> {
    let (request, _context) = event.into_parts();

    let response = match request {
        ContourRequest::Coordinator(req) => handle_coordinator(req).await?,
        ContourRequest::Worker(req) => handle_worker(req).await?,
    };

    // Always gzip the response
    let compressed = gzip_compress(&response)?;

    Ok(compressed)
}

fn gzip_compress(data: &[u8]) -> Result<Vec<u8>, Error> {
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use std::io::Write;

    let mut encoder = GzEncoder::new(Vec::new(), Compression::fast());
    encoder.write_all(data)?;
    Ok(encoder.finish()?)
}
```

### 3.2 Coordinator Handler

```rust
// src/coordinator.rs

use aws_sdk_lambda::Client as LambdaClient;
use aws_sdk_lambda::primitives::Blob;
use futures::future::join_all;
use std::time::Instant;

pub async fn handle_coordinator(req: CoordinatorRequest) -> Result<Vec<u8>, Error> {
    let start = Instant::now();

    // Initialize AWS Lambda client
    let config = aws_config::load_from_env().await;
    let lambda_client = LambdaClient::new(&config);

    // Generate band intervals from thresholds
    let bands: Vec<(f32, f32)> = req.thresholds
        .windows(2)
        .map(|w| (w[0], w[1]))
        .collect();

    let num_bands = bands.len();

    // Create worker requests
    let worker_requests: Vec<WorkerRequest> = bands
        .iter()
        .map(|(lower, upper)| WorkerRequest {
            zarr_bucket: req.zarr_bucket.clone(),
            zarr_key: req.zarr_key.clone(),
            width: req.width,
            height: req.height,
            x_origin: req.x_origin,
            y_origin: req.y_origin,
            x_step: req.x_step,
            y_step: req.y_step,
            lower: *lower,
            upper: *upper,
            precision: req.precision,
        })
        .collect();

    // Invoke all workers in parallel
    let worker_futures: Vec<_> = worker_requests
        .into_iter()
        .map(|worker_req| {
            invoke_worker(&lambda_client, &req.function_name, worker_req)
        })
        .collect();

    let worker_results = join_all(worker_futures).await;

    // Collect results, fail if any worker failed
    let mut band_results = Vec::with_capacity(num_bands);
    for result in worker_results {
        let worker_response = result?;
        band_results.push(worker_response.band);
    }

    // Build coordinator response
    let response = CoordinatorResponse {
        stats: CoordinatorStats {
            total_ms: start.elapsed().as_secs_f64() * 1000.0,
            num_bands,
            num_workers: num_bands,
        },
        bands: band_results,
    };

    // Serialize to JSON
    let json = serde_json::to_vec(&response)?;

    Ok(json)
}

async fn invoke_worker(
    client: &LambdaClient,
    function_name: &str,
    request: WorkerRequest,
) -> Result<WorkerResponse, Error> {
    // Wrap in worker mode envelope
    let payload = serde_json::json!({
        "mode": "worker",
        "zarr_bucket": request.zarr_bucket,
        "zarr_key": request.zarr_key,
        "width": request.width,
        "height": request.height,
        "x_origin": request.x_origin,
        "y_origin": request.y_origin,
        "x_step": request.x_step,
        "y_step": request.y_step,
        "lower": request.lower,
        "upper": request.upper,
        "precision": request.precision,
    });

    let payload_bytes = serde_json::to_vec(&payload)?;

    // Invoke synchronously with 30 second timeout
    let response = tokio::time::timeout(
        std::time::Duration::from_secs(30),
        client
            .invoke()
            .function_name(function_name)
            .payload(Blob::new(payload_bytes))
            .send()
    ).await??;

    // Check for function error
    if let Some(error) = response.function_error() {
        return Err(format!("Worker error: {}", error).into());
    }

    // Decompress and parse response
    let response_payload = response.payload()
        .ok_or("No payload in worker response")?;

    let decompressed = gzip_decompress(response_payload.as_ref())?;
    let worker_response: WorkerResponse = serde_json::from_slice(&decompressed)?;

    Ok(worker_response)
}

fn gzip_decompress(data: &[u8]) -> Result<Vec<u8>, Error> {
    use flate2::read::GzDecoder;
    use std::io::Read;

    let mut decoder = GzDecoder::new(data);
    let mut decompressed = Vec::new();
    decoder.read_to_end(&mut decompressed)?;
    Ok(decompressed)
}
```

### 3.3 Worker Handler

```rust
// src/worker.rs

use std::time::Instant;

pub async fn handle_worker(req: WorkerRequest) -> Result<Vec<u8>, Error> {
    let start = Instant::now();

    // Phase 1: Fetch Zarr data
    let zarr_start = Instant::now();
    let grid_data = fetch_zarr_grid(&req.zarr_bucket, &req.zarr_key).await?;
    let zarr_fetch_ms = zarr_start.elapsed().as_secs_f64() * 1000.0;

    // Validate dimensions
    if grid_data.len() != req.width * req.height {
        return Err(format!(
            "Grid size mismatch: expected {}x{}={}, got {}",
            req.width, req.height, req.width * req.height, grid_data.len()
        ).into());
    }

    // Phase 2: Generate contour
    let contour_start = Instant::now();
    let polygons = generate_band_contour(
        &grid_data,
        req.width,
        req.height,
        req.x_origin,
        req.y_origin,
        req.x_step,
        req.y_step,
        req.lower,
        req.upper,
        req.precision,
    )?;
    let contour_ms = contour_start.elapsed().as_secs_f64() * 1000.0;

    // Build response
    let response = WorkerResponse {
        stats: WorkerStats {
            zarr_fetch_ms,
            contour_ms,
            total_ms: start.elapsed().as_secs_f64() * 1000.0,
        },
        band: BandResult {
            lower: req.lower,
            upper: req.upper,
            polygons,
        },
    };

    // Serialize to JSON
    let json = serde_json::to_vec(&response)?;

    Ok(json)
}
```

---

## Part 4: Core Algorithm Changes

### 4.1 New Flat-Array Contour Function

Add to `src/marching_squares.rs`:

```rust
/// Process a single isoband from a flat f32 array
///
/// Returns polygons in grid-space coordinates (f32).
/// Coordinates represent positions in the grid where:
/// - x = column index (0 to width-1, with sub-pixel interpolation)
/// - y = row index (0 to height-1, with sub-pixel interpolation)
///
/// # Arguments
/// * `values` - Flat row-major array of grid values (f32 from Zarr)
/// * `width` - Number of columns in the grid
/// * `height` - Number of rows in the grid
/// * `lower` - Lower threshold (inclusive)
/// * `upper` - Upper threshold (exclusive)
/// * `precision` - Decimal places for coordinate rounding
///
/// # Returns
/// Vector of polygons, each with exterior ring and optional holes.
/// Winding order follows RFC 7946 (exterior CCW, holes CW).
pub fn process_band_flat(
    values: &[f32],
    width: usize,
    height: usize,
    lower: f32,
    upper: f32,
    precision: u32,
) -> Vec<Polygon> {
    assert_eq!(values.len(), width * height, "Grid size mismatch");

    let rows = height;
    let cols = width;

    // Create cell grid (rows-1 × cols-1 cells)
    let mut cells: Vec<Vec<Option<Shape>>> = Vec::with_capacity(rows - 1);

    for r in 0..rows - 1 {
        let mut row_cells = Vec::with_capacity(cols - 1);
        for c in 0..cols - 1 {
            // Get corner values
            let tl = values[r * cols + c];
            let tr = values[r * cols + c + 1];
            let bl = values[(r + 1) * cols + c];
            let br = values[(r + 1) * cols + c + 1];

            // Skip cells with NaN values
            if tl.is_nan() || tr.is_nan() || bl.is_nan() || br.is_nan() {
                row_cells.push(None);
                continue;
            }

            // Create shape using grid indices as coordinates
            // Top-left corner is at (c, r), bottom-right at (c+1, r+1)
            let shape = Shape::create_from_corners(
                c as f64, r as f64,           // top-left (x, y)
                (c + 1) as f64, r as f64,     // top-right
                (c + 1) as f64, (r + 1) as f64, // bottom-right
                c as f64, (r + 1) as f64,     // bottom-left
                tl as f64, tr as f64, bl as f64, br as f64, // corner values
                lower as f64,
                upper as f64,
                c, r,
                r == 0,           // top_edge
                c == cols - 2,    // right_edge
                r == rows - 2,    // bottom_edge
                c == 0,           // left_edge
            );

            row_cells.push(shape);
        }
        cells.push(row_cells);
    }

    // Walk polygons using existing algorithm
    let raw_polygons = walk_all_polygons(&mut cells, precision);

    // Resolve nesting (exterior vs holes)
    let nested_polygons = resolve_polygon_nesting(raw_polygons);

    // Orient polygons (RFC 7946: exterior CCW, holes CW)
    let oriented = orient_polygons_grid_space(nested_polygons);

    // Convert to output format with f32 coordinates
    convert_to_output_polygons(oriented, precision)
}

/// Convert internal polygon format to output Polygon structs with f32 coords
fn convert_to_output_polygons(
    polygons: Vec<Vec<Vec<Position>>>,
    precision: u32,
) -> Vec<Polygon> {
    let factor = 10_f32.powi(precision as i32);

    polygons
        .into_iter()
        .map(|rings| {
            let exterior = rings[0]
                .iter()
                .map(|pos| {
                    [
                        (pos[0] as f32 * factor).round() / factor,
                        (pos[1] as f32 * factor).round() / factor,
                    ]
                })
                .collect();

            let holes = rings[1..]
                .iter()
                .map(|ring| {
                    ring.iter()
                        .map(|pos| {
                            [
                                (pos[0] as f32 * factor).round() / factor,
                                (pos[1] as f32 * factor).round() / factor,
                            ]
                        })
                        .collect()
                })
                .collect();

            Polygon { exterior, holes }
        })
        .collect()
}
```

### 4.2 New Shape Constructor for Grid Coordinates

Add to `src/shape.rs`:

```rust
/// Create a shape from explicit corner coordinates and values
/// Used for flat-array processing where coordinates are grid indices
#[allow(clippy::too_many_arguments)]
pub fn create_from_corners(
    tl_x: f64, tl_y: f64,
    tr_x: f64, tr_y: f64,
    br_x: f64, br_y: f64,
    bl_x: f64, bl_y: f64,
    tl_val: f64, tr_val: f64, bl_val: f64, br_val: f64,
    lower: f64,
    upper: f64,
    x: usize,
    y: usize,
    top_edge: bool,
    right_edge: bool,
    bottom_edge: bool,
    left_edge: bool,
) -> Option<Shape> {
    // Compute ternary classification
    let classify = |v: f64| -> u8 {
        if v < lower { 0 }
        else if v >= upper { 2 }
        else { 1 }
    };

    let tl_class = classify(tl_val);
    let tr_class = classify(tr_val);
    let br_class = classify(br_val);
    let bl_class = classify(bl_val);

    // Compute shape value (ternary encoding)
    let value = tl_class * 64 + tr_class * 16 + br_class * 4 + bl_class;

    // Skip empty cases (all below or all above)
    if value == 0 || value == 170 {
        return None;
    }

    // Create corner points
    let top_left = Point::new(tl_x, tl_y);
    let top_right = Point::new(tr_x, tr_y);
    let bottom_right = Point::new(br_x, br_y);
    let bottom_left = Point::new(bl_x, bl_y);

    // Build shape with edges
    let mut shape = Shape::new_internal(
        top_left,
        top_right,
        bottom_right,
        bottom_left,
        value,
        lower,
        upper,
        x,
        y,
        top_edge,
        right_edge,
        bottom_edge,
        left_edge,
        tl_val,
        tr_val,
        bl_val,
        br_val,
    );

    shape.build_edges();

    Some(shape)
}
```

---

## Part 5: Zarr Integration

### 5.1 Zarr Fetching (copied from grib-inspector)

```rust
// src/zarr.rs

use aws_sdk_s3::Client as S3Client;
use bytes::Bytes;

/// Fetch a Zarr array from S3 and return as flat f32 vector
pub async fn fetch_zarr_grid(bucket: &str, key: &str) -> Result<Vec<f32>, Error> {
    let config = aws_config::load_from_env().await;
    let s3_client = S3Client::new(&config);

    // Read .zarray metadata to understand chunking
    let zarray_key = format!("{}/.zarray", key);
    let zarray_meta = fetch_s3_object(&s3_client, bucket, &zarray_key).await?;
    let zarray: ZarrayMetadata = serde_json::from_slice(&zarray_meta)?;

    // Calculate chunks needed
    let shape = &zarray.shape;
    let chunks = &zarray.chunks;
    let dtype = &zarray.dtype;

    let num_chunk_rows = (shape[0] + chunks[0] - 1) / chunks[0];
    let num_chunk_cols = (shape[1] + chunks[1] - 1) / chunks[1];

    // Fetch all chunks in parallel
    let mut chunk_futures = Vec::new();
    for chunk_row in 0..num_chunk_rows {
        for chunk_col in 0..num_chunk_cols {
            let chunk_key = format!("{}/{}.{}", key, chunk_row, chunk_col);
            chunk_futures.push(fetch_and_decode_chunk(
                &s3_client,
                bucket,
                &chunk_key,
                &zarray.compressor,
                dtype,
                chunks[0] * chunks[1],
            ));
        }
    }

    let chunk_results = futures::future::join_all(chunk_futures).await;

    // Assemble chunks into full array
    let mut full_array = vec![0.0f32; shape[0] * shape[1]];

    let mut chunk_idx = 0;
    for chunk_row in 0..num_chunk_rows {
        for chunk_col in 0..num_chunk_cols {
            let chunk_data = chunk_results[chunk_idx].as_ref()?;

            // Copy chunk data into full array
            let row_start = chunk_row * chunks[0];
            let col_start = chunk_col * chunks[1];

            for local_row in 0..chunks[0] {
                let global_row = row_start + local_row;
                if global_row >= shape[0] {
                    break;
                }

                for local_col in 0..chunks[1] {
                    let global_col = col_start + local_col;
                    if global_col >= shape[1] {
                        break;
                    }

                    let local_idx = local_row * chunks[1] + local_col;
                    let global_idx = global_row * shape[1] + global_col;

                    full_array[global_idx] = chunk_data[local_idx];
                }
            }

            chunk_idx += 1;
        }
    }

    Ok(full_array)
}

#[derive(Deserialize)]
struct ZarrayMetadata {
    shape: Vec<usize>,
    chunks: Vec<usize>,
    dtype: String,
    compressor: Option<CompressorConfig>,
}

#[derive(Deserialize)]
struct CompressorConfig {
    id: String,
    // Other compressor-specific fields
}

async fn fetch_s3_object(
    client: &S3Client,
    bucket: &str,
    key: &str,
) -> Result<Bytes, Error> {
    let response = client
        .get_object()
        .bucket(bucket)
        .key(key)
        .send()
        .await?;

    let bytes = response.body.collect().await?.into_bytes();
    Ok(bytes)
}

async fn fetch_and_decode_chunk(
    client: &S3Client,
    bucket: &str,
    key: &str,
    compressor: &Option<CompressorConfig>,
    dtype: &str,
    expected_elements: usize,
) -> Result<Vec<f32>, Error> {
    let compressed = fetch_s3_object(client, bucket, key).await?;

    // Decompress based on compressor type
    let decompressed = match compressor {
        Some(c) if c.id == "blosc" => decompress_blosc(&compressed)?,
        Some(c) if c.id == "zlib" => decompress_zlib(&compressed)?,
        None => compressed.to_vec(),
        _ => return Err("Unsupported compressor".into()),
    };

    // Convert bytes to f32 based on dtype
    let values = match dtype {
        "<f4" | "f4" => bytes_to_f32_le(&decompressed),
        ">f4" => bytes_to_f32_be(&decompressed),
        _ => return Err(format!("Unsupported dtype: {}", dtype).into()),
    };

    Ok(values)
}

fn bytes_to_f32_le(bytes: &[u8]) -> Vec<f32> {
    bytes
        .chunks_exact(4)
        .map(|chunk| f32::from_le_bytes(chunk.try_into().unwrap()))
        .collect()
}

fn decompress_blosc(data: &[u8]) -> Result<Vec<u8>, Error> {
    // Use blosc crate for decompression
    let decompressed = blosc::decompress_bytes(data)?;
    Ok(decompressed)
}
```

---

## Part 6: Project Structure

### 6.1 Cargo.toml

```toml
[package]
name = "zarr-contour"
version = "0.1.0"
edition = "2021"

[[bin]]
name = "bootstrap"
path = "src/bin/lambda.rs"

[dependencies]
# Core
geo-marching-squares = { path = "../geo-marching-squares" }
serde = { version = "1.0", features = ["derive"] }
serde_json = "1.0"
tokio = { version = "1", features = ["full"] }
futures = "0.3"
tracing = "0.1"
tracing-subscriber = { version = "0.3", features = ["json", "env-filter"] }

# AWS
aws-config = "1.0"
aws-sdk-s3 = "1.0"
aws-sdk-lambda = "1.0"
lambda_runtime = "0.8"

# Compression
flate2 = "1.0"
blosc = "0.2"

# Zarr
bytes = "1.0"
```

### 6.2 File Structure

```
zarr-contour/
├── Cargo.toml
├── src/
│   ├── lib.rs              # Public API
│   ├── bin/
│   │   └── lambda.rs       # Lambda entry point
│   ├── coordinator.rs      # Coordinator mode handler
│   ├── worker.rs           # Worker mode handler
│   ├── zarr.rs             # Zarr S3 fetching
│   ├── types.rs            # Request/response types
│   └── contour.rs          # Contour generation wrapper
└── tests/
    ├── integration_tests.rs
    └── test_data/
```

---

## Part 7: grib-inspector Integration

### 7.1 Calling the Lambda from grib-inspector

```rust
// In grib-inspector/src/api/contour_service.rs

use aws_sdk_lambda::Client as LambdaClient;

impl ContourService {
    pub async fn generate_isobands_via_lambda(
        &self,
        zarr_bucket: &str,
        zarr_key: &str,
        grid_info: &GridInfo,
        thresholds: &[f32],
        precision: u32,
    ) -> Result<Vec<BandResult>> {
        let lambda_client = LambdaClient::new(&self.aws_config);

        // Build coordinator request
        let request = serde_json::json!({
            "mode": "coordinator",
            "zarr_bucket": zarr_bucket,
            "zarr_key": zarr_key,
            "width": grid_info.width,
            "height": grid_info.height,
            "x_origin": grid_info.x_origin,
            "y_origin": grid_info.y_origin,
            "x_step": grid_info.x_step,
            "y_step": grid_info.y_step,
            "thresholds": thresholds,
            "precision": precision,
            "function_name": "zarr-contour",
        });

        // Invoke Lambda
        let response = lambda_client
            .invoke()
            .function_name("zarr-contour")
            .payload(Blob::new(serde_json::to_vec(&request)?))
            .send()
            .await?;

        // Decompress and parse response
        let payload = response.payload().ok_or("No payload")?;
        let decompressed = gzip_decompress(payload.as_ref())?;
        let coord_response: CoordinatorResponse = serde_json::from_slice(&decompressed)?;

        // Transform grid-space polygons to lat/lon
        let geo_bands = self.transform_bands_to_geographic(
            &coord_response.bands,
            grid_info,
        ).await?;

        Ok(geo_bands)
    }

    async fn transform_bands_to_geographic(
        &self,
        bands: &[BandResult],
        grid_info: &GridInfo,
    ) -> Result<Vec<GeoBandResult>> {
        let proj = self.create_proj_transformer(grid_info)?;

        let geo_bands: Vec<GeoBandResult> = bands
            .iter()
            .map(|band| {
                let geo_polygons = band.polygons
                    .iter()
                    .map(|poly| {
                        // Transform exterior ring
                        let geo_exterior = self.transform_ring(&poly.exterior, &proj)?;

                        // Transform holes
                        let geo_holes: Vec<_> = poly.holes
                            .iter()
                            .map(|hole| self.transform_ring(hole, &proj))
                            .collect::<Result<Vec<_>>>()?;

                        Ok(GeoPolygon {
                            exterior: geo_exterior,
                            holes: geo_holes,
                        })
                    })
                    .collect::<Result<Vec<_>>>()?;

                // Re-orient after PROJ transform (may have flipped winding)
                let oriented = orient_polygons_geographic(&geo_polygons);

                Ok(GeoBandResult {
                    lower: band.lower,
                    upper: band.upper,
                    polygons: oriented,
                })
            })
            .collect::<Result<Vec<_>>>()?;

        Ok(geo_bands)
    }

    fn transform_ring(
        &self,
        ring: &[[f32; 2]],
        proj: &Proj,
    ) -> Result<Vec<[f64; 2]>> {
        ring.iter()
            .map(|[grid_x, grid_y]| {
                // Convert grid index to projected coordinates
                let proj_x = self.grid_info.x_origin + (*grid_x as f64) * self.grid_info.x_step;
                let proj_y = self.grid_info.y_origin + (*grid_y as f64) * self.grid_info.y_step;

                // Transform to lat/lon
                let (lon, lat) = proj.transform(proj_x, proj_y)?;

                Ok([lon, lat])
            })
            .collect()
    }
}
```

---

## Part 8: Deployment

### 8.1 Build Script

```bash
#!/bin/bash
# build-lambda.sh

set -e

# Build for Lambda (Amazon Linux 2)
cargo lambda build --release --output-format zip

# Output location
echo "Lambda package: target/lambda/zarr-contour/bootstrap.zip"
```

### 8.2 Terraform/CloudFormation Configuration

```hcl
# terraform/lambda.tf

resource "aws_lambda_function" "zarr_contour" {
  function_name = "zarr-contour"
  role          = aws_iam_role.zarr_contour_role.arn
  handler       = "bootstrap"
  runtime       = "provided.al2"

  filename         = "../target/lambda/zarr-contour/bootstrap.zip"
  source_code_hash = filebase64sha256("../target/lambda/zarr-contour/bootstrap.zip")

  memory_size = 1024  # 1 GB
  timeout     = 60    # 60 seconds (coordinator needs time for workers)

  environment {
    variables = {
      RUST_LOG = "info"
    }
  }
}

resource "aws_iam_role" "zarr_contour_role" {
  name = "zarr-contour-role"

  assume_role_policy = jsonencode({
    Version = "2012-10-17"
    Statement = [{
      Action = "sts:AssumeRole"
      Effect = "Allow"
      Principal = {
        Service = "lambda.amazonaws.com"
      }
    }]
  })
}

resource "aws_iam_role_policy" "zarr_contour_policy" {
  name = "zarr-contour-policy"
  role = aws_iam_role.zarr_contour_role.id

  policy = jsonencode({
    Version = "2012-10-17"
    Statement = [
      {
        Effect = "Allow"
        Action = [
          "s3:GetObject",
        ]
        Resource = "arn:aws:s3:::${var.zarr_bucket}/*"
      },
      {
        Effect = "Allow"
        Action = [
          "lambda:InvokeFunction",
        ]
        Resource = aws_lambda_function.zarr_contour.arn
      },
      {
        Effect = "Allow"
        Action = [
          "logs:CreateLogGroup",
          "logs:CreateLogStream",
          "logs:PutLogEvents",
        ]
        Resource = "arn:aws:logs:*:*:*"
      }
    ]
  })
}
```

### 8.3 Lambda Configuration Notes

| Setting | Value | Rationale |
|---------|-------|-----------|
| Memory | 1024 MB | Balance of cost and performance |
| Timeout | 60s | Coordinator waits for workers |
| Reserved Concurrency | None | Allow scaling |
| Provisioned Concurrency | Optional | For consistent cold start times |

---

## Part 9: Testing Plan

### 9.1 Unit Tests

```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_process_band_flat_simple() {
        // 3x3 grid with values:
        // 5  5  5
        // 5 15 15
        // 5 15 15
        let values: Vec<f32> = vec![
            5.0, 5.0, 5.0,
            5.0, 15.0, 15.0,
            5.0, 15.0, 15.0,
        ];

        let polygons = process_band_flat(&values, 3, 3, 10.0, 20.0, 5);

        assert!(!polygons.is_empty());
        assert!(!polygons[0].exterior.is_empty());
    }

    #[test]
    fn test_coordinator_request_parsing() {
        let json = r#"{
            "mode": "coordinator",
            "zarr_bucket": "my-bucket",
            "zarr_key": "data/grid",
            "width": 1799,
            "height": 1059,
            "x_origin": -2699020.0,
            "y_origin": -1588806.0,
            "x_step": 3000.0,
            "y_step": 3000.0,
            "thresholds": [273.15, 283.15, 293.15],
            "function_name": "zarr-contour"
        }"#;

        let request: ContourRequest = serde_json::from_str(json).unwrap();

        match request {
            ContourRequest::Coordinator(r) => {
                assert_eq!(r.thresholds.len(), 3);
            }
            _ => panic!("Expected coordinator request"),
        }
    }

    #[test]
    fn test_worker_request_parsing() {
        let json = r#"{
            "mode": "worker",
            "zarr_bucket": "my-bucket",
            "zarr_key": "data/grid",
            "width": 1799,
            "height": 1059,
            "x_origin": -2699020.0,
            "y_origin": -1588806.0,
            "x_step": 3000.0,
            "y_step": 3000.0,
            "lower": 273.15,
            "upper": 283.15,
            "precision": 5
        }"#;

        let request: ContourRequest = serde_json::from_str(json).unwrap();

        match request {
            ContourRequest::Worker(r) => {
                assert_eq!(r.lower, 273.15);
                assert_eq!(r.upper, 283.15);
            }
            _ => panic!("Expected worker request"),
        }
    }
}
```

### 9.2 Integration Tests

```rust
// tests/integration_tests.rs

#[tokio::test]
async fn test_worker_with_real_zarr() {
    // Requires AWS credentials and real Zarr data
    if std::env::var("RUN_INTEGRATION_TESTS").is_err() {
        return;
    }

    let request = WorkerRequest {
        zarr_bucket: "my-zarr-bucket".to_string(),
        zarr_key: "hrrr/tmp/2025/01/01/00/00".to_string(),
        width: 1799,
        height: 1059,
        x_origin: -2699020.0,
        y_origin: -1588806.0,
        x_step: 3000.0,
        y_step: 3000.0,
        lower: 273.15,
        upper: 283.15,
        precision: 5,
    };

    let response = handle_worker(request).await.unwrap();
    let worker_response: WorkerResponse = serde_json::from_slice(&response).unwrap();

    assert!(worker_response.stats.total_ms > 0.0);
    println!("Zarr fetch: {}ms", worker_response.stats.zarr_fetch_ms);
    println!("Contour: {}ms", worker_response.stats.contour_ms);
}
```

---

## Part 10: Implementation Checklist

### Phase 1: Project Setup
- [ ] Create `zarr-contour` project directory
- [ ] Initialize Cargo.toml with dependencies
- [ ] Set up project structure
- [ ] Add geo-marching-squares as path dependency

### Phase 2: Core Types
- [ ] Implement request/response types in `types.rs`
- [ ] Add serde serialization tests
- [ ] Verify gzip compression/decompression

### Phase 3: Algorithm Changes
- [ ] Add `process_band_flat()` to geo-marching-squares
- [ ] Add `Shape::create_from_corners()` constructor
- [ ] Add f32 coordinate output support
- [ ] Test with simple grids

### Phase 4: Zarr Integration
- [ ] Copy Zarr fetching code from grib-inspector
- [ ] Adapt for standalone use
- [ ] Test with real S3 Zarr data
- [ ] Handle blosc decompression

### Phase 5: Worker Implementation
- [ ] Implement `handle_worker()`
- [ ] Test locally with mock data
- [ ] Test with real Zarr data
- [ ] Verify timing statistics

### Phase 6: Coordinator Implementation
- [ ] Implement `handle_coordinator()`
- [ ] Implement parallel worker invocation
- [ ] Handle worker failures
- [ ] Test with multiple bands

### Phase 7: Lambda Handler
- [ ] Implement main handler with mode dispatch
- [ ] Add gzip response compression
- [ ] Test locally with SAM/cargo-lambda
- [ ] Verify payload sizes

### Phase 8: Deployment
- [ ] Create build script
- [ ] Create Terraform/CloudFormation config
- [ ] Deploy to dev environment
- [ ] Test end-to-end

### Phase 9: grib-inspector Integration
- [ ] Add Lambda client to contour_service
- [ ] Implement `generate_isobands_via_lambda()`
- [ ] Add PROJ transformation after Lambda response
- [ ] Test full flow
- [ ] Compare results with current implementation

### Phase 10: Validation
- [ ] Compare output with current contour-rs output
- [ ] Verify winding order after PROJ transform
- [ ] Performance benchmarks
- [ ] Load testing

---

## Part 11: Success Criteria

### Functional Requirements
- [ ] Coordinator correctly fans out to workers
- [ ] Workers generate valid grid-space polygons
- [ ] Polygon winding order is correct (RFC 7946)
- [ ] gzip compression reduces response size
- [ ] PROJ transformation in grib-inspector produces correct lat/lon
- [ ] Output matches current contour-rs output (within precision)

### Performance Requirements
- [ ] Single band: < 1,000ms (currently ~580ms + 440ms Zarr)
- [ ] 13 bands: < 1,500ms (currently ~8,000ms)
- [ ] Speedup: > 5x compared to current sequential
- [ ] Cold start: < 2,000ms additional

### Quality Requirements
- [ ] All tests passing
- [ ] No clippy warnings
- [ ] Comprehensive error handling
- [ ] Structured logging with tracing

---

## Appendix A: Winding Order Considerations

### Grid Space vs Geographic Space

```
Grid Space (Lambda output):
- Row 0 = top of grid (lowest Y in screen coords, but could be highest lat)
- Y increases downward (row index increases)
- CCW in grid space = positive signed area

Geographic Space (after PROJ):
- Latitude increases upward
- PROJ may flip Y axis depending on projection
- Must re-verify winding after transformation
```

### Safe Winding Handling

```rust
// In grib-inspector, after PROJ transform:
fn ensure_correct_winding(polygons: &mut [GeoPolygon]) {
    for poly in polygons.iter_mut() {
        // Exterior should be CCW (positive signed area in lat/lon)
        if signed_area_geo(&poly.exterior) < 0.0 {
            poly.exterior.reverse();
        }

        // Holes should be CW (negative signed area in lat/lon)
        for hole in poly.holes.iter_mut() {
            if signed_area_geo(hole) > 0.0 {
                hole.reverse();
            }
        }
    }
}
```

---

## Appendix B: Payload Size Estimates

### Worker Response Size

```
Per polygon:
- Average ring: 500 points × 2 × 4 bytes = 4 KB
- Average polygon: 2 rings × 4 KB = 8 KB
- Per band: ~10 polygons × 8 KB = 80 KB

Coordinator response (13 bands):
- 13 × 80 KB = 1 MB uncompressed
- With gzip (~5x): ~200 KB compressed

Lambda limits:
- Sync response: 6 MB ✓
- With compression: easily within limits
```

### Worst Case

```
Complex weather pattern:
- 50 polygons per band
- 2000 points per ring
- 3 rings per polygon

Per band: 50 × 3 × 2000 × 8 = 2.4 MB
13 bands: 31 MB uncompressed

With gzip: ~6 MB (borderline)

Mitigation: Simplify polygons if response > 4 MB
```

---

**End of Implementation Plan**
