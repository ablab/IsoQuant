# Barcode Detection for Single-Cell and Spatial Transcriptomics

Platform-specific barcode calling using k-mer indexing and sequence alignment. Detects barcodes, UMIs, and structural features (primers, linkers, polyT tails) in long reads.

## Supported Platforms

### Stereo-seq
- **Standard mode**: `StereoBarcodeDetector`
- **Splitting mode**: `StereoSplttingBarcodeDetector` (detects multiple barcodes per read)
- **Barcode**: 25bp
- **UMI**: 10bp
- **Features**: PC1 primer, linker, polyT tail, TSO (splitting mode only)

### 10x Genomics
- **v3**: `TenXBarcodeDetector`
- **Visium HD**: `VisiumHDBarcodeDetector` (dual barcode)
- **Barcode**: 16bp (v3) or 16bp+15bp (Visium HD)
- **UMI**: 12bp (v3), 9bp (Visium HD)
- **Features**: R1 primer, polyT tail

### Curio
- **Standard**: `DoubleBarcodeDetector`
- **Illumina-compatible**: `IlluminaDoubleBarcodeDetector`
- **Brute force**: `BruteForceDoubleBarcodeDetector` (exact matching only)
- **Barcode**: 8bp + 6bp (split by linker)
- **UMI**: 9bp
- **Features**: PCR primer, linker, polyT tail

## Detection Algorithm

1. **PolyT Detection**: Find polyT/polyA tail using sliding window
2. **Anchor Search**: Locate structural features (primer, linker, R1) using k-mer indexing
3. **Barcode Extraction**: Extract sequence between structural features
4. **Barcode Matching**: K-mer index lookup + Smith-Waterman alignment against whitelist
5. **UMI Extraction**: Extract UMI sequence between barcode and polyT
6. **Strand Detection**: Check forward and reverse complement, keep best match

## Result Classes

### `BarcodeDetectionResult`
Base class storing:
- Read ID
- Detected barcode and UMI sequences
- Barcode alignment score
- UMI validity flag
- Detected strand

### `DoubleBarcodeDetectionResult`
Extends base with positions of:
- PolyT tail start
- Primer end
- Linker start/end

### `StereoBarcodeDetectionResult`
Extends double barcode with:
- TSO 5' position

### `TenXBarcodeDetectionResult`
Extends base with positions of:
- PolyT tail start
- R1 primer end

### `SplittingBarcodeDetectionResult`
Container for multiple detection patterns per read (e.g., stereo splitting mode)

## Statistics Tracking

**`ReadStats`** accumulates:
- Total reads processed
- Barcodes detected
- Valid UMIs
- Feature detection counts (primer, linker, polyT, TSO, R1)

**Usage**:
```python
stats = ReadStats()
for result in detection_results:
    stats.add_read(result)
print(stats)  # Human-readable summary
```

## Memory Optimization

For large barcode whitelists (e.g., spatial transcriptomics):

### `SharedMemoryStereoBarcodeDetector`
- Uses shared memory for barcode index across parallel workers
- Efficient for whitelists > 1M barcodes
- Reduces memory footprint in multiprocessing

### `SharedMemoryWrapper`
- Generic wrapper for any detector class
- Serializes only shared memory references (not full index)
- Enables efficient process forking

## Usage Example

```python
from src.barcode_calling.barcode_callers import StereoBarcodeDetector, ReadStats

# Load barcode whitelist
with open("whitelist.txt") as f:
    barcodes = [line.strip() for line in f]

# Initialize detector
detector = StereoBarcodeDetector(barcodes, min_score=21)

# Process reads
stats = ReadStats()
results = []

for read_id, sequence in reads:
    result = detector.find_barcode_umi(read_id, sequence)
    results.append(result)
    stats.add_read(result)

# Print statistics
print(stats)
```

## Performance Considerations

- **K-mer size**: Larger k = fewer false positives, but less tolerance for errors
  - Stereo-seq uses k=14 (25bp barcodes)
  - 10x uses k=6-7 (16bp barcodes)
  - Curio uses k=5-6 (14bp total barcode)

- **Min score**: Minimum alignment score for barcode acceptance
  - Higher scores = stricter matching
  - Adjust based on sequencing error rate

- **Shared memory**: Use for whitelists > 1M barcodes to reduce memory usage in parallel processing

## Detector Selection

```python
from src.barcode_calling.barcode_callers import (
    StereoBarcodeDetector,          # Standard stereo-seq
    StereoSplttingBarcodeDetector,  # Stereo-seq with read splitting
    TenXBarcodeDetector,            # 10x v3
    VisiumHDBarcodeDetector,        # Visium HD (dual barcode)
    DoubleBarcodeDetector,          # Curio
)

# Initialize based on platform
detector = TenXBarcodeDetector(barcode_whitelist)

# Detect barcodes
result = detector.find_barcode_umi(read_id, sequence)

# Check result
if result.is_valid():
    print(f"Barcode: {result.barcode}, UMI: {result.UMI}")
```

## Testing

See `tests/test_barcode_callers.py` for comprehensive unit tests covering:
- Result class initialization and methods
- Coordinate updating and shifting
- Result comparison logic
- Statistics accumulation
- Feature detection
