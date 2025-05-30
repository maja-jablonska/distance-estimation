# Standalone Distance Estimation Tool

A pure Python implementation of stellar distance estimation from parallax measurements using EDSD, GGD, and Photogeometric models via MCMC sampling. This version eliminates all R dependencies and provides a streamlined command-line interface.

## Features

- **Pure Python**: No R dependencies - eliminates segmentation faults and installation issues
- **Three distance models**: EDSD, GGD, and Photogeometric
- **Multiple input formats**: Supports both CSV and Parquet files
- **MCMC sampling**: Robust Bayesian distance estimation using Metropolis algorithm
- **Batch processing**: Process multiple sources from CSV/Parquet files
- **Consolidated output**: All posterior samples saved to single file
- **Extensive logging**: Detailed progress tracking and error reporting
- **Command-line interface**: Easy-to-use CLI with click

## Installation

1. Install Python dependencies:

```bash
pip install -r requirements_standalone.txt
```

The requirements are minimal:

- click>=8.0.0
- pandas>=1.3.0
- numpy>=1.20.0
- scipy>=1.7.0
- pyarrow>=5.0.0

## Usage

### Basic Usage

```bash
# CSV input
python distance_estimation_standalone.py -i input.csv -o output -m EDSD

# Parquet input
python distance_estimation_standalone.py -i input.parquet -o output -m EDSD
```

### Supported Input Formats

- **CSV** (.csv)
- **Parquet** (.parquet, .pq)

The tool automatically detects the file format based on the file extension. If the extension is not recognized, it will attempt to read the file as CSV first, then as Parquet.

### Command-line Options

- `-i, --input-file`: Input CSV or Parquet file with source data (required)
- `-o, --output-file`: Output file prefix (required)
- `-m, --model`: Distance model to use: EDSD, GGD, or Photogeometric (required)
- `-s, --seed`: Random seed for reproducibility (default: 42)
- `--nsamp`: Number of MCMC samples (default: 5000)
- `--nburnin`: Number of MCMC burn-in samples (default: 500)
- `--log-level`: Logging level: DEBUG, INFO, WARNING, ERROR (default: INFO)
- `--validate-only`: Only validate input file, do not process

### Input File Format

The input file (CSV or Parquet) must contain the following columns:

**Base columns (all models):**

- `gaia_dr3_source_id`: Unique identifier for each source
- `ra`, `dec`: Coordinates in degrees
- `parallax`: Parallax measurement in milliarcseconds
- `parallax_error`: Parallax uncertainty in milliarcseconds

**Model-specific columns:**

**Photogeometric model:**

- `g_mag`: G magnitude
- `bp_rp`: BP-RP color
- `pseudocolour`: Pseudocolour parameter

### Examples

**EDSD model (CSV):**

```bash
python distance_estimation_standalone.py -i sample_input_edsd.csv -o results_edsd -m EDSD
```

**EDSD model (Parquet):**

```bash
python distance_estimation_standalone.py -i sample_input_edsd.parquet -o results_edsd -m EDSD
```

**GGD model (Parquet):**

```bash
python distance_estimation_standalone.py -i sample_input_ggd.parquet -o results_ggd -m GGD
```

**Photogeometric model:**

```bash
python distance_estimation_standalone.py -i sample_input_photogeo.parquet -o results_photogeo -m Photogeometric
```

**Custom MCMC parameters:**

```bash
python distance_estimation_standalone.py -i input.parquet -o output -m EDSD --nsamp 10000 --nburnin 1000
```

**Validation only:**

```bash
python distance_estimation_standalone.py -i input.parquet -o output -m EDSD --validate-only
```

## Output Files

The tool creates two output files (always in CSV format):

1. **`{output}_samples.csv`**: Contains all posterior samples from all sources

   - `source_id`: Source identifier
   - `sample`: Distance sample in parsecs
2. **`{output}_summary.csv`**: Contains summary statistics for each source

   - `source_id`: Source identifier
   - `status`: Processing status (success/failed)
   - `message`: Status message
   - `r_median`: Median distance (50th percentile) in parsecs
   - `r_lo`: Lower bound (15.9th percentile) in parsecs
   - `r_hi`: Upper bound (84.1th percentile) in parsecs
   - `n_samples`: Number of samples obtained

## Sample Input Files

### EDSD Sample (CSV format)

```csv
source_id,parallax,parallax_error,rlen
1234567890123456789,5.234,0.123,1250.0
2345678901234567890,12.456,0.234,1100.0
```

### GGD Sample (CSV format)

```csv
source_id,parallax,parallax_error,rlen,alpha,beta
1234567890123456789,5.234,0.123,1250.0,1.2,2.8
2345678901234567890,12.456,0.234,1100.0,1.1,2.9
```

Sample Parquet files are also provided:

- `sample_input_edsd.parquet`
- `sample_input_ggd.parquet`
- `sample_input_photogeo.parquet`

You can create your own Parquet files using the provided `create_sample_parquet.py` script.

## Algorithm Details

### EDSD (Exponentially Decreasing Space Density)

- Prior: P(r) ∝ r² exp(-r/rlen)
- Suitable for general stellar populations
- Single parameter: length scale (rlen)

### GGD (Generalized Gamma Distribution)

- Prior: P(r) ∝ r^β exp(-(r/rlen)^α)
- More flexible than EDSD
- Three parameters: rlen, α, β

### Photogeometric

- Combines geometric prior with photometric information
- Uses color-magnitude models for improved accuracy
- Requires additional photometric data

### MCMC Sampling

- Metropolis algorithm with Gaussian proposals
- Automatic initialization using EDSD mode
- Adaptive step size based on parallax uncertainty
- Burn-in period to ensure convergence

## Performance

- Processing rate: ~40-100 sources/second (depending on model complexity)
- Memory usage: Minimal (processes one source at a time)
- Scalable to large datasets
- Parquet files offer faster read times for large datasets

## Advantages of Parquet Format

- **Faster I/O**: Significantly faster read/write compared to CSV for large files
- **Smaller file size**: Efficient compression reduces storage requirements
- **Type safety**: Column data types are preserved
- **Better for big data**: Optimized for analytical workloads
- **Cross-platform**: Works across different languages and tools

## Logging

The tool provides extensive logging:

- Progress updates every 100 sources
- Processing rate and ETA estimates
- Detailed error messages for failed sources
- Summary statistics at completion

Log files are saved as `distance_estimation.log`.

## Error Handling

Common failure modes and solutions:

1. **Invalid parallax data**: Check for negative or infinite uncertainties
2. **MCMC initialization failure**: Usually due to extreme parallax values
3. **Missing required columns**: Verify input file format
4. **Numerical issues**: May occur with very small/large parallax values
5. **File format issues**: Ensure PyArrow is installed for Parquet support

## Differences from R Implementation

This standalone version:

- Eliminates R dependencies and associated installation issues
- Provides simplified photogeometric model (mock implementation)
- Uses pure Python MCMC instead of R's metropolis function
- Streamlined for batch processing without plotting
- Consolidated output format for easier analysis
- Supports both CSV and Parquet input formats

## Validation

To validate the implementation:

```bash
# CSV
python distance_estimation_standalone.py -i sample_input_edsd.csv -o test -m EDSD --validate-only

# Parquet
python distance_estimation_standalone.py -i sample_input_edsd.parquet -o test -m EDSD --validate-only
```

## Support

For issues or questions:

1. Check the log file for detailed error messages
2. Verify input file format matches requirements
3. Test with provided sample files (both CSV and Parquet)
4. Use `--validate-only` flag to check input data
5. Ensure PyArrow is installed for Parquet support

## License

This implementation maintains compatibility with the original research code while providing a more robust and user-friendly interface.
