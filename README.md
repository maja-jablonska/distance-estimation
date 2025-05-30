# Distance Estimation CLI Tool

A command-line tool for estimating stellar distances from parallax measurements using EDSD, GGD, and Photogeometric models via MCMC sampling.

Rewritten [this repository](https://github.com/ElisaHaas25/Interactive-Distance-Estimation) using Claude Opus 4.

## Features

- **Three distance models**: EDSD, GGD, and Photogeometric
- **MCMC sampling**: Robust Bayesian distance estimation
- **Batch processing**: Process multiple sources from CSV files
- **Consolidated output**: All posterior samples saved to single file
- **Extensive logging**: Detailed progress tracking and error reporting
- **No plotting dependencies**: Streamlined for batch processing

## Installation

1. Install Python dependencies:

```bash
pip install -r requirements.txt
```

2. Ensure you have the original `packages.py` and `functions.py` files in the same directory.

## Usage

### Basic Usage

```bash
python distance_estimation_cli.py -i input.csv -o output -m EDSD
```

### Command Line Options

- `--input-file, -i`: Input CSV file with source data (required)
- `--output-file, -o`: Output file prefix (required)
- `--model, -m`: Distance model - EDSD, GGD, or Photogeometric (required)
- `--seed, -s`: Random seed for reproducibility (default: 42)
- `--nsamp`: Number of MCMC samples (default: 5000)
- `--nburnin`: Number of MCMC burn-in samples (default: 500)
- `--log-level`: Logging level - DEBUG, INFO, WARNING, ERROR (default: INFO)
- `--validate-only`: Only validate input file, don't process

### Input File Format

#### EDSD Model

Required columns:

- `source_id`: Unique identifier
- `parallax`: Parallax in milliarcseconds
- `parallax_error`: Parallax uncertainty in milliarcseconds
- `rlen`: Length scale parameter in parsecs

Example:

```csv
source_id,parallax,parallax_error,rlen
1234567890123456789,5.234,0.123,1250.0
2345678901234567890,12.456,0.234,1100.0
```

#### GGD Model

Additional required columns:

- `alpha`: Shape parameter
- `beta`: Scale parameter

Example:

```csv
source_id,parallax,parallax_error,rlen,alpha,beta
1234567890123456789,5.234,0.123,1250.0,1.2,2.8
2345678901234567890,12.456,0.234,1100.0,1.1,2.9
```

#### Photogeometric Model

Additional required columns:

- `hp`: HEALPix pixel number
- `phot_g_mean_mag`: Gaia G magnitude
- `bp_rp`: Gaia BP-RP color
- `pseudocolour`: Gaia pseudocolour

### Output Files

Two output files are created:

1. **`{output}_samples.csv`**: All posterior samples from all sources

   - Columns: `source_id`, `sample`
   - Each row is one posterior sample
2. **`{output}_summary.csv`**: Summary statistics for each source

   - Columns: `source_id`, `status`, `message`, `r_median`, `r_lo`, `r_hi`, `n_samples`
   - One row per input source

### Examples

#### Process EDSD model with custom MCMC parameters:

```bash
python distance_estimation_cli.py \
    -i data/my_sources.csv \
    -o results/edsd_results \
    -m EDSD \
    --nsamp 10000 \
    --nburnin 1000 \
    --seed 123
```

#### Validate input file only:

```bash
python distance_estimation_cli.py \
    -i data/my_sources.csv \
    -o dummy \
    -m GGD \
    --validate-only
```

#### Enable debug logging:

```bash
python distance_estimation_cli.py \
    -i data/my_sources.csv \
    -o results/debug_run \
    -m Photogeometric \
    --log-level DEBUG
```

## Logging

The tool provides extensive logging:

- **File logging**: All messages saved to `distance_estimation.log`
- **Console logging**: Progress and important messages displayed
- **Progress tracking**: Regular updates on processing status
- **Error handling**: Detailed error messages and stack traces

Log levels:

- `DEBUG`: Detailed per-source processing information
- `INFO`: Progress updates and summary statistics (default)
- `WARNING`: Non-fatal issues
- `ERROR`: Fatal errors and failures

## Performance

- **Processing rate**: Typically 1-10 sources per second depending on model complexity
- **Memory usage**: Scales linearly with number of sources
- **Progress tracking**: Regular updates every 100 sources

## Error Handling

The tool handles various error conditions gracefully:

- Invalid parallax values (infinite, negative uncertainties)
- MCMC convergence failures
- Missing or malformed input data
- R interface errors (for Photogeometric model)

Failed sources are logged but processing continues for remaining sources.

## Dependencies

- Python 3.7+
- NumPy, SciPy, Pandas
- Astropy, AstroQuery
- HEALPy
- Click (command-line interface)
- rpy2 (for Photogeometric model)

See `requirements.txt` for specific versions.

## Troubleshooting

### Common Issues

1. **R interface errors**: Ensure R and required R packages are installed for Photogeometric model
2. **Memory issues**: Reduce `--nsamp` for large datasets
3. **Slow processing**: Check MCMC convergence, reduce sample count if needed
4. **Input validation failures**: Check column names and data types match requirements

### Getting Help

Use the `--help` flag for quick reference:

```bash
python distance_estimation_cli.py --help
```

Enable debug logging to diagnose issues:

```bash
python distance_estimation_cli.py --log-level DEBUG ...
```
