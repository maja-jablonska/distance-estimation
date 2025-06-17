#!/usr/bin/env python3
"""
Utility script to merge chunk results from distance estimation processing.
"""

import click
import pandas as pd
import numpy as np
from pathlib import Path
import logging
import glob

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def merge_chunk_files(output_prefix: str, cleanup: bool = False) -> None:
    """
    Merge chunk files into consolidated results.
    
    Args:
        output_prefix: Output file prefix used in original processing
        cleanup: Whether to remove chunk files after merging
    """
    output_path = Path(output_prefix)
    output_dir = output_path.parent
    
    # Find all chunk files
    chunk_samples_pattern = str(output_dir / f"{output_path.stem}_chunk_*_samples.parquet")
    chunk_summary_pattern = str(output_dir / f"{output_path.stem}_chunk_*_summary.csv")
    
    chunk_samples_files = sorted(glob.glob(chunk_samples_pattern))
    chunk_summary_files = sorted(glob.glob(chunk_summary_pattern))
    
    logger.info(f"Found {len(chunk_samples_files)} sample chunk files")
    logger.info(f"Found {len(chunk_summary_files)} summary chunk files")
    
    # Merge samples
    if chunk_samples_files:
        logger.info("Merging sample chunks...")
        all_samples = {}
        
        for chunk_file in chunk_samples_files:
            logger.info(f"Processing {chunk_file}")
            chunk_df = pd.read_parquet(chunk_file)
            
            for _, row in chunk_df.iterrows():
                source_id = row['source_id']
                samples = row['samples']
                all_samples[source_id] = samples
        
        # Save consolidated samples
        samples_file = output_dir / f"{output_path.stem}_samples.parquet"
        samples_df = pd.DataFrame([
            {'source_id': source_id, 'samples': samples}
            for source_id, samples in all_samples.items()
        ])
        samples_df.to_parquet(samples_file, index=False)
        logger.info(f"Saved consolidated samples to {samples_file}")
        
        total_samples = sum(len(samples) for samples in all_samples.values())
        logger.info(f"Total samples: {total_samples} from {len(all_samples)} sources")
    
    # Merge summaries
    if chunk_summary_files:
        logger.info("Merging summary chunks...")
        all_summaries = []
        
        for chunk_file in chunk_summary_files:
            logger.info(f"Processing {chunk_file}")
            chunk_df = pd.read_csv(chunk_file)
            all_summaries.append(chunk_df)
        
        # Save consolidated summary
        summary_file = output_dir / f"{output_path.stem}_summary.csv"
        summary_df = pd.concat(all_summaries, ignore_index=True)
        summary_df.to_csv(summary_file, index=False)
        logger.info(f"Saved consolidated summary to {summary_file}")
        logger.info(f"Total sources: {len(summary_df)}")
    
    # Cleanup chunk files if requested
    if cleanup:
        logger.info("Cleaning up chunk files...")
        for chunk_file in chunk_samples_files + chunk_summary_files:
            Path(chunk_file).unlink()
            logger.debug(f"Removed {chunk_file}")
        logger.info(f"Removed {len(chunk_samples_files + chunk_summary_files)} chunk files")


@click.command()
@click.argument('output_prefix', type=str)
@click.option('--cleanup', is_flag=True, 
              help='Remove chunk files after merging')
def main(output_prefix: str, cleanup: bool):
    """
    Merge chunk results from distance estimation processing.
    
    OUTPUT_PREFIX should be the same prefix used in the original distance estimation run.
    
    This will look for files matching:
    - {output_prefix}_chunk_XXXX_samples.parquet
    - {output_prefix}_chunk_XXXX_summary.csv
    
    And create consolidated files:
    - {output_prefix}_samples.parquet
    - {output_prefix}_summary.csv
    
    Example:
        python merge_chunks.py results --cleanup
    """
    try:
        merge_chunk_files(output_prefix, cleanup)
        logger.info("Chunk merging completed successfully")
    except Exception as e:
        logger.error(f"Error merging chunks: {e}")
        raise


if __name__ == '__main__':
    main() 