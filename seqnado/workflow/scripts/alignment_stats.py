import pandas as pd
from pathlib import Path
from loguru import logger

logger.add(snakemake.log[0], level="INFO")

processed_data = {}

for f in snakemake.input:
    sample = Path(f).stem.replace('alignment_stats_', '')
    try:
        df = pd.read_csv(f, sep='\t')

        # Drop unchanged steps
        df = df.loc[df['Reads After'].shift() != df['Reads After']].reset_index(drop=True)

        steps = df['Step'].tolist()
        counts = df['Reads After'].tolist()
        deltas = []

        for i in range(len(counts) - 1):
            step_name = f"Lost at {steps[i+1]}"
            delta = counts[i] - counts[i+1]
            deltas.append((step_name, delta))

        # Final retained reads
        deltas.append(("Retained", counts[-1]))

        # Store in dict
        processed_data[sample] = dict(deltas)

    except Exception as e:
        logger.error(f"Error reading {f}: {e}")
        continue

# Convert to DataFrame
if processed_data:
    stacked_df = pd.DataFrame(processed_data).T.fillna(0).reset_index()
    stacked_df.rename(columns={'index': 'sample'}, inplace=True)
    stacked_df.to_csv(snakemake.output[0], sep='\t', index=False)
    logger.info(f"Alignment data written to {snakemake.output[0]}")
else:
    logger.warning("No valid input files found; no output written.")
