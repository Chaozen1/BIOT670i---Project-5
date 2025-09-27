
# Capstone Microbiome Results (Kraken2 + Bracken)

## Database
- Built with Kraken2 using **SILVA 138.1 16S rRNA reference**.

## Pipeline
1. Classified reads with Kraken2.
2. Re-estimated abundances with Bracken.
3. Compiled genus-level matrices across all samples.
4. Visualized top genera with barplots & heatmaps.

## Contents
- `bracken_outputs/` → per-sample Bracken abundance tables.
- `matrices/` → merged count tables (`.tsv`).
- `plots/` → visualization of top 10 genera.
- `README.md` → this file.

## Notes
- Bracken outputs are relative to the SILVA 16S database.
- Results can be re-generated using Kraken2 + Bracken with the same database.
