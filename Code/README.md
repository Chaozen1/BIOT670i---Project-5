# Code Directory

This folder contains helper scripts for the BIOT670i Project 5 pipeline.

## Running Kraken2 + Bracken

We provide a helper script at [`Code/scripts/run_kraken_bracken.sh`](scripts/run_kraken_bracken.sh) to process FASTQ files through **Kraken2** and **Bracken**.

### Requirements
- Kraken2 and Bracken installed and available in your `$PATH`
- Access to the Kraken2 database (e.g., SILVA, RefSeq, or standard)
- Input FASTQ file(s)

### Usage
From the repo root:

```bash
./Code/scripts/run_kraken_bracken.sh input.fastq output_prefix

### Inputs
- **FASTQ file** (can be uncompressed `.fastq` or gzipped `.fastq.gz`)
  - Example: `Data/fastq/SRR8594636.fastq.gz`
- **Output prefix** (basename for output files)
  - Example: `SRR8594636`

---

### Outputs
For each run, the script produces:

- `Results/Bracken/<prefix>.kraken` → raw Kraken2 classification
- `Results/Bracken/<prefix>.kreport` → Kraken2 report (hierarchical breakdown)
- `Results/Bracken/<prefix>.bracken.genus` → Bracken re-estimated abundances at the genus level  
  *(you can adapt the script for species, family, phylum if needed)*

Example command:
```bash
./Code/scripts/run_kraken_bracken.sh Data/fastq/SRR8594636.fastq.gz SRR8594636
```bash
./Code/scripts/run_kraken_bracken.sh Data/fastq/SRR8594636.fastq.gz SRR859463

