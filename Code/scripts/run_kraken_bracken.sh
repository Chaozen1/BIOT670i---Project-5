#!/bin/bash
set -euo pipefail

# ====== CONFIG ======
WORKDIR="$HOME/capstone"
DB="$WORKDIR/db/k2_standard_8_20250714"   # change if your DB path is different
THREADS=4                                  # bump if your Mac can handle it

# ====== PREP ======
mkdir -p "$WORKDIR/reports" "$WORKDIR/class" "$WORKDIR/bracken"

# choose the right decompressor (mac uses gzcat)
if command -v gzcat >/dev/null 2>&1; then
  ZCAT=gzcat
else
  ZCAT=zcat
fi

# ====== 1) KRAKEN2 on every FASTQ (single-end) ======
for fq in "$WORKDIR"/fastq/*.fastq.gz; do
  base="$(basename "$fq" .fastq.gz)"
  echo ">>> Kraken2: $base"
  kraken2 \
    --db "$DB" \
    --threads "$THREADS" \
    --memory-mapping \
    --gzip-compressed \
    --use-names \
    --report "$WORKDIR/reports/${base}.kreport" \
    --output "$WORKDIR/class/${base}.kraken2" \
    "$fq"
done

# ====== 2) READ LENGTH ONCE (for Bracken) ======
FIRST="$(ls -1 "$WORKDIR"/fastq/*.fastq.gz | head -n 1)"
READLEN="$($ZCAT "$FIRST" | awk 'NR%4==2{print length($0); exit}')"
echo ">>> Estimated read length = $READLEN"

# ====== 3) BRACKEN on each report (species level) ======
for rep in "$WORKDIR"/reports/*.kreport; do
  base="$(basename "$rep" .kreport)"
  echo ">>> Bracken: $base"
  bracken \
    -d "$DB" \
    -i "$rep" \
    -o "$WORKDIR/bracken/${base}.bracken" \
    -r "$READLEN" \
    -l S \
    -t 10
done

# ====== 4) COMBINE outputs ======
echo ">>> Combine Kraken2 reports"
python3 -m KrakenTools.combine_kreports \
  -r "$WORKDIR"/reports/*.kreport \
  -o "$WORKDIR/reports/kraken_combined.kreport"

echo ">>> Combine Bracken outputs"
python3 -m KrakenTools.combine_bracken_outputs \
  -i "$WORKDIR"/bracken/*.bracken \
  -o "$WORKDIR/bracken/bracken_combined.txt"

echo ">>> DONE. Results in:"
echo "    $WORKDIR/reports  (Kraken2 + combined)"
echo "    $WORKDIR/bracken  (Bracken + combined)"
