#!/usr/bin/env bash
# Create (or refresh) the lab-wide shared directory holding the reference
# data and the genotype database:
#
#   $SCRNA_SHARED_ROOT/reference/hg38/   read-only for everyone
#   $SCRNA_SHARED_ROOT/db/variants.db    group-writable, guarded by variants.db.lock
#   $SCRNA_SHARED_ROOT/db/backups/
#
# Usage (from the repo root):
#   bash scripts/setup_shared.sh [SRC_REFERENCE_DIR] [SRC_SNP_PANEL_VCF] [SRC_DB]
#
# An existing shared database is never overwritten.
set -euo pipefail

SHARED_ROOT="${SCRNA_SHARED_ROOT:-/groups/nils/data/scrna_pipeline}"
SRC_REF="${1:-data/reference}"
SRC_SNP_SOURCE="${2:-validation/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf}"
SRC_DB="${3:-results/variants.db}"

REF="$SHARED_ROOT/reference/hg38"
DB_DIR="$SHARED_ROOT/db"

mkdir -p "$REF" "$DB_DIR/backups"
chmod 2775 "$SHARED_ROOT" "$SHARED_ROOT/reference" "$DB_DIR" "$DB_DIR/backups"

echo "Copying reference data to $REF ..."
chmod -R u+w "$REF"
rsync -a --info=progress2 \
    "$SRC_REF/genome.fa" \
    "$SRC_REF/genome.fa.fai" \
    "$SRC_REF/genes.gtf" \
    "$SRC_REF/genome1K.hg38.common_snps.vcf.gz" \
    "$SRC_REF/genome1K.hg38.common_snps.vcf.gz.tbi" \
    "$SRC_REF/star_index_hg38_oh99" \
    "$SRC_REF/bwa_index_hg38" \
    "$SRC_SNP_SOURCE" \
    "$REF/"
# Read-only, so a pipeline run can never rebuild or overwrite a shared index.
chmod -R a-w,a+rX "$REF"

if [ -e "$DB_DIR/variants.db" ]; then
    echo "Shared database already exists — leaving it untouched."
else
    echo "Copying database $SRC_DB -> $DB_DIR/variants.db ..."
    cp -p "$SRC_DB" "$DB_DIR/variants.db"
fi
touch "$DB_DIR/variants.db.lock"
chmod 664 "$DB_DIR/variants.db" "$DB_DIR/variants.db.lock"

echo "Done."
du -sh "$REF" "$DB_DIR/variants.db"
