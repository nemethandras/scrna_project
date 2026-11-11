"""Default locations of the lab-wide shared reference data and genotype database.

Override with SCRNA_SHARED_ROOT (whole shared directory) or SCRNA_DB (database only).
"""
import os

SHARED_ROOT   = os.environ.get("SCRNA_SHARED_ROOT", "/groups/nils/data/scrna_pipeline")
REFERENCE_DIR = os.path.join(SHARED_ROOT, "reference", "hg38")
DEFAULT_DB    = os.environ.get("SCRNA_DB", os.path.join(SHARED_ROOT, "db", "variants.db"))

REFERENCE_GENOME = os.path.join(REFERENCE_DIR, "genome.fa")
ANNOTATION_GTF   = os.path.join(REFERENCE_DIR, "genes.gtf")
STAR_INDEX       = os.path.join(REFERENCE_DIR, "star_index_hg38_oh99")
BWA_INDEX        = os.path.join(REFERENCE_DIR, "bwa_index_hg38")
SNPS_VCF         = os.path.join(REFERENCE_DIR, "genome1K.hg38.common_snps.vcf.gz")
SNPS_SOURCE      = os.path.join(REFERENCE_DIR, "genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf")
