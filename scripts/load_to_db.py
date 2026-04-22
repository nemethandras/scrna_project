import os
import sqlite3
from grist_api import GristDocAPI
from cyvcf2 import VCF

GRIST_SERVER = os.environ["GRIST_SERVER"]
GRIST_API_KEY = os.environ["GRIST_API_KEY"]
GRIST_DOC_ID  = os.environ["GRIST_DOC_ID"]

grist = GristDocAPI(GRIST_DOC_ID, server=GRIST_SERVER, api_key=GRIST_API_KEY)

def load_run_to_grist(sample_name, vcf_path, flagstat_path, run_metadata):
    """Push run summary to Grist — metadata only, not raw variants."""

    # Count variants
    vcf = VCF(vcf_path)
    variant_count = sum(1 for _ in vcf)

    # Parse mapping rate from flagstat
    with open(flagstat_path) as f:
        for line in f:
            if "mapped (" in line:
                mapping_rate = float(line.split("(")[1].split("%")[0])
                break

    # Push to Grist Runs table
    grist.add_records("Runs", [{
        "sample_name":       sample_name,
        "pipeline_version":  run_metadata["pipeline_version"],
        "reference_genome":  run_metadata["reference_genome"],
        "run_date":          run_metadata["run_date"],
        "total_variants":    variant_count,
        "mapping_rate":      mapping_rate,
        "vcf_file_path":     vcf_path,
    }])

    print(f"Pushed run summary for {sample_name} to Grist")


def load_variants_to_sqlite(vcf_path, db_path="variants.db"):
    """Push full variant data to SQLite — for proper querying at scale."""

    conn = sqlite3.connect(db_path)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS variants (
            chromosome TEXT,
            position   INTEGER,
            ref_allele TEXT,
            alt_allele TEXT,
            quality    REAL,
            depth      INTEGER,
            genotype   TEXT,
            sample     TEXT
        )
    """)

    vcf = VCF(vcf_path)
    sample_name = vcf_path.split("/")[-1].replace(".filtered.vcf", "")

    rows = []
    for variant in vcf:
        rows.append((
            variant.CHROM,
            variant.POS,
            variant.REF,
            str(variant.ALT[0]),
            variant.QUAL,
            variant.INFO.get("DP"),
            variant.genotypes[0],
            sample_name
        ))

    conn.executemany(
        "INSERT INTO variants VALUES (?,?,?,?,?,?,?,?)", rows
    )
    conn.commit()
    conn.close()
    print(f"Loaded {len(rows)} variants to SQLite for {sample_name}")