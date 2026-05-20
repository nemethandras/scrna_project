# scripts/validate_vcf.py
"""
Validates an external (non-pipeline) VCF file before it is used for matching
or loaded into the database.

Checks performed:
  1. Reference genome detection (hg38 vs hg19/GRCh37) using contig lengths
     and chromosome naming conventions.
  2. Chromosome naming style (chr-prefixed vs bare).
  3. Rejects hg19 VCFs with a clear message; does not auto-liftover.

Usage (standalone):
    python scripts/validate_vcf.py path/to/sample.vcf[.gz]

Returns exit code 0 if the VCF passes, non-zero otherwise.
Can also be imported:
    from scripts.validate_vcf import validate_external_vcf
"""

import sys
import argparse
from cyvcf2 import VCF


# Chromosome 1 length differs between builds — reliable discriminator.
# GRCh37/hg19: 249,250,621   GRCh38/hg38: 248,956,422
_CHR1_LEN_HG19 = 249_250_621
_CHR1_LEN_HG38 = 248_956_422

# Accepted chr1 aliases in VCF headers
_CHR1_ALIASES = {"1", "chr1"}


def _detect_build_from_contigs(header_contigs):
    """
    Return 'hg38', 'hg19', or 'unknown' based on contig length metadata
    in the VCF header.  header_contigs is a dict {name: length_or_None}.
    """
    for alias in _CHR1_ALIASES:
        length = header_contigs.get(alias)
        if length is None:
            continue
        if length == _CHR1_LEN_HG38:
            return "hg38"
        if length == _CHR1_LEN_HG19:
            return "hg19"
    return "unknown"


def _detect_build_from_reference_line(header_str):
    """
    Scan the raw header string for a ##reference= line and look for
    build keywords.  Less reliable than contig lengths but useful as
    a fallback.
    """
    for line in header_str.splitlines():
        if not line.startswith("##reference="):
            continue
        low = line.lower()
        if "grch38" in low or "hg38" in low:
            return "hg38"
        if "grch37" in low or "hg19" in low or "b37" in low or "hs37" in low:
            return "hg19"
    return "unknown"


def _chrom_naming_style(vcf):
    """
    Return 'chr-prefixed', 'bare', or 'mixed' based on the first few
    records or contig names in the VCF.
    """
    names = set()
    # Prefer header contigs if present
    for c in vcf.seqnames:
        names.add(c)
        if len(names) >= 5:
            break

    has_chr = any(n.startswith("chr") for n in names)
    has_bare = any(not n.startswith("chr") for n in names)

    if has_chr and not has_bare:
        return "chr-prefixed"
    if has_bare and not has_chr:
        return "bare"
    return "mixed"


def validate_external_vcf(vcf_path, require_hg38=True):
    """
    Validate an external VCF file.

    Parameters
    ----------
    vcf_path : str
        Path to the VCF or VCF.gz file.
    require_hg38 : bool
        If True, raise ValueError for confirmed hg19 VCFs.

    Returns
    -------
    dict with keys: build, naming_style, contig_count, warnings
    """
    vcf = VCF(vcf_path)
    warnings = []

    # --- Build detection via contig lengths ---
    header_contigs = {}
    for contig in vcf.seqlens:
        pass  # cyvcf2 exposes seqlens as a list parallel to seqnames
    # cyvcf2 API: seqnames and seqlens are parallel lists
    header_contigs = dict(zip(vcf.seqnames, vcf.seqlens))

    build = _detect_build_from_contigs(header_contigs)

    if build == "unknown":
        # Fall back to ##reference= line
        build = _detect_build_from_reference_line(vcf.raw_header)

    if build == "unknown":
        warnings.append(
            "Could not determine reference build from contig lengths or ##reference= header. "
            "Verify manually that this VCF was called against GRCh38/hg38."
        )
    elif build == "hg19" and require_hg38:
        raise ValueError(
            f"VCF '{vcf_path}' appears to be aligned to hg19/GRCh37 "
            f"(chr1 length = {header_contigs.get('1') or header_contigs.get('chr1'):,}).\n"
            f"This pipeline uses GRCh38/hg38. Please liftover the VCF first.\n"
            f"Suggested tool: CrossMap (crossmap.sourceforge.net) with the hg19→hg38 chain file."
        )

    # --- Chromosome naming style ---
    naming = _chrom_naming_style(vcf)
    if naming == "bare":
        warnings.append(
            "Chromosome names are bare (e.g. '1', 'X') rather than chr-prefixed. "
            "The pipeline reference uses 'chr1', 'chrX' etc. "
            "Rename with: bcftools annotate --rename-chrs workflow/chr_name_conv.txt"
        )
    elif naming == "mixed":
        warnings.append(
            "VCF contains a mix of chr-prefixed and bare chromosome names — this is unusual "
            "and may indicate a malformed file."
        )

    result = {
        "build":         build,
        "naming_style":  naming,
        "contig_count":  len(header_contigs),
        "warnings":      warnings,
    }
    vcf.close()
    return result


def main():
    p = argparse.ArgumentParser(description="Validate an external VCF for use with this pipeline")
    p.add_argument("vcf", help="Path to VCF or VCF.gz file")
    p.add_argument("--allow-hg19", action="store_true",
                   help="Warn about hg19 instead of failing (useful for inspection only)")
    args = p.parse_args()

    print(f"\nValidating: {args.vcf}")
    try:
        result = validate_external_vcf(args.vcf, require_hg38=not args.allow_hg19)
    except ValueError as e:
        print(f"\nFAIL: {e}\n", file=sys.stderr)
        sys.exit(1)

    print(f"  Reference build : {result['build']}")
    print(f"  Chrom naming    : {result['naming_style']}")
    print(f"  Contigs in header: {result['contig_count']}")

    if result["warnings"]:
        print("\nWarnings:")
        for w in result["warnings"]:
            print(f"  - {w}")
    else:
        print("\nOK — no issues found.")

    sys.exit(0)


if __name__ == "__main__":
    main()
