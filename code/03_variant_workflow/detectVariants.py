#!/usr/bin/env python3
"""
detectVariants.py
-----------------
UMI-aware variant detection from minimap2-remapped BAMs. By default the script
attempts UMI-collapsed counting for Cell Ranger-derived BAMs that have preserved
CB/UB tags. Use --umi-off to force raw read-count mode for BAMs without UMIs.

Main outputs per BAM:
  <sample>.vector.variant.tsv
  <sample>.control_exons.variant.tsv
  <sample>.vector.variant_plot.png/.pdf
  <sample>.control_exons.variant_plot.png/.pdf
  <sample>.vector.significant_variants.vcf
  <sample>.control_exons.significant_variants.vcf
  <sample>.vector.cb_mutation.tsv   # only when --mut is supplied and --umi-off is not used

Statistical model:
  For each tested alternative allele/event A,C,G,T,Ins,Del:
      pval = P[X >= observed_count | X ~ Binomial(depth, error_rate)]
      sig  = pval < base_alpha / (n_bases * n_tests)

  n_bases is calculated from the number of positions written to the corresponding
  variant.tsv file. Reference-matching base counts are reported but not tested:
      <ref>_pval = NA, <ref>_sig = FALSE

UMI-collapsed allele assignment:
  A unique molecule is defined by CB + UB at a position. For each molecule/position,
  read-level observations vote for one of A,C,G,T,Ins,Del. The majority allele wins;
  ties are marked ambiguous and excluded. Reads missing CB or UB are excluded by
  default and summarized in the console preflight report and run log. Use --umi-off
  to disable UMI collapsing and treat each aligned read as unique.

Control exons:
  By default the script uses hard-coded hg38/GRCh38, UCSC-style, 1-based inclusive
  CDS-overlapping exon intervals for TRAC, CD3D, CD3E, CD3G, and CD247. These
  coordinates were downloaded from Ensembl REST using the helper script
  download_control_exons_v2.py and reviewed before hard-coding. You can override
  them with --control-exons-tsv.

Requirements:
  Python >= 3.6, pysam. Rscript + ggplot2 are required for plots unless --skip-plots. scipy is optional but faster for binomial p-values.

Example:
  python3 detectVariants.py \
    --bam PATIENT_ID_TIMEPOINT__Patient070_Late.vs.lentiviral_vector_ciltacel__hg38.sorted.bam \
    --reference lentiviral_vector_ciltacel__hg38.fasta \
    --out detectVariants_out

Preflight only:
  python3 detectVariants.py --bam sample.sorted.bam --reference ref.fasta --out out --preflight-only

  Preflight results are printed directly to the console.
"""

import argparse
import csv
import datetime as dt
import gzip
import html
import importlib
import json
import math
import os
import re
import platform
import shutil
import subprocess
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

BASES = ("A", "C", "G", "T")
ALLELES = ("A", "C", "G", "T", "Ins", "Del")
DEFAULT_VECTOR_KEYWORDS = (
    "vector", "idecel", "ciltacel", "axicel", "tisacel",
    "lisocel", "kymriah", "obecel", "car", "CAR", "leucel", "agene", "GC33"
)
CONTROL_GENES = ("TRAC", "CD3D", "CD3E", "CD3G", "CD247")

# Hard-coded control exon coordinates generated with download_control_exons_v2.py.
# Coordinates are GRCh38/hg38, UCSC-style chromosome names, 1-based inclusive.
# These are CDS-overlapping exon intervals from the selected Ensembl canonical
# transcript for each gene. TRAC is annotated as TR_C_gene rather than ordinary
# protein_coding, but Ensembl still returned CDS-overlapping exon intervals.
CONTROL_EXONS_HG38 = [
    {"gene": "CD247", "chrom": "chr1",  "start": 167431681, "end": 167431746, "strand": -1, "transcript_id": "ENST00000362089", "transcript_name": "CD247-201", "coordinate_type": "cds_overlapping_exon", "source": "Ensembl REST lookup/symbol expand=1, GRCh38/hg38"},
    {"gene": "CD247", "chrom": "chr1",  "start": 167433024, "end": 167433059, "strand": -1, "transcript_id": "ENST00000362089", "transcript_name": "CD247-201", "coordinate_type": "cds_overlapping_exon", "source": "Ensembl REST lookup/symbol expand=1, GRCh38/hg38"},
    {"gene": "CD247", "chrom": "chr1",  "start": 167434020, "end": 167434076, "strand": -1, "transcript_id": "ENST00000362089", "transcript_name": "CD247-201", "coordinate_type": "cds_overlapping_exon", "source": "Ensembl REST lookup/symbol expand=1, GRCh38/hg38"},
    {"gene": "CD247", "chrom": "chr1",  "start": 167435399, "end": 167435434, "strand": -1, "transcript_id": "ENST00000362089", "transcript_name": "CD247-201", "coordinate_type": "cds_overlapping_exon", "source": "Ensembl REST lookup/symbol expand=1, GRCh38/hg38"},
    {"gene": "CD247", "chrom": "chr1",  "start": 167438570, "end": 167438650, "strand": -1, "transcript_id": "ENST00000362089", "transcript_name": "CD247-201", "coordinate_type": "cds_overlapping_exon", "source": "Ensembl REST lookup/symbol expand=1, GRCh38/hg38"},
    {"gene": "CD247", "chrom": "chr1",  "start": 167439344, "end": 167439400, "strand": -1, "transcript_id": "ENST00000362089", "transcript_name": "CD247-201", "coordinate_type": "cds_overlapping_exon", "source": "Ensembl REST lookup/symbol expand=1, GRCh38/hg38"},
    {"gene": "CD247", "chrom": "chr1",  "start": 167440664, "end": 167440767, "strand": -1, "transcript_id": "ENST00000362089", "transcript_name": "CD247-201", "coordinate_type": "cds_overlapping_exon", "source": "Ensembl REST lookup/symbol expand=1, GRCh38/hg38"},
    {"gene": "CD247", "chrom": "chr1",  "start": 167518408, "end": 167518465, "strand": -1, "transcript_id": "ENST00000362089", "transcript_name": "CD247-201", "coordinate_type": "cds_overlapping_exon", "source": "Ensembl REST lookup/symbol expand=1, GRCh38/hg38"},
    {"gene": "CD3D",  "chrom": "chr11", "start": 118339162, "end": 118339227, "strand": -1, "transcript_id": "ENST00000300692", "transcript_name": "CD3D-201",  "coordinate_type": "cds_overlapping_exon", "source": "Ensembl REST lookup/symbol expand=1, GRCh38/hg38"},
    {"gene": "CD3D",  "chrom": "chr11", "start": 118339451, "end": 118339494, "strand": -1, "transcript_id": "ENST00000300692", "transcript_name": "CD3D-201",  "coordinate_type": "cds_overlapping_exon", "source": "Ensembl REST lookup/symbol expand=1, GRCh38/hg38"},
    {"gene": "CD3D",  "chrom": "chr11", "start": 118339775, "end": 118339906, "strand": -1, "transcript_id": "ENST00000300692", "transcript_name": "CD3D-201",  "coordinate_type": "cds_overlapping_exon", "source": "Ensembl REST lookup/symbol expand=1, GRCh38/hg38"},
    {"gene": "CD3D",  "chrom": "chr11", "start": 118340375, "end": 118340593, "strand": -1, "transcript_id": "ENST00000300692", "transcript_name": "CD3D-201",  "coordinate_type": "cds_overlapping_exon", "source": "Ensembl REST lookup/symbol expand=1, GRCh38/hg38"},
    {"gene": "CD3D",  "chrom": "chr11", "start": 118342553, "end": 118342607, "strand": -1, "transcript_id": "ENST00000300692", "transcript_name": "CD3D-201",  "coordinate_type": "cds_overlapping_exon", "source": "Ensembl REST lookup/symbol expand=1, GRCh38/hg38"},
    {"gene": "CD3E",  "chrom": "chr11", "start": 118304953, "end": 118305001, "strand": 1,  "transcript_id": "ENST00000361763", "transcript_name": "CD3E-201",  "coordinate_type": "cds_overlapping_exon", "source": "Ensembl REST lookup/symbol expand=1, GRCh38/hg38"},
    {"gene": "CD3E",  "chrom": "chr11", "start": 118307288, "end": 118307308, "strand": 1,  "transcript_id": "ENST00000361763", "transcript_name": "CD3E-201",  "coordinate_type": "cds_overlapping_exon", "source": "Ensembl REST lookup/symbol expand=1, GRCh38/hg38"},
    {"gene": "CD3E",  "chrom": "chr11", "start": 118308427, "end": 118308441, "strand": 1,  "transcript_id": "ENST00000361763", "transcript_name": "CD3E-201",  "coordinate_type": "cds_overlapping_exon", "source": "Ensembl REST lookup/symbol expand=1, GRCh38/hg38"},
    {"gene": "CD3E",  "chrom": "chr11", "start": 118312153, "end": 118312170, "strand": 1,  "transcript_id": "ENST00000361763", "transcript_name": "CD3E-201",  "coordinate_type": "cds_overlapping_exon", "source": "Ensembl REST lookup/symbol expand=1, GRCh38/hg38"},
    {"gene": "CD3E",  "chrom": "chr11", "start": 118312618, "end": 118312866, "strand": 1,  "transcript_id": "ENST00000361763", "transcript_name": "CD3E-201",  "coordinate_type": "cds_overlapping_exon", "source": "Ensembl REST lookup/symbol expand=1, GRCh38/hg38"},
    {"gene": "CD3E",  "chrom": "chr11", "start": 118313707, "end": 118313874, "strand": 1,  "transcript_id": "ENST00000361763", "transcript_name": "CD3E-201",  "coordinate_type": "cds_overlapping_exon", "source": "Ensembl REST lookup/symbol expand=1, GRCh38/hg38"},
    {"gene": "CD3E",  "chrom": "chr11", "start": 118314448, "end": 118314494, "strand": 1,  "transcript_id": "ENST00000361763", "transcript_name": "CD3E-201",  "coordinate_type": "cds_overlapping_exon", "source": "Ensembl REST lookup/symbol expand=1, GRCh38/hg38"},
    {"gene": "CD3E",  "chrom": "chr11", "start": 118315486, "end": 118315542, "strand": 1,  "transcript_id": "ENST00000361763", "transcript_name": "CD3E-201",  "coordinate_type": "cds_overlapping_exon", "source": "Ensembl REST lookup/symbol expand=1, GRCh38/hg38"},
    {"gene": "CD3G",  "chrom": "chr11", "start": 118344424, "end": 118344478, "strand": 1,  "transcript_id": "ENST00000532917", "transcript_name": "CD3G-206",  "coordinate_type": "cds_overlapping_exon", "source": "Ensembl REST lookup/symbol expand=1, GRCh38/hg38"},
    {"gene": "CD3G",  "chrom": "chr11", "start": 118349027, "end": 118349050, "strand": 1,  "transcript_id": "ENST00000532917", "transcript_name": "CD3G-206",  "coordinate_type": "cds_overlapping_exon", "source": "Ensembl REST lookup/symbol expand=1, GRCh38/hg38"},
    {"gene": "CD3G",  "chrom": "chr11", "start": 118349743, "end": 118349970, "strand": 1,  "transcript_id": "ENST00000532917", "transcript_name": "CD3G-206",  "coordinate_type": "cds_overlapping_exon", "source": "Ensembl REST lookup/symbol expand=1, GRCh38/hg38"},
    {"gene": "CD3G",  "chrom": "chr11", "start": 118350552, "end": 118350683, "strand": 1,  "transcript_id": "ENST00000532917", "transcript_name": "CD3G-206",  "coordinate_type": "cds_overlapping_exon", "source": "Ensembl REST lookup/symbol expand=1, GRCh38/hg38"},
    {"gene": "CD3G",  "chrom": "chr11", "start": 118351628, "end": 118351671, "strand": 1,  "transcript_id": "ENST00000532917", "transcript_name": "CD3G-206",  "coordinate_type": "cds_overlapping_exon", "source": "Ensembl REST lookup/symbol expand=1, GRCh38/hg38"},
    {"gene": "CD3G",  "chrom": "chr11", "start": 118352404, "end": 118352469, "strand": 1,  "transcript_id": "ENST00000532917", "transcript_name": "CD3G-206",  "coordinate_type": "cds_overlapping_exon", "source": "Ensembl REST lookup/symbol expand=1, GRCh38/hg38"},
    {"gene": "TRAC",  "chrom": "chr14", "start": 22547506,  "end": 22547778,  "strand": 1,  "transcript_id": "ENST00000611116", "transcript_name": "TRAC-201",  "coordinate_type": "cds_overlapping_exon", "source": "Ensembl REST lookup/symbol expand=1, GRCh38/hg38"},
    {"gene": "TRAC",  "chrom": "chr14", "start": 22549638,  "end": 22549682,  "strand": 1,  "transcript_id": "ENST00000611116", "transcript_name": "TRAC-201",  "coordinate_type": "cds_overlapping_exon", "source": "Ensembl REST lookup/symbol expand=1, GRCh38/hg38"},
    {"gene": "TRAC",  "chrom": "chr14", "start": 22550557,  "end": 22550663,  "strand": 1,  "transcript_id": "ENST00000611116", "transcript_name": "TRAC-201",  "coordinate_type": "cds_overlapping_exon", "source": "Ensembl REST lookup/symbol expand=1, GRCh38/hg38"},
]

try:
    import pysam  # type: ignore
except Exception:  # pragma: no cover
    pysam = None

try:
    from scipy.stats import binom as scipy_binom  # type: ignore
except Exception:  # pragma: no cover
    scipy_binom = None

class Region:
    def __init__(self, gene, chrom, start, end, transcript_id="", source=""):
        self.gene = gene
        self.chrom = chrom
        self.start = int(start)  # 1-based inclusive
        self.end = int(end)      # 1-based inclusive
        self.transcript_id = transcript_id or ""
        self.source = source or ""


class PileupStats:
    def __init__(self):
        self.missing_cbub_reads = 0
        self.ambiguous_umis = 0
        self.usable_umis = 0
        self.reads_seen = 0


class MutationSpec:
    """Requested mutation using vector/reference coordinates: ref + 1-based position + alt."""
    def __init__(self, ref_base: str, pos1: int, alt: str):
        self.ref_base = ref_base.upper()
        self.pos1 = int(pos1)
        self.alt = alt
        self.label = f"{self.ref_base}{self.pos1}{self.alt}"

    def as_dict(self) -> dict:
        return {"ref_base": self.ref_base, "pos1": self.pos1, "alt": self.alt, "label": self.label}


def canonical_alt_label(alt: str) -> Optional[str]:
    """Return canonical allele/event label for --mut alternatives."""
    if alt is None:
        return None
    s = str(alt).strip()
    if len(s) == 1 and s.upper() in BASES:
        return s.upper()
    if s.lower() == "ins":
        return "Ins"
    if s.lower() == "del":
        return "Del"
    return None


def parse_mutation_specs(raw_values: Optional[Sequence[str]]) -> List[MutationSpec]:
    """
    Parse --mut/-mut values such as G222A or G222A,G314A.

    Format: [reference base][1-based position][mutated base/event], where
    reference base is A/C/G/T and mutated base/event is A/C/G/T/Ins/Del.
    """
    specs: List[MutationSpec] = []
    seen = set()
    for raw_group in raw_values or []:
        for raw in str(raw_group).split(","):
            s = raw.strip()
            if not s:
                continue
            m = re.fullmatch(r"([ACGTacgt])([1-9][0-9]*)([A-Za-z]+)", s)
            if not m:
                raise SystemExit(
                    f"ERROR: invalid --mut value '{s}'. Expected format like G222A or G222Ins. "
                    "Reference base must be A/C/G/T; position must be 1-based; alternate must be A/C/G/T/Ins/Del."
                )
            ref_base = m.group(1).upper()
            pos1 = int(m.group(2))
            alt = canonical_alt_label(m.group(3))
            if alt is None:
                raise SystemExit(
                    f"ERROR: invalid --mut alternate in '{s}'. Alternate must be A/C/G/T/Ins/Del."
                )
            label = f"{ref_base}{pos1}{alt}"
            if label in seen:
                continue
            seen.add(label)
            specs.append(MutationSpec(ref_base, pos1, alt))
    return specs


def timestamp() -> str:
    return dt.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")


def now_human() -> str:
    return dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def format_int(x) -> str:
    try:
        return f"{int(x):,}"
    except Exception:
        return str(x)


def file_size_human(path: Path) -> str:
    size = path.stat().st_size
    units = ["B", "KB", "MB", "GB", "TB"]
    value = float(size)
    for unit in units:
        if value < 1024 or unit == units[-1]:
            return f"{value:.2f} {unit}"
        value /= 1024
    return f"{value:.2f} TB"


class Logger:
    def __init__(self, log_path: Path):
        self.log_path = log_path
        self.log_path.parent.mkdir(parents=True, exist_ok=True)

    def write(self, message: str, also_stderr: bool = True) -> None:
        line = f"[{now_human()}] {message}"
        with open(self.log_path, "a", encoding="utf-8") as handle:
            handle.write(line + "\n")
        if also_stderr:
            print(line, file=sys.stderr)

    def section(self, title: str) -> None:
        self.write("")
        self.write("=" * 80)
        self.write(title)
        self.write("=" * 80)


GLOBAL_LOGGER = None


def log(msg: str) -> None:
    if GLOBAL_LOGGER is not None:
        GLOBAL_LOGGER.write(msg)
    else:
        sys.stderr.write(f"[{now_human()}] {msg}\n")
        sys.stderr.flush()


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="UMI-collapsed vector/control-exon variant detection from remapped BAMs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--bam", nargs="*", default=[], help="One or more sorted/indexed BAM files")
    p.add_argument("--bam-list", help="Text file with one BAM path per line")
    p.add_argument("--reference", required=True, help="Reference FASTA used for minimap2 alignment")
    p.add_argument("--out", required=True, help="Output directory")
    p.add_argument("--preflight-only", action="store_true", help="Run dependency/input checks and exit")
    p.add_argument("--error-rate", type=float, default=0.0024, help="Background error rate")
    p.add_argument("--base-alpha", type=float, default=0.001, help="Family-wise alpha before Bonferroni adjustment")
    p.add_argument("--n-tests", type=int, default=6, help="Allele/event tests per position")
    p.add_argument("--min-mapq", type=int, default=0, help="Minimum read MAPQ")
    p.add_argument("--min-base-quality", type=int, default=0, help="Minimum base quality for base observations")
    p.add_argument(
        "--umi-off",
        action="store_true",
        help=(
            "Disable UMI collapsing and treat each aligned read as unique. "
            "Use this for bulk/SRA/non-Cell-Ranger BAMs or BAMs without CB/UB tags."
        ),
    )
    p.add_argument(
        "-mut", "--mut",
        action="append",
        default=[],
        help=(
            "Comma-separated requested vector mutation(s), e.g. G222A,G314A or G222Ins. "
            "Format is [reference base][1-based position][mutated base/event], where alternate is A/C/G/T/Ins/Del. "
            "Only valid when UMI collapsing is enabled; do not combine with --umi-off."
        ),
    )
    p.add_argument("--vector-keywords", default=",".join(DEFAULT_VECTOR_KEYWORDS), help="Comma-separated vector contig keywords")
    p.add_argument("--control-exons-tsv", help="Optional 1-based inclusive exon BED-like TSV: gene,chrom,start,end,transcript_id,source")
    p.add_argument("--skip-control", action="store_true", help="Do not generate control-exon TSV/plots")
    p.add_argument("--skip-plots", action="store_true", help="Do not generate static PNG plots")
    p.add_argument("--skip-html", action="store_true", help="Do not generate interactive HTML variant browsers")
    p.add_argument("--html-flank", type=int, default=5, help="Number of reference bases on each side for interactive HTML sequence context")
    p.add_argument("--html-px-per-position", type=float, default=2.2, help="Intrinsic horizontal pixels per variant-table row for the interactive HTML SVG before browser scaling")
    p.add_argument("--html-max-width", type=int, default=60000, help="Maximum intrinsic SVG width for interactive HTML before browser scaling")
    p.add_argument("--skip-vcf", action="store_true", help="Do not write simple VCF files for significant variant calls")
    p.add_argument(
        "--write-region-bams",
        action="store_true",
        help=(
            "Write indexed BAM subsets for reads overlapping the detected vector region(s) "
            "and, unless --skip-control is used, the control exon regions"
        ),
    )
    p.add_argument("--keep-uncovered-vector-positions", action="store_true", default=True,
                   help="Write all vector reference positions, including zero-depth positions")
    p.add_argument("--version", action="version", version="detectVariants.py 0.3.5")
    return p.parse_args()


def collect_bams(args: argparse.Namespace) -> List[Path]:
    bams = [Path(x) for x in args.bam]
    if args.bam_list:
        with open(args.bam_list) as fh:
            for line in fh:
                s = line.strip()
                if s and not s.startswith("#"):
                    bams.append(Path(s))
    seen, out = set(), []
    for b in bams:
        bp = b.expanduser()
        if bp not in seen:
            out.append(bp)
            seen.add(bp)
    if not out:
        raise SystemExit("ERROR: provide at least one BAM via --bam or --bam-list")
    return out


def sample_name_from_bam(bam: Path) -> Tuple[str, str, str, str]:
    name = bam.name
    for suffix in (".sorted.bam", ".bam"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    sample_id = name.split(".vs.", 1)[0]
    patient_id = "NA"
    timepoint = "NA"
    # Pattern used by demultiplexed outputs: PREFIX__Patient070_Late.vs....bam
    if "__" in sample_id:
        right = sample_id.split("__", 1)[1]
    else:
        right = sample_id
    m = re.match(r"([^_]+)_(.+)$", right)
    if m:
        patient_id, timepoint = m.group(1), m.group(2)
    return sample_id, patient_id, timepoint, bam.name


def open_reference(ref_path: Path):
    if pysam is None:
        raise RuntimeError("pysam is not available")
    try:
        return pysam.FastaFile(str(ref_path))
    except Exception as e:
        raise RuntimeError(f"could not open reference FASTA {ref_path}: {e}")


def reference_base(ref, chrom: str, pos1: int) -> str:
    try:
        b = ref.fetch(chrom, pos1 - 1, pos1).upper()
        return b if b in BASES else "N"
    except Exception:
        return "N"


def matches_vector(chrom: str, keywords: Sequence[str]) -> bool:
    c = chrom.lower()
    return any(k.lower() in c for k in keywords if k)


def bam_index_exists(bam: Path) -> bool:
    return Path(str(bam) + ".bai").exists() or bam.with_suffix(".bai").exists()


def record_preflight(rows: List[Tuple[str, str, str]], name: str, status: str, details: str = "") -> None:
    rows.append((name, status, details))


def status_label(status: str) -> str:
    if status == "PASS":
        return "[ok]"
    if status == "WARN":
        return "[warning]"
    return "[error]"


def check_python_package(rows: List[Tuple[str, str, str]], logger: Logger, package_name: str, required: bool = True) -> None:
    try:
        module = importlib.import_module(package_name)
        version = getattr(module, "__version__", "unknown")
        path = getattr(module, "__file__", "")
        logger.write(f"[ok] Python package available: {package_name} {version}")
        if path:
            logger.write(f"[info] {package_name} path: {path}", also_stderr=False)
        record_preflight(rows, f"python_package:{package_name}", "PASS", f"{version}; {path}")
    except Exception as e:
        status = "FAIL" if required else "WARN"
        logger.write(f"[{ 'error' if required else 'warning' }] Python package missing: {package_name}; {e}")
        record_preflight(rows, f"python_package:{package_name}", status, str(e))


def check_executable(rows: List[Tuple[str, str, str]], logger: Logger, executable: str, required: bool = False) -> None:
    exe = shutil.which(executable)
    if exe is None:
        status = "FAIL" if required else "WARN"
        logger.write(f"[{ 'error' if required else 'warning' }] Executable not found on PATH: {executable}")
        record_preflight(rows, f"executable:{executable}", status, "not found on PATH")
        return
    logger.write(f"[ok] Executable available: {executable} -> {exe}")
    record_preflight(rows, f"executable:{executable}", "PASS", exe)
    if executable == "samtools":
        try:
            result = subprocess.run(["samtools", "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            if result.returncode == 0:
                first = result.stdout.splitlines()[0] if result.stdout else "unknown"
                logger.write(f"[ok] samtools version: {first}")
                record_preflight(rows, "samtools_version", "PASS", first)
            else:
                logger.write("[warning] Could not determine samtools version")
                record_preflight(rows, "samtools_version", "WARN", result.stderr.strip())
        except Exception as e:
            logger.write(f"[warning] Could not determine samtools version: {e}")
            record_preflight(rows, "samtools_version", "WARN", str(e))


def check_r_plotting_dependency(rows, logger):
    logger.section("R plotting dependency preflight")

    rscript = shutil.which("Rscript")
    if rscript is None:
        logger.write("[error] Rscript not found on PATH; required for plots unless --skip-plots")
        record_preflight(rows, "executable:Rscript", "FAIL", "required for plots unless --skip-plots")
        return

    logger.write(f"[ok] Executable available: Rscript -> {rscript}")
    record_preflight(rows, "executable:Rscript", "PASS", rscript)

    try:
        version = subprocess.run(["Rscript", "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        version_text = (version.stdout + version.stderr).strip()
        if version_text:
            logger.write(f"[ok] {version_text}")
            record_preflight(rows, "Rscript_version", "PASS", version_text)
    except Exception as e:
        logger.write(f"[warning] Could not determine Rscript version: {e}")
        record_preflight(rows, "Rscript_version", "WARN", str(e))

    r_code = """
required <- c("ggplot2")
for (pkg in required) {
  ok <- requireNamespace(pkg, quietly = TRUE)
  if (!ok) {
    cat(paste0(pkg, "\tMISSING\t\tPackage not installed or not loadable\n"))
  } else {
    ver <- as.character(utils::packageVersion(pkg))
    path <- find.package(pkg)
    cat(paste0(pkg, "\tOK\t", ver, "\t", path, "\n"))
  }
}
"""
    try:
        result = subprocess.run(["Rscript", "-e", r_code], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    except Exception as e:
        logger.write(f"[error] Could not run Rscript package check: {e}")
        record_preflight(rows, "R_package_check", "FAIL", str(e))
        return

    if result.stderr.strip():
        logger.write("[R plotting dependency stderr]\n" + result.stderr.strip(), also_stderr=False)

    if result.returncode != 0:
        logger.write("[error] R package check failed while checking ggplot2")
        record_preflight(rows, "R_package_check", "FAIL", result.stderr.strip())
        return

    found_ggplot2 = False
    for line in result.stdout.strip().splitlines():
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        pkg = parts[0]
        status = parts[1]
        version = parts[2] if len(parts) > 2 else ""
        path = parts[3] if len(parts) > 3 else ""
        if status == "OK":
            logger.write(f"[ok] R package available: {pkg} {version}")
            record_preflight(rows, f"R_package:{pkg}", "PASS", f"{version}; {path}")
            if pkg == "ggplot2":
                found_ggplot2 = True
        else:
            logger.write(f"[error] R package missing: {pkg}")
            record_preflight(rows, f"R_package:{pkg}", "FAIL", "Package not installed or not loadable")

    if not found_ggplot2:
        logger.write('[error] Missing required R package: ggplot2. Install in R with install.packages("ggplot2")')
        record_preflight(rows, "R_package:ggplot2", "FAIL", "Install in R with install.packages('ggplot2')")


def check_preflight(bams: List[Path], ref_path: Path, outdir: Path, args: argparse.Namespace, logger: Logger) -> List[Tuple[str, str, str]]:
    rows: List[Tuple[str, str, str]] = []

    logger.section("Forced preflight")
    logger.write(f"[info] Output directory: {outdir.resolve()}")
    logger.write(f"[info] Number of BAM inputs: {format_int(len(bams))}")
    logger.write(f"[info] --preflight-only: {args.preflight_only}")
    logger.write(f"[info] --skip-plots: {args.skip_plots}")
    logger.write(f"[info] --skip-html: {getattr(args, 'skip_html', False)}")
    logger.write(f"[info] --skip-control: {args.skip_control}")
    logger.write(f"[info] --umi-off: {args.umi_off}")
    if args.umi_off:
        logger.write("[info] --umi-off enabled: UMI tag absence will not block variant counting; raw read counts will be used")
    mut_specs = getattr(args, "mut_specs", [])
    if mut_specs:
        labels = ",".join(m.label for m in mut_specs)
        logger.write(f"[info] --mut requested: {labels}")
        record_preflight(rows, "mut_specs", "PASS", labels)
        if args.umi_off:
            logger.write("[error] --mut requires CB/UB-tagged UMI mode; remove --umi-off to write the CB mutation table")
            record_preflight(rows, "mut_requires_umi_mode", "FAIL", "--mut cannot be combined with --umi-off")

    logger.section("Python dependency preflight")
    py_version = platform.python_version()
    py_path = sys.executable
    logger.write(f"[ok] Python executable: {py_path}")
    logger.write(f"[ok] Python version: {py_version}")
    record_preflight(rows, "python", "PASS" if sys.version_info >= (3, 6) else "FAIL", f"{py_version}; {py_path}")
    if sys.version_info < (3, 8):
        logger.write(f"[warning] Python >= 3.8 is recommended; found {py_version}")
        record_preflight(rows, "python_version_recommendation", "WARN", f"Python >= 3.8 recommended; found {py_version}")

    if pysam is None:
        logger.write("[error] Python package missing: pysam")
        record_preflight(rows, "python_package:pysam", "FAIL", "required")
    else:
        check_python_package(rows, logger, "pysam", required=True)

    if args.skip_plots:
        logger.write("[info] Plotting preflight skipped because --skip-plots was requested")
        record_preflight(rows, "plotting:Rscript_ggplot2", "WARN", "skipped because --skip-plots")
    else:
        check_r_plotting_dependency(rows, logger)

    if scipy_binom is None:
        logger.write("[warning] Python package scipy unavailable; exact Python fallback will be used for binomial p-values")
        record_preflight(rows, "python_package:scipy", "WARN", "optional; exact Python fallback will be used")
    else:
        check_python_package(rows, logger, "scipy", required=False)

    logger.section("Executable dependency preflight")
    check_executable(rows, logger, "samtools", required=False)

    logger.section("Reference preflight")
    if not ref_path.exists():
        logger.write(f"[error] Reference FASTA not found: {ref_path}")
        record_preflight(rows, "reference_exists", "FAIL", str(ref_path))
    else:
        logger.write(f"[ok] Reference FASTA exists: {ref_path}")
        logger.write(f"[info] Reference FASTA size: {file_size_human(ref_path)}")
        record_preflight(rows, "reference_exists", "PASS", str(ref_path))

    if pysam is not None and ref_path.exists():
        try:
            ref = open_reference(ref_path)
            refs = list(ref.references)
            logger.write("[ok] Reference FASTA opened successfully")
            logger.write(f"[info] Number of reference sequences: {format_int(len(refs))}")
            if refs:
                preview = ", ".join(refs[:10])
                logger.write(f"[info] First reference sequences: {preview}")
            record_preflight(rows, "reference_readable", "PASS", ",".join(refs[:5]) + ("..." if len(refs) > 5 else ""))
            ref.close()
        except Exception as e:
            logger.write(f"[error] Could not open reference FASTA with pysam: {e}")
            record_preflight(rows, "reference_readable", "FAIL", str(e))

    keywords = [x.strip() for x in args.vector_keywords.split(",") if x.strip()]
    logger.write(f"[info] Vector keyword list: {', '.join(keywords)}")

    for bam in bams:
        logger.section(f"BAM preflight: {bam.name}")
        if not bam.exists():
            logger.write(f"[error] BAM not found: {bam}")
            record_preflight(rows, f"bam_exists:{bam.name}", "FAIL", str(bam))
            continue

        logger.write(f"[ok] BAM exists: {bam}")
        logger.write(f"[info] BAM size: {file_size_human(bam)}")
        record_preflight(rows, f"bam_exists:{bam.name}", "PASS", str(bam))

        if bam_index_exists(bam):
            logger.write("[ok] BAM index found")
            record_preflight(rows, f"bam_index:{bam.name}", "PASS", "required for region pileups")
        else:
            logger.write("[error] BAM index not found. Run: samtools index <bam>")
            record_preflight(rows, f"bam_index:{bam.name}", "FAIL", "required for region pileups")

        if pysam is None:
            continue

        try:
            with pysam.AlignmentFile(str(bam), "rb") as bf:
                refs = list(bf.references)
                lengths = list(bf.lengths)
                logger.write("[ok] BAM opened successfully")
                logger.write(f"[info] Number of reference sequences in BAM header: {format_int(len(refs))}")
                if refs:
                    preview = []
                    for ref_name, length in list(zip(refs, lengths))[:10]:
                        preview.append(f"{ref_name}:{length}")
                    logger.write(f"[info] First reference sequences: {', '.join(preview)}")
                record_preflight(rows, f"bam_readable:{bam.name}", "PASS", f"{len(refs)} references")

                vector_contigs = [r for r in refs if matches_vector(r, keywords)]
                if vector_contigs:
                    logger.write(f"[ok] Vector-like contigs matched: {', '.join(vector_contigs[:10])}")
                    record_preflight(rows, f"vector_contig_match:{bam.name}", "PASS", ",".join(vector_contigs[:10]))
                else:
                    logger.write("[error] No BAM contigs matched the vector keyword list")
                    logger.write(f"[info] Keywords: {', '.join(keywords)}")
                    record_preflight(rows, f"vector_contig_match:{bam.name}", "FAIL", "no contigs matched keywords")

                logger.section(f"BAM CB/UB tag sampling: {bam.name}")
                cb = ub = both = sampled = 0
                top_cb: Counter = Counter()
                for rec in bf.fetch(until_eof=True):
                    sampled += 1
                    has_cb = rec.has_tag("CB")
                    has_ub = rec.has_tag("UB")
                    cb += int(has_cb)
                    ub += int(has_ub)
                    both += int(has_cb and has_ub)
                    if has_cb:
                        try:
                            top_cb.update([rec.get_tag("CB")])
                        except Exception:
                            pass
                    if sampled >= 10000:
                        break
                logger.write(f"[info] Sampled records: {format_int(sampled)}")
                logger.write(f"[info] Sampled records with CB: {format_int(cb)}")
                logger.write(f"[info] Sampled records with UB: {format_int(ub)}")
                logger.write(f"[info] Sampled records with both CB and UB: {format_int(both)}")
                logger.write(f"[info] Unique CB values in sample: {format_int(len(top_cb))}")
                if top_cb:
                    logger.write("[info] Top sampled CB tags:")
                    for tag, count in top_cb.most_common(10):
                        logger.write(f"  {tag}\t{count}", also_stderr=False)
                if both > 0:
                    logger.write("[ok] UMI tags detected: sampled records include both CB and UB tags; UMI-collapsed counting is available")
                    record_preflight(rows, f"CB_UB_tags_sample:{bam.name}", "PASS", f"sampled={sampled}; CB={cb}; UB={ub}; both={both}")
                else:
                    if args.umi_off:
                        logger.write("[warning] No sampled records had both CB and UB tags; proceeding is allowed because --umi-off was requested")
                        record_preflight(rows, f"CB_UB_tags_sample:{bam.name}", "WARN", f"sampled={sampled}; CB={cb}; UB={ub}; both={both}; --umi-off enabled")
                    else:
                        logger.write("[warning] No sampled records had both CB and UB tags; use --umi-off if you want to proceed with raw read counts instead of UMI-collapsed counts")
                        record_preflight(rows, f"CB_UB_tags_sample:{bam.name}", "WARN", f"sampled={sampled}; CB={cb}; UB={ub}; both={both}; consider --umi-off")
        except Exception as e:
            logger.write(f"[error] Could not open/read BAM with pysam: {e}")
            record_preflight(rows, f"bam_readable:{bam.name}", "FAIL", str(e))

    logger.section("Preflight summary")
    counts = Counter(status for _check, status, _details in rows)
    logger.write(f"[info] PASS: {format_int(counts.get('PASS', 0))}")
    logger.write(f"[info] WARN: {format_int(counts.get('WARN', 0))}")
    logger.write(f"[info] FAIL: {format_int(counts.get('FAIL', 0))}")

    if counts.get("FAIL", 0):
        logger.write("[error] Preflight failed")
        for check, status, details in rows:
            if status == "FAIL":
                logger.write(f"[error] {check} | {details}")
    elif counts.get("WARN", 0):
        logger.write("[ok] Preflight complete with warnings")
    else:
        logger.write("[ok] Preflight complete")

    return rows


def hardcoded_control_exons() -> List[Region]:
    """Return built-in hg38 control exon intervals as Region objects."""
    return [
        Region(
            gene=r["gene"],
            chrom=r["chrom"],
            start=int(r["start"]),
            end=int(r["end"]),
            transcript_id=r.get("transcript_id", ""),
            source=r.get("source", "hard-coded CONTROL_EXONS_HG38"),
        )
        for r in CONTROL_EXONS_HG38
    ]


def read_control_exons_tsv(path: Path) -> List[Region]:
    """
    Read optional control exon override TSV.

    Accepted coordinate column names:
      gene, chrom, start, end
    or the downloader output:
      gene, chrom, exon_start, exon_end

    Coordinates must be 1-based inclusive.
    """
    regions: List[Region] = []
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            raise RuntimeError(f"{path} has no header")
        fields = set(reader.fieldnames)
        if {"gene", "chrom", "start", "end"}.issubset(fields):
            start_col, end_col = "start", "end"
        elif {"gene", "chrom", "exon_start", "exon_end"}.issubset(fields):
            start_col, end_col = "exon_start", "exon_end"
        else:
            raise RuntimeError(
                f"{path} must contain either gene/chrom/start/end or "
                "gene/chrom/exon_start/exon_end"
            )
        for row in reader:
            tx = row.get("transcript_id", "")
            source = row.get("source", row.get("coordinate_type", str(path)))
            regions.append(Region(
                gene=row["gene"],
                chrom=row["chrom"],
                start=int(row[start_col]),
                end=int(row[end_col]),
                transcript_id=tx,
                source=source,
            ))
    return regions

def write_control_regions(regions: Sequence[Region], out_path: Path) -> None:
    with open(out_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["gene", "chrom", "start", "end", "transcript_id", "source"])
        for r in regions:
            w.writerow([r.gene, r.chrom, r.start, r.end, r.transcript_id, r.source])


def allele_from_pread(pread, min_base_quality: int) -> Optional[str]:
    if pread.is_refskip:
        return None
    if pread.is_del:
        return "Del"
    # Treat insertion as a mutually exclusive event for depth accounting.
    if pread.indel and pread.indel > 0:
        return "Ins"
    qpos = pread.query_position
    if qpos is None:
        return None
    aln = pread.alignment
    if aln.query_qualities is not None and qpos < len(aln.query_qualities):
        if aln.query_qualities[qpos] < min_base_quality:
            return None
    b = aln.query_sequence[qpos].upper()
    return b if b in BASES else None


def get_cb_ub(rec) -> Optional[Tuple[str, str]]:
    try:
        cb = rec.get_tag("CB")
        ub = rec.get_tag("UB")
    except Exception:
        return None
    if not cb or not ub:
        return None
    return str(cb), str(ub)


def collapse_votes_to_counts(votes: Dict[Tuple[str, str], Counter], stats: PileupStats) -> Dict[str, int]:
    counts = {a: 0 for a in ALLELES}
    for _umi, counter in votes.items():
        if not counter:
            continue
        max_count = max(counter.values())
        winners = [a for a, c in counter.items() if c == max_count]
        if len(winners) != 1:
            stats.ambiguous_umis += 1
            continue
        allele = winners[0]
        if allele in counts:
            counts[allele] += 1
            stats.usable_umis += 1
    return counts



def _read_unique_key(rec) -> Tuple[str, int, int, int, str, int, int, int]:
    """Conservative key used to avoid duplicate writes from overlapping regions."""
    return (
        rec.query_name,
        int(rec.flag),
        int(rec.reference_id),
        int(rec.reference_start),
        rec.cigarstring or "",
        int(rec.next_reference_id),
        int(rec.next_reference_start),
        int(rec.template_length),
    )


def subset_bam_to_regions(
    bam_path: Path,
    out_bam_path: Path,
    regions: Sequence[Region],
    label: str,
) -> int:
    """
    Write a BAM containing reads overlapping one or more 1-based inclusive regions.

    Region.start/end are 1-based inclusive. pysam.fetch expects 0-based half-open
    intervals, so start is shifted by -1 and end is passed unchanged.

    Returns the number of records written.
    """
    if pysam is None:
        raise RuntimeError("pysam is required")

    out_bam_path.parent.mkdir(parents=True, exist_ok=True)
    n_written = 0
    n_skipped_missing_contig = 0
    seen = set()

    with pysam.AlignmentFile(str(bam_path), "rb") as in_bam:
        bam_refs = set(in_bam.references)
        with pysam.AlignmentFile(str(out_bam_path), "wb", template=in_bam) as out_bam:
            for r in regions:
                if r.chrom not in bam_refs:
                    n_skipped_missing_contig += 1
                    continue
                try:
                    iterator = in_bam.fetch(r.chrom, r.start - 1, r.end)
                except Exception as e:
                    log(f"WARNING: could not subset {label} BAM for {r.gene} {r.chrom}:{r.start}-{r.end}: {e}")
                    continue
                for rec in iterator:
                    key = _read_unique_key(rec)
                    if key in seen:
                        continue
                    seen.add(key)
                    out_bam.write(rec)
                    n_written += 1

    try:
        pysam.index(str(out_bam_path))
        log(f"Wrote indexed {label} BAM subset: {out_bam_path} ({format_int(n_written)} records)")
    except Exception as e:
        log(f"WARNING: wrote {label} BAM subset but could not index it: {out_bam_path}; {e}")

    if n_skipped_missing_contig:
        log(f"WARNING: skipped {format_int(n_skipped_missing_contig)} {label} region(s) absent from BAM header")

    return n_written

def count_region_positions(
    bam,
    region: Region,
    min_mapq: int,
    min_base_quality: int,
    use_umi: bool = True,
) -> Tuple[Dict[int, Dict[str, int]], PileupStats]:
    stats = PileupStats()
    position_counts: Dict[int, Dict[str, int]] = {}
    try:
        iterator = bam.pileup(
            contig=region.chrom,
            start=region.start - 1,
            stop=region.end,
            truncate=True,
            stepper="samtools",
            ignore_overlaps=False,
            ignore_orphans=False,
            min_base_quality=0,
            min_mapping_quality=min_mapq,
            max_depth=1000000,
        )
    except Exception as e:
        log(f"WARNING: could not pileup {region.gene} {region.chrom}:{region.start}-{region.end}: {e}")
        return position_counts, stats

    for col in iterator:
        pos1 = col.reference_pos + 1
        if use_umi:
            votes: Dict[Tuple[str, str], Counter] = defaultdict(Counter)
        else:
            read_counts = {a: 0 for a in ALLELES}

        for pread in col.pileups:
            rec = pread.alignment
            stats.reads_seen += 1
            if rec.is_unmapped or rec.mapping_quality < min_mapq:
                continue

            allele = allele_from_pread(pread, min_base_quality)
            if allele is None:
                continue

            if use_umi:
                key = get_cb_ub(rec)
                if key is None:
                    stats.missing_cbub_reads += 1
                    continue
                votes[key][allele] += 1
            else:
                if allele in read_counts:
                    read_counts[allele] += 1
                    stats.usable_umis += 1  # retained for summary compatibility; represents usable read observations in read-count mode

        if use_umi:
            position_counts[pos1] = collapse_votes_to_counts(votes, stats)
        else:
            position_counts[pos1] = read_counts
    return position_counts, stats


def count_cb_mutation_rows(
    bam,
    ref,
    vector_regions: Sequence[Region],
    mut_specs: Sequence[MutationSpec],
    min_mapq: int,
    min_base_quality: int,
) -> List[dict]:
    """
    Count raw read support per cell barcode for requested vector mutations.

    This is intentionally read-level, not UMI-collapsed, because the requested
    output columns are Nbr_reads and Nbr_reads_mut. Reads must still carry both
    CB and UB tags so this output is only produced in the normal CB/UMI mode.
    """
    totals: Dict[Tuple[str, str], Dict[str, int]] = defaultdict(lambda: {"Nbr_reads": 0, "Nbr_reads_mut": 0})
    covered_labels = set()

    for spec in mut_specs:
        regions_for_spec = [r for r in vector_regions if r.start <= spec.pos1 <= r.end]
        if not regions_for_spec:
            log(f"WARNING: requested mutation {spec.label} is outside all matched vector contig ranges; no CB rows will be written for it")
            continue

        for r in regions_for_spec:
            actual_ref = reference_base(ref, r.chrom, spec.pos1)
            if actual_ref != "N" and actual_ref != spec.ref_base:
                log(
                    f"WARNING: requested mutation {spec.label} has reference base {spec.ref_base}, "
                    f"but {r.chrom}:{spec.pos1} in the reference FASTA is {actual_ref}; counting {spec.alt} observations anyway"
                )

            try:
                iterator = bam.pileup(
                    contig=r.chrom,
                    start=spec.pos1 - 1,
                    stop=spec.pos1,
                    truncate=True,
                    stepper="samtools",
                    ignore_overlaps=False,
                    ignore_orphans=False,
                    min_base_quality=0,
                    min_mapping_quality=min_mapq,
                    max_depth=1000000,
                )
            except Exception as e:
                log(f"WARNING: could not pileup requested mutation {spec.label} at {r.chrom}:{spec.pos1}: {e}")
                continue

            for col in iterator:
                pos1 = col.reference_pos + 1
                if pos1 != spec.pos1:
                    continue
                covered_labels.add(spec.label)
                for pread in col.pileups:
                    rec = pread.alignment
                    if rec.is_unmapped or rec.mapping_quality < min_mapq:
                        continue
                    key = get_cb_ub(rec)
                    if key is None:
                        continue
                    allele = allele_from_pread(pread, min_base_quality)
                    if allele is None:
                        continue
                    cb, _ub = key
                    row_key = (cb, spec.label)
                    totals[row_key]["Nbr_reads"] += 1
                    if allele == spec.alt:
                        totals[row_key]["Nbr_reads_mut"] += 1

    rows: List[dict] = []
    for spec in mut_specs:
        spec_keys = sorted(k for k in totals if k[1] == spec.label)
        for cb, mut_label in spec_keys:
            n_total = int(totals[(cb, mut_label)]["Nbr_reads"])
            n_mut = int(totals[(cb, mut_label)]["Nbr_reads_mut"])
            if n_total <= 0:
                continue
            rows.append({
                "CB": cb,
                "mut": mut_label,
                "Mutation": "yes" if n_mut > 0 else "no",
                "Nbr_reads": n_total,
                "Nbr_reads_mut": n_mut,
            })
        if spec.label not in covered_labels:
            log(f"WARNING: no pileup column was observed for requested mutation {spec.label}; output may contain no rows for this mutation")
    return rows


def write_cb_mutation_tsv(rows: Sequence[dict], path: Path) -> None:
    fields = ["CB", "mut", "Mutation", "Nbr_reads", "Nbr_reads_mut"]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for row in rows:
            w.writerow(row)


def binom_sf(k: int, n: int, p: float) -> float:
    if k <= 0:
        return 1.0
    if n <= 0 or k > n:
        return 0.0
    if p <= 0:
        return 0.0 if k > 0 else 1.0
    if p >= 1:
        return 1.0
    if scipy_binom is not None:
        return float(scipy_binom.sf(k - 1, n, p))
    # Exact fallback: sum PMF from k..n using recurrence.
    log_term = (math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)
                + k * math.log(p) + (n - k) * math.log1p(-p))
    if log_term < -745:
        term = 0.0
    else:
        term = math.exp(log_term)
    total = term
    i = k
    odds = p / (1.0 - p)
    while i < n and term > 0.0:
        term *= ((n - i) / (i + 1)) * odds
        total += term
        i += 1
        if term < total * 1e-15:
            break
    return min(1.0, total)


def min_detectable_count(depth: int, error_rate: float, alpha_adj: float) -> Optional[int]:
    if depth <= 0:
        return None
    lo, hi = 1, depth
    answer = None
    while lo <= hi:
        mid = (lo + hi) // 2
        if binom_sf(mid, depth, error_rate) < alpha_adj:
            answer = mid
            hi = mid - 1
        else:
            lo = mid + 1
    return answer


def make_rows_for_regions(
    bam_path: Path,
    sample_id: str,
    patient_id: str,
    timepoint: str,
    ref,
    regions: Sequence[Region],
    position_counts_by_key: Dict[Tuple[str, int], Dict[str, int]],
    error_rate: float,
    base_alpha: float,
    n_tests: int,
    count_mode: str = "umi",
) -> List[dict]:
    rows: List[dict] = []
    # n_bases = number of positions in the output table, including zero-depth positions.
    n_bases = sum(max(0, r.end - r.start + 1) for r in regions)
    alpha_adj = base_alpha / (max(1, n_bases) * n_tests)
    minfreq_cache: Dict[int, Tuple[Optional[int], Optional[float]]] = {}

    for r in regions:
        for pos1 in range(r.start, r.end + 1):
            counts = position_counts_by_key.get((r.chrom, pos1), {a: 0 for a in ALLELES})
            depth = sum(int(counts.get(a, 0)) for a in ALLELES)
            ref_base = reference_base(ref, r.chrom, pos1)
            row = {
                "bam_file": bam_path.name,
                "sample_id": sample_id,
                "patient_id": patient_id,
                "timepoint": timepoint,
                "gene": r.gene,
                "chrom_ref": r.chrom,
                "chrom_pos": pos1,
                "reference": ref_base,
                "depth": depth,
                "count_mode": count_mode,
                "umi_depth": depth if count_mode == "umi" else "NA",
                "n_bases": n_bases,
                "n_tests": n_tests,
                "error_rate": error_rate,
                "base_alpha": base_alpha,
                "alpha_adjusted": alpha_adj,
            }
            for a in ALLELES:
                row[a] = int(counts.get(a, 0))
            for a in ALLELES:
                row[f"{a}_perc"] = (int(counts.get(a, 0)) / depth) if depth > 0 else 0.0
            if depth not in minfreq_cache:
                mdc = min_detectable_count(depth, error_rate, alpha_adj)
                minfreq_cache[depth] = (mdc, None if mdc is None else mdc / depth)
            row["min_sig_count"] = "NA" if minfreq_cache[depth][0] is None else minfreq_cache[depth][0]
            row["min_sig_alt_freq"] = "NA" if minfreq_cache[depth][1] is None else minfreq_cache[depth][1]
            for a in ALLELES:
                if a == ref_base:
                    row[f"{a}_pval"] = "NA"
                    row[f"{a}_sig"] = "FALSE"
                else:
                    k = int(counts.get(a, 0))
                    pv = binom_sf(k, depth, error_rate) if depth > 0 else 1.0
                    row[f"{a}_pval"] = pv
                    row[f"{a}_sig"] = "TRUE" if pv < alpha_adj else "FALSE"
            rows.append(row)
    return rows


def write_variant_tsv(rows: Sequence[dict], path: Path) -> None:
    fields = [
        "bam_file", "sample_id", "patient_id", "timepoint", "gene",
        "chrom_ref", "chrom_pos", "reference", "depth", "count_mode", "umi_depth",
        *ALLELES,
        *[f"{a}_perc" for a in ALLELES],
        *[f"{a}_pval" for a in ALLELES],
        *[f"{a}_sig" for a in ALLELES],
        "min_sig_count", "min_sig_alt_freq",
        "n_bases", "n_tests", "error_rate", "base_alpha", "alpha_adjusted",
    ]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for row in rows:
            w.writerow(row)


def vcf_escape_info_value(value) -> str:
    """Return a conservative VCF INFO-safe string value."""
    s = str(value)
    s = s.replace(" ", "_")
    # VCF INFO strings cannot safely contain semicolons, equals signs, commas, tabs, or newlines.
    return re.sub(r"[;=,\t\r\n]", "_", s)


def write_simple_vcf(rows: Sequence[dict], path: Path, sample_name: str, source_tsv: str) -> int:
    """
    Write a simple per-sample VCF containing significant alternative calls.

    This caller records SNV base counts and aggregate insertion/deletion event counts.
    Therefore:
      - A/C/G/T calls are written as standard single-base REF/ALT records.
      - Ins/Del calls are written as symbolic ALT records (<INS>/<DEL>) because the
        exact inserted/deleted sequence is not represented in the variant TSV.

    Returns the number of VCF records written.
    """
    records = []
    for row in rows:
        ref = str(row.get("reference", "N")).upper()
        if ref not in BASES:
            continue
        depth = int(float(row.get("depth", row.get("umi_depth", 0)) or 0))
        chrom = str(row.get("chrom_ref", row.get("gene", ".")))
        gene = str(row.get("gene", "."))
        try:
            pos = int(float(row.get("chrom_pos", 0)))
        except Exception:
            continue
        if pos <= 0:
            continue

        for allele in ALLELES:
            if allele == ref:
                continue
            if str(row.get(f"{allele}_sig", "FALSE")) != "TRUE":
                continue
            try:
                count = int(float(row.get(allele, 0) or 0))
            except Exception:
                count = 0
            try:
                freq = float(row.get(f"{allele}_perc", 0.0) or 0.0)
            except Exception:
                freq = 0.0
            pval = row.get(f"{allele}_pval", "NA")
            pval_str = "." if str(pval) in ("", "NA", "nan", "None") else str(pval)

            if allele in BASES:
                alt = allele
                event = "SNV"
            elif allele == "Ins":
                alt = "<INS>"
                event = "INS"
            elif allele == "Del":
                alt = "<DEL>"
                event = "DEL"
            else:
                continue

            records.append({
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "gene": gene,
                "event": event,
                "allele": allele,
                "count": count,
                "depth": depth,
                "freq": freq,
                "pval": pval_str,
                "patient_id": row.get("patient_id", "NA"),
                "timepoint": row.get("timepoint", "NA"),
                "sample_id": row.get("sample_id", sample_name),
            })

    records.sort(key=lambda x: (x["chrom"], x["pos"], x["alt"]))

    with open(path, "w", encoding="utf-8", newline="") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write("##source=detectVariants.py\n")
        out.write(f"##source_tsv={vcf_escape_info_value(source_tsv)}\n")
        out.write('##INFO=<ID=GENE,Number=1,Type=String,Description="Gene or region label">\n')
        out.write('##INFO=<ID=EVENT,Number=1,Type=String,Description="Event type from pileup table: SNV, INS, or DEL">\n')
        out.write('##INFO=<ID=ALLELE,Number=1,Type=String,Description="Original allele/event column from variant TSV">\n')
        out.write('##INFO=<ID=AO,Number=1,Type=Integer,Description="Alternate allele/event count used for calling">\n')
        out.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth used for calling; raw read depth when --umi-off, otherwise UMI depth">\n')
        out.write('##INFO=<ID=AF,Number=1,Type=Float,Description="Alternate allele/event frequency">\n')
        out.write('##INFO=<ID=PVAL,Number=1,Type=Float,Description="One-sided binomial p-value">\n')
        out.write('##INFO=<ID=PATIENT,Number=1,Type=String,Description="Patient ID parsed from BAM name">\n')
        out.write('##INFO=<ID=TIMEPOINT,Number=1,Type=String,Description="Timepoint parsed from BAM name">\n')
        out.write('##ALT=<ID=INS,Description="Insertion event; exact inserted sequence is not represented in the pileup table">\n')
        out.write('##ALT=<ID=DEL,Description="Deletion event; exact deleted sequence is not represented in the pileup table">\n')
        out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype placeholder for presence of significant event">\n')
        out.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth used for calling; raw read depth when --umi-off, otherwise UMI depth">\n')
        out.write('##FORMAT=<ID=AO,Number=1,Type=Integer,Description="Alternate allele/event count used for calling">\n')
        out.write('##FORMAT=<ID=AF,Number=1,Type=Float,Description="Alternate allele/event frequency">\n')
        out.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{vcf_escape_info_value(sample_name)}\n")

        for rec in records:
            info = ";".join([
                f"GENE={vcf_escape_info_value(rec['gene'])}",
                f"EVENT={rec['event']}",
                f"ALLELE={rec['allele']}",
                f"AO={rec['count']}",
                f"DP={rec['depth']}",
                f"AF={rec['freq']:.8g}",
                f"PVAL={rec['pval']}",
                f"PATIENT={vcf_escape_info_value(rec['patient_id'])}",
                f"TIMEPOINT={vcf_escape_info_value(rec['timepoint'])}",
            ])
            fmt = "GT:DP:AO:AF"
            sample_field = f"0/1:{rec['depth']}:{rec['count']}:{rec['freq']:.8g}"
            out.write(
                f"{rec['chrom']}\t{rec['pos']}\t.\t{rec['ref']}\t{rec['alt']}\t.\tPASS\t"
                f"{info}\t{fmt}\t{sample_field}\n"
            )

    return len(records)



def any_alt_sig(row: dict) -> bool:
    ref = str(row["reference"])
    for a in ALLELES:
        if a != ref and str(row.get(f"{a}_sig", "FALSE")) == "TRUE":
            return True
    return False


def max_alt_perc(row: dict) -> Tuple[float, str]:
    ref = str(row["reference"])
    best_a, best_v = "", 0.0
    for a in ALLELES:
        if a == ref:
            continue
        v = float(row.get(f"{a}_perc", 0.0) or 0.0)
        if v > best_v:
            best_v, best_a = v, a
    return best_v, best_a


def read_variant_rows(path: Path) -> List[dict]:
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_r_plot_helper(script_path: Path) -> None:
    r_code = r'''
suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: plot_detectVariants.R <variant.tsv> <out_prefix> <title> <control TRUE/FALSE>")
}

tsv <- args[[1]]
out_prefix <- args[[2]]
plot_title <- args[[3]]
control <- tolower(args[[4]]) %in% c("true", "1", "yes", "y")

alleles <- c("A", "C", "G", "T", "Ins", "Del")
base_colors <- c(
  "A"   = "#4DAF4A",
  "C"   = "#377EB8",
  "Del" = "#808080",
  "G"   = "#000000",
  "Ins" = "#F781BF",
  "T"   = "#E41A1C",
  "None" = "#D9D9D9"
)

d <- read.delim(tsv, stringsAsFactors = FALSE, check.names = FALSE)
if (nrow(d) == 0) {
  stop(paste("No rows in", tsv))
}

boolish <- function(x) {
  tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y")
}

get_ref <- function(row) {
  ref <- as.character(row[["reference"]])
  if (is.na(ref) || ref == "") ref <- "."
  ref
}

max_alt <- numeric(nrow(d))
max_alt_name <- character(nrow(d))
any_sig <- logical(nrow(d))

for (i in seq_len(nrow(d))) {
  ref <- get_ref(d[i, , drop = FALSE])
  best_v <- 0
  best_a <- "None"
  sig <- FALSE

  for (a in alleles) {
    if (a == ref) next
    pcol <- paste0(a, "_perc")
    scol <- paste0(a, "_sig")
    v <- 0
    if (pcol %in% names(d)) {
      v <- suppressWarnings(as.numeric(d[[pcol]][i]))
      if (is.na(v)) v <- 0
    }
    if (v > best_v) {
      best_v <- v
      best_a <- a
    }
    if (scol %in% names(d) && boolish(d[[scol]][i])) {
      sig <- TRUE
    }
  }

  max_alt[i] <- best_v * 100
  max_alt_name[i] <- best_a
  any_sig[i] <- sig
}

d$max_alt_percent <- max_alt
d$max_alt_name <- factor(max_alt_name, levels = c("A", "C", "G", "T", "Ins", "Del", "None"))
d$any_alt_sig <- any_sig
if ("depth" %in% names(d)) {
  d$depth_num <- suppressWarnings(as.numeric(d$depth))
} else {
  d$depth_num <- suppressWarnings(as.numeric(d$umi_depth))
}
count_mode_label <- if ("count_mode" %in% names(d)) unique(as.character(d$count_mode))[1] else "umi"
depth_axis_label <- ifelse(count_mode_label == "read", "Read depth", "UMI depth")
d$min_sig_alt_percent <- suppressWarnings(as.numeric(d$min_sig_alt_freq)) * 100
d$chrom_pos_num <- suppressWarnings(as.numeric(d$chrom_pos))
if (all(is.na(d$chrom_pos_num))) d$chrom_pos_num <- seq_len(nrow(d))

if (control && "gene" %in% names(d)) {
  d$block <- as.character(d$gene)
} else {
  d$block <- as.character(d$chrom_ref)
}

d$x_visual <- NA_real_
centers <- data.frame(block = character(), center = numeric(), stringsAsFactors = FALSE)
exon_rects <- data.frame(block = character(), exon = integer(), xmin = numeric(), xmax = numeric(), start = numeric(), end = numeric(), stringsAsFactors = FALSE)
axis_breaks <- numeric()
axis_labels <- character()
boundaries <- numeric()

if (control) {
  offset <- 0
  for (b in unique(d$block)) {
    idx_block <- which(d$block == b)
    idx_block <- idx_block[order(d$chrom_pos_num[idx_block])]
    block_x <- numeric()
    exon_start <- 1
    exon_id <- 0

    while (exon_start <= length(idx_block)) {
      exon_end <- exon_start
      while (exon_end < length(idx_block)) {
        p1 <- d$chrom_pos_num[idx_block[exon_end]]
        p2 <- d$chrom_pos_num[idx_block[exon_end + 1]]
        if (is.na(p1) || is.na(p2) || p2 != p1 + 1) break
        exon_end <- exon_end + 1
      }

      exon_idx <- idx_block[exon_start:exon_end]
      exon_id <- exon_id + 1
      x <- seq_len(length(exon_idx)) + offset
      d$x_visual[exon_idx] <- x
      block_x <- c(block_x, x)

      exon_rects <- rbind(exon_rects, data.frame(
        block = b,
        exon = exon_id,
        xmin = min(x) - 0.5,
        xmax = max(x) + 0.5,
        start = min(d$chrom_pos_num[exon_idx], na.rm = TRUE),
        end = max(d$chrom_pos_num[exon_idx], na.rm = TRUE),
        stringsAsFactors = FALSE
      ))
      axis_breaks <- c(axis_breaks, min(x), max(x))
      axis_labels <- c(axis_labels, as.character(min(d$chrom_pos_num[exon_idx], na.rm = TRUE)), as.character(max(d$chrom_pos_num[exon_idx], na.rm = TRUE)))

      offset <- max(x) + 6
      exon_start <- exon_end + 1
    }

    if (length(block_x) > 0) {
      centers <- rbind(centers, data.frame(block = b, center = mean(range(block_x)), stringsAsFactors = FALSE))
      boundaries <- c(boundaries, max(block_x) + 12)
      offset <- max(block_x) + 25
    }
  }
} else {
  d <- d[order(d$chrom_pos_num), , drop = FALSE]
  d$x_visual <- d$chrom_pos_num
  axis_breaks <- pretty(d$chrom_pos_num, n = 12)
  axis_breaks <- axis_breaks[axis_breaks >= min(d$chrom_pos_num, na.rm = TRUE) & axis_breaks <= max(d$chrom_pos_num, na.rm = TRUE)]
  axis_labels <- as.character(axis_breaks)
  centers <- data.frame(block = unique(d$block), center = mean(range(d$x_visual, na.rm = TRUE)), stringsAsFactors = FALSE)
}

star_df <- d[d$any_alt_sig, , drop = FALSE]
if (nrow(star_df) > 0) {
  star_df$star_y <- 98
}

base_theme <- theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

p1 <- ggplot(d, aes(x = x_visual, y = max_alt_percent, fill = max_alt_name)) +
  geom_col(width = 1) +
  scale_fill_manual(values = base_colors, drop = FALSE) +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.02))) +
  labs(title = plot_title, x = NULL, y = "Max alt allele (%)", fill = "Best alt") +
  base_theme +
  theme(legend.position = "top")

if (nrow(star_df) > 0) {
  p1 <- p1 + geom_text(data = star_df, aes(x = x_visual, y = star_y, color = max_alt_name), label = "*", size = 2.7, inherit.aes = FALSE, show.legend = FALSE) +
    scale_color_manual(values = base_colors, drop = FALSE)
}

p2 <- ggplot(d, aes(x = x_visual, y = min_sig_alt_percent)) +
  geom_line(linewidth = 0.35, na.rm = TRUE) +
  geom_point(size = 0.35, na.rm = TRUE) +
  labs(x = NULL, y = "Min significant alt (%)") +
  base_theme

p3_x_label <- if (control) NULL else "Reference position"
p3_bottom_margin <- if (control) 28 else 8

p3 <- ggplot(d, aes(x = x_visual, y = depth_num)) +
  geom_col(width = 1, fill = "grey35") +
  labs(x = p3_x_label, y = depth_axis_label) +
  base_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks.x = element_line(),
    plot.margin = margin(t = 5.5, r = 5.5, b = p3_bottom_margin, l = 5.5)
  ) +
  scale_x_continuous(breaks = axis_breaks, labels = axis_labels)

if (length(boundaries) > 1) {
  for (cx in boundaries[-length(boundaries)]) {
    p1 <- p1 + geom_vline(xintercept = cx, alpha = 0.2, linewidth = 0.25, linetype = "dashed")
    p2 <- p2 + geom_vline(xintercept = cx, alpha = 0.2, linewidth = 0.25, linetype = "dashed")
    p3 <- p3 + geom_vline(xintercept = cx, alpha = 0.2, linewidth = 0.25, linetype = "dashed")
  }
}

p4 <- NULL
if (control && nrow(exon_rects) > 0) {
  p4 <- ggplot() +
    geom_rect(data = exon_rects, aes(xmin = xmin, xmax = xmax, ymin = 0.25, ymax = 0.75), fill = "grey75", color = "grey25", linewidth = 0.25) +
    geom_text(data = centers, aes(x = center, y = 1.05, label = block), fontface = "bold", size = 3.2) +
    scale_x_continuous(breaks = axis_breaks, labels = axis_labels) +
    coord_cartesian(ylim = c(0, 1.25), clip = "off") +
    labs(x = "Exon diagram; tick labels are actual genomic positions", y = NULL) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.margin = margin(t = 28, r = 5.5, b = 5.5, l = 5.5)
    )

  if (length(boundaries) > 1) {
    for (cx in boundaries[-length(boundaries)]) {
      p4 <- p4 + geom_vline(xintercept = cx, alpha = 0.2, linewidth = 0.25, linetype = "dashed")
    }
  }
}

draw_panels <- function() {
  if (control && !is.null(p4)) {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(5, 1, heights = unit(c(2.1, 1.55, 1.55, 0.35, 1.0), "null"))))
    print(p1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
    print(p2, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
    print(p3, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
    print(p4, vp = viewport(layout.pos.row = 5, layout.pos.col = 1))
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(3, 1, heights = unit(c(2.1, 1.55, 1.55), "null"))))
    print(p1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
    print(p2, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
    print(p3, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
  }
}

png_file <- paste0(out_prefix, ".png")
png_opened <- FALSE
plot_height <- ifelse(control && !is.null(p4), 10, 8)

if (requireNamespace("ragg", quietly = TRUE)) {
  ragg::agg_png(filename = png_file, width = 18, height = plot_height, units = "in", res = 200)
  png_opened <- TRUE
} else if (isTRUE(capabilities("cairo"))) {
  png(filename = png_file, width = 18, height = plot_height, units = "in", res = 200, type = "cairo")
  png_opened <- TRUE
} else {
  warning("No headless-safe PNG device available: install R package ragg or use an R build with cairo support. Skipping PNG output.")
}

if (png_opened) {
  draw_panels()
  dev.off()
}
'''
    script_path.write_text(r_code)


def plot_variant_tsv(tsv: Path, out_prefix: Path, title: str, control: bool = False) -> None:
    rows = read_variant_rows(tsv)
    if not rows:
        log(f"WARNING: no rows in {tsv}; skipping plot")
        return

    rscript = shutil.which("Rscript")
    if rscript is None:
        log("WARNING: Rscript is unavailable; skipping plot")
        return

    helper = out_prefix.parent / "plot_detectVariants.R"
    write_r_plot_helper(helper)

    cmd = [
        rscript,
        str(helper),
        str(tsv),
        str(out_prefix),
        title,
        "TRUE" if control else "FALSE",
    ]
    log("[cmd] " + " ".join(cmd))
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.stdout.strip():
        log("[R plot stdout]\n" + result.stdout.strip())
    if result.stderr.strip():
        log("[R plot stderr]\n" + result.stderr.strip())
    if result.returncode != 0:
        raise RuntimeError("R plotting failed for %s\n%s" % (tsv, result.stderr))
    log(f"Wrote plot file: {out_prefix}.png")



BASE_COLORS_HTML = {
    "A": "#4DAF4A",
    "C": "#377EB8",
    "Del": "#808080",
    "G": "#000000",
    "Ins": "#F781BF",
    "T": "#E41A1C",
    "None": "#D9D9D9",
}


class ReferenceContextAccessor:
    def __init__(self, ref_path: Path):
        if pysam is None:
            raise RuntimeError("pysam is required for interactive HTML sequence context")
        self.ref = pysam.FastaFile(str(ref_path))

    def fetch_context(self, chrom: str, pos1: int, ref_base: str, flank: int) -> str:
        ref_base = (ref_base or "N").upper()
        if not chrom or pos1 <= 0:
            return "N" * flank + f"[{ref_base}]" + "N" * flank
        try:
            start0 = max(0, pos1 - 1 - flank)
            end0 = pos1 + flank
            seq = self.ref.fetch(chrom, start0, end0).upper()
            left_pad = max(0, flank - (pos1 - 1 - start0))
            seq = "N" * left_pad + seq
            need = flank * 2 + 1
            seq = seq + "N" * max(0, need - len(seq))
            seq = seq[:need]
            left = seq[:flank]
            center = seq[flank] if len(seq) > flank and seq[flank] in BASES else ref_base
            right = seq[flank + 1:]
            return f"{left}[{center}]{right}"
        except Exception:
            return "N" * flank + f"[{ref_base}]" + "N" * flank

    def close(self) -> None:
        try:
            self.ref.close()
        except Exception:
            pass


def context_plain(ctx: str) -> str:
    return re.sub(r"[^ACGTNacgtn]", "", ctx).upper()


def float_or_zero(x) -> float:
    try:
        if x is None or str(x) in ("", "NA", "nan", "None"):
            return 0.0
        v = float(x)
        if math.isnan(v):
            return 0.0
        return v
    except Exception:
        return 0.0


def int_or_zero(x) -> int:
    try:
        return int(float(x))
    except Exception:
        return 0


def read_variant_tsv_for_html(tsv: Path, ref_accessor: ReferenceContextAccessor, flank: int) -> List[dict]:
    rows: List[dict] = []
    with open(tsv, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for idx, row in enumerate(reader, start=1):
            ref_base = str(row.get("reference", "N")).upper()
            chrom = str(row.get("chrom_ref", ""))
            pos = int_or_zero(row.get("chrom_pos")) or idx
            depth = int_or_zero(row.get("depth", row.get("umi_depth", 0)))
            counts = {a: int_or_zero(row.get(a, 0)) for a in ALLELES}
            percs = {a: float_or_zero(row.get(f"{a}_perc", 0.0)) * 100.0 for a in ALLELES}
            pvals = {a: row.get(f"{a}_pval", "NA") for a in ALLELES}
            sigs = {a: str(row.get(f"{a}_sig", "FALSE")).strip().lower() in {"true", "t", "1", "yes", "y"} for a in ALLELES}
            best_a = "None"
            best_v = 0.0
            for a in ALLELES:
                if a == ref_base:
                    continue
                if percs[a] > best_v:
                    best_a = a
                    best_v = percs[a]
            sig_alleles = [a for a in ALLELES if a != ref_base and sigs[a]]
            ctx = ref_accessor.fetch_context(chrom, pos, ref_base, flank)
            rows.append({
                "gene": row.get("gene", ""), "chrom": chrom, "pos": pos, "reference": ref_base,
                "depth": depth, "count_mode": row.get("count_mode", ""), "counts": counts,
                "percs": percs, "pvals": pvals, "sigs": sigs, "best_alt": best_a,
                "best_alt_percent": best_v, "any_sig": bool(sig_alleles), "sig_alleles": sig_alleles,
                "min_sig_alt_percent": float_or_zero(row.get("min_sig_alt_freq", 0.0)) * 100.0,
                "min_sig_count": row.get("min_sig_count", "NA"), "alpha_adjusted": row.get("alpha_adjusted", ""),
                "error_rate": row.get("error_rate", ""), "context": ctx, "context_plain": context_plain(ctx),
            })
    return rows


def build_interactive_svg(rows: List[dict], px: float, max_width: int) -> str:
    n = len(rows)
    width = min(max_width, max(1200, int(n * px) + 120))
    scale_x = (width - 100) / max(1, n - 1)
    left = 60
    panel_h = 160
    gap = 32
    p1_y = 30
    p2_y = p1_y + panel_h + gap
    p3_y = p2_y + panel_h + gap
    total_h = p3_y + panel_h + 50
    max_depth = max([r["depth"] for r in rows] + [1])

    def x_for(i0):
        return left + i0 * scale_x

    def y_percent(base_y, val):
        val = max(0.0, min(100.0, val))
        return base_y + panel_h - (val / 100.0) * panel_h

    def y_depth(val):
        return p3_y + panel_h - (val / max_depth) * panel_h

    elems = [f'<svg id="variantSvg" viewBox="0 0 {width} {total_h}" data-base-height="{total_h}" preserveAspectRatio="none" xmlns="http://www.w3.org/2000/svg">']
    elems.append('<style>.axis{stroke:#333;stroke-width:1}.grid{stroke:#ddd;stroke-width:0.7}.bar{shape-rendering:crispEdges}.sigmark{stroke-width:1.2}.search-hit .halo{display:block}.halo{display:none;fill:none;stroke:#ffbf00;stroke-width:2}.sigflag{font-size:10px;fill:#111}.hiddenByFilter{display:none}.sig-hit{cursor:pointer}.pos{cursor:pointer}</style>')
    for label, y in [("Max alt allele (%)", p1_y), ("Min significant alt/event (%)", p2_y), ("Depth", p3_y)]:
        elems.append(f'<text x="8" y="{y + 14}" font-size="12" font-family="Arial" font-weight="bold">{html.escape(label)}</text>')
        elems.append(f'<line class="axis" x1="{left}" y1="{y + panel_h}" x2="{width - 20}" y2="{y + panel_h}"/>')
        elems.append(f'<line class="axis" x1="{left}" y1="{y}" x2="{left}" y2="{y + panel_h}"/>')
        for pct in [0, 25, 50, 75, 100]:
            yy = y_percent(y, pct)
            elems.append(f'<line class="grid" x1="{left}" y1="{yy:.2f}" x2="{width - 20}" y2="{yy:.2f}"/>')
            if label != "Depth":
                elems.append(f'<text x="30" y="{yy + 4:.2f}" font-size="10" font-family="Arial">{pct}</text>')
        if label == "Depth":
            for frac in [0, 0.25, 0.5, 0.75, 1.0]:
                yy = p3_y + panel_h - frac * panel_h
                elems.append(f'<text x="22" y="{yy + 4:.2f}" font-size="10" font-family="Arial">{int(max_depth * frac)}</text>')
    bar_w = max(1.0, min(3.0, scale_x * 0.85))
    for j, r in enumerate(rows):
        x = x_for(j)
        cls = "pos sig-pos" if r["any_sig"] else "pos"
        data_ctx = html.escape(r["context_plain"])
        data_sig = "1" if r["any_sig"] else "0"
        elems.append(f'<g class="{cls}" data-row-index="{j}" data-context="{data_ctx}" data-sig="{data_sig}" data-pos="{r["pos"]}">')
        elems.append(f'<rect class="halo" x="{x - bar_w:.2f}" y="{p1_y - 4}" width="{bar_w * 2:.2f}" height="{p3_y + panel_h - p1_y + 8}"/>')
        color = BASE_COLORS_HTML.get(r["best_alt"], BASE_COLORS_HTML["None"])
        y = y_percent(p1_y, r["best_alt_percent"])
        elems.append(f'<rect class="bar" x="{x - bar_w/2:.2f}" y="{y:.2f}" width="{bar_w:.2f}" height="{p1_y + panel_h - y:.2f}" fill="{color}"/>')
        for k, a in enumerate(r["sig_alleles"]):
            av = r["percs"].get(a, 0.0)
            ay = y_percent(p1_y, av)
            ax = x + (k - (len(r["sig_alleles"]) - 1)/2.0) * max(1.3, bar_w)
            acol = BASE_COLORS_HTML.get(a, "#111")
            elems.append(f'<g class="sig-hit" data-row-index="{j}" data-allele="{a}" data-mutation="{r["reference"]}>{a}">')
            elems.append(f'<line class="sigmark" x1="{ax:.2f}" y1="{p1_y + panel_h:.2f}" x2="{ax:.2f}" y2="{ay:.2f}" stroke="{acol}"/>')
            elems.append(f'<circle cx="{ax:.2f}" cy="{ay:.2f}" r="1.7" fill="{acol}"/>')
            sy = max(p1_y + 10, ay - 5)
            elems.append(f'<text class="sigflag" x="{ax - 2:.2f}" y="{sy:.2f}">*</text>')
            elems.append(f'<rect x="{ax - 4:.2f}" y="{min(ay, p1_y + panel_h) - 2:.2f}" width="8" height="{abs((p1_y + panel_h) - ay) + 8:.2f}" fill="transparent"/>')
            elems.append('</g>')
        mv = r["min_sig_alt_percent"]
        if mv >= 0:
            y2 = y_percent(p2_y, min(mv, 100.0))
            elems.append(f'<rect class="bar" x="{x - bar_w/2:.2f}" y="{y2:.2f}" width="{bar_w:.2f}" height="{p2_y + panel_h - y2:.2f}" fill="#555" opacity="0.7"/>')
        yd = y_depth(r["depth"])
        elems.append(f'<rect class="bar" x="{x - bar_w/2:.2f}" y="{yd:.2f}" width="{bar_w:.2f}" height="{p3_y + panel_h - yd:.2f}" fill="#555"/>')
        elems.append('</g>')
    if rows:
        step = max(1, len(rows) // 12)
        for j in range(0, len(rows), step):
            x = x_for(j)
            elems.append(f'<text x="{x:.2f}" y="{p3_y + panel_h + 18}" font-size="10" font-family="Arial" transform="rotate(45 {x:.2f},{p3_y + panel_h + 18})">{html.escape(str(rows[j]["pos"]))}</text>')
    elems.append('</svg>')
    return "\n".join(elems)


def significant_table_html(rows: List[dict]) -> str:
    sigs = []
    for r in rows:
        for a in r["sig_alleles"]:
            sigs.append({"gene": r["gene"], "chrom": r["chrom"], "pos": r["pos"], "ref": r["reference"], "allele": a, "mutation": f"{r['reference']}>{a}", "count": r["counts"].get(a, 0), "percent": r["percs"].get(a, 0.0), "depth": r["depth"], "pval": r["pvals"].get(a, "NA"), "context": r["context"]})
    headers = ["gene", "chrom", "pos", "ref", "allele", "mutation", "count", "percent", "depth", "pval", "context"]
    parts = ['<table id="sigTable"><thead><tr>']
    for h in headers:
        parts.append(f'<th onclick="sortTable(\'{h}\')">{html.escape(h)}</th>')
    parts.append('</tr></thead><tbody>')
    for r in sigs:
        parts.append('<tr>')
        for h in headers:
            val = f"{r[h]:.5f}" if h == "percent" else r[h]
            parts.append(f'<td>{html.escape(str(val))}</td>')
        parts.append('</tr>')
    parts.append('</tbody></table>')
    return "\n".join(parts)


def write_interactive_html(rows: List[dict], out: Path, title: str, svg: str, source_tsv: Path, reference: Path) -> None:
    n_sig_pos = sum(1 for r in rows if r["any_sig"])
    n_sig_calls = sum(len(r["sig_alleles"]) for r in rows)
    legend = " ".join(f'<span class="swatch" style="background:{c}"></span>{html.escape(a)}' for a, c in BASE_COLORS_HTML.items() if a != "None")
    hover_rows = []
    for r in rows:
        hover_rows.append({"gene": r["gene"], "chrom": r["chrom"], "pos": r["pos"], "reference": r["reference"], "depth": r["depth"], "count_mode": r["count_mode"], "best_alt": r["best_alt"], "best_alt_percent": r["best_alt_percent"], "sig_alleles": r["sig_alleles"], "counts": r["counts"], "percs": r["percs"], "pvals": r["pvals"], "sigs": r["sigs"], "context": r["context"], "context_plain": r["context_plain"], "min_sig_count": r["min_sig_count"], "min_sig_alt_percent": r["min_sig_alt_percent"], "alpha_adjusted": r["alpha_adjusted"], "error_rate": r["error_rate"]})
    data_json = json.dumps({"n": len(rows), "rows": hover_rows}, separators=(",", ":"))
    table_html = significant_table_html(rows)
    template = """<!doctype html>
<html lang=\"en\"><head><meta charset=\"utf-8\"><title>__TITLE__</title>
<style>
body{font-family:Arial,sans-serif;margin:18px;color:#222} h1{font-size:20px;margin-bottom:4px}.meta{color:#555;font-size:13px;margin-bottom:12px}.controls{position:sticky;top:0;background:white;z-index:10;padding:10px 0;border-bottom:1px solid #ddd}.controls input[type=text]{width:360px;padding:5px}.legend{margin-top:8px;font-size:13px}.swatch{display:inline-block;width:12px;height:12px;margin-left:12px;margin-right:4px;border:1px solid #aaa;vertical-align:-1px}.plot-wrap{overflow-x:auto;border:1px solid #ccc;padding:8px;margin-top:12px}.note{font-size:13px;color:#555}.zoom-controls{margin-top:8px;display:flex;align-items:center;gap:8px;flex-wrap:wrap}.zoom-controls button{padding:2px 8px}#zoomPct{min-width:56px;display:inline-block;font-variant-numeric:tabular-nums}#hoverPanel{margin-top:12px;border:1px solid #ddd;background:#fafafa;padding:10px 12px;font-size:12px}#hoverTitle{font-weight:bold;margin-bottom:4px}#hoverMutation{color:#444;margin-bottom:6px}.hover-grid{display:grid;grid-template-columns:repeat(2,minmax(220px,1fr));gap:8px 18px}.hover-mono{font-family:Menlo,Consolas,monospace}.hover-call-table{border-collapse:collapse;width:100%;margin-top:8px;font-size:12px}.hover-call-table th,.hover-call-table td{border:1px solid #ddd;padding:3px 5px;text-align:left}.hover-call-table th{background:#f2f2f2}table{border-collapse:collapse;font-size:12px;margin-top:12px;width:100%}th,td{border:1px solid #ddd;padding:4px 6px}th{background:#f2f2f2;cursor:pointer;position:sticky;top:48px}tr:nth-child(even){background:#fafafa}.hiddenByFilter{display:none}
</style></head><body>
<h1>__TITLE__</h1><div class=\"meta\">Source TSV: __SOURCE_TSV__<br>Reference FASTA: __REFERENCE__<br>Rows: __N_ROWS__; significant positions: __N_SIG_POS__; significant allele/event calls: __N_SIG_CALLS__</div>
<div class=\"controls\"><label>Search 11-mer/context: <input id=\"seqSearch\" type=\"text\" placeholder=\"e.g. ACACGCTGCAC or CACG[C]TG\" oninput=\"applyFilters()\"></label><label style=\"margin-left:12px\"><input id=\"sigOnly\" type=\"checkbox\" onchange=\"applyFilters()\"> show only significant positions</label><button onclick=\"clearFilters()\">Clear</button><div class=\"zoom-controls\"><strong>Horizontal zoom:</strong><button type=\"button\" onclick=\"zoomOut()\">-</button><button type=\"button\" onclick=\"zoomReset()\">Reset</button><button type=\"button\" onclick=\"zoomIn()\">+</button><input id=\"zoomSlider\" type=\"range\" min=\"1\" max=\"30\" step=\"0.1\" value=\"1\" oninput=\"setZoom(parseFloat(this.value))\"><span id=\"zoomPct\">100%</span><span class=\"note\">Default view shows the full region. Zoom changes width only; height and vertical alignment stay fixed.</span></div><div class=\"legend\">__LEGEND__</div><div class=\"note\">Click any bar or significance marker to pin the exact mutation, all six call counts, p-values, significance, and the 11-mer context. Yellow outline indicates sequence search match.</div></div>
<div class=\"plot-wrap\" id=\"plotWrap\">__SVG__</div><div id=\"hoverPanel\"><div id=\"hoverTitle\">Click a bar to inspect it</div><div id=\"hoverMutation\">Focused mutation: none</div><div id=\"hoverBody\" class=\"note\">Click a position bar or significance marker to show the full 6-call support at that position plus sequence context.</div></div>
<h2>Significant allele/event calls</h2><div class=\"note\">Click column headers to sort. Context can be copied and searched above.</div>__SIG_TABLE__
<script>
const META=__DATA_JSON__;
function normSeq(s){return(s||'').toUpperCase().replace(/[^ACGTN]/g,'');}
function mutationLabel(ref,allele){return allele?`${ref}>${allele}`:'none';}
function escapeHtml(s){return String(s).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;').replace(/\"/g,'&quot;').replace(/'/g,'&#39;');}
function callTable(row,focusAllele){const order=['A','C','G','T','Ins','Del'];let out='<table class=\"hover-call-table\"><thead><tr><th>Call</th><th>Mutation</th><th>Count</th><th>%</th><th>p-value</th><th>Significant</th></tr></thead><tbody>';order.forEach(a=>{const mutation=mutationLabel(row.reference,a);const isFocus=focusAllele&&a===focusAllele;const sig=row.sigs&&row.sigs[a]?'yes':'no';const count=row.counts&&row.counts[a]!=null?row.counts[a]:0;const perc=row.percs&&row.percs[a]!=null?Number(row.percs[a]).toFixed(4):'0.0000';const pval=row.pvals&&row.pvals[a]!=null?row.pvals[a]:'NA';out+=`<tr${isFocus?' style=\"background:#fff8d6\"':''}><td>${escapeHtml(a)}</td><td>${escapeHtml(mutation)}</td><td>${escapeHtml(count)}</td><td>${escapeHtml(perc)}</td><td>${escapeHtml(pval)}</td><td>${escapeHtml(sig)}</td></tr>`;});out+='</tbody></table>';return out;}
function renderInfo(rowIndex,focusAllele){const row=META.rows[rowIndex];if(!row)return;const focusMutation=focusAllele?mutationLabel(row.reference,focusAllele):(row.sig_alleles&&row.sig_alleles.length?row.sig_alleles.map(a=>mutationLabel(row.reference,a)).join(', '):mutationLabel(row.reference,row.best_alt));document.getElementById('hoverTitle').textContent=`${row.gene}  ${row.chrom}:${row.pos}`;document.getElementById('hoverMutation').textContent=`Focused mutation: ${focusMutation}`;const sigList=row.sig_alleles&&row.sig_alleles.length?row.sig_alleles.map(a=>mutationLabel(row.reference,a)).join(', '):'none';const body=`<div class=\"hover-grid\"><div><strong>Reference:</strong> ${escapeHtml(row.reference)}<br><strong>Depth:</strong> ${escapeHtml(row.depth)} (${escapeHtml(row.count_mode)})<br><strong>Best alt:</strong> ${escapeHtml(row.best_alt)} (${Number(row.best_alt_percent||0).toFixed(4)}%)</div><div><strong>Context:</strong> <span class=\"hover-mono\">${escapeHtml(row.context)}</span><br><strong>Significant mutation(s):</strong> ${escapeHtml(sigList)}<br><strong>Min significant alt/event:</strong> ${Number(row.min_sig_alt_percent||0).toFixed(4)}%</div></div><div class=\"note\" style=\"margin-top:8px\">All six calls shown below simultaneously for this position.</div>${callTable(row,focusAllele)}`;document.getElementById('hoverBody').innerHTML=body;}
function applyFilters(){const q=normSeq(document.getElementById('seqSearch').value);const sigOnly=document.getElementById('sigOnly').checked;document.querySelectorAll('#variantSvg .pos').forEach(g=>{const ctx=g.getAttribute('data-context')||'';const sig=g.getAttribute('data-sig')==='1';const match=q.length===0||ctx.includes(q);const visible=match&&(!sigOnly||sig);g.classList.toggle('hiddenByFilter',!visible);g.classList.toggle('search-hit',q.length>0&&match);});}
function clearFilters(){document.getElementById('seqSearch').value='';document.getElementById('sigOnly').checked=false;applyFilters();}
let zoomFactor=1.0;function applyZoom(){const svg=document.getElementById('variantSvg');if(!svg)return;const baseHeight=parseFloat(svg.getAttribute('data-base-height'))||svg.viewBox.baseVal.height||svg.height.baseVal.value;svg.style.width=Math.round(zoomFactor*100)+'%';svg.style.height=baseHeight+'px';document.getElementById('zoomPct').textContent=Math.round(zoomFactor*100)+'%';document.getElementById('zoomSlider').value=zoomFactor;}
function setZoom(v){zoomFactor=Math.min(30,Math.max(1,v||1));applyZoom();}function zoomIn(){setZoom(zoomFactor*1.25);}function zoomOut(){setZoom(zoomFactor/1.25);}function zoomReset(){setZoom(1.0);}
function initClick(){const svg=document.getElementById('variantSvg');if(!svg)return;svg.addEventListener('click',ev=>{const sig=ev.target.closest('.sig-hit');if(sig){const idx=parseInt(sig.getAttribute('data-row-index')||'-1',10);const allele=sig.getAttribute('data-allele')||'';if(idx>=0)renderInfo(idx,allele);return;}const pos=ev.target.closest('.pos');if(pos){const idx=parseInt(pos.getAttribute('data-row-index')||'-1',10);if(idx>=0)renderInfo(idx,null);}});}
function sortTable(key){const table=document.getElementById('sigTable');const headers=Array.from(table.querySelectorAll('th')).map(th=>th.textContent);const idx=headers.indexOf(key);const tbody=table.querySelector('tbody');const rows=Array.from(tbody.querySelectorAll('tr'));const numeric=['pos','count','percent','depth'].includes(key);const asc=table.getAttribute('data-sort-key')!==key||table.getAttribute('data-sort-dir')==='desc';rows.sort((a,b)=>{let av=a.children[idx].textContent;let bv=b.children[idx].textContent;if(numeric){av=parseFloat(av)||0;bv=parseFloat(bv)||0;}if(av<bv)return asc?-1:1;if(av>bv)return asc?1:-1;return 0;});rows.forEach(r=>tbody.appendChild(r));table.setAttribute('data-sort-key',key);table.setAttribute('data-sort-dir',asc?'asc':'desc');}
document.addEventListener('DOMContentLoaded',()=>{applyFilters();initClick();applyZoom();});
</script></body></html>"""
    replacements = {"__TITLE__": html.escape(title), "__SOURCE_TSV__": html.escape(str(source_tsv)), "__REFERENCE__": html.escape(str(reference)), "__N_ROWS__": f"{len(rows):,}", "__N_SIG_POS__": f"{n_sig_pos:,}", "__N_SIG_CALLS__": f"{n_sig_calls:,}", "__LEGEND__": legend, "__SVG__": svg, "__SIG_TABLE__": table_html, "__DATA_JSON__": data_json}
    html_text = template
    for key, value in replacements.items():
        html_text = html_text.replace(key, value)
    out.write_text(html_text, encoding="utf-8")


def generate_interactive_html(tsv: Path, ref_path: Path, out_html: Path, title: str, flank: int, px_per_position: float, max_width: int) -> None:
    if pysam is None:
        log("WARNING: pysam is unavailable; skipping interactive HTML")
        return
    accessor = ReferenceContextAccessor(ref_path)
    try:
        rows = read_variant_tsv_for_html(tsv, accessor, flank)
    finally:
        accessor.close()
    if not rows:
        log(f"WARNING: no rows in {tsv}; skipping interactive HTML")
        return
    svg = build_interactive_svg(rows, px_per_position, max_width)
    write_interactive_html(rows, out_html, title, svg, tsv, ref_path)
    log(f"Wrote interactive HTML: {out_html}")


def process_bam(bam_path: Path, ref_path: Path, outdir: Path, args: argparse.Namespace, control_regions: Sequence[Region]) -> None:
    if pysam is None:
        raise RuntimeError("pysam is required")
    sample_id, patient_id, timepoint, bam_file = sample_name_from_bam(bam_path)
    log(f"Processing {bam_file} as sample_id={sample_id}, patient_id={patient_id}, timepoint={timepoint}")
    keywords = [x.strip() for x in args.vector_keywords.split(",") if x.strip()]

    ref = open_reference(ref_path)
    stats_summary = []
    use_umi = not args.umi_off
    count_mode = "umi" if use_umi else "read"
    if args.umi_off:
        log("[info] --umi-off enabled for this BAM: using raw read counts; each aligned read is treated as unique")
    else:
        log("[info] UMI collapsing enabled: reads without CB/UB tags will be excluded from allele counts")
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        vector_regions = []
        for chrom, length in zip(bam.references, bam.lengths):
            if matches_vector(chrom, keywords):
                vector_regions.append(Region(gene="vector", chrom=chrom, start=1, end=int(length), source="BAM header vector keyword"))
        if not vector_regions:
            raise RuntimeError(f"No vector-like contigs matched keywords in {bam_path}: {keywords}")

        if getattr(args, "write_region_bams", False):
            vector_bam = outdir / f"{sample_id}.vector_region.bam"
            subset_bam_to_regions(bam_path, vector_bam, vector_regions, "vector-region")

        # Count vector regions.
        vector_counts: Dict[Tuple[str, int], Dict[str, int]] = {}
        vector_stats = PileupStats()
        for r in vector_regions:
            log(f"Pileup vector {r.chrom}:{r.start}-{r.end}")
            pc, st = count_region_positions(bam, r, args.min_mapq, args.min_base_quality, use_umi=use_umi)
            vector_counts.update({(r.chrom, pos): counts for pos, counts in pc.items()})
            vector_stats.missing_cbub_reads += st.missing_cbub_reads
            vector_stats.ambiguous_umis += st.ambiguous_umis
            vector_stats.usable_umis += st.usable_umis
            vector_stats.reads_seen += st.reads_seen
        vector_rows = make_rows_for_regions(
            bam_path, sample_id, patient_id, timepoint, ref, vector_regions, vector_counts,
            args.error_rate, args.base_alpha, args.n_tests, count_mode=count_mode
        )
        vector_tsv = outdir / f"{sample_id}.vector.variant.tsv"
        write_variant_tsv(vector_rows, vector_tsv)
        log(f"Wrote {vector_tsv}")

        mut_specs = getattr(args, "mut_specs", [])
        if mut_specs and use_umi:
            cb_mut_rows = count_cb_mutation_rows(
                bam,
                ref,
                vector_regions,
                mut_specs,
                args.min_mapq,
                args.min_base_quality,
            )
            cb_mut_tsv = outdir / f"{sample_id}.vector.cb_mutation.tsv"
            write_cb_mutation_tsv(cb_mut_rows, cb_mut_tsv)
            log(f"Wrote {cb_mut_tsv} ({format_int(len(cb_mut_rows))} CB-mut rows)")
        elif mut_specs and not use_umi:
            log("WARNING: --mut was requested but --umi-off is enabled; skipping CB mutation table")

        if not args.skip_vcf:
            vector_vcf = outdir / f"{sample_id}.vector.significant_variants.vcf"
            n_vcf = write_simple_vcf(vector_rows, vector_vcf, sample_id, vector_tsv.name)
            log(f"Wrote {vector_vcf} ({format_int(n_vcf)} records)")
        stats_summary.append(("vector", vector_stats))

        if not args.skip_control:
            if getattr(args, "write_region_bams", False):
                control_bam = outdir / f"{sample_id}.control_exons_region.bam"
                subset_bam_to_regions(bam_path, control_bam, control_regions, "control-exons-region")

            control_counts: Dict[Tuple[str, int], Dict[str, int]] = {}
            control_stats = PileupStats()
            bam_refs = set(bam.references)
            for r in control_regions:
                if r.chrom not in bam_refs:
                    log(f"WARNING: skipping {r.gene} {r.chrom}:{r.start}-{r.end}; contig absent from BAM")
                    continue
                log(f"Pileup control {r.gene} {r.chrom}:{r.start}-{r.end}")
                pc, st = count_region_positions(bam, r, args.min_mapq, args.min_base_quality, use_umi=use_umi)
                control_counts.update({(r.chrom, pos): counts for pos, counts in pc.items()})
                control_stats.missing_cbub_reads += st.missing_cbub_reads
                control_stats.ambiguous_umis += st.ambiguous_umis
                control_stats.usable_umis += st.usable_umis
                control_stats.reads_seen += st.reads_seen
            control_rows = make_rows_for_regions(
                bam_path, sample_id, patient_id, timepoint, ref, control_regions, control_counts,
                args.error_rate, args.base_alpha, args.n_tests, count_mode=count_mode
            )
            control_tsv = outdir / f"{sample_id}.control_exons.variant.tsv"
            write_variant_tsv(control_rows, control_tsv)
            log(f"Wrote {control_tsv}")
            if not args.skip_vcf:
                control_vcf = outdir / f"{sample_id}.control_exons.significant_variants.vcf"
                n_vcf = write_simple_vcf(control_rows, control_vcf, sample_id, control_tsv.name)
                log(f"Wrote {control_vcf} ({format_int(n_vcf)} records)")
            stats_summary.append(("control_exons", control_stats))

    if not args.skip_plots:
        plot_variant_tsv(vector_tsv, outdir / f"{sample_id}.vector.variant_plot", f"{sample_id} vector variants", control=False)
        if not args.skip_control:
            plot_variant_tsv(control_tsv, outdir / f"{sample_id}.control_exons.variant_plot", f"{sample_id} control exon variants", control=True)

    if not getattr(args, "skip_html", False):
        generate_interactive_html(
            vector_tsv,
            ref_path,
            outdir / f"{sample_id}.vector.variant_interactive.html",
            f"{sample_id} vector interactive variants",
            getattr(args, "html_flank", 5),
            getattr(args, "html_px_per_position", 2.2),
            getattr(args, "html_max_width", 60000),
        )
        if not args.skip_control:
            generate_interactive_html(
                control_tsv,
                ref_path,
                outdir / f"{sample_id}.control_exons.variant_interactive.html",
                f"{sample_id} control exon interactive variants",
                getattr(args, "html_flank", 5),
                getattr(args, "html_px_per_position", 2.2),
                getattr(args, "html_max_width", 60000),
            )

    summary_path = outdir / f"{sample_id}.umi_collapse_summary.tsv"
    with open(summary_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample_id", "region_set", "count_mode", "reads_seen", "missing_cbub_reads", "usable_observations", "ambiguous_umis_excluded"])
        for region_set, st in stats_summary:
            w.writerow([sample_id, region_set, count_mode, st.reads_seen, st.missing_cbub_reads, st.usable_umis, st.ambiguous_umis])
    log(f"Wrote {summary_path}")


def main() -> None:
    global GLOBAL_LOGGER

    args = parse_args()
    args.mut_specs = parse_mutation_specs(getattr(args, "mut", []))
    outdir = Path(args.out).expanduser()
    outdir.mkdir(parents=True, exist_ok=True)
    ref_path = Path(args.reference).expanduser()
    bams = collect_bams(args)

    log_path = outdir / f"preflight_and_run_{timestamp()}.log"
    logger = Logger(log_path)
    GLOBAL_LOGGER = logger

    logger.section("detectVariants.py started")
    logger.write(f"[info] Command: {' '.join(sys.argv)}")
    logger.write(f"[info] Output directory: {outdir.resolve()}")
    logger.write(f"[info] Log path: {log_path}")

    with open(outdir / "run.command.txt", "w") as fh:
        fh.write(" ".join(map(str, sys.argv)) + "\n")
    logger.write(f"[ok] Wrote command log: {outdir / 'run.command.txt'}")

    parameters = vars(args).copy()
    parameters["mut_specs"] = [m.as_dict() for m in getattr(args, "mut_specs", [])]
    with open(outdir / "parameters.json", "w") as fh:
        json.dump(parameters, fh, indent=2, sort_keys=True)
    logger.write(f"[ok] Wrote parameter JSON: {outdir / 'parameters.json'}")

    rows = check_preflight(bams, ref_path, outdir, args, logger)
    has_fail = any(status == "FAIL" for _check, status, _details in rows)
    if args.preflight_only:
        if has_fail:
            raise SystemExit(2)
        return
    if has_fail:
        raise SystemExit("ERROR: preflight failed; see console output above")

    control_regions: List[Region] = []
    if not args.skip_control:
        logger.section("Control exon configuration")
        if args.control_exons_tsv:
            logger.write(f"[info] Loading control exons from TSV override: {args.control_exons_tsv}")
            control_regions = read_control_exons_tsv(Path(args.control_exons_tsv).expanduser())
        else:
            logger.write("[info] Using hard-coded hg38 control CDS-overlapping exon coordinates")
            control_regions = hardcoded_control_exons()
        write_control_regions(control_regions, outdir / "control_exons.used_regions.tsv")
        logger.write(f"[ok] Using {format_int(len(control_regions))} control exon intervals")
        logger.write(f"[ok] Wrote control exon manifest: {outdir / 'control_exons.used_regions.tsv'}")

    logger.section("Variant detection")
    for bam in bams:
        process_bam(bam, ref_path, outdir, args, control_regions)

    logger.section("detectVariants.py complete")
    logger.write(f"[done] Output directory: {outdir.resolve()}")
    logger.write(f"[done] Log file: {log_path}")


if __name__ == "__main__":
    main()
