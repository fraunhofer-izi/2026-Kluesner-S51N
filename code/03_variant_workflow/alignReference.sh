#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# alignReference.sh
#
# Usage:
#   Single-end/default FASTQ mode:
#     alignReference.sh READS.fastq[.gz] [MORE_READS.fastq[.gz] ...] reference.fasta
#     alignReference.sh READS.fastq[.gz] [MORE_READS.fastq[.gz] ...] -cat ref1.fasta ref2.fa [...]
#
#   Paired FASTQ mode, strict exactly two FASTQ inputs:
#     alignReference.sh --paired R1.fastq[.gz] R2.fastq[.gz] reference.fasta
#     alignReference.sh --paired R1.fastq[.gz] R2.fastq[.gz] -cat ref1.fasta ref2.fa [...]
#
#   Cell Ranger / alignment BAM mode, one alignment input per run:
#     alignReference.sh READS.bam|READS.cram|READS.sam reference.fasta
#     alignReference.sh READS.bam|READS.cram|READS.sam -cat ref1.fasta ref2.fa [...]
#
#   Paired alignment mode, one alignment input per run:
#     alignReference.sh --paired READS.bam|READS.cram|READS.sam reference.fasta
#     alignReference.sh --paired READS.bam|READS.cram|READS.sam -cat ref1.fasta ref2.fa [...]
#
# Backward-compatible form for a single reads file plus concatenated references:
#     alignReference.sh -cat READS.fastq[.gz]|READS.bam ref1.fasta ref2.fa [...]
#
# Behavior:
#   • If not running inside Slurm: submits itself via sbatch and exits
#   • Inside Slurm job:
#       - (Optional) concatenates FASTA references
#       - Indexes reference FASTA
#       - FASTQ input, default mode:
#           one or more FASTQs are concatenated and aligned as single-end reads
#       - FASTQ input, --paired mode:
#           exactly two FASTQs are aligned as paired-end reads
#       - BAM/CRAM/SAM input, default Cell Ranger mode:
#           name-collates, converts to four FASTQ bins (R1, R2, singleton,
#           other), concatenates every nonempty bin, and aligns as single-end
#       - BAM/CRAM/SAM input, --paired mode:
#           name-collates, converts to four FASTQ bins, then aligns only R1/R2
#           as paired-end reads; singleton/other bins are ignored
#       - Preserves selected BAM tags in FASTQ comments and back into minimap2
#         SAM/BAM output via minimap2 -y for BAM/CRAM/SAM inputs
#       - Sorts & indexes with samtools
#       - Writes QC summaries: flagstat, idxstats, minimap2 log, command log
#
# Important FASTQ-mode warnings:
#   • Default multi-FASTQ mode assumes all provided FASTQs are alignable reads
#     from the same sample, such as multiple lanes of R2 or true single-end data.
#   • Do not provide raw 10x I1/I2/R1/R2 files together unless you intentionally
#     want all of them aligned. For standard 10x GEX, R1 is usually cell barcode
#     + UMI, not biological insert sequence.
#   • FASTQ input mode does not create CB/UB BAM tags from raw 10x R1 barcode/UMI
#     reads. For barcode/UMI-aware remapping, prefer Cell Ranger BAM input.
#   • Use --paired only when both FASTQs contain biologically alignable insert
#     sequence. In --paired FASTQ mode, exactly two FASTQ files are required.
#
# Output directory:
#   <READS_DIR>/<reads_stem>__vs__<ref1>__<ref2>__.../
#
# Final output:
#   One sorted BAM per run:
#     <reads_stem>.vs.<ref_stem>.sorted.bam
###############################################################################

# -----------------------------
# Defaults (override via env)
# -----------------------------
PRESET="${PRESET:-sr}"          # sr | map-ont | map-hifi | map-pb
CPUS="${CPUS:-32}"
MEM="${MEM:-250G}"
TIME="${TIME:-200:00:00}"
PARTITION="${PARTITION:-}"
ACCOUNT="${ACCOUNT:-}"

# Tags copied from input BAM/CRAM/SAM into FASTQ comments, then back into
# minimap2 SAM/BAM output via minimap2 -y. For Cell Ranger BAMs, CB/UB are the
# key corrected barcode/UMI tags; CR/UR are raw; CY/UY are qualities.
PRESERVE_BAM_TAGS="${PRESERVE_BAM_TAGS:-CB,CR,CY,UB,UR,UY,RG,BC,QT}"

usage() {
  cat >&2 <<'USAGE'
Usage:
  # Default FASTQ mode: one or more alignable FASTQs, concatenated and aligned single-end
  alignReference.sh READS.fastq[.gz] [MORE_READS.fastq[.gz] ...] reference.fasta
  alignReference.sh READS.fastq[.gz] [MORE_READS.fastq[.gz] ...] -cat ref1.fasta ref2.fa [...]

  # Paired FASTQ mode: strict exactly two FASTQ inputs, aligned paired-end
  alignReference.sh --paired R1.fastq[.gz] R2.fastq[.gz] reference.fasta
  alignReference.sh --paired R1.fastq[.gz] R2.fastq[.gz] -cat ref1.fasta ref2.fa [...]

  # Cell Ranger / alignment BAM mode: one BAM/CRAM/SAM input per run
  alignReference.sh READS.bam|READS.cram|READS.sam reference.fasta
  alignReference.sh READS.bam|READS.cram|READS.sam -cat ref1.fasta ref2.fa [...]

  # Paired alignment mode: use only R1/R2 bins from the BAM-derived FASTQs
  alignReference.sh --paired READS.bam|READS.cram|READS.sam reference.fasta
  alignReference.sh --paired READS.bam|READS.cram|READS.sam -cat ref1.fasta ref2.fa [...]

  # Backward-compatible single-read-file reference concatenation form
  alignReference.sh -cat READS.fastq[.gz]|READS.bam ref1.fasta ref2.fa [...]

Examples:
  alignReference.sh sample.fastq.gz vector.fasta
  alignReference.sh sample_L001_R2.fastq.gz sample_L002_R2.fastq.gz -cat vector.fasta hg38.fa
  alignReference.sh --paired sample_R1.fastq.gz sample_R2.fastq.gz vector.fasta
  alignReference.sh cellranger.bam vector.fasta
  alignReference.sh --paired cellranger.bam -cat vector.fasta hg38.fa

Important notes:
  * Default multi-FASTQ mode concatenates all provided FASTQs and aligns them as single-end.
  * Only provide FASTQs that should be aligned. Do not blindly provide raw 10x I1/I2/R1/R2 files.
  * Standard 10x GEX R1 is usually cell barcode + UMI, not biological insert sequence.
  * FASTQ input mode does not convert raw 10x R1 barcode/UMI sequence into CB/UB BAM tags.
  * For Cell Ranger BAM input, CB/UB/CR/UR tags are preserved through samtools fastq -T and minimap2 -y.
  * Multiple BAM/CRAM/SAM inputs are not supported; run one alignment file per job.
USAGE
}

strip_read_ext() {
  local b="$1"
  for ext in .fastq.gz .fq.gz .fastq .fq .bam .cram .sam; do
    b="${b%$ext}"
  done
  printf '%s' "$b"
}

read_kind_for() {
  local b
  b="$(basename "$1")"
  case "$b" in
    *.fastq|*.fq|*.fastq.gz|*.fq.gz)
      printf 'fastq'
      ;;
    *.bam|*.cram|*.sam)
      printf 'alignment'
      ;;
    *)
      printf 'unsupported'
      ;;
  esac
}

cat_fastqs_to() {
  local out="$1"
  shift
  : > "$out"
  local fq
  for fq in "$@"; do
    case "$fq" in
      *.gz)
        gzip -dc "$fq" >> "$out"
        ;;
      *)
        cat "$fq" >> "$out"
        ;;
    esac
  done
}

append_nonempty_fastq() {
  local src="$1"
  local dest="$2"
  if [[ -s "$src" ]]; then
    cat "$src" >> "$dest"
    return 0
  fi
  return 1
}

# -----------------------------
# Parse args
# -----------------------------
PAIRED_MODE=0
CAT_MODE=0
READ_INPUTS=()
REF_INPUTS=()
POSITIONAL=()

if [[ $# -eq 0 ]]; then
  usage
  exit 1
fi

while [[ $# -gt 0 ]]; do
  case "$1" in
    --paired)
      PAIRED_MODE=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    -cat)
      CAT_MODE=1
      shift
      if [[ ${#POSITIONAL[@]} -eq 0 ]]; then
        # Backward-compatible form:
        #   alignReference.sh -cat READS ref1 ref2 [...]
        [[ $# -ge 2 ]] || { usage; exit 1; }
        READ_INPUTS=("$1")
        shift
        REF_INPUTS=("$@")
        set --
      else
        # Preferred form for multiple read files:
        #   alignReference.sh READ1 [READ2 ...] -cat ref1 ref2 [...]
        [[ $# -ge 1 ]] || { usage; exit 1; }
        READ_INPUTS=("${POSITIONAL[@]}")
        REF_INPUTS=("$@")
        set --
      fi
      ;;
    --)
      shift
      while [[ $# -gt 0 ]]; do
        POSITIONAL+=("$1")
        shift
      done
      ;;
    -*)
      echo "ERROR: unknown option: $1" >&2
      usage
      exit 1
      ;;
    *)
      POSITIONAL+=("$1")
      shift
      ;;
  esac
done

if [[ ${#READ_INPUTS[@]} -eq 0 && ${#REF_INPUTS[@]} -eq 0 ]]; then
  [[ ${#POSITIONAL[@]} -ge 2 ]] || { usage; exit 1; }
  last_idx=$((${#POSITIONAL[@]} - 1))
  READ_INPUTS=("${POSITIONAL[@]:0:$last_idx}")
  REF_INPUTS=("${POSITIONAL[$last_idx]}")
fi

[[ ${#READ_INPUTS[@]} -ge 1 ]] || { echo "ERROR: no reads input provided" >&2; usage; exit 1; }
[[ ${#REF_INPUTS[@]} -ge 1 ]] || { echo "ERROR: no reference FASTA provided" >&2; usage; exit 1; }

READ_INPUTS_ABS=()
for r in "${READ_INPUTS[@]}"; do
  [[ -f "$r" ]] || { echo "ERROR: reads file not found: $r" >&2; exit 2; }
  READ_INPUTS_ABS+=("$(readlink -f "$r")")
done

REF_INPUTS_ABS=()
for f in "${REF_INPUTS[@]}"; do
  [[ -f "$f" ]] || { echo "ERROR: reference FASTA not found: $f" >&2; exit 2; }
  REF_INPUTS_ABS+=("$(readlink -f "$f")")
done

# -----------------------------
# Validate read input kind
# -----------------------------
READS_KIND="$(read_kind_for "${READ_INPUTS_ABS[0]}")"
if [[ "$READS_KIND" == "unsupported" ]]; then
  echo "ERROR: unsupported reads input extension: ${READ_INPUTS_ABS[0]}" >&2
  echo "Supported: .fastq, .fq, .fastq.gz, .fq.gz, .bam, .cram, .sam" >&2
  exit 2
fi

for r in "${READ_INPUTS_ABS[@]}"; do
  this_kind="$(read_kind_for "$r")"
  if [[ "$this_kind" == "unsupported" ]]; then
    echo "ERROR: unsupported reads input extension: $r" >&2
    echo "Supported: .fastq, .fq, .fastq.gz, .fq.gz, .bam, .cram, .sam" >&2
    exit 2
  fi
  if [[ "$this_kind" != "$READS_KIND" ]]; then
    echo "ERROR: mixed read input types are not supported" >&2
    echo "All reads inputs must be FASTQ-like, or a single BAM/CRAM/SAM." >&2
    exit 2
  fi
done

if [[ "$READS_KIND" == "alignment" && ${#READ_INPUTS_ABS[@]} -ne 1 ]]; then
  echo "ERROR: multiple BAM/CRAM/SAM inputs are not supported; run one alignment file per job." >&2
  exit 2
fi

if [[ "$READS_KIND" == "fastq" && $PAIRED_MODE -eq 1 && ${#READ_INPUTS_ABS[@]} -ne 2 ]]; then
  echo "ERROR: --paired FASTQ mode requires exactly two FASTQ inputs: R1 and R2." >&2
  exit 2
fi

# -----------------------------
# Derive paths
# -----------------------------
READS_DIR="$(dirname "${READ_INPUTS_ABS[0]}")"
FIRST_READS_BASE="$(basename "${READ_INPUTS_ABS[0]}")"
READS_STEM="$(strip_read_ext "$FIRST_READS_BASE")"
if [[ "$READS_KIND" == "fastq" && ${#READ_INPUTS_ABS[@]} -gt 1 && $PAIRED_MODE -eq 0 ]]; then
  READS_STEM="${READS_STEM}__${#READ_INPUTS_ABS[@]}fastqs"
elif [[ "$READS_KIND" == "fastq" && ${#READ_INPUTS_ABS[@]} -eq 2 && $PAIRED_MODE -eq 1 ]]; then
  READS_STEM="${READS_STEM}__paired"
fi

# Build reference stem from all FASTA names (order-preserving)
REF_STEM="$(
  for f in "${REF_INPUTS_ABS[@]}"; do
    b="$(basename "$f")"
    b="${b%.fasta}"
    b="${b%.fa}"
    printf "%s__" "$b"
  done
)"
REF_STEM="${REF_STEM%__}"

OUTDIR="${READS_DIR}/${READS_STEM}__vs__${REF_STEM}"
mkdir -p "$OUTDIR" "$OUTDIR/logs" "$OUTDIR/tmp"

REF_ABS="$OUTDIR/${REF_STEM}.fasta"

BAM="$OUTDIR/${READS_STEM}.vs.${REF_STEM}.sorted.bam"
MM2_LOG="$OUTDIR/${READS_STEM}.vs.${REF_STEM}.minimap2.log"
FLAGSTAT="$OUTDIR/${READS_STEM}.vs.${REF_STEM}.flagstat.txt"
IDXSTATS="$OUTDIR/${READS_STEM}.vs.${REF_STEM}.idxstats.txt"
CMDLOG="$OUTDIR/run.command.txt"
CONVERT_LOG="$OUTDIR/${READS_STEM}.samtools_fastq.log"

MERGED_FASTQ="$OUTDIR/tmp/${READS_STEM}.merged.single_end.fastq"
FASTQ_R1="$OUTDIR/tmp/${READS_STEM}.R1.from_alignment.fastq"
FASTQ_R2="$OUTDIR/tmp/${READS_STEM}.R2.from_alignment.fastq"
FASTQ_SINGLETON="$OUTDIR/tmp/${READS_STEM}.singleton.from_alignment.fastq"
FASTQ_OTHER="$OUTDIR/tmp/${READS_STEM}.other.from_alignment.fastq"
FASTQ_ALL="$OUTDIR/tmp/${READS_STEM}.all.from_alignment.fastq"

# -----------------------------
# Submit to Slurm if needed
# -----------------------------
if [[ -z "${SLURM_JOB_ID:-}" ]]; then
  sbatch_args=(
    --job-name="mm2_${READS_STEM}"
    --cpus-per-task="$CPUS"
    --mem="$MEM"
    --time="$TIME"
    --output="$OUTDIR/logs/slurm_%x_%j.out"
    --error="$OUTDIR/logs/slurm_%x_%j.err"
    --nodelist=ribnode[012]
    # --nodelist=ribnode[003-012,016]
    --nodes=1
    --ntasks=1
  )
  [[ -n "$PARTITION" ]] && sbatch_args+=(--partition="$PARTITION")
  [[ -n "$ACCOUNT"   ]] && sbatch_args+=(--account="$ACCOUNT")

  submit_args=()
  [[ $PAIRED_MODE -eq 1 ]] && submit_args+=(--paired)
  submit_args+=("${READ_INPUTS_ABS[@]}")
  if [[ ${#REF_INPUTS_ABS[@]} -gt 1 || $CAT_MODE -eq 1 ]]; then
    submit_args+=(-cat "${REF_INPUTS_ABS[@]}")
  else
    submit_args+=("${REF_INPUTS_ABS[0]}")
  fi

  sbatch "${sbatch_args[@]}" "$0" "${submit_args[@]}"
  exit 0
fi

# -----------------------------
# Inside Slurm job
# -----------------------------
{
  echo "Date: $(date)"
  echo "Host: $(hostname)"
  echo "Job:  ${SLURM_JOB_ID}"
  printf 'Cmd:  %q' "$0"
  [[ $PAIRED_MODE -eq 1 ]] && printf ' %q' "--paired"
  printf ' %q' "${READ_INPUTS_ABS[@]}"
  if [[ ${#REF_INPUTS_ABS[@]} -gt 1 || $CAT_MODE -eq 1 ]]; then
    printf ' %q' "-cat"
    printf ' %q' "${REF_INPUTS_ABS[@]}"
  else
    printf ' %q' "${REF_INPUTS_ABS[0]}"
  fi
  printf '\n'
  echo "Input kind: $READS_KIND"
  echo "Paired mode: $PAIRED_MODE"
  echo "Preset: $PRESET"
  echo "Preserved BAM tags: $PRESERVE_BAM_TAGS"
  echo "Read inputs:"
  printf '  %s\n' "${READ_INPUTS_ABS[@]}"
  echo "Reference inputs:"
  printf '  %s\n' "${REF_INPUTS_ABS[@]}"
} > "$CMDLOG"

# -----------------------------
# Load modules (explicit & robust)
# -----------------------------
if command -v module >/dev/null 2>&1; then
  module purge 2>/dev/null || true

  module load GCC/14.3.0>/dev/null || module load GCCcore/14.3.0 2>/dev/null || true
  module load zlib/1.3.1 2>/dev/null || true

  module load SAMtools/1.22.1-GCC-14.3.0 2>/dev/null \
    || { echo "ERROR: SAMtools module not found" >&2; exit 3; }

  module load minimap2/2.30-GCCcore-14.3.0 2>/dev/null || true
fi

need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: missing $1" >&2; exit 3; }; }
need minimap2
need samtools
if [[ "$READS_KIND" == "fastq" && $PAIRED_MODE -eq 0 && ${#READ_INPUTS_ABS[@]} -gt 1 ]]; then
  need gzip
fi

THREADS="${SLURM_CPUS_PER_TASK:-$CPUS}"

# -----------------------------
# Build reference FASTA
# -----------------------------
echo "[1/4] Building reference FASTA"
: > "$REF_ABS"
for f in "${REF_INPUTS_ABS[@]}"; do
  grep -q '^>' "$f" || { echo "ERROR: FASTA missing header: $f" >&2; exit 2; }
  cat "$f" >> "$REF_ABS"
  printf "\n" >> "$REF_ABS"
done
samtools faidx "$REF_ABS"

# -----------------------------
# Align → sorted BAM
# -----------------------------
echo "[2/4] Aligning reads"
: > "$MM2_LOG"

if [[ "$READS_KIND" == "alignment" ]]; then
  echo "[2/4] Name-collating and converting alignment input to FASTQ bins"
  : > "$CONVERT_LOG"
  rm -f "$FASTQ_R1" "$FASTQ_R2" "$FASTQ_SINGLETON" "$FASTQ_OTHER" "$FASTQ_ALL"

  # Cell Ranger BAMs are usually coordinate-sorted. Name-collation groups mates
  # when they exist, while samtools fastq -0 captures reads that are neither
  # READ1 nor READ2. Selected CB/UB/CR/UR tags are written to FASTQ comments.
  samtools collate -@ "$THREADS" -u -O "${READ_INPUTS_ABS[0]}" 2>> "$CONVERT_LOG" \
    | samtools fastq \
        -@ "$THREADS" \
        -n \
        -T "$PRESERVE_BAM_TAGS" \
        -1 "$FASTQ_R1" \
        -2 "$FASTQ_R2" \
        -s "$FASTQ_SINGLETON" \
        -0 "$FASTQ_OTHER" \
        - 2>> "$CONVERT_LOG"

  if [[ $PAIRED_MODE -eq 1 ]]; then
    echo "[2/4] --paired requested for alignment input; using only R1/R2 FASTQ bins"
    if [[ ! -s "$FASTQ_R1" || ! -s "$FASTQ_R2" ]]; then
      echo "ERROR: --paired was requested, but paired R1/R2 FASTQs were not produced from ${READ_INPUTS_ABS[0]}" >&2
      echo "R1: $FASTQ_R1" >&2
      echo "R2: $FASTQ_R2" >&2
      echo "See log: $CONVERT_LOG" >&2
      exit 4
    fi

    minimap2 -y -t "$THREADS" -x "$PRESET" --MD -a "$REF_ABS" "$FASTQ_R1" "$FASTQ_R2" 2>> "$MM2_LOG" \
      | samtools sort -@ "$THREADS" -o "$BAM" -
  else
    echo "[2/4] Default alignment-input mode; concatenating all nonempty FASTQ bins and aligning single-end"
    : > "$FASTQ_ALL"
    appended=0
    append_nonempty_fastq "$FASTQ_R1" "$FASTQ_ALL" && appended=1 || true
    append_nonempty_fastq "$FASTQ_R2" "$FASTQ_ALL" && appended=1 || true
    append_nonempty_fastq "$FASTQ_SINGLETON" "$FASTQ_ALL" && appended=1 || true
    append_nonempty_fastq "$FASTQ_OTHER" "$FASTQ_ALL" && appended=1 || true

    if [[ $appended -eq 0 || ! -s "$FASTQ_ALL" ]]; then
      echo "ERROR: no FASTQ records were produced from ${READ_INPUTS_ABS[0]}" >&2
      echo "Expected at least one nonempty bin among:" >&2
      echo "  $FASTQ_R1" >&2
      echo "  $FASTQ_R2" >&2
      echo "  $FASTQ_SINGLETON" >&2
      echo "  $FASTQ_OTHER" >&2
      echo "See log: $CONVERT_LOG" >&2
      exit 4
    fi

    minimap2 -y -t "$THREADS" -x "$PRESET" --MD -a "$REF_ABS" "$FASTQ_ALL" 2>> "$MM2_LOG" \
      | samtools sort -@ "$THREADS" -o "$BAM" -
  fi

else
  if [[ $PAIRED_MODE -eq 1 ]]; then
    echo "[2/4] Aligning exactly two FASTQ inputs in paired-end mode"
    minimap2 -t "$THREADS" -x "$PRESET" --MD -a "$REF_ABS" "${READ_INPUTS_ABS[0]}" "${READ_INPUTS_ABS[1]}" 2> "$MM2_LOG" \
      | samtools sort -@ "$THREADS" -o "$BAM" -
  else
    if [[ ${#READ_INPUTS_ABS[@]} -eq 1 ]]; then
      echo "[2/4] Aligning one FASTQ input as single-end"
      minimap2 -t "$THREADS" -x "$PRESET" --MD -a "$REF_ABS" "${READ_INPUTS_ABS[0]}" 2> "$MM2_LOG" \
        | samtools sort -@ "$THREADS" -o "$BAM" -
    else
      echo "[2/4] Concatenating ${#READ_INPUTS_ABS[@]} FASTQ inputs, then aligning as single-end"
      cat_fastqs_to "$MERGED_FASTQ" "${READ_INPUTS_ABS[@]}"
      minimap2 -t "$THREADS" -x "$PRESET" --MD -a "$REF_ABS" "$MERGED_FASTQ" 2> "$MM2_LOG" \
        | samtools sort -@ "$THREADS" -o "$BAM" -
    fi
  fi
fi

# -----------------------------
# QC
# -----------------------------
echo "[3/4] Indexing + QC"
samtools index "$BAM"
samtools flagstat "$BAM" > "$FLAGSTAT"
samtools idxstats "$BAM" > "$IDXSTATS"

# -----------------------------
# Done
# -----------------------------
echo "[4/4] Done"
echo "Output directory:"
echo "  $OUTDIR"
echo "Reference:"
echo "  $REF_ABS"
echo "BAM:"
echo "  $BAM"
if [[ "$READS_KIND" == "alignment" ]]; then
  echo "Converted FASTQ bins:"
  echo "  R1:        $FASTQ_R1"
  echo "  R2:        $FASTQ_R2"
  echo "  singleton: $FASTQ_SINGLETON"
  echo "  other:     $FASTQ_OTHER"
  if [[ $PAIRED_MODE -eq 0 ]]; then
    echo "Combined single-end FASTQ:"
    echo "  $FASTQ_ALL"
  fi
elif [[ "$READS_KIND" == "fastq" && $PAIRED_MODE -eq 0 && ${#READ_INPUTS_ABS[@]} -gt 1 ]]; then
  echo "Merged single-end FASTQ:"
  echo "  $MERGED_FASTQ"
fi
