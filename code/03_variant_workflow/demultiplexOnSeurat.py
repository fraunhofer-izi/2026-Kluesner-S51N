#!/usr/bin/env python3
from __future__ import annotations

"""
demultiplexOnSeurat.py

Demultiplex a Cell Ranger BAM using CB tags and grouping information from
obj@meta.data in a Seurat/SeuratObject .Rds/.rda/.RData file.

Key features:
  - R dependency preflight with live logging
  - Seurat metadata filtering via --filter-on R expression
  - barcode normalization for merged objects, e.g. minded_T3-C_AAAC...-1 -> AAAC...-1
  - duplicate normalized barcode policies for multi-reference Seurat objects
  - BAM-only output by default
  - unassigned BAM named from --filter-on or --unassigned-label
"""

import argparse
import csv
import datetime as dt
import importlib
import platform
import re
import shutil
import subprocess
import sys
import tempfile
import time
from collections import defaultdict
from pathlib import Path

try:
    import pysam  # type: ignore
except Exception as e:
    pysam = None
    PYSAM_IMPORT_ERROR = e
else:
    PYSAM_IMPORT_ERROR = None

DEMUX_SEP = "__"

R_HELPER = r'''
suppressPackageStartupMessages({
  library(SeuratObject)
})

args <- commandArgs(trailingOnly = TRUE)
object_path <- args[[1]]
filter_on <- args[[2]]
demux_on <- args[[3]]
out_tsv <- args[[4]]
preflight_tsv <- args[[5]]
demux_sep <- "__"

message("[R] Starting Seurat metadata extraction")
message("[R] object_path: ", object_path)
message("[R] filter_on: ", filter_on)
message("[R] demux_on: ", demux_on)
message("[R] demux_sep: ", demux_sep)

load_seurat_object <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "rds") {
    message("[R] readRDS() starting")
    obj <- readRDS(path)
    message("[R] readRDS() finished")
    return(obj)
  }
  if (ext %in% c("rda", "rdata")) {
    message("[R] load() starting")
    env <- new.env(parent = emptyenv())
    loaded <- load(path, envir = env)
    message("[R] load() finished; objects loaded: ", paste(loaded, collapse = ","))
    if (length(loaded) == 0) stop("No objects found in .rda/.RData file")
    for (nm in loaded) {
      candidate <- env[[nm]]
      if (inherits(candidate, "Seurat")) return(candidate)
    }
    stop("No Seurat object found in .rda/.RData file")
  }
  stop("Unsupported object file extension. Use .rds, .rda, or .RData")
}

normalize_filter <- function(expr) {
  if (is.na(expr) || expr == "" || expr == "NULL") return(NULL)
  # Convenience support for simple unquoted equality/inequality only.
  if (!grepl("[\"']", expr)) {
    expr <- gsub("==\\s*([A-Za-z0-9_.:-]+)", "== '\\1'", expr)
    expr <- gsub("!=\\s*([A-Za-z0-9_.:-]+)", "!= '\\1'", expr)
  }
  return(expr)
}

obj <- load_seurat_object(object_path)
if (!inherits(obj, "Seurat")) stop("Loaded object is not a Seurat object")

message("[R] Extracting obj@meta.data")
md <- obj@meta.data
message("[R] meta.data rows before filter: ", nrow(md), " columns: ", ncol(md))

demux_cols <- trimws(unlist(strsplit(demux_on, ",", fixed = TRUE)))
demux_cols <- demux_cols[nzchar(demux_cols)]
if (length(demux_cols) == 0) stop("--demux-on did not contain any metadata column names")
missing_demux_cols <- setdiff(demux_cols, colnames(md))
if (length(missing_demux_cols) > 0) {
  stop(paste0("--demux-on column(s) not found in meta.data: ", paste(missing_demux_cols, collapse = ","),
              ". Available columns: ", paste(colnames(md), collapse = ",")))
}
message("[R] demux columns: ", paste(demux_cols, collapse = ","))

n_cells_before_filter <- nrow(md)
n_columns <- ncol(md)
filter_expr <- normalize_filter(filter_on)
message("[R] normalized filter expression: ", ifelse(is.null(filter_expr), "NULL", filter_expr))

if (!is.null(filter_expr)) {
  message("[R] Applying metadata filter expression to obj@meta.data")
  message("[R] filter expression: ", filter_expr)
  keep <- eval(parse(text = filter_expr), envir = md, enclos = parent.frame())
  if (!is.logical(keep)) stop("--filter-on did not evaluate to a logical vector")
  if (length(keep) != nrow(md)) stop("--filter-on logical vector length does not match meta.data rows")
  keep[is.na(keep)] <- FALSE
  md <- md[keep, , drop = FALSE]
  message("[R] rows after filter: ", nrow(md))
  message("[R] rows removed by filter: ", n_cells_before_filter - nrow(md))
}

n_cells_after_filter <- nrow(md)
message("[R] Removing rows with NA/empty values in any demux column")
demux_values <- as.data.frame(lapply(md[, demux_cols, drop = FALSE], as.character), stringsAsFactors = FALSE)
keep_demux <- stats::complete.cases(demux_values)
for (col in demux_cols) keep_demux <- keep_demux & !is.na(demux_values[[col]]) & nzchar(demux_values[[col]])
md <- md[keep_demux, , drop = FALSE]
demux_values <- demux_values[keep_demux, , drop = FALSE]
n_cells_after_non_na_demux <- nrow(md)
message("[R] Rows after removing NA/empty demux values: ", n_cells_after_non_na_demux)

if (length(demux_cols) == 1) {
  group <- demux_values[[demux_cols[[1]]]]
} else {
  group <- apply(demux_values[, demux_cols, drop = FALSE], 1, paste, collapse = demux_sep)
}

out <- data.frame(cell_id = rownames(md), group = as.character(group), stringsAsFactors = FALSE)
message("[R] Writing raw demultiplex map: ", out_tsv)
write.table(out, file = out_tsv, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

group_counts <- as.data.frame(table(out$group), stringsAsFactors = FALSE)
colnames(group_counts) <- c("group", "n_cells")

preflight <- data.frame(
  metric = c("object_path", "object_class", "metadata_columns", "cells_before_filter",
             "filter_on_original", "filter_on_evaluated", "cells_after_filter",
             "cells_removed_by_filter", "demux_on", "demux_columns", "demux_sep",
             "cells_after_removing_na_demux", "n_demux_groups"),
  value = c(object_path, paste(class(obj), collapse = ","), n_columns, n_cells_before_filter,
            filter_on, ifelse(is.null(filter_expr), "NULL", filter_expr), n_cells_after_filter,
            n_cells_before_filter - n_cells_after_filter, demux_on, paste(demux_cols, collapse = ","), demux_sep,
            n_cells_after_non_na_demux, length(unique(out$group))),
  stringsAsFactors = FALSE
)
message("[R] Writing Seurat preflight TSV: ", preflight_tsv)
write.table(preflight, file = preflight_tsv, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

group_count_path <- sub("\\.tsv$", "_group_counts.tsv", preflight_tsv)
message("[R] Writing group-count TSV: ", group_count_path)
write.table(group_counts, file = group_count_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
message("[R] Finished Seurat metadata extraction")
'''


def timestamp() -> str:
    return dt.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")


def now_human() -> str:
    return dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def format_int(x) -> str:
    return f"{int(x):,}"


def format_elapsed(seconds: float) -> str:
    seconds = int(seconds)
    return f"{seconds//3600:02d}:{(seconds%3600)//60:02d}:{seconds%60:02d}"


def sanitize_filename(value: str) -> str:
    value = str(value)
    value = re.sub(r"[^\w.-]+", "_", value)
    value = re.sub(r"_+", "_", value).strip("_")
    return value or "NA"


def filter_to_filename_label(filter_on: str | None, fallback: str = "all") -> str:
    """Convert an R-style --filter-on expression into a safe filename label."""
    if filter_on is None or str(filter_on).strip() == "":
        return fallback
    label = str(filter_on).strip()
    label = label.replace("&&", "&").replace("||", "|")
    label = label.replace('"', "").replace("'", "")
    label = re.sub(r"\s*==\s*", "_", label)
    label = re.sub(r"\s*!=\s*", "_not_", label)
    label = re.sub(r"\s*>=\s*", "_gte_", label)
    label = re.sub(r"\s*<=\s*", "_lte_", label)
    label = re.sub(r"\s*>\s*", "_gt_", label)
    label = re.sub(r"\s*<\s*", "_lt_", label)
    label = re.sub(r"\s*&\s*", "__", label)
    label = re.sub(r"\s*\|\s*", "__OR__", label)
    label = label.replace("grepl(", "grepl_").replace(")", "")
    label = re.sub(r"[^\w.\-]+", "_", label)
    label = re.sub(r"_+", "_", label).strip("_")
    return label or fallback


class Logger:
    def __init__(self, log_path: Path):
        self.log_path = log_path
        self.log_path.parent.mkdir(parents=True, exist_ok=True)

    def write(self, message: str, also_stderr: bool = True):
        line = f"[{now_human()}] {message}"
        with open(self.log_path, "a") as handle:
            handle.write(line + "\n")
            handle.flush()
        if also_stderr:
            print(line, file=sys.stderr, flush=True)

    def section(self, title: str):
        self.write("")
        self.write("=" * 80)
        self.write(title)
        self.write("=" * 80)


def get_pysam():
    if pysam is None:
        raise RuntimeError(f"Could not import pysam. Install with: python3 -m pip install pysam\nOriginal error: {PYSAM_IMPORT_ERROR}")
    return pysam


def require_executable(name: str, logger: Logger | None = None, required: bool = True):
    exe = shutil.which(name)
    if exe is None:
        if required:
            raise RuntimeError(f"Required executable not found on PATH: {name}")
        if logger:
            logger.write(f"[warning] Executable not found on PATH: {name}")
        return None
    if logger:
        logger.write(f"[ok] Executable available: {name} -> {exe}")
    return exe


def run_cmd(cmd: list[str], logger: Logger | None = None):
    cmd_str = " ".join(map(str, cmd))
    if logger:
        logger.write(f"[cmd] {cmd_str}")
        logger.write("[info] Subprocess started; streaming output...")
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)
    output_lines: list[str] = []
    assert process.stdout is not None
    for line in process.stdout:
        line = line.rstrip("\n")
        output_lines.append(line)
        if logger:
            logger.write(f"[subprocess] {line}")
    return_code = process.wait()
    if logger:
        logger.write(f"[info] Subprocess finished with exit code {return_code}")
    output = "\n".join(output_lines)
    if return_code != 0:
        raise RuntimeError("Command failed:\n" + cmd_str + "\n\nOUTPUT:\n" + output)
    return output


def file_size_human(path: Path) -> str:
    value = float(path.stat().st_size)
    for unit in ["B", "KB", "MB", "GB", "TB"]:
        if value < 1024 or unit == "TB":
            return f"{value:.2f} {unit}"
        value /= 1024
    return f"{value:.2f} TB"


def run_dependency_preflight(args, out_dir: Path, logger: Logger):
    logger.section("Dependency preflight")
    logger.write(f"[ok] Python executable: {sys.executable}")
    logger.write(f"[ok] Python version: {platform.python_version()}")
    try:
        p = importlib.import_module("pysam")
        logger.write(f"[ok] Python package available: pysam {getattr(p, '__version__', 'unknown')}")
    except Exception as e:
        raise RuntimeError(f"Required Python package is not available: pysam\nInstall with: python3 -m pip install pysam\nError: {e}")

    rscript_path = require_executable("Rscript", logger=logger, required=True)
    rver = subprocess.run(["Rscript", "--version"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    logger.write(f"[ok] {rver.stdout.strip() if rver.stdout.strip() else 'Rscript version unknown'}")

    r_check = r'''
cat("[R dependency check] .libPaths():\n")
cat(paste(.libPaths(), collapse=";"), "\n")
required <- c("sp", "future", "SeuratObject")
for (pkg in required) {
  cat(paste0("[checking]\t", pkg, "\n"))
  ok <- requireNamespace(pkg, quietly = TRUE)
  if (!ok) { cat(paste0(pkg, "\tMISSING\t\tPackage not installed or namespace could not be loaded\n")); next }
  load_ok <- TRUE; load_msg <- ""
  tryCatch({ suppressPackageStartupMessages(library(pkg, character.only = TRUE)) }, error=function(e){ load_ok <<- FALSE; load_msg <<- conditionMessage(e) })
  if (!load_ok) cat(paste0(pkg, "\tLOAD_FAILED\t\t", load_msg, "\n")) else cat(paste0(pkg, "\tOK\t", as.character(utils::packageVersion(pkg)), "\t", find.package(pkg), "\n"))
}
stack_ok <- TRUE; stack_msg <- ""
tryCatch({ suppressPackageStartupMessages(library(SeuratObject)) }, error=function(e){ stack_ok <<- FALSE; stack_msg <<- conditionMessage(e) })
if (!stack_ok) cat(paste0("SeuratObject_stack\tLOAD_FAILED\t\t", stack_msg, "\n")) else cat("SeuratObject_stack\tOK\t\tSeuratObject loaded successfully\n")
'''
    res = subprocess.run(["Rscript", "-e", r_check], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    missing: list[str] = []
    for line in res.stdout.strip().splitlines():
        logger.write(f"[R dependency] {line}")
        if line.startswith("[") or not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) >= 2 and parts[1] != "OK":
            missing.append(parts[0])
    if res.returncode != 0 or missing:
        raise RuntimeError(
            "R dependency preflight failed. Missing/unloadable: " + ", ".join(missing) +
            "\nSuggested fix:\nmodule load R\nmkdir -p ~/R/libs\nexport R_LIBS_USER=~/R/libs\n" +
            "Rscript -e 'install.packages(c(\"sp\", \"future\", \"SeuratObject\"), lib=\"~/R/libs\", repos=\"https://cloud.r-project.org\")'"
        )
    logger.write("[ok] Dependency preflight complete")


def check_bam_readable(bam_path: Path, logger: Logger):
    logger.section("BAM preflight")
    ps = get_pysam()
    if not bam_path.exists():
        raise FileNotFoundError(f"Input BAM not found: {bam_path}")
    logger.write(f"[ok] BAM exists: {bam_path}")
    logger.write(f"[info] BAM size: {file_size_human(bam_path)}")
    with ps.AlignmentFile(str(bam_path), "rb") as bam:
        logger.write("[ok] BAM opened successfully")
        logger.write(f"[info] Number of reference sequences in header: {len(bam.references)}")
        if bam.references:
            logger.write("[info] First reference sequences: " + ", ".join([f"{r}:{l}" for r, l in list(zip(bam.references, bam.lengths))[:10]]))
        logger.write(f"[info] BAM index detected by pysam: {bam.has_index()}")
        if not bam.has_index():
            logger.write("[warning] BAM does not appear indexed. Streaming will still work; percent denominator may be unavailable.")


def sample_bam_cb_tags(bam_path: Path, max_records: int, logger: Logger) -> dict[str, int]:
    logger.section("BAM CB tag sampling")
    ps = get_pysam()
    counts: dict[str, int] = {}
    total = with_cb = missing_cb = 0
    with ps.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam.fetch(until_eof=True):
            total += 1
            try:
                cb = read.get_tag("CB")
                with_cb += 1
                counts[cb] = counts.get(cb, 0) + 1
            except KeyError:
                missing_cb += 1
            if total >= max_records:
                break
    logger.write(f"[info] Sampled records: {format_int(total)}")
    logger.write(f"[info] Sampled records with CB: {format_int(with_cb)}")
    logger.write(f"[info] Sampled records missing CB: {format_int(missing_cb)}")
    logger.write(f"[info] Unique CB values in sample: {format_int(len(counts))}")
    logger.write("[info] Top sampled CB tags:")
    for cb, count in sorted(counts.items(), key=lambda x: x[1], reverse=True)[:10]:
        logger.write(f"  {cb}\t{count}", also_stderr=False)
    return counts


def extract_metadata_map(object_path: Path, filter_on: str | None, demux_on: str, out_tsv: Path, preflight_tsv: Path, logger: Logger):
    require_executable("Rscript", logger=logger)
    logger.section("Seurat metadata preflight")
    logger.write(f"[info] Seurat object path: {object_path}")
    logger.write(f"[info] Seurat object size: {file_size_human(object_path)}")
    logger.write(f"[info] --filter-on: {filter_on if filter_on else 'NULL'}")
    logger.write(f"[info] --demux-on: {demux_on}")
    with tempfile.NamedTemporaryFile("w", suffix=".R", delete=False) as handle:
        helper_path = Path(handle.name)
        handle.write(R_HELPER)
    try:
        run_cmd(["Rscript", str(helper_path), str(object_path), filter_on if filter_on else "NULL", demux_on, str(out_tsv), str(preflight_tsv)], logger=logger)
    finally:
        helper_path.unlink(missing_ok=True)
    logger.write(f"[ok] Wrote temporary raw metadata barcode map: {out_tsv}")
    logger.write(f"[ok] Wrote temporary Seurat preflight TSV: {preflight_tsv}")


def read_raw_map(path: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            cell_id = row.get("cell_id") or row.get("barcode")
            group = row.get("group")
            if cell_id is None or group is None:
                raise RuntimeError(f"Raw map must contain cell_id/barcode and group columns: {path}")
            rows.append({"cell_id": cell_id, "group": group})
    return rows


def infer_prefix(cell_id: str, regex: str) -> str:
    m = re.search(regex, cell_id)
    if m:
        return cell_id[:m.start()].rstrip("_:-|") or "<none>"
    if "_" in cell_id:
        return cell_id.rsplit("_", 1)[0]
    return "<none>"


def normalize_barcode(cell_id: str, strategy: str, regex: str) -> str | None:
    if strategy == "exact":
        return cell_id
    if strategy == "regex":
        m = re.search(regex, cell_id)
        return m.group(1) if m else None
    if strategy == "strip-prefix":
        return cell_id.rsplit("_", 1)[-1]
    raise ValueError(f"Unknown barcode strategy: {strategy}")


def compare_barcode_strategies(raw_rows: list[dict[str, str]], sampled_cb_counts: dict[str, int], args, logger: Logger) -> str:
    sampled = set(sampled_cb_counts)
    strategies = ["exact", "regex", "strip-prefix"] if args.barcode_match_strategy == "auto" else [args.barcode_match_strategy]
    report: list[dict] = []
    for strategy in strategies:
        normalized = [normalize_barcode(r["cell_id"], strategy, args.barcode_regex) for r in raw_rows]
        normalized_nonnull = [x for x in normalized if x]
        unique = set(normalized_nonnull)
        overlap = unique & sampled
        report.append({
            "strategy": strategy,
            "unique_normalized_barcodes": len(unique),
            "sampled_bam_unique_cb": len(sampled),
            "sampled_overlap": len(overlap),
            "sampled_overlap_frac": (len(overlap) / len(sampled)) if sampled else 0,
        })
    logger.write("[info] Barcode strategy comparison:")
    for r in report:
        logger.write(f"  {r['strategy']} overlap={r['sampled_overlap']}/{r['sampled_bam_unique_cb']} ({100*r['sampled_overlap_frac']:.2f}%) unique_norm={format_int(r['unique_normalized_barcodes'])}")
    if args.barcode_match_strategy == "auto":
        best = max(report, key=lambda r: (r["sampled_overlap"], r["unique_normalized_barcodes"]))
        logger.write(f"[info] Auto-selected barcode strategy: {best['strategy']}")
        return str(best["strategy"])
    logger.write(f"[info] User-selected barcode strategy: {args.barcode_match_strategy}")
    return args.barcode_match_strategy


def build_barcode_map(raw_rows: list[dict[str, str]], sampled_cb_counts: dict[str, int], args, logger: Logger) -> dict[str, str]:
    logger.section("Barcode normalization preflight")
    logger.write(f"[info] Raw metadata rows before barcode normalization: {format_int(len(raw_rows))}")
    prefixes = defaultdict(int)
    for row in raw_rows:
        prefixes[infer_prefix(row["cell_id"], args.barcode_regex)] += 1
    logger.write("[info] Top inferred Seurat cell prefixes before barcode normalization:")
    for prefix, count in sorted(prefixes.items(), key=lambda x: x[1], reverse=True)[:20]:
        logger.write(f"  {prefix}\t{count}", also_stderr=False)
    logger.write("[info] Example raw Seurat cell IDs after R metadata filter:")
    for row in raw_rows[:10]:
        logger.write(f"  {row['cell_id']} -> {row['group']}", also_stderr=False)
    logger.write("[info] Example sampled BAM CB tags:")
    for cb in list(sampled_cb_counts.keys())[:10]:
        logger.write(f"  {cb}", also_stderr=False)

    strategy = compare_barcode_strategies(raw_rows, sampled_cb_counts, args, logger)
    buckets: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in raw_rows:
        nb = normalize_barcode(row["cell_id"], strategy, args.barcode_regex)
        if nb:
            buckets[nb].append({"cell_id": row["cell_id"], "group": row["group"]})

    barcode_to_group: dict[str, str] = {}
    duplicate_count = same_group_dups = conflicting_dups = ambiguous_excluded = 0
    ambiguous_examples: list[str] = []
    usable_examples: list[str] = []

    for nb, entries in buckets.items():
        groups = sorted({e["group"] for e in entries})
        is_dup = len(entries) > 1
        if is_dup:
            duplicate_count += 1
        if len(groups) == 1:
            if is_dup:
                same_group_dups += 1
            if args.duplicate_normalized_barcode_policy == "error" and is_dup:
                continue
            barcode_to_group[nb] = groups[0]
            if len(usable_examples) < 10:
                usable_examples.append(f"{nb}\tgroup={groups[0]}\tn_source_cells={len(entries)}\tsource_cell_ids={';'.join([e['cell_id'] for e in entries[:8]])}")
        else:
            conflicting_dups += 1
            if args.duplicate_normalized_barcode_policy == "first":
                barcode_to_group[nb] = entries[0]["group"]
            else:
                ambiguous_excluded += 1
                if len(ambiguous_examples) < 10:
                    ambiguous_examples.append(f"{nb}\tgroups={';'.join(groups)}\tn_source_cells={len(entries)}\tsource_cell_ids={';'.join([e['cell_id'] for e in entries[:8]])}")

    if duplicate_count and args.duplicate_normalized_barcode_policy == "error":
        examples = []
        for nb, entries in buckets.items():
            if len(entries) > 1:
                examples.append(f"  {nb}: " + "; ".join([f"{e['cell_id']}->{e['group']}" for e in entries[:8]]))
            if len(examples) >= 10:
                break
        raise RuntimeError(
            "Barcode normalization created duplicate normalized barcodes.\n"
            "This usually means the Seurat object contains multiple reference-alignment rows or multiple libraries/samples.\n"
            f"Strategy: {strategy}\nDuplicate normalized barcodes: {format_int(duplicate_count)}\nExamples:\n" + "\n".join(examples) +
            "\nUse --duplicate-normalized-barcode-policy same-group-ok if duplicates agree on --demux-on, or unassigned-ambiguous to exclude conflicting duplicates."
        )

    if args.duplicate_normalized_barcode_policy == "same-group-ok" and conflicting_dups:
        raise RuntimeError(
            f"Found {format_int(conflicting_dups)} normalized barcodes mapping to multiple {args.demux_on} groups. "
            "Use --duplicate-normalized-barcode-policy unassigned-ambiguous to send these to unassigned, or filter more strictly."
        )

    if ambiguous_examples:
        logger.write("[info] Example ambiguous normalized barcodes excluded:", also_stderr=False)
        for x in ambiguous_examples:
            logger.write(f"[ambiguous] {x}", also_stderr=False)
    if usable_examples:
        logger.write("[info] Example usable normalized barcode mappings:", also_stderr=False)
        for x in usable_examples:
            logger.write(f"[usable_map] {x}", also_stderr=False)

    sampled = set(sampled_cb_counts)
    overlap = sampled & set(barcode_to_group)
    overlap_frac = len(overlap) / len(sampled) if sampled else 0
    logger.write(f"[info] Duplicate normalized barcodes: {format_int(duplicate_count)}")
    logger.write(f"[info] Duplicate normalized barcodes with same {args.demux_on}: {format_int(same_group_dups)}")
    logger.write(f"[info] Duplicate normalized barcodes with conflicting {args.demux_on}: {format_int(conflicting_dups)}")
    logger.write(f"[info] Usable normalized barcodes: {format_int(len(barcode_to_group))}")
    logger.write(f"[info] Ambiguous normalized barcodes excluded: {format_int(ambiguous_excluded)}")
    logger.write(f"[info] Final sampled CB overlap: {format_int(len(overlap))}/{format_int(len(sampled))} ({100*overlap_frac:.2f}%)")
    return barcode_to_group


def estimate_bam_total_records_for_progress(bam_path: Path, logger: Logger) -> int | None:
    ps = get_pysam()
    try:
        with ps.AlignmentFile(str(bam_path), "rb") as bam:
            if not bam.has_index():
                logger.write("[warning] No BAM index available; progress bar will show records processed without percent denominator.")
                return None
        stats = ps.idxstats(str(bam_path))
        total = 0
        for line in stats.strip().splitlines():
            parts = line.split("\t")
            if len(parts) >= 4:
                total += int(parts[2]) + int(parts[3])
        if total > 0:
            logger.write(f"[info] Progress denominator from BAM index: {format_int(total)} records")
            return total
    except Exception as e:
        logger.write(f"[warning] Could not estimate total records from BAM index: {e}")
    return None


def print_progress(summary: dict, start_time: float, logger: Logger, total_estimate: int | None = None, final: bool = False, width: int = 36):
    elapsed = time.time() - start_time
    processed = summary["total_records"]
    assigned = summary["assigned_records"]
    with_cb = summary["records_with_CB"]
    rate = processed / elapsed if elapsed > 0 else 0
    assigned_cb_pct = 100 * assigned / with_cb if with_cb else 0
    if total_estimate:
        frac = min(processed / total_estimate, 1.0)
        filled = int(frac * width)
        bar = "#" * filled + "-" * (width - filled)
        msg = (f"\r[demux] |{bar}| {100*frac:6.2f}% {format_int(processed)}/{format_int(total_estimate)} records "
               f"assigned={format_int(assigned)} unassigned={format_int(summary['unassigned_records'])} "
               f"missing_CB={format_int(summary['missing_CB_records'])} assigned_CB={assigned_cb_pct:.2f}% "
               f"elapsed={format_elapsed(elapsed)} rate={rate:,.0f}/s")
    else:
        msg = (f"\r[demux] processed={format_int(processed)} assigned={format_int(assigned)} "
               f"unassigned={format_int(summary['unassigned_records'])} missing_CB={format_int(summary['missing_CB_records'])} "
               f"assigned_CB={assigned_cb_pct:.2f}% elapsed={format_elapsed(elapsed)} rate={rate:,.0f}/s")
    print(msg, file=sys.stderr, end="\n" if final else "", flush=True)
    logger.write(("[progress] processed={processed} total_estimate={total} with_CB={with_cb} assigned={assigned} "
                  "unassigned={unassigned} missing_CB={missing} assigned_CB={pct:.2f}% elapsed={elapsed} rate={rate:.0f} reads/sec").format(
        processed=format_int(processed), total=format_int(total_estimate) if total_estimate else "NA", with_cb=format_int(with_cb),
        assigned=format_int(assigned), unassigned=format_int(summary['unassigned_records']), missing=format_int(summary['missing_CB_records']),
        pct=assigned_cb_pct, elapsed=format_elapsed(elapsed), rate=rate), also_stderr=False)


def split_bam_by_barcode_group(bam_path: Path, barcode_to_group: dict[str, str], out_dir: Path, demux_on: str, logger: Logger, threads: int, progress_every: int, unassigned_label: str):
    logger.section("BAM demultiplexing")
    ps = get_pysam()
    summary = {"total_records": 0, "records_with_CB": 0, "assigned_records": 0, "unassigned_records": 0, "missing_CB_records": 0}
    group_counts: dict[str, int] = {}
    writers: dict[str, object] = {}
    total_estimate = estimate_bam_total_records_for_progress(bam_path, logger)

    in_bam = ps.AlignmentFile(str(bam_path), "rb")
    header = in_bam.header.to_dict()

    def get_writer(group_name: str):
        if group_name not in writers:
            out_bam = out_dir / f"{sanitize_filename(demux_on.replace(',', '_'))}__{sanitize_filename(group_name)}.bam"
            logger.write(f"[open] Creating BAM for group {group_name}: {out_bam}")
            writers[group_name] = ps.AlignmentFile(str(out_bam), "wb", header=header, threads=threads)
        return writers[group_name]

    unassigned_bam = out_dir / f"unassigned__{sanitize_filename(unassigned_label)}.bam"
    logger.write(f"[open] Creating unassigned BAM: {unassigned_bam}")
    unassigned_writer = ps.AlignmentFile(str(unassigned_bam), "wb", header=header, threads=threads)

    start_time = time.time()
    try:
        for read in in_bam.fetch(until_eof=True):
            summary["total_records"] += 1
            try:
                cb = read.get_tag("CB")
            except KeyError:
                summary["missing_CB_records"] += 1
                unassigned_writer.write(read)
            else:
                summary["records_with_CB"] += 1
                group = barcode_to_group.get(cb)
                if group is None:
                    summary["unassigned_records"] += 1
                    unassigned_writer.write(read)
                else:
                    get_writer(group).write(read)
                    summary["assigned_records"] += 1
                    group_counts[group] = group_counts.get(group, 0) + 1
            if summary["total_records"] % progress_every == 0:
                print_progress(summary, start_time, logger, total_estimate=total_estimate)
    finally:
        in_bam.close()
        for writer in writers.values():
            writer.close()
        unassigned_writer.close()
    print_progress(summary, start_time, logger, total_estimate=total_estimate, final=True)

    logger.section("Indexing output BAMs")
    for bam in sorted(out_dir.glob("*.bam")):
        try:
            logger.write(f"[index] {bam}")
            ps.index(str(bam))
        except Exception as e:
            logger.write(f"[warning] Could not index {bam}: {e}")

    logger.section("Demultiplex summary")
    for key, value in summary.items():
        logger.write(f"[summary] {key}\t{value}")
    assigned_pct = 100 * summary["assigned_records"] / summary["total_records"] if summary["total_records"] else 0
    unassigned_pct = 100 * summary["unassigned_records"] / summary["total_records"] if summary["total_records"] else 0
    logger.write(f"[summary] assigned_pct\t{assigned_pct:.2f}")
    logger.write(f"[summary] unassigned_pct\t{unassigned_pct:.2f}")

    logger.section("Demultiplex group counts")
    if group_counts:
        for group, count in sorted(group_counts.items()):
            logger.write(f"[group_count] {group}\t{count}")
    else:
        logger.write("[group_count] No assigned groups")
    logger.write("[ok] Demultiplex summary and group counts recorded in log")
    return summary


def write_run_manifest(args, logger: Logger):
    logger.section("Run manifest")
    for key, value in sorted(vars(args).items()):
        logger.write(f"[manifest] {key}\t{value}", also_stderr=False)
    logger.write(f"[manifest] timestamp\t{now_human()}", also_stderr=False)
    logger.write(f"[manifest] command\t{' '.join(sys.argv)}", also_stderr=False)
    logger.write("[ok] Run manifest recorded in log")


def run_preflight(args, out_dir: Path, logger: Logger) -> dict[str, str]:
    logger.section("Forced preflight")
    run_dependency_preflight(args, out_dir, logger)
    bam_path = Path(args.bam).resolve()
    object_path = Path(args.object).resolve()
    if not bam_path.exists():
        raise FileNotFoundError(f"Input BAM not found: {bam_path}")
    if not object_path.exists():
        raise FileNotFoundError(f"Seurat object not found: {object_path}")
    check_bam_readable(bam_path, logger)
    with tempfile.TemporaryDirectory(prefix="demultiplexOnSeurat_") as tmp:
        tmp_dir = Path(tmp)
        logger.write(f"[info] Using temporary directory for intermediate audit tables: {tmp_dir}")
        raw_map = tmp_dir / "demultiplex_map.raw.tsv"
        preflight_seurat = tmp_dir / "preflight_seurat.tsv"
        extract_metadata_map(object_path, args.filter_on, args.demux_on, raw_map, preflight_seurat, logger)
        raw_rows = read_raw_map(raw_map)
        if not raw_rows:
            raise RuntimeError("No metadata rows left after --filter-on and demux NA removal")
        logger.write("[info] Temporary metadata handoff files were read successfully; they will be deleted automatically")
    sampled_cb_counts = sample_bam_cb_tags(bam_path, args.sample_records, logger)
    barcode_to_group = build_barcode_map(raw_rows, sampled_cb_counts, args, logger)
    if not barcode_to_group:
        raise RuntimeError("No usable normalized barcodes remain after duplicate/ambiguity handling")
    logger.write("[ok] Preflight complete")
    return barcode_to_group


def parse_args():
    p = argparse.ArgumentParser(prog="demultiplexOnSeurat", description="Demultiplex a Cell Ranger BAM using Seurat metadata and CB tags.")
    p.add_argument("--bam", required=True, help="Input Cell Ranger BAM")
    p.add_argument("--object", required=True, help="Seurat object path: .rds, .Rds, .rda, or .RData")
    p.add_argument("--filter-on", default=None, help="Optional R expression evaluated in obj@meta.data")
    p.add_argument("--demux-on", required=True, help="Metadata column(s) used to split BAM. Comma-separated for combined groups, e.g. PATIENT_ID,TIMEPOINT")
    p.add_argument("--out", default="demultiplexOnSeurat_out", help="Output directory")
    p.add_argument("--threads", type=int, default=4, help="Threads to use for BAM writing/indexing where supported. Default: 4")
    p.add_argument("--barcode-match-strategy", choices=["auto", "exact", "regex", "strip-prefix"], default="auto")
    p.add_argument("--barcode-regex", default=r"([ACGT]{16}-[0-9]+)$", help="Regex with one capture group extracting the 10x barcode from Seurat rownames")
    p.add_argument("--duplicate-normalized-barcode-policy", choices=["error", "same-group-ok", "unassigned-ambiguous", "first"], default="error", help="How to handle duplicated normalized barcodes after barcode normalization")
    p.add_argument("--preflight-only", action="store_true", help="Run preflight and exit before demultiplexing")
    p.add_argument("--unassigned-label", default=None, help="Optional custom label for the unassigned BAM filename. Default is derived from --filter-on.")
    p.add_argument("--sample-records", type=int, default=200_000, help="Number of BAM records to sample during preflight. Default: 200,000")
    p.add_argument("--progress-every", type=int, default=1_000_000, help="Print progress every N BAM records. Default: 1,000,000")
    # Backward-compatible accepted options from older versions.
    p.add_argument("--output-format", choices=["bam", "fastq", "both"], default="bam", help="Accepted for compatibility; this version writes BAM only.")
    p.add_argument("--strip-gem-suffix", action="store_true", help=argparse.SUPPRESS)
    p.add_argument("--no-unassigned", action="store_true", help=argparse.SUPPRESS)
    p.add_argument("--fastq-tags", default="CB,CR,UB,UR", help=argparse.SUPPRESS)
    return p.parse_args()


def main():
    args = parse_args()
    if args.output_format != "bam":
        print("[warning] This version writes BAM only; --output-format is accepted only for compatibility.", file=sys.stderr, flush=True)
    if args.no_unassigned:
        print("[warning] This version always writes an unassigned BAM; --no-unassigned is ignored.", file=sys.stderr, flush=True)

    out_dir = Path(args.out).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    logger = Logger(out_dir / f"preflight_and_run_{timestamp()}.log")
    logger.section("demultiplexOnSeurat started")
    logger.write(f"[info] Command: {' '.join(sys.argv)}")
    logger.write(f"[info] Output directory: {out_dir}")
    logger.write(f"[info] Log path: {logger.log_path}")
    write_run_manifest(args, logger)

    barcode_to_group = run_preflight(args, out_dir, logger)
    if args.preflight_only:
        logger.section("Preflight-only mode complete")
        logger.write("[done] --preflight-only was requested; exiting before BAM demultiplexing")
        return

    label = sanitize_filename(args.unassigned_label) if args.unassigned_label else filter_to_filename_label(args.filter_on, fallback="all")
    logger.write(f"[info] Unassigned BAM label: {label}")

    logger.section("Run configuration")
    logger.write(f"[info] Loaded {format_int(len(barcode_to_group))} usable metadata barcodes")
    logger.write("[info] Output format: bam")
    logger.write(f"[info] Threads: {args.threads}")
    logger.write(f"[info] Progress interval: every {format_int(args.progress_every)} BAM records")
    summary = split_bam_by_barcode_group(Path(args.bam).resolve(), barcode_to_group, out_dir, args.demux_on, logger, threads=args.threads, progress_every=args.progress_every, unassigned_label=label)

    logger.section("demultiplexOnSeurat complete")
    logger.write(f"[done] Output directory: {out_dir}")
    logger.write(f"[done] Log file: {logger.log_path}")
    for key, value in summary.items():
        print(f"{key}\t{value}")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"[fatal] {e}", file=sys.stderr, flush=True)
        sys.exit(1)
