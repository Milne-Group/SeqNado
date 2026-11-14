#!/usr/bin/env bash
# Split a BAM into reference (endogenous) and spike-in (exogenous) BAMs by contig name.
#
# Contigs are classified from the BAM header (@SQ lines):
#   1. spike-in  - name matches --exo-regex (takes precedence)
#   2. reference - name matches --ref-regex
#   3. unassigned - neither; reads on these contigs are dropped
# Reads are then extracted with `samtools view -L <contigs.bed>`.
#
# Regexes are POSIX ERE, built by seqnado.workflow.helpers.spikein.
#
# Usage (CLI):
#   split_bam_by_contig.sh --bam in.bam --ref-regex RE --exo-regex RE \
#       [--ref-out ref.bam] [--exo-out exo.bam] \
#       [--on-empty-ref error|warn] [--on-empty-exo error|warn] \
#       [--index] [--threads N] [--log file] [--label sample]
#
# Called from Snakemake rules via `shell:` with the path from
# `workflow.source_path(...)`, which Snakemake bind-mounts into containers.
# Bash-3 compatible (macOS) so it can be unit-tested directly.

set -euo pipefail

MAX_LISTED=20

usage() {
    sed -n '2,/^$/p' "$0" | sed 's/^# \{0,1\}//'
}

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
BAM="" REF_RE="" EXO_RE="" REF_OUT="" EXO_OUT="" LOG="" LABEL=""
ON_EMPTY_REF="error" ON_EMPTY_EXO="error" INDEX=0 THREADS=1

while [[ $# -gt 0 ]]; do
    case "$1" in
        --bam) BAM="$2"; shift 2 ;;
        --ref-regex) REF_RE="$2"; shift 2 ;;
        --exo-regex) EXO_RE="$2"; shift 2 ;;
        --ref-out) REF_OUT="$2"; shift 2 ;;
        --exo-out) EXO_OUT="$2"; shift 2 ;;
        --on-empty-ref) ON_EMPTY_REF="$2"; shift 2 ;;
        --on-empty-exo) ON_EMPTY_EXO="$2"; shift 2 ;;
        --index) INDEX=1; shift ;;
        --threads) THREADS="$2"; shift 2 ;;
        --log) LOG="$2"; shift 2 ;;
        --label) LABEL="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

if [[ -n "$LOG" ]]; then
    mkdir -p "$(dirname "$LOG")"
    : > "$LOG"
fi

log() {
    local line
    line="[$(date '+%Y-%m-%d %H:%M:%S')] split_bam_by_contig${LABEL:+ ($LABEL)}: $*"
    if [[ -n "$LOG" ]]; then echo "$line" >> "$LOG"; else echo "$line" >&2; fi
}

die() {
    log "ERROR: $*"
    # Always surface errors on stderr so they appear in the Snakemake output too.
    if [[ -n "$LOG" ]]; then echo "ERROR: $*" >&2; fi
    exit 1
}

[[ -n "$BAM" ]] || die "--bam is required"
[[ -f "$BAM" ]] || die "input BAM not found: $BAM"
[[ -n "$REF_OUT" || -n "$EXO_OUT" ]] || die "at least one of --ref-out / --exo-out is required"
for mode in "$ON_EMPTY_REF" "$ON_EMPTY_EXO"; do
    [[ "$mode" == "error" || "$mode" == "warn" ]] || die "--on-empty-* must be 'error' or 'warn', got '$mode'"
done
command -v samtools >/dev/null 2>&1 || die "samtools not found on PATH"

# Matches only an empty name (contig names are never empty): "match nothing"
# when a class is not configured. Portable across gawk/mawk/BSD awk.
NEVER='^$'
REF_RE="${REF_RE:-$NEVER}"
EXO_RE="${EXO_RE:-$NEVER}"

log "input BAM: $BAM"
log "reference contig pattern: $REF_RE"
log "spike-in contig pattern:  $EXO_RE"

# ---------------------------------------------------------------------------
# Classify contigs from the header
# ---------------------------------------------------------------------------
WORKDIR="$(mktemp -d "${TMPDIR:-/tmp}/split_bam_by_contig.XXXXXX")"
trap 'rm -rf "$WORKDIR"' EXIT
REF_BED="$WORKDIR/ref.bed"
EXO_BED="$WORKDIR/exo.bed"
UNASSIGNED="$WORKDIR/unassigned.txt"
: > "$REF_BED"; : > "$EXO_BED"; : > "$UNASSIGNED"

samtools view -H "$BAM" > "$WORKDIR/header.sam" || die "could not read header of $BAM"

REF_RE="$REF_RE" EXO_RE="$EXO_RE" awk -F'\t' \
    -v ref_bed="$REF_BED" -v exo_bed="$EXO_BED" -v unassigned="$UNASSIGNED" '
    $1 == "@SQ" {
        name = ""; len = "";
        for (i = 2; i <= NF; i++) {
            if ($i ~ /^SN:/) name = substr($i, 4);
            else if ($i ~ /^LN:/) len = substr($i, 4);
        }
        if (name == "") next;
        if (name ~ ENVIRON["EXO_RE"]) print name "\t0\t" len > exo_bed;
        else if (name ~ ENVIRON["REF_RE"]) print name "\t0\t" len > ref_bed;
        else print name > unassigned;
    }' "$WORKDIR/header.sam" \
    || die "failed to classify contigs (invalid regex? ref='$REF_RE' exo='$EXO_RE')"

n_total=$(grep -c '^@SQ' "$WORKDIR/header.sam" || true)
n_ref=$(wc -l < "$REF_BED" | tr -d ' ')
n_exo=$(wc -l < "$EXO_BED" | tr -d ' ')
n_un=$(wc -l < "$UNASSIGNED" | tr -d ' ')

list_names() {
    # $1: file whose first column holds contig names
    local n
    n=$(wc -l < "$1" | tr -d ' ')
    if [[ "$n" -eq 0 ]]; then echo "-"; return; fi
    cut -f1 "$1" | head -n "$MAX_LISTED" | paste -sd ' ' -
    if [[ "$n" -gt "$MAX_LISTED" ]]; then echo "... (+$((n - MAX_LISTED)) more)"; fi
}

log "header contigs: $n_total"
log "  reference:  $n_ref  [$(list_names "$REF_BED" | paste -sd ' ' -)]"
log "  spike-in:   $n_exo  [$(list_names "$EXO_BED" | paste -sd ' ' -)]"
log "  unassigned: $n_un  [$(list_names "$UNASSIGNED" | paste -sd ' ' -)] (reads on these are dropped)"

check_empty() {
    # $1: class label, $2: count, $3: mode, $4: config hint
    if [[ "$2" -eq 0 ]]; then
        local msg="no $1 contigs found in $BAM. Contigs in header: $(grep '^@SQ' "$WORKDIR/header.sam" | sed -n 's/.*SN:\([^[:space:]]*\).*/\1/p' | head -n 50 | paste -sd ' ' -). $4"
        if [[ "$3" == "error" ]]; then die "$msg"; else log "WARNING: $msg"; fi
    fi
    return 0
}

if [[ -n "$REF_OUT" ]]; then
    check_empty "reference" "$n_ref" "$ON_EMPTY_REF" \
        "Set 'endogenous_contig_prefix' in your config to the reference contig prefix(es)."
fi
if [[ -n "$EXO_OUT" ]]; then
    check_empty "spike-in" "$n_exo" "$ON_EMPTY_EXO" \
        "Set 'exogenous_contig_prefix' in your config to the spike-in contig name(s)/prefix(es)."
fi

# ---------------------------------------------------------------------------
# Extract reads
# ---------------------------------------------------------------------------
extract() {
    # $1: class label, $2: bed, $3: output BAM
    local label="$1" bed="$2" out="$3"
    mkdir -p "$(dirname "$out")"
    if [[ -s "$bed" ]]; then
        samtools view -b -@ "$THREADS" -L "$bed" -o "$out" "$BAM" \
            || die "samtools view failed while writing $label BAM $out"
    else
        # No contigs: emit a valid, header-only BAM.
        samtools view -b -H -o "$out" "$BAM" || die "could not write empty $label BAM $out"
    fi
    if [[ "$INDEX" -eq 1 ]]; then
        samtools index "$out" || die "samtools index failed for $out"
    fi
    log "$label BAM: $out ($(samtools view -c -@ "$THREADS" "$out") alignments)"
}

log "input alignments: $(samtools view -c -@ "$THREADS" "$BAM")"
if [[ -n "$REF_OUT" ]]; then extract "reference" "$REF_BED" "$REF_OUT"; fi
if [[ -n "$EXO_OUT" ]]; then extract "spike-in" "$EXO_BED" "$EXO_OUT"; fi
log "done"
