#!/usr/bin/env bash
# run_benchmarks.sh -- T71: the Rust and R baselines on the same data.
#
#   usage: scripts/run_benchmarks.sh <data_dir> <out_dir> [repeats]
#
# Runs each implementation three times by default and reports the median wall time
# and peak RSS, so a single noisy run cannot become a headline number. Also diffs
# the two outputs, because a benchmark of two programs that disagree is worthless.
#
# Time and memory come from /usr/bin/time -v, which reports the kernel's own
# high-water mark for RSS rather than anything sampled.

set -euo pipefail

DATA_DIR=${1:?usage: run_benchmarks.sh <data_dir> <out_dir> [repeats]}
OUT_DIR=${2:?usage: run_benchmarks.sh <data_dir> <out_dir> [repeats]}
REPEATS=${3:-3}

# The Tier A R interpreter (docs/reference-version.md). Override to use another one.
R_BIN=${METHYTFR_R:-/scratch/mdra00001/micromamba/envs/methyltfr-ref/bin/Rscript}
# Threads for both sides. R's baseline is single-threaded here because
# computeDeviation is called one motif at a time; the parallel R path is
# methylTFR::methyltfr_core, measured separately in docs/benchmarks.md.
THREADS=${METHYTFR_THREADS:-8}

REPO=$(cd "$(dirname "$0")/.." && pwd)
mkdir -p "$OUT_DIR"

measure() {
    # measure <label> <csv> <command...>
    local label=$1 csv=$2
    shift 2
    local i
    for i in $(seq 1 "$REPEATS"); do
        /usr/bin/time -v -o "$OUT_DIR/$label.$i.time" "$@" >"$csv.$i" 2>"$OUT_DIR/$label.$i.log"
    done
}

summarise() {
    local label=$1 csv=$2
    local times=() rss=()
    local i
    for i in $(seq 1 "$REPEATS"); do
        times+=("$(awk -F': ' '/Elapsed \(wall clock\) time/ {print $2}' "$OUT_DIR/$label.$i.time")")
        rss+=("$(awk -F': ' '/Maximum resident set size/ {print $2}' "$OUT_DIR/$label.$i.time")")
        cp "$csv.$i" "$OUT_DIR/$label.run$i.csv"
    done
    # /usr/bin/time -v prints h:mm:ss or m:ss.ss; convert to seconds for the median.
    local median_time median_rss
    median_time=$(printf '%s\n' "${times[@]}" | awk -F: '{
        s = $NF
        for (i = 1; i < NF; i++) s = s + 60 * $i
        print s
    }' | sort -n | awk '{a[NR] = $0} END {print a[int((NR+1)/2)]}')
    median_rss=$(printf '%s\n' "${rss[@]}" | sort -n | awk '{a[NR] = $0} END {print a[int((NR+1)/2)]}')
    printf '%-10s median wall %8.3f s   median peak RSS %8.1f MiB\n' \
        "$label" "$median_time" "$(awk -v k="$median_rss" 'BEGIN {print k / 1024}')"
}

echo "data:    $DATA_DIR"
echo "threads: $THREADS"
echo "repeats: $REPEATS"
echo

cargo build --release --manifest-path "$REPO/Cargo.toml"

echo "--- Rust"
RAYON_NUM_THREADS=$THREADS measure rust "$OUT_DIR/rust.csv" \
    "$REPO/target/release/methyltfr" run \
    --format bismarkcov \
    --annotation "$DATA_DIR/annotation" \
    "$DATA_DIR/sample0.bismarkCov"
summarise rust "$OUT_DIR/rust.csv"

echo
echo "--- R (methylTFR 0.99.9, single-threaded: computeDeviation per motif)"
measure r "$OUT_DIR/r.csv" \
    "$R_BIN" "$REPO/scripts/bench_r.R" "$DATA_DIR" "$OUT_DIR/r.csv"
summarise r "$OUT_DIR/r.csv"

echo
echo "--- output comparison"
if diff -q "$OUT_DIR/rust.run1.csv" "$OUT_DIR/r.run1.csv" >/dev/null; then
    echo "byte-identical"
else
    # Compare numerically: R's mean() accumulates in long double and its %*%
    # goes through BLAS, so the last one or two digits can differ.
    "$R_BIN" "$REPO/scripts/compare_csv.R" "$OUT_DIR/r.run1.csv" "$OUT_DIR/rust.run1.csv"
fi