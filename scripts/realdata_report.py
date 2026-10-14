#!/usr/bin/env python3
"""Turn a real-data run into the numbers the README figures and prose quote.

    usage: scripts/realdata_report.py <run_dir>

    <run_dir> holds the artefacts of one R-against-Rust comparison on real data:

        rust-full.csv          the Rust CLI output
        r-full.csv             the R reference output, same shape
        rust-full.time         /usr/bin/time -v for the Rust run
        r-full.time            /usr/bin/time -v for the R run
        r-full.csv.timings.tsv per-motif R wall time and motif size
        manifest.tsv           what was run: accessions, sizes, checksums

Writes docs/data/realdata_parity.tsv (one row per motif) and
docs/data/realdata_summary.json (the handful of scalars quoted in prose), so
the figures read measurements rather than transcriptions.

No dependencies beyond the standard library.
"""

import csv
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "docs" / "data"


def read_time_v(path):
    """Pull elapsed seconds and peak RSS from /usr/bin/time -v output."""
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    out = {}
    m = re.search(r"Elapsed \(wall clock\) time.*?:\s*([0-9:.]+)", text)
    if m:
        parts = [float(p) for p in m.group(1).split(":")]
        secs = 0.0
        for p in parts:
            secs = secs * 60 + p
        out["wall_s"] = secs
    m = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", text)
    if m:
        out["peak_rss_mib"] = int(m.group(1)) / 1024
    return out


def read_csv(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def main(run_dir):
    run_dir = Path(run_dir)
    rust = read_csv(run_dir / "rust-full.csv")
    r = read_csv(run_dir / "r-full.csv")

    if [x["motif"] for x in rust] != [x["motif"] for x in r]:
        raise SystemExit("motif order differs between the two outputs")

    rows = []
    worst = 0.0
    identical = 0
    for a, b in zip(rust, r):
        d = abs(float(a["deviation"]) - float(b["deviation"]))
        e = abs(float(a["expected_deviation"]) - float(b["expected_deviation"]))
        err = max(d, e)
        rows.append(
            {
                "motif": a["motif"],
                "r_deviation": b["deviation"],
                "rust_deviation": a["deviation"],
                "r_expected": b["expected_deviation"],
                "rust_expected": a["expected_deviation"],
                "abs_dev": f"{d:.6e}",
                "abs_exp": f"{e:.6e}",
                "max_abs": f"{err:.6e}",
            }
        )
        worst = max(worst, err)
        if err == 0.0:
            identical += 1

    DATA.mkdir(parents=True, exist_ok=True)
    parity = DATA / "realdata_parity.tsv"
    with open(parity, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    rt = read_time_v(run_dir / "rust-full.time")
    pt = read_time_v(run_dir / "r-full.time")

    timings = run_dir / "r-full.csv.timings.tsv"
    tsum = {}
    if timings.exists():
        with open(timings, newline="", encoding="utf-8") as f:
            tr = list(csv.DictReader(f, delimiter="\t"))
        # `seconds` is written as elapsed-since-start, so take the differences to
        # get per-motif cost. Done here rather than re-running the 38-minute R
        # job, and it keeps the report correct for timings files from either
        # version of the writer.
        elapsed = [float(x["seconds"]) for x in tr]
        delta = [b - a for a, b in zip([0.0] + elapsed, elapsed)]
        slowest = max(range(len(tr)), key=lambda i: delta[i]) if tr else 0
        tsum = {
            "motifs": len(tr),
            "tfbs_total": sum(int(x["n_tfbs"]) for x in tr),
            "loop_seconds": elapsed[-1] if elapsed else 0.0,
            "slowest_motif": tr[slowest]["motif"] if tr else "",
            "slowest_motif_seconds": delta[slowest] if tr else 0.0,
            "slowest_motif_sites": int(tr[slowest]["n_tfbs"]) if tr else 0,
        }

    devs = [float(x["rust_deviation"]) for x in rows]
    summary = {
        "motifs": len(rows),
        "max_abs_diff": worst,
        "identical_motifs": identical,
        "rust": rt,
        "r": pt,
        "r_timings": tsum,
        "speedup_wall": (pt.get("wall_s", 0) / rt["wall_s"]) if rt.get("wall_s") else 0,
        "memory_ratio": (pt.get("peak_rss_mib", 0) / rt["peak_rss_mib"])
        if rt.get("peak_rss_mib")
        else 0,
        "deviation_min": min(devs),
        "deviation_max": max(devs),
        "deviation_mean": sum(devs) / len(devs),
    }
    (DATA / "realdata_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    print(f"motifs compared      {summary['motifs']}")
    print(f"identical to R       {summary['identical_motifs']}")
    print(f"max abs difference   {worst:.3e}")
    print(f"Rust   {rt.get('wall_s', 0):.1f} s  {rt.get('peak_rss_mib', 0):.0f} MiB")
    print(f"R      {pt.get('wall_s', 0):.1f} s  {pt.get('peak_rss_mib', 0):.0f} MiB")
    print(f"wall speed-up        {summary['speedup_wall']:.1f}x")
    print(f"wrote {parity.relative_to(ROOT)}")
    print(f"wrote {(DATA / 'realdata_summary.json').relative_to(ROOT)}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    main(sys.argv[1])
