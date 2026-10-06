#!/usr/bin/env python3
"""Draw the README figures as SVGs, one light and one dark version of each.

    usage: scripts/plot_readme_figures.py

Benchmark numbers are the medians reported in docs/benchmarks; edit them here
when a benchmark is re-run. The parity figure reads docs/data/parity_cases.tsv,
which the differential suite writes:

    METHYLTFR_PARITY_REPORT=docs/data/parity_cases.tsv \\
        cargo test --release --test differential_vs_r

No dependencies beyond the standard library.
"""

import csv
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "docs" / "img"

# (row label, implementation, value)
TIME = [("R, 1 thread", "r", 28.5), ("Rust, 1 thread", "rust", 3.25), ("Rust, 8 threads", "rust", 1.37)]
MEMORY = [("R", "r", 953), ("Rust", "rust", 207)]
# threads -> wall seconds, 28M sites x 500 motifs
SCALING = [(1, 498.7), (2, 264.8), (4, 150.9), (8, 92.9), (16, 77.6)]
TOLERANCE_EXP = -10

THEMES = {
    "light": {"r": "#2a78d6", "rust": "#eb6834", "text": "#0b0b0b", "muted": "#52514e", "axis": "#d9d8d2"},
    "dark": {"r": "#3987e5", "rust": "#d95926", "text": "#ffffff", "muted": "#c3c2b7", "axis": "#3a3a37"},
}

WIDTH = 680
FONT = "-apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, Arial, sans-serif"
STYLE = (
    f"text{{font-family:{FONT};font-size:13px}}.h{{font-weight:600;font-size:14px}}"
    ".v{font-weight:600}.c{font-size:11.5px}"
)


def text(x, y, s, fill, cls="", anchor="start"):
    return f'<text x="{x:.1f}" y="{y:.1f}" class="{cls}" fill="{fill}" text-anchor="{anchor}">{s}</text>'


def hbar(x, y, w, h, fill):
    """Horizontal bar: square at the baseline, 4px rounded at the data end."""
    r = min(4, w / 2)
    return (
        f'<path d="M{x},{y} h{w - r:.1f} a{r},{r} 0 0 1 {r},{r} v{h - 2 * r} '
        f'a{r},{r} 0 0 1 -{r},{r} h-{w - r:.1f} z" fill="{fill}"/>'
    )


def column(x, base, w, h, fill):
    """Vertical bar: square at the baseline, 4px rounded at the top."""
    r = min(4, h / 2, w / 2)
    return (
        f'<path d="M{x:.1f},{base} v-{h - r:.1f} a{r},{r} 0 0 1 {r},-{r} h{w - 2 * r:.1f} '
        f'a{r},{r} 0 0 1 {r},{r} v{h - r:.1f} z" fill="{fill}"/>'
    )


def svg(height, label, body):
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {WIDTH} {height}" width="{WIDTH}" '
        f'height="{height}" role="img" aria-label="{label}"><style>{STYLE}</style>' + "".join(body) + "</svg>\n"
    )


# --- Rust against R -----------------------------------------------------------


def benchmark(t):
    x0, bar_max, bar_h, row = 150, 360, 22, 36

    def panel(rows, top, title, fmt, note):
        out = [text(x0, top, title, t["text"], "h")]
        peak, base = max(v for _, _, v in rows), rows[0][2]
        y = top + 14
        for label, impl, value in rows:
            w = bar_max * value / peak
            mid = y + bar_h / 2 + 4.5
            out.append(text(x0 - 12, mid, label, t["muted"], anchor="end"))
            out.append(hbar(x0, y, w, bar_h, t[impl]))
            shown = fmt(value)
            out.append(text(x0 + w + 10, mid, shown, t["text"], "v"))
            if value != base:
                # A second <text>, placed by an estimated width: tspan offsets are
                # not honoured by every SVG renderer.
                out.append(text(x0 + w + 16 + 7.6 * len(shown), mid, note(base / value), t["muted"]))
            y += row
        out.append(f'<line x1="{x0}" y1="{top + 8}" x2="{x0}" y2="{y - 8}" stroke="{t["axis"]}"/>')
        return out, y

    body = [
        f'<rect x="{x0}" y="6" width="12" height="12" rx="2" fill="{t["r"]}"/>',
        text(x0 + 18, 16.5, "methylTFR (R)", t["muted"]),
        f'<rect x="{x0 + 130}" y="6" width="12" height="12" rx="2" fill="{t["rust"]}"/>',
        text(x0 + 148, 16.5, "methylTFR-rs (Rust)", t["muted"]),
    ]
    rows, y = panel(TIME, 52, "Wall time", lambda v: f"{v:g} s", lambda k: f"{k:.1f}× faster")
    body += rows
    rows, y = panel(MEMORY, y + 22, "Peak memory", lambda v: f"{v:g} MiB", lambda k: f"{k:.1f}× less")
    body += rows
    height = y + 22
    body.append(
        text(x0, height - 6, "2M methylation sites, 50 motifs, synthetic data. "
             "Median of three runs, AMD Ryzen 7 3700X.", t["muted"], "c")
    )
    return svg(height, "Wall time 28.5 s in R, 3.25 s in Rust on one thread, 1.37 s on eight. "
               "Peak memory 953 MiB in R, 207 MiB in Rust.", body)


# --- Thread scaling -----------------------------------------------------------


def scaling(t):
    x0, base, plot_h, slot, bar_w = 150, 230, 170, 84, 24
    peak, one = SCALING[0][1], SCALING[0][1]
    body = [text(x0, 22, "Wall time by thread count", t["text"], "h")]
    for i, (threads, secs) in enumerate(SCALING):
        cx = x0 + slot * i + slot / 2
        h = plot_h * secs / peak
        body.append(column(cx - bar_w / 2, base, bar_w, h, t["rust"]))
        body.append(text(cx, base - h - 8, f"{secs:.0f} s", t["text"], "v", "middle"))
        body.append(text(cx, base + 19, f"{threads}", t["muted"], anchor="middle"))
        if threads > 1:
            body.append(text(cx, base + 37, f"{one / secs:.1f}× faster", t["muted"], "c", "middle"))
    body.append(f'<line x1="{x0}" y1="{base}" x2="{x0 + slot * len(SCALING)}" y2="{base}" stroke="{t["axis"]}"/>')
    body.append(text(x0 - 12, base + 19, "threads", t["muted"], anchor="end"))
    height = base + 70
    body.append(
        text(x0, height - 6, "28M methylation sites, 7.7M GC windows, 500 motifs × 250k binding sites, synthetic. "
             "Median of three runs.", t["muted"], "c")
    )
    label = "Wall time by threads: " + ", ".join(f"{n} threads {s:.0f} s" for n, s in SCALING) + "."
    return svg(height, label, body)


# --- Parity -------------------------------------------------------------------


def parity_counts():
    """Cases per decade of their largest absolute difference from R."""
    counts, n = {}, 0
    with open(ROOT / "docs" / "data" / "parity_cases.tsv", newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["upstream_error"] == "true":
                continue
            n += 1
            err = float(row["abs"])
            key = None if err == 0 else math.floor(math.log10(err))
            counts[key] = counts.get(key, 0) + 1
    return counts, n


def parity(t):
    counts, n = parity_counts()
    lowest = min(k for k in counts if k is not None)
    decades = list(range(lowest, TOLERANCE_EXP))
    slots = [None] + decades
    x0, base, plot_h, bar_w = 150, 210, 150, 24
    slot = 440 / len(slots)
    peak = max(counts.values())
    body = [text(x0, 22, f"Largest difference from R, per test case ({n} cases)", t["text"], "h")]
    for i, key in enumerate(slots):
        cx = x0 + slot * i + slot / 2
        c = counts.get(key, 0)
        if c:
            h = max(plot_h * c / peak, 3)
            body.append(column(cx - bar_w / 2, base, bar_w, h, t["rust"]))
            body.append(text(cx, base - h - 8, str(c), t["text"], "v", "middle"))
        label = "identical" if key is None else f"1e{key}"
        body.append(text(cx, base + 19, label, t["muted"], "c", "middle"))
    right = x0 + slot * len(slots)
    body.append(f'<line x1="{x0}" y1="{base}" x2="{right}" y2="{base}" stroke="{t["axis"]}"/>')
    body.append(f'<line x1="{right}" y1="{base - plot_h}" x2="{right}" y2="{base + 6}" stroke="{t["muted"]}"/>')
    body.append(text(right + 8, base - plot_h + 10, "test tolerance", t["muted"], "c"))
    body.append(text(right + 8, base - plot_h + 25, "1e-10", t["text"], "v"))
    body.append(text(x0 - 12, base + 19, "difference", t["muted"], "c", "end"))
    height = base + 52
    body.append(
        text(x0, height - 6, "Randomized cases generated and scored by methylTFR 0.99.9. "
             "Each bar is one decade of absolute error.", t["muted"], "c")
    )
    worst = max(k for k in counts if k is not None)
    label = (f"Of {n} randomized cases, {counts.get(None, 0)} match R exactly; the largest difference in any case "
             f"is of order 1e{worst}, against a tolerance of 1e{TOLERANCE_EXP}.")
    return svg(height, label, body)


# --- Real data ----------------------------------------------------------------

REAL_SAMPLE = "ENCODE ENCSR663MXB, primary T cell, 58.6M CpGs"
# Two short lines: at 680px wide, 11.5px text holds about 105 characters, and a
# single 190-character caption runs off the right edge.
REAL_NOTE = [
    "ENCODE ENCSR663MXB primary T cell, the paper's own validation sample. 40.4M sites at coverage >= 5,",
    "632 JASPAR2020 motifs, 263.4M binding sites. AMD Ryzen 7 3700X.",
]


def caption(t, x, y, lines):
    out = []
    for i, line in enumerate(lines):
        out.append(text(x, y + i * 15, line, t["muted"], "c"))
    return out, y + (len(lines) - 1) * 15


def load_realdata():
    """The measured real-data run, as written by scripts/realdata_report.py."""
    summary = json.loads(
        (ROOT / "docs" / "data" / "realdata_summary.json").read_text(encoding="utf-8")
    )
    with open(ROOT / "docs" / "data" / "realdata_parity.tsv", newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    return summary, rows


def decade_histogram(t, counts, n, title, note):
    """Bars per decade of absolute error, with the tolerance marked."""
    lowest = min(k for k in counts if k is not None)
    slots = [None] + list(range(lowest, TOLERANCE_EXP))
    x0, base, plot_h, bar_w = 150, 210, 150, 24
    slot = 440 / len(slots)
    peak = max(counts.values())
    body = [text(x0, 22, f"{title} ({n} motifs)", t["text"], "h")]
    for i, key in enumerate(slots):
        cx = x0 + slot * i + slot / 2
        c = counts.get(key, 0)
        if c:
            h = max(plot_h * c / peak, 3)
            body.append(column(cx - bar_w / 2, base, bar_w, h, t["rust"]))
            body.append(text(cx, base - h - 8, str(c), t["text"], "v", "middle"))
        label = "identical" if key is None else f"1e{key}"
        body.append(text(cx, base + 19, label, t["muted"], "c", "middle"))
    right = x0 + slot * len(slots)
    body.append(f'<line x1="{x0}" y1="{base}" x2="{right}" y2="{base}" stroke="{t["axis"]}"/>')
    body.append(f'<line x1="{right}" y1="{base - plot_h}" x2="{right}" y2="{base + 6}" stroke="{t["muted"]}"/>')
    body.append(text(right + 8, base - plot_h + 10, "test tolerance", t["muted"], "c"))
    body.append(text(right + 8, base - plot_h + 25, "1e-10", t["text"], "v"))
    body.append(text(x0 - 12, base + 19, "difference", t["muted"], "c", "end"))
    height = base + 52 + (len(note) - 1) * 15
    cap, _ = caption(t, x0, height - 6, note)
    body += cap
    worst = max(k for k in counts if k is not None)
    label = (
        f"Of {n} motifs on real whole-genome data, {counts.get(None, 0)} match R exactly; "
        f"the largest difference in any motif is of order 1e{worst}, against a tolerance of 1e{TOLERANCE_EXP}."
    )
    return svg(height, label, body)


def realparity(t):
    summary, rows = load_realdata()
    counts = {}
    for r in rows:
        err = float(r["max_abs"])
        key = None if err == 0 else math.floor(math.log10(err))
        counts[key] = counts.get(key, 0) + 1
    return decade_histogram(
        t, counts, len(rows), "Largest difference from R on real data", REAL_NOTE
    )


def realbench(t):
    """R against R on the real sample: wall time and peak memory."""
    s, _ = load_realdata()
    rows_time = [
        ("R, 1 thread", "r", s["r"]["wall_s"]),
        ("Rust, 16 threads", "rust", s["rust"]["wall_s"]),
    ]
    rows_mem = [
        ("R", "r", s["r"]["peak_rss_mib"]),
        ("Rust", "rust", s["rust"]["peak_rss_mib"]),
    ]
    x0, bar_max, bar_h, row = 150, 330, 22, 36

    def panel(rows, top, title, fmt, note):
        out = [text(x0, top, title, t["text"], "h")]
        peak, base = max(v for _, _, v in rows), rows[0][2]
        y = top + 14
        for label, impl, value in rows:
            w = bar_max * value / peak
            mid = y + bar_h / 2 + 4.5
            out.append(text(x0 - 12, mid, label, t["muted"], anchor="end"))
            out.append(hbar(x0, y, w, bar_h, t[impl]))
            shown = fmt(value)
            out.append(text(x0 + w + 10, mid, shown, t["text"], "v"))
            if value != base:
                out.append(text(x0 + w + 16 + 7.6 * len(shown), mid,
                                note(base / value), t["muted"]))
            y += row
        out.append(f'<line x1="{x0}" y1="{top + 8}" x2="{x0}" y2="{y - 8}" stroke="{t["axis"]}"/>')
        return out, y

    body = [
        f'<rect x="{x0}" y="6" width="12" height="12" rx="2" fill="{t["r"]}"/>',
        text(x0 + 18, 16.5, "methylTFR (R)", t["muted"]),
        f'<rect x="{x0 + 130}" y="6" width="12" height="12" rx="2" fill="{t["rust"]}"/>',
        text(x0 + 148, 16.5, "methylTFR-rs (Rust)", t["muted"]),
    ]
    rows, y = panel(rows_time, 52, "Wall time", lambda v: f"{v:.0f} s",
                    lambda k: f"{k:.0f}x faster")
    body += rows
    rows, y = panel(rows_mem, y + 22, "Peak memory", lambda v: f"{v / 1024:.1f} GiB",
                    lambda k: f"{k:.1f}x less")
    body += rows
    height = y + 22 + (len(REAL_NOTE) - 1) * 15
    cap, _ = caption(t, x0, height - 6, REAL_NOTE)
    body += cap
    label = (
        f"On real whole-genome data, wall time {s['r']['wall_s'] / 60:.0f} min in R against "
        f"{s['rust']['wall_s']:.0f} s in Rust, a factor of {s['speedup_wall']:.0f}. Peak memory "
        f"{s['r']['peak_rss_mib'] / 1024:.1f} GiB against {s['rust']['peak_rss_mib'] / 1024:.1f} GiB."
    )
    return svg(height, label, body)


def realresult(t):
    """The deviation landscape across all 632 motifs, ranked.

    This is the point of running on the paper's own sample: not a speed number
    but 632 TF activity scores from a real methylome, ordered, so the shape of
    the distribution is visible.
    """
    s, rows = load_realdata()
    devs = sorted((float(r["rust_deviation"]) for r in rows), reverse=True)
    n = len(devs)
    x0, y0, plot_w, plot_h = 60, 40, 560, 170
    lo, hi = min(devs), max(devs)
    span = hi - lo or 1.0

    def px(i):
        return x0 + plot_w * i / (n - 1)

    def py(v):
        return y0 + plot_h * (hi - v) / span

    body = [text(60, 22, f"Transcription factor activity across {n} motifs, ranked", t["text"], "h")]
    zero_y = py(0.0)
    body.append(f'<line x1="{x0}" y1="{zero_y:.1f}" x2="{x0 + plot_w}" y2="{zero_y:.1f}" '
                f'stroke="{t["muted"]}" stroke-dasharray="3 3"/>')
    pts = " ".join(f"{px(i):.1f},{py(v):.1f}" for i, v in enumerate(devs))
    area = f"{x0},{zero_y:.1f} {pts} {x0 + plot_w},{zero_y:.1f}"
    body.append(f'<polygon points="{area}" fill="{t["rust"]}" opacity="0.16"/>')
    body.append(f'<polyline points="{pts}" fill="none" stroke="{t["rust"]}" stroke-width="1.8"/>')

    for v in (hi, 0.0, lo):
        y = py(v)
        body.append(text(x0 - 8, y + 4, f"{v:+.2f}", t["muted"], anchor="end"))
    body.append(text(x0, y0 + plot_h + 20, "motifs, ranked from the highest positive deviation to the lowest",
                     t["muted"], "c"))
    body.append(text(x0 - 8, y0 - 12, "deviation", t["muted"], "c", "end"))
    height = y0 + plot_h + 40 + (len(REAL_NOTE) - 1) * 15
    cap, _ = caption(t, 60, height - 4, REAL_NOTE)
    body += cap
    label = (
        f"Deviation scores for {n} motifs on a real T-cell methylome range from {lo:+.2f} to {hi:+.2f}, "
        f"mean {s['deviation_mean']:+.2f}."
    )
    return svg(height, label, body)


def column_down(x, top, w, h, fill):
    """Vertical bar hanging from `top`: square at the top, rounded at the value end."""
    r = min(4, h / 2, w / 2)
    return (
        f'<path d="M{x:.1f},{top:.1f} v{h - r:.1f} a{r},{r} 0 0 1 {r},{r} '
        f'h{w - 2 * r:.1f} a{r},{r} 0 0 1 {r},-{r} v-{h - r:.1f} z" fill="{fill}"/>'
    )


def footprint(t):
    """This run against the values the methylTFR paper published for T cells.

    A bar per motif around a zero line: negative is methylation depleted at the
    motif centre relative to the flanks, the sign convention of the numbers in
    the paper's Figure 2E. FOSL1::JUND has no published single value -- the
    paper only states it is depleted in all four populations -- so it gets one
    bar and no comparison.
    """
    with open(ROOT / "docs" / "data" / "realdata_footprint.tsv", newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    rows = [r for r in rows if r["motif"] in ("CEBPB", "SPI1", "FOSL1::JUND")]
    LIMIT = 0.10
    x0, zero, plot_h, slot = 118, 172, 108, 168
    scale = plot_h / LIMIT

    def py(v):
        return zero - v * scale

    body = [text(40, 22, "Methylation at the motif centre, against the paper's T-cell values", t["text"], "h")]
    f'<rect x="40" y="32" width="12" height="12" rx="2" fill="{t["rust"]}"/>'
    body.append(f'<rect x="40" y="32" width="12" height="12" rx="2" fill="{t["rust"]}"/>')
    body.append(text(58, 42.5, "this run", t["muted"]))
    body.append(f'<rect x="140" y="32" width="12" height="12" rx="2" fill="{t["r"]}"/>')
    body.append(text(158, 42.5, "paper, T cells", t["muted"]))

    for i, r in enumerate(rows):
        base = x0 + slot * i
        ours = float(r["depletion"])
        theirs = float(r["paper_t_cell"]) if r["paper_t_cell"] else None
        for j, v in enumerate((ours, theirs)):
            if v is None:
                body.append(text(base + 34 + j * 56 + 15, zero + 18, "not", t["muted"], "c", "middle"))
                body.append(text(base + 34 + j * 56 + 15, zero + 32, "published", t["muted"], "c", "middle"))
                continue
            cx = base + 34 + j * 56
            h = max(abs(v) * scale, 3)
            if v < 0:
                body.append(column_down(cx, zero, 30, h, t["rust" if j == 0 else "r"]))
                body.append(text(cx + 15, zero + h + 15, f"{v:+.2f}", t["text"], "v", "middle"))
            else:
                body.append(column(cx, zero, 30, h, t["rust" if j == 0 else "r"]))
                body.append(text(cx + 15, zero - h - 8, f"+{v:.2f}", t["text"], "v", "middle"))
        body.append(text(base + 78, zero + plot_h + 34, r["motif"], t["muted"], anchor="middle"))

    body.append(f'<line x1="{x0 - 14}" y1="{zero}" x2="{x0 + slot * len(rows) - 10}" '
                f'y2="{zero}" stroke="{t["axis"]}"/>')
    for v in (LIMIT, 0.0, -LIMIT):
        y = py(v)
        label = "0" if v == 0 else f"{v:+.2f}"
        body.append(text(x0 - 20, y + 4, label, t["muted"], "c", "end"))
        if v != 0:
            body.append(f'<line x1="{x0 - 14}" y1="{y:.1f}" x2="{x0 + slot * len(rows) - 10}" '
                        f'y2="{y:.1f}" stroke="{t["axis"]}" opacity="0.45"/>')

    lines = [
        "Centre against flank methylation, whole-genome, one ENCODE T cell.",
        "The paper's summary statistic is in its analysis code, not the package, so the statistic",
        "is defined in scripts/footprint_check.R rather than copied. Same scale and sign.",
    ]
    cap_top = zero + plot_h + 54
    cap, _ = caption(t, 40, cap_top, lines)
    body += cap
    height = cap_top + (len(lines) - 1) * 15 + 10
    vals = {r["motif"]: float(r["depletion"]) for r in rows}
    label = (
        f"Methylation depletion at motif centres on a real T cell: SPI1 {vals['SPI1']:+.3f} against "
        f"the paper's -0.09, CEBPB {vals['CEBPB']:+.3f} against the paper's 0.01, and FOSL1::JUND "
        f"{vals['FOSL1::JUND']:+.3f}, which the paper reports as depleted without a single value."
    )
    return svg(height, label, body)


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    figures = [benchmark, scaling, parity]
    if (ROOT / "docs" / "data" / "realdata_summary.json").exists():
        figures += [realparity, realbench, realresult, footprint]
    for figure in figures:
        for name, theme in THEMES.items():
            path = OUT / f"{figure.__name__}-{name}.svg"
            path.write_text(figure(theme), encoding="utf-8")
            print(path.relative_to(ROOT))