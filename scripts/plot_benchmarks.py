#!/usr/bin/env python3
"""Draw the README benchmark figure as two SVGs (light and dark).

    usage: scripts/plot_benchmarks.py

Numbers are the medians reported in README.md ("Performance"); edit them here
when the benchmark is re-run. No dependencies beyond the standard library.
"""

from pathlib import Path

OUT = Path(__file__).resolve().parent.parent / "docs" / "img"

# (row label, implementation, value)
TIME = [("R, 1 thread", "r", 28.5), ("Rust, 1 thread", "rust", 3.25), ("Rust, 8 threads", "rust", 1.37)]
MEMORY = [("R", "r", 953), ("Rust", "rust", 207)]

THEMES = {
    "light": {"r": "#2a78d6", "rust": "#eb6834", "text": "#0b0b0b", "muted": "#52514e", "axis": "#d9d8d2"},
    "dark": {"r": "#3987e5", "rust": "#d95926", "text": "#ffffff", "muted": "#c3c2b7", "axis": "#3a3a37"},
}

WIDTH, X0, BAR_MAX, BAR_H, ROW = 680, 150, 360, 22, 36
FONT = "-apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, Arial, sans-serif"


def bar(x, y, w, h, fill):
    """Square at the baseline, 4px rounded at the data end."""
    r = min(4, w / 2)
    return (
        f'<path d="M{x},{y} h{w - r:.1f} a{r},{r} 0 0 1 {r},{r} v{h - 2 * r} '
        f'a{r},{r} 0 0 1 -{r},{r} h-{w - r:.1f} z" fill="{fill}"/>'
    )


def panel(rows, top, title, fmt, note, t):
    out = [f'<text x="{X0}" y="{top}" class="h" fill="{t["text"]}">{title}</text>']
    peak = max(v for _, _, v in rows)
    base = rows[0][2]
    y = top + 14
    for label, impl, value in rows:
        w = BAR_MAX * value / peak
        mid = y + BAR_H / 2 + 4.5
        out.append(f'<text x="{X0 - 12}" y="{mid}" class="l" fill="{t["muted"]}" text-anchor="end">{label}</text>')
        out.append(bar(X0, y, w, BAR_H, t[impl]))
        text = fmt(value)
        tx = X0 + w + 10
        out.append(f'<text x="{tx:.1f}" y="{mid}" class="v" fill="{t["text"]}">{text}</text>')
        if value != base:
            # A second <text>, placed by an estimated width: tspan offsets are not
            # honoured by every SVG renderer.
            out.append(f'<text x="{tx + 7.6 * len(text) + 6:.1f}" y="{mid}" class="l" fill="{t["muted"]}">{note(base / value)}</text>')
        y += ROW
    out.append(f'<line x1="{X0}" y1="{top + 8}" x2="{X0}" y2="{y - 8}" stroke="{t["axis"]}" stroke-width="1"/>')
    return out, y


def figure(t):
    body = [
        f'<rect x="{X0}" y="6" width="12" height="12" rx="2" fill="{t["r"]}"/>',
        f'<text x="{X0 + 18}" y="16.5" class="l" fill="{t["muted"]}">methylTFR (R)</text>',
        f'<rect x="{X0 + 130}" y="6" width="12" height="12" rx="2" fill="{t["rust"]}"/>',
        f'<text x="{X0 + 148}" y="16.5" class="l" fill="{t["muted"]}">methylTFR-rs (Rust)</text>',
    ]
    rows, y = panel(TIME, 52, "Wall time", lambda v: f"{v:g} s", lambda k: f"{k:.1f}× faster", t)
    body += rows
    rows, y = panel(MEMORY, y + 22, "Peak memory", lambda v: f"{v:g} MiB", lambda k: f"{k:.1f}× less", t)
    body += rows
    height = y + 22
    body.append(
        f'<text x="{X0}" y="{height - 6}" class="c" fill="{t["muted"]}">2M methylation sites, 50 motifs, synthetic data. '
        "Median of three runs, AMD Ryzen 7 3700X.</text>"
    )
    style = f"text{{font-family:{FONT};font-size:13px}}.h{{font-weight:600;font-size:14px}}.v{{font-weight:600}}.c{{font-size:11.5px}}"
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {WIDTH} {height}" width="{WIDTH}" height="{height}" '
        'role="img" aria-label="Benchmark: methylTFR-rs against methylTFR. Wall time 28.5 s in R, 3.25 s in Rust on '
        'one thread, 1.37 s on eight. Peak memory 953 MiB in R, 207 MiB in Rust.">'
        f"<style>{style}</style>" + "".join(body) + "</svg>\n"
    )


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    for name, theme in THEMES.items():
        path = OUT / f"benchmark-{name}.svg"
        path.write_text(figure(theme), encoding="utf-8")
        print(path)
