//! T10: every R numeric table in `tests/fixtures/rmath/` is compared bit for bit.
//!
//! Not a tolerance. `AGENT_PLAN.md` section 3 says `rmath` tables "must be
//! bit-identical", because the difference between R's `round(x, 6)` and the
//! naive `(x * 1e6).round() / 1e6` is exactly the kind of last-digit difference
//! a tolerance would hide.

use std::path::{Path, PathBuf};

use methyltfr::rmath::{cut_index, midpoint, profile_grid, round_half_even, round6, seq_len_out};

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/rmath")
        .join(name)
}

fn rows(path: &Path) -> Vec<String> {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    text.lines()
        .skip(1)
        .filter(|l| !l.is_empty())
        .map(|l| l.to_string())
        .collect()
}

#[test]
fn round6_is_bit_identical_to_r() {
    let path = fixture("round6.tsv");
    let mut checked = 0usize;
    let mut mismatches: Vec<String> = Vec::new();
    for line in rows(&path) {
        let (x, r) = line.split_once('\t').expect("two columns");
        let x: f64 = x.parse().expect("x");
        let r: f64 = r.parse().expect("r");
        let got = round6(x);
        checked += 1;
        if got.to_bits() != r.to_bits() && mismatches.len() < 10 {
            mismatches.push(format!("round6({x:?}) = {got:?}, R says {r:?}"));
        }
    }
    assert!(
        checked >= 200_000,
        "only {checked} rows, the plan wants >= 200000"
    );
    assert!(
        mismatches.is_empty(),
        "{} of {checked} mismatches, first few:\n{}",
        mismatches.len(),
        mismatches.join("\n")
    );
    eprintln!("round6.tsv: {checked} rows, 0 mismatches");
}

#[test]
fn round6_contains_the_one_over_128_tie() {
    let text = std::fs::read_to_string(fixture("round6.tsv")).unwrap();
    assert!(
        text.contains("0.0078125\t0.0078120000000000004\n"),
        "the plan's acceptance row 0.0078125 -> 0.007812 is missing"
    );
}

#[test]
fn round_half_even_is_bit_identical_to_r() {
    let path = fixture("round0.tsv");
    let mut checked = 0usize;
    let mut mismatches: Vec<String> = Vec::new();
    for line in rows(&path) {
        let (x, r) = line.split_once('\t').expect("two columns");
        let x: f64 = x.parse().expect("x");
        let r: f64 = r.parse().expect("r");
        let got = round_half_even(x);
        checked += 1;
        if got.to_bits() != r.to_bits() && mismatches.len() < 10 {
            mismatches.push(format!("round({x:?}) = {got:?}, R says {r:?}"));
        }
    }
    assert!(checked >= 20_000, "only {checked} rows");
    assert!(
        mismatches.is_empty(),
        "{} of {checked} mismatches, first few:\n{}",
        mismatches.len(),
        mismatches.join("\n")
    );
    eprintln!("round0.tsv: {checked} rows, 0 mismatches");
}

#[test]
fn seq_grid_is_bit_identical_to_r() {
    let path = fixture("seq.tsv");
    let mut checked = 0usize;
    let mut mismatches: Vec<String> = Vec::new();
    for line in rows(&path) {
        let (l, positions) = line.split_once('\t').expect("two columns");
        let l: usize = l.parse().expect("L");
        let got = profile_grid(l);
        checked += 1;
        if got.len() != l {
            mismatches.push(format!("L={l}: produced {} values", got.len()));
            continue;
        }
        for (i, want) in positions.split(',').enumerate() {
            let want: f64 = want.parse().expect("position");
            if got[i].to_bits() != want.to_bits() {
                if mismatches.len() < 10 {
                    mismatches.push(format!(
                        "L={l} index {i}: got {:?}, R says {:?}",
                        got[i], want
                    ));
                }
                break;
            }
        }
    }
    assert_eq!(checked, 1200, "L = 1..1200");
    assert!(
        mismatches.is_empty(),
        "{} of {checked} grids differ, first few:\n{}",
        mismatches.len(),
        mismatches.join("\n")
    );
    eprintln!("seq.tsv: {checked} grids (L = 1..1200), 0 mismatches");
}

#[test]
fn seq_len_out_edge_cases() {
    assert_eq!(seq_len_out(1.0, 2.0, 0), Vec::<f64>::new());
    assert_eq!(seq_len_out(-0.0, 0.0, 1)[0].to_bits(), (-0.0f64).to_bits());
    assert_eq!(seq_len_out(3.0, 3.0, 4), vec![3.0; 4]);
    // The last element is exactly `to`, not `from + (n - 1) * by`.
    let g = seq_len_out(-256.0, 256.0, 512);
    assert_eq!(g[511], 256.0);
}

#[test]
fn cut_index_is_identical_to_r() {
    let path = fixture("cut.tsv");
    let mut checked = 0usize;
    let mut mismatches: Vec<String> = Vec::new();
    for line in rows(&path) {
        let (x, want) = line.split_once('\t').expect("two columns");
        let x: f64 = x.parse().expect("x");
        let want: Option<usize> = if want == "NA" {
            None
        } else {
            Some(want.parse::<usize>().expect("index") - 1)
        };
        checked += 1;
        if cut_index(x) != want && mismatches.len() < 10 {
            mismatches.push(format!(
                "cut_index({x}) = {:?}, R says {want:?}",
                cut_index(x)
            ));
        }
    }
    assert!(checked >= 1500, "only {checked} rows");
    assert!(
        mismatches.is_empty(),
        "{} of {checked} mismatches, first few:\n{}",
        mismatches.len(),
        mismatches.join("\n")
    );
    eprintln!("cut.tsv: {checked} rows, 0 mismatches");
}

#[test]
fn midpoint_is_identical_to_r() {
    let path = fixture("midpoint.tsv");
    let mut checked = 0usize;
    let mut mismatches: Vec<String> = Vec::new();
    for line in rows(&path) {
        let mut f = line.split('\t');
        let start: i64 = f.next().expect("start").parse().expect("start");
        let end: i64 = f.next().expect("end").parse().expect("end");
        let want: f64 = f.next().expect("mid").parse().expect("mid");
        let got = midpoint(start, end);
        checked += 1;
        if (got as f64).to_bits() != want.to_bits() && mismatches.len() < 10 {
            mismatches.push(format!(
                "midpoint({start}, {end}) = {got} ({got:?}), R says {want:?}"
            ));
        }
    }
    assert!(checked >= 8000, "only {checked} rows");
    assert!(
        mismatches.is_empty(),
        "{} of {checked} mismatches, first few:\n{}",
        mismatches.len(),
        mismatches.join("\n")
    );
    eprintln!("midpoint.tsv: {checked} rows, 0 mismatches");
}
