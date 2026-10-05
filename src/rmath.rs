//! R numeric semantics (`docs/AGENT_PLAN.md` section 2.1).
//!
//! Every function in this module exists because R's base numerics differ from
//! Rust's in a way that changes the last digit of a result, or in a way that
//! changes which side of a boundary a value lands on. The oracles in
//! `tests/fixtures/rmath/` were produced by R 4.3.3 and are compared bit for
//! bit by `tests/rmath.rs`.
//!
//! Two rules from the plan are load bearing here:
//!
//! * `f64::round()` rounds halves away from zero, R's `round()` rounds halves
//!   to even. Never use `f64::round()` for anything that mirrors R.
//! * R's `round(x, d)` for `d > 0` is not `(x * 1e6).round() / 1e6`. Doing that
//!   naive thing gives `0.007813` where R gives `0.007812` for `x = 1/128`.

/// R's `round(x)`: nearest integer, ties to even.
///
/// `round(2.5) == 2`, `round(0.5) == 0`, `round(-0.5) == -0`, `round(3.5) == 4`.
#[inline]
pub fn round_half_even(x: f64) -> f64 {
    // `round_ties_even` is exactly R's rule, including the sign of zero.
    x.round_ties_even()
}

/// R's `round(x, d)` for `d > 0`, with all arithmetic in `f64`.
///
/// This is a transcription of the rule in `AGENT_PLAN.md` section 2.1, which
/// was verified against R 4.3.3 on 7.5 M values with `d = 6` and zero
/// mismatches. `d` must be positive; callers that need `d <= 0` should use
/// [`round_half_even`] on a scaled value instead.
pub fn round_digits(x: f64, d: i32) -> f64 {
    // NaN, infinity and zero are returned unchanged: the algorithm below turns
    // them into nonsense, and R returns them unchanged too.
    if !x.is_finite() || x == 0.0 {
        return x;
    }
    if x < 0.0 {
        return -round_digits(-x, d);
    }
    let p = 10f64.powi(d);
    let x10 = p * x;
    let i10 = x10.floor();
    let xd = i10 / p;
    let xu = x10.ceil() / p;
    let dd = x - xd;
    let du = xu - x;
    if dd < du || (dd == du && i10 % 2.0 == 0.0) {
        xd
    } else {
        xu
    }
}

/// R's `round(x, 6)` -- the rounding every methylation score goes through
/// (`AGENT_PLAN.md` section 2.2).
#[inline]
pub fn round6(x: f64) -> f64 {
    round_digits(x, 6)
}

/// R's `seq(from, to, length.out = n)`.
///
/// The last element is exactly `to`, not `from + (n - 1) * by`, and for `n > 2`
/// with `from == to` the result is `from` repeated. For even `n` the resulting
/// grid skips zero, which is what makes the 512-column BATF grid run
/// `-256 .. -1, 1 .. 256` rather than `-256 .. 255`.
pub fn seq_len_out(from: f64, to: f64, n: usize) -> Vec<f64> {
    if n == 0 {
        return Vec::new();
    }
    if n == 1 {
        return vec![from];
    }
    if n > 2 && from == to {
        return vec![from; n];
    }
    let by = (to - from) / (n - 1) as f64;
    let mut out = Vec::with_capacity(n);
    for i in 0..n - 1 {
        out.push(from + i as f64 * by);
    }
    out.push(to);
    out
}

/// The breaks of R's `cut(x, c(-250, -200, -25, 25, 200, 250))`.
pub const CUT_BREAKS: [f64; 6] = [-250.0, -200.0, -25.0, 25.0, 200.0, 250.0];

/// R's `cut(x, c(-250, -200, -25, 25, 200, 250))` as a 0-based interval index.
///
/// Intervals are open on the left and closed on the right, so a value equal to
/// a break belongs to the interval above it:
///
/// | index | interval |
/// |---|---|
/// | 0 | (-250, -200] |
/// | 1 | (-200, -25] |
/// | 2 | (-25, 25] |
/// | 3 | (25, 200] |
/// | 4 | (200, 250] |
///
/// `None` for anything at or below `-250`, above `250`, or NaN.
#[inline]
pub fn cut_index(x: f64) -> Option<usize> {
    // The index is the number of breaks strictly below x, which puts each break
    // itself in the interval above it.
    let mut idx = 0usize;
    for b in CUT_BREAKS {
        if x > b {
            idx += 1;
        }
    }
    if idx == 0 || idx > CUT_BREAKS.len() - 1 {
        None
    } else {
        Some(idx - 1)
    }
}

/// R's `a %/% b` for `b > 0`: floor division, so `-131 %/% 2 == -66`.
#[inline]
pub fn floor_div2(a: i64) -> i64 {
    a.div_euclid(2)
}

/// The profile grid `x` for a profile of length `L`:
/// `round(seq(-floor(L/2), floor(L/2), length.out = L))`.
///
/// Odd `L` gives `-h .. h`. Even `L` skips zero.
#[inline]
pub fn profile_grid(l: usize) -> Vec<f64> {
    let h = (l / 2) as f64;
    seq_len_out(-h, h, l)
        .into_iter()
        .map(round_half_even)
        .collect()
}

/// R's `mid_point`: `round(end + (start - end) / 2)` with the half-integer
/// resolved by half-to-even.
///
/// Deliberately computed in `f64` first, as upstream does. Replacing it with
/// integer arithmetic changes the answer for every even-width TFBS.
#[inline]
pub fn midpoint(start: i64, end: i64) -> i64 {
    let v = end as f64 + ((start - end) as f64) / 2.0;
    round_half_even(v) as i64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_half_even_matches_r() {
        assert_eq!(round_half_even(2.5), 2.0);
        assert_eq!(round_half_even(0.5), 0.0);
        assert_eq!(round_half_even(-0.5).to_bits(), (-0.0f64).to_bits());
        assert_eq!(round_half_even(3.5), 4.0);
        assert_eq!(round_half_even(-2.5), -2.0);
        assert_eq!(round_half_even(2.4), 2.0);
        assert_eq!(round_half_even(2.6), 3.0);
    }

    #[test]
    fn round6_is_not_the_naive_formula() {
        // 1/128 is the case AGENT_PLAN.md section 1 calls out.
        assert_eq!(round6(1.0 / 128.0), 0.007812);
        assert_ne!(round6(1.0 / 128.0), (1.0f64 / 128.0 * 1e6).round() / 1e6);
        assert_eq!(round6(1.0 / 128.0).to_bits(), 0.007812f64.to_bits());
    }

    #[test]
    fn cut_index_is_right_closed() {
        assert_eq!(cut_index(-250.0), None);
        assert_eq!(cut_index(-249.0), Some(0));
        assert_eq!(cut_index(-200.0), Some(0));
        assert_eq!(cut_index(-199.0), Some(1));
        assert_eq!(cut_index(-25.0), Some(1));
        assert_eq!(cut_index(-24.0), Some(2));
        assert_eq!(cut_index(0.0), Some(2));
        assert_eq!(cut_index(25.0), Some(2));
        assert_eq!(cut_index(26.0), Some(3));
        assert_eq!(cut_index(200.0), Some(3));
        assert_eq!(cut_index(201.0), Some(4));
        assert_eq!(cut_index(250.0), Some(4));
        assert_eq!(cut_index(250.5), None);
        assert_eq!(cut_index(-251.0), None);
        assert_eq!(cut_index(f64::NAN), None);
    }

    #[test]
    fn seq_even_length_skips_zero() {
        let g = profile_grid(512);
        assert_eq!(g.len(), 512);
        assert_eq!(g[0], -256.0);
        assert_eq!(g[255], -1.0);
        assert_eq!(g[256], 1.0);
        assert_eq!(g[511], 256.0);
        assert!(!g.contains(&0.0));
    }

    #[test]
    fn seq_odd_length_has_zero() {
        let g = profile_grid(5);
        assert_eq!(g, vec![-2.0, -1.0, 0.0, 1.0, 2.0]);
    }

    #[test]
    fn seq_length_one_keeps_negative_zero() {
        let g = profile_grid(1);
        assert_eq!(g.len(), 1);
        assert_eq!(g[0].to_bits(), (-0.0f64).to_bits());
    }

    #[test]
    fn floor_div2_floors() {
        assert_eq!(floor_div2(-131), -66);
        assert_eq!(floor_div2(-130), -65);
        assert_eq!(floor_div2(4), 2);
    }

    #[test]
    fn midpoint_uses_half_to_even() {
        // Even width 411: the exact half is on an integer, no tie.
        assert_eq!(midpoint(1, 411), 206);
        // Even width 410: the exact result is 205.5, which rounds to 206.
        assert_eq!(midpoint(1, 410), 206);
        // Even width 412: 206.5 rounds to 206.
        assert_eq!(midpoint(1, 412), 206);
    }
}
