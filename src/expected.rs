//! The expected profile (`docs/AGENT_PLAN.md` section 2.7, T21).
//!
//! `expected[j] = sum_i gcfreq[i][j] * mean[i]`, summed in ascending `i`, over
//! the `x[j] = round(seq(-floor(L/2), floor(L/2), length.out = L))` grid.
//!
//! Upstream computes the inner product with `t(gcfreq) %*% binMsites[, 2]`,
//! which goes through BLAS. With five terms the accumulation order is the only
//! difference from a sequential `f64` sum, worth one or two ulp; the plan's
//! tolerance is 1e-10, and `tests/expected_profile.rs` prints the measured
//! maximum.

use crate::error::{Error, Result};
use crate::gc::GcBins;
use crate::model::GcFreq;
use crate::rmath::profile_grid;

/// An `(x, value)` profile: positions and the value at each.
#[derive(Clone, Debug, PartialEq)]
pub struct Profile {
    pub x: Vec<f64>,
    pub value: Vec<f64>,
}

impl Profile {
    pub fn len(&self) -> usize {
        self.x.len()
    }

    pub fn is_empty(&self) -> bool {
        self.x.is_empty()
    }
}

/// `computeExpectations(binMsites, gcfreq)`.
///
/// `bins` must have exactly as many populated bins as `gcfreq` has rows, and
/// those bins must be `1..=rows`. Upstream only checks that both arguments are
/// matrices and lets `%*%` fail with "non-conformable arguments"; the plan
/// makes that an explicit error and adds the `1..=rows` check, which is recorded
/// in `docs/divergences.md` under "Not divergences".
pub fn compute_expectations(gcfreq: &GcFreq, bins: &GcBins) -> Result<Profile> {
    if bins.len() != gcfreq.rows {
        return Err(Error::run(format!(
            "gcfreq has {} rows but {} GC bins are populated; \
             the multiplication is non-conformable",
            gcfreq.rows,
            bins.len()
        )));
    }
    for (i, &bin) in bins.bin.iter().enumerate() {
        if bin as usize != i + 1 {
            return Err(Error::run(format!(
                "populated GC bins must be 1..={}, but position {i} is bin {bin}",
                gcfreq.rows
            )));
        }
    }

    let cols = gcfreq.cols;
    let mut value = vec![0.0f64; cols];
    for (i, &m) in bins.mean.iter().enumerate() {
        let row = gcfreq.row(i);
        for (j, v) in value.iter_mut().enumerate() {
            *v += row[j] * m;
        }
    }
    Ok(Profile {
        x: profile_grid(cols),
        value,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn bins(labels: &[u8], means: &[f64]) -> GcBins {
        GcBins {
            bin: labels.to_vec(),
            mean: means.to_vec(),
            n_hits: means.iter().map(|_| 1u64).collect(),
        }
    }

    #[test]
    fn grid_skips_zero_when_even() {
        let f = GcFreq::new(1, 4, vec![0.25; 4]).unwrap();
        let p = compute_expectations(&f, &bins(&[1], &[1.0])).unwrap();
        assert_eq!(p.x, vec![-2.0, -1.0, 1.0, 2.0]);
        assert_eq!(p.value, vec![0.25; 4]);
    }

    #[test]
    fn grid_lengths_1_2_4_5() {
        for (l, expect) in [
            (1usize, vec![-0.0]),
            (2, vec![-1.0, 1.0]),
            (4, vec![-2.0, -1.0, 1.0, 2.0]),
            (5, vec![-2.0, -1.0, 0.0, 1.0, 2.0]),
        ] {
            let f = GcFreq::new(1, l, vec![1.0; l]).unwrap();
            let p = compute_expectations(&f, &bins(&[1], &[1.0])).unwrap();
            assert_eq!(p.x, expect, "L = {l}");
            assert_eq!(p.value, vec![1.0; l]);
        }
    }

    #[test]
    fn weights_are_applied_per_bin() {
        // Two bins, L = 2: expected = gcfreq[0][j] * m0 + gcfreq[1][j] * m1.
        let f = GcFreq::new(2, 2, vec![0.5, 0.25, 0.25, 0.75]).unwrap();
        let p = compute_expectations(&f, &bins(&[1, 2], &[0.4, 0.8])).unwrap();
        assert_eq!(p.value[0], 0.5 * 0.4 + 0.25 * 0.8);
        assert_eq!(p.value[1], 0.25 * 0.4 + 0.75 * 0.8);
    }

    #[test]
    fn missing_bin_is_an_error() {
        let f = GcFreq::new(5, 3, vec![1.0 / 15.0; 15]).unwrap();
        let e = compute_expectations(&f, &bins(&[1, 2, 3, 4], &[1.0; 4])).unwrap_err();
        assert!(e.to_string().contains("non-conformable"));
    }

    #[test]
    fn bins_must_be_one_based_and_contiguous() {
        let f = GcFreq::new(3, 2, vec![0.5; 6]).unwrap();
        let e = compute_expectations(&f, &bins(&[1, 3, 5], &[1.0; 3])).unwrap_err();
        assert!(e.to_string().contains("must be 1..=3"));
    }

    #[test]
    fn odd_length_grid_keeps_negative_zero_at_one_column() {
        // L = 1: R's round(seq(-0, 0, length.out = 1)) is -0, and the bit pattern
        // differs from +0, so the test compares bits.
        let f = GcFreq::new(1, 1, vec![1.0]).unwrap();
        let p = compute_expectations(&f, &bins(&[1], &[1.0])).unwrap();
        assert_eq!(p.x[0].to_bits(), (-0.0f64).to_bits());
    }
}
