//! The observed profile and `dev_helper` (`docs/AGENT_PLAN.md` sections 2.4
//! step 4-7 and 2.8, T22).
//!
//! `dev_helper` is five lines of R with three traps:
//!
//! * `cut(x, c(-250, -200, -25, 25, 200, 250))` is right-closed, so a value
//!   exactly on a break belongs to the interval *above* it.
//! * The numerator is `mean[(k + 1) %/% 2]` in R, 1-based. For `k = 1..5` that is
//!   the 3rd, 2nd, 2nd, 1st and 1st entry. Do not "simplify" it to `k / 2`.
//! * `k == 0` gives `NA`, not `NaN`, and a zero denominator gives `Inf` or
//!   `NaN`, which are passed through rather than caught.

use crate::error::{Error, Result};
use crate::expected::Profile;
use crate::intervals::{Query, StartIndex, motif_width, resize_center};
use crate::model::{GcFreq, Range, Site};
use crate::rmath::{cut_index, round_half_even};

/// The observed `(x, value)` pairs of one motif, before `dev_helper`.
pub type ObservedProfile = Profile;

/// Step 1-3 and 5-7 of `AGENT_PLAN.md` section 2.4.
///
/// * `width` is the motif's `W`, from the **first** TFBS before any filtering.
/// * The binding sites are resized about their centres; `new_start` may be zero
///   or negative and is kept.
/// * With an enhancer, resized TFBS overlapping at least one enhancer region
///   survive.
/// * A site is *within* a TFBS when the same chromosome, `tfbs.start <= s &&
///   tfbs.end <= e` reversed -- concretely `tfbs.start <= site.start &&
///   tfbs.end >= site.end` -- and the strands are compatible. Every
///   (site, TFBS) pair is one hit, so a site inside `k` TFBS yields `k` hits.
/// * Each hit contributes `x = site.start - midpoint` and `value = site.score`.
///   `site.start`, never `site.end`; never strand-flipped.
pub fn compute_observed(
    tfbs: &[Range],
    width: i64,
    enhancer: Option<&[Range]>,
    sites: &[Site],
    ignore_strand: bool,
    motif: &str,
) -> Result<ObservedProfile> {
    let resized = resize_center(tfbs, width);
    let filtered: Vec<Range> = match enhancer {
        None => resized,
        Some(e) => {
            // The enhancer obeys the strand rule too (section 2.5), and the
            // ranges being filtered are the *resized* TFBS, not the originals.
            let keep = crate::intervals::subset_by_overlaps(&resized, e, ignore_strand);
            keep.iter()
                .zip(resized.iter())
                .filter(|(k, _)| **k)
                .map(|(_, r)| *r)
                .collect()
        }
    };

    let index = StartIndex::build(&filtered);
    let mut x = Vec::new();
    let mut value = Vec::new();
    for s in sites {
        for t in index.within_uniform(&filtered, Query::from_site(s, ignore_strand), width) {
            let tfbs = filtered[t];
            let mid = crate::rmath::midpoint(tfbs.start, tfbs.end);
            x.push((s.start - mid) as f64);
            value.push(s.score);
        }
    }

    if x.is_empty() {
        return Err(Error::run(format!(
            "No methylation sites found in the {motif} binding sites"
        )));
    }
    Ok(Profile { x, value })
}

/// `W` for a motif: the width of its first TFBS, plus 130, taken before any
/// filtering.
pub fn motif_width_of(tfbs: &[Range]) -> Result<i64> {
    tfbs.first()
        .map(motif_width)
        .ok_or_else(|| Error::run("binding sites are empty"))
}

/// `dev_helper(data)` over a profile.
///
/// Returns `None` for `NA`, which is what `k == 0` produces. `Inf` and `NaN` from
/// a zero denominator are returned as they are, which is why the return type is
/// `Option<f64>` and not a validated type.
pub fn dev_helper(profile: &Profile) -> Option<f64> {
    let mut sum = [0.0f64; 5];
    let mut n = [0u64; 5];
    for (&x, &v) in profile.x.iter().zip(profile.value.iter()) {
        if let Some(i) = cut_index(x) {
            sum[i] += v;
            n[i] += 1;
        }
    }

    // Populated intervals in ascending index order. With five fixed intervals
    // there is nothing to sort, which is why this is two fixed-size arrays.
    let mut means: Vec<f64> = Vec::with_capacity(5);
    for i in 0..5 {
        if n[i] > 0 {
            means.push(sum[i] / n[i] as f64);
        }
    }
    let k = means.len();
    if k == 0 {
        return None;
    }
    // R's `mean[(k + 1) %/% 2]`, 1-based.
    let mid = (k + 1).div_euclid(2) - 1;
    let denominator = (means[0] + means[k - 1]) / 2.0;
    Some(means[mid] / denominator)
}

/// `dev` for one sample and one motif: the observed deviation, the expected
/// deviation of the observed profile, and their difference.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Deviation {
    pub observed: Option<f64>,
    pub expected: Option<f64>,
}

impl Deviation {
    /// `deviation = D(observed) - D(expected)`, and `NA` when either side is.
    pub fn deviation(&self) -> Option<f64> {
        match (self.observed, self.expected) {
            (Some(o), Some(e)) => Some(o - e),
            _ => None,
        }
    }
}

/// The expected deviation of the expected profile: `D(expected profile)`.
///
/// Naming note: upstream's `dev_helper` is applied to `exp_meth`, the *expected*
/// profile, and that value is what it returns as `exp_dev`. So
/// `Deviation::expected` here is `D(expected profile)` and
/// `Deviation::observed` is `D(observed profile)`; `deviation()` is their
/// difference, exactly as `computeDeviation` returns it.
pub fn deviations(observed: &Profile, expected: &Profile) -> Deviation {
    Deviation {
        observed: dev_helper(observed),
        expected: dev_helper(expected),
    }
}

/// Build the expected profile for a motif from its bins and frequency matrix.
pub fn expected_profile(gcfreq: &GcFreq, bins: &crate::gc::GcBins) -> Result<Profile> {
    crate::expected::compute_expectations(gcfreq, bins)
}

/// The profile grid midpoint rule, exposed for the table-driven tests: the
/// `round` here is R's, not `f64::round`.
#[inline]
pub fn midpoint_of(start: i64, end: i64) -> f64 {
    round_half_even(end as f64 + ((start - end) as f64) / 2.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::Strand;

    fn p(x: &[f64], v: &[f64]) -> Profile {
        Profile {
            x: x.to_vec(),
            value: v.to_vec(),
        }
    }

    #[test]
    fn k_zero_is_na() {
        // Everything outside the cut range.
        let pr = p(&[-300.0, -251.0, 251.0, 300.0], &[1.0, 1.0, 1.0, 1.0]);
        assert_eq!(dev_helper(&pr), None);
    }

    #[test]
    fn k_one_is_one() {
        let pr = p(&[0.0, 1.0, -1.0], &[0.5, 0.5, 0.5]);
        let d = dev_helper(&pr).unwrap();
        assert_eq!(d, 1.0, "k = 1 divides the only mean by itself");
    }

    #[test]
    fn k_two_uses_the_first_mean() {
        // (-250,-200]: 0.5, (200,250]: 1.5 -> numerator is mean[1] -> 0.5/1.0
        let pr = p(&[-225.0, 225.0], &[0.5, 1.5]);
        assert_eq!(dev_helper(&pr).unwrap(), 0.5 / ((0.5 + 1.5) / 2.0));
    }

    #[test]
    fn k_three_and_four_use_the_second_mean() {
        // k = 3 and k = 4 both pick index 1 in 0-based terms.
        let pr3 = p(&[-225.0, 0.0, 225.0], &[0.2, 0.4, 0.6]);
        assert_eq!(dev_helper(&pr3).unwrap(), 0.4 / ((0.2 + 0.6) / 2.0));
        let pr4 = p(&[-225.0, 0.0, 100.0, 225.0], &[0.2, 0.4, 0.8, 0.6]);
        assert_eq!(dev_helper(&pr4).unwrap(), 0.4 / ((0.2 + 0.6) / 2.0));
    }

    #[test]
    fn k_five_uses_the_middle_mean() {
        let xs = [-225.0, -100.0, 0.0, 100.0, 225.0];
        let vs = [0.2, 0.4, 0.6, 0.8, 1.0];
        let pr = p(&xs, &vs);
        assert_eq!(dev_helper(&pr).unwrap(), 0.6 / ((0.2 + 1.0) / 2.0));
    }

    #[test]
    fn break_values_belong_to_the_interval_above() {
        // -200 is in (-250,-200]; -25 in (-200,-25]; 25 in (-25,25];
        // 200 in (25,200]; 250 in (200,250].
        let pr = p(
            &[-200.0, -25.0, 25.0, 200.0, 250.0],
            &[1.0, 2.0, 3.0, 4.0, 5.0],
        );
        // k = 5, numerator is the middle: 3.0; denominator (1.0 + 5.0)/2.
        assert_eq!(dev_helper(&pr).unwrap(), 3.0 / ((1.0 + 5.0) / 2.0));
    }

    #[test]
    fn zero_denominator_is_passed_through() {
        // Only interval 0 populated, and its mean is 0 -> 0 / 0 = NaN.
        let pr = p(&[-225.0], &[0.0]);
        let d = dev_helper(&pr).unwrap();
        assert!(d.is_nan(), "0/0 must be NaN, not an error");
        // A negative denominator gives -Inf.
        let pr = p(&[-225.0], &[-1.0]);
        assert_eq!(dev_helper(&pr).unwrap(), 1.0);
    }

    #[test]
    fn observed_x_is_site_start_minus_midpoint() {
        let tfbs = vec![Range {
            chr: 0,
            start: 1000,
            end: 1410,
            strand: Strand::Star,
        }];
        let width = motif_width_of(&tfbs).unwrap();
        let sites = vec![Site {
            chr: 0,
            start: 1200,
            end: 1201,
            strand: Strand::Star,
            score: 0.75,
            coverage: 1.0,
        }];
        let obs = compute_observed(&tfbs, width, None, &sites, true, "T").unwrap();
        // (411 - 541) %/% 2 == -65, so the TFBS becomes 935..1475, midpoint
        // 1205, and x = site.start - midpoint = -5. Never site.end.
        assert_eq!(obs.x, vec![-5.0]);
        assert_eq!(obs.value, vec![0.75]);
    }

    #[test]
    fn no_hits_is_an_error_naming_the_motif() {
        let tfbs = vec![Range {
            chr: 0,
            start: 1_000_000,
            end: 1_000_410,
            strand: Strand::Star,
        }];
        let width = motif_width_of(&tfbs).unwrap();
        let sites = vec![Site {
            chr: 0,
            start: 5,
            end: 6,
            strand: Strand::Star,
            score: 1.0,
            coverage: 1.0,
        }];
        let e = compute_observed(&tfbs, width, None, &sites, true, "BATF").unwrap_err();
        assert_eq!(
            e.to_string(),
            "No methylation sites found in the BATF binding sites"
        );
    }

    #[test]
    fn a_site_inside_k_tfbs_yields_k_hits() {
        let tfbs = vec![
            Range {
                chr: 0,
                start: 100,
                end: 510,
                strand: Strand::Star,
            },
            Range {
                chr: 0,
                start: 300,
                end: 710,
                strand: Strand::Star,
            },
        ];
        let width = motif_width_of(&tfbs).unwrap();
        let sites = vec![Site {
            chr: 0,
            start: 400,
            end: 401,
            strand: Strand::Star,
            score: 0.5,
            coverage: 1.0,
        }];
        let obs = compute_observed(&tfbs, width, None, &sites, true, "T").unwrap();
        assert_eq!(obs.x.len(), 2);
        assert_eq!(obs.value, vec![0.5, 0.5]);
    }

    #[test]
    fn deviations_difference_is_na_when_either_side_is() {
        let obs = p(&[300.0], &[1.0]);
        let exp = p(&[0.0], &[1.0]);
        let d = deviations(&obs, &exp);
        assert_eq!(d.observed, None);
        assert_eq!(d.expected, Some(1.0));
        assert_eq!(d.deviation(), None);
    }
}
