//! GC bin means (`docs/AGENT_PLAN.md` section 2.3, T20).
//!
//! The one rule that is easy to get wrong: the mean is over **overlap hits**,
//! not over sites. GC windows are abutting 30 bp tiles and a BED-style site has
//! width 2, so a site whose second cytosine falls in the next window is counted
//! twice -- 44 of the 1000 BATF fixture sites are like that. A site inside `k`
//! windows contributes `k` terms to the sums.

use crate::error::{Error, Result};
use crate::intervals::{Query, StartIndex};
use crate::model::{GcWindow, Methylome, Site};

/// The per-bin result of [`bin_means`], in ascending bin order, containing only
/// populated bins -- which is what `addGCBintoMethylome` returns after its
/// `order(gcbin)`.
#[derive(Clone, Debug, PartialEq)]
pub struct GcBins {
    pub bin: Vec<u8>,
    pub mean: Vec<f64>,
    pub n_hits: Vec<u64>,
}

impl GcBins {
    pub fn len(&self) -> usize {
        self.bin.len()
    }

    pub fn is_empty(&self) -> bool {
        self.bin.is_empty()
    }

    /// The bin's index in [`Self::bin`], for a given bin label.
    pub fn position_of(&self, bin: u8) -> Option<usize> {
        self.bin.iter().position(|&b| b == bin)
    }
}

/// `addGCBintoMethylome(msites, gcdist, ignoreStrand)`.
///
/// Accumulation order is fixed by the plan: sites in input order, then windows
/// by ascending start, each as a plain sequential `f64` sum. Upstream's `mean()`
/// accumulates in long double instead, which is a difference of one or two ulp
/// (see `docs/reference-version.md`, "T04 oracle agreement"); it is inside the
/// 1e-10 the plan allows and far below anything a downstream deviation can see.
///
/// The accumulator is a fixed 256-entry table indexed by bin label, so the hot
/// loop allocates nothing and touches no `HashMap`. Bin labels above 255 cannot
/// be represented in `f64`'s exact integer range *and* in a `u8`, and
/// `read_gc_windows` rejects them at parse time.
pub fn bin_means(
    sites: &[Site],
    windows: &[GcWindow],
    index: &StartIndex,
    ignore_strand: bool,
) -> Result<GcBins> {
    let mut sum = [0.0f64; 256];
    let mut hits = [0u64; 256];
    let mut any = false;

    for site in sites {
        for w in index.any_overlaps(windows, Query::from_site(site, ignore_strand)) {
            let bin = windows[w].gc_bin as usize;
            sum[bin] += site.score;
            hits[bin] += 1;
            any = true;
        }
    }

    if !any {
        return Err(Error::run(
            "No methylation sites found in the GC distribution",
        ));
    }

    // Ascending bin order, populated bins only.
    let mut bin = Vec::new();
    let mut mean = Vec::new();
    let mut n_hits = Vec::new();
    for (i, &n) in hits.iter().enumerate() {
        if n > 0 {
            bin.push(i as u8);
            mean.push(sum[i] / n as f64);
            n_hits.push(n);
        }
    }
    Ok(GcBins { bin, mean, n_hits })
}

/// Reduce the GC windows to those overlapping at least one enhancer region
/// (`AGENT_PLAN.md` section 2.6). Done once per run, before any bin mean, so the
/// bin means themselves change with the enhancer.
pub fn reduce_windows_by_enhancer<E: crate::intervals::Interval>(
    windows: &[GcWindow],
    enhancer: &[E],
    ignore_strand: bool,
) -> Vec<GcWindow> {
    let keep = crate::intervals::subset_by_overlaps(windows, enhancer, ignore_strand);
    keep.iter()
        .zip(windows)
        .filter(|(k, _)| **k)
        .map(|(_, w)| *w)
        .collect()
}

/// Convenience wrapper: build the index and compute the means.
pub fn bin_means_indexed(
    sites: &[Site],
    windows: &[GcWindow],
    ignore_strand: bool,
) -> Result<GcBins> {
    let index = StartIndex::build(windows);
    bin_means(sites, windows, &index, ignore_strand)
}

/// Convenience wrapper taking a [`Methylome`].
pub fn bin_means_for(
    methylome: &Methylome,
    windows: &[GcWindow],
    ignore_strand: bool,
) -> Result<GcBins> {
    bin_means_indexed(&methylome.sites, windows, ignore_strand)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::Strand;

    fn w(start: i64, end: i64, bin: u8) -> GcWindow {
        GcWindow {
            chr: 0,
            start,
            end,
            strand: Strand::Star,
            gc_bin: bin,
        }
    }

    fn site(start: i64, end: i64, score: f64) -> Site {
        Site {
            chr: 0,
            start,
            end,
            strand: Strand::Star,
            score,
            coverage: 1.0,
        }
    }

    #[test]
    fn mean_is_over_hits_not_sites() {
        // Two abutting windows, one site straddling both, scores 1.0 and 0.0.
        let g = vec![w(10, 20, 1), w(21, 30, 2)];
        let sites = vec![site(20, 21, 0.5)];
        let b = bin_means_indexed(&sites, &g, true).unwrap();
        assert_eq!(b.bin, vec![1, 2]);
        assert_eq!(b.n_hits, vec![1, 1]);
        assert_eq!(b.mean, vec![0.5, 0.5]);
    }

    #[test]
    fn one_site_in_two_windows_of_the_same_bin_counts_twice() {
        let g = vec![w(10, 20, 1), w(21, 30, 1)];
        let sites = vec![site(20, 21, 0.5)];
        let b = bin_means_indexed(&sites, &g, true).unwrap();
        assert_eq!(b.bin, vec![1]);
        assert_eq!(b.n_hits, vec![2]);
        assert_eq!(b.mean, vec![0.5]);
    }

    #[test]
    fn unpopulated_bins_are_absent() {
        let g = vec![w(10, 20, 2), w(21, 30, 4)];
        let sites = vec![site(15, 16, 1.0)];
        let b = bin_means_indexed(&sites, &g, true).unwrap();
        assert_eq!(b.bin, vec![2]);
        assert_eq!(b.position_of(2), Some(0));
        assert_eq!(b.position_of(4), None);
    }

    #[test]
    fn no_hits_is_an_error() {
        let g = vec![w(10, 20, 1)];
        let sites = vec![site(100, 101, 1.0)];
        let e = bin_means_indexed(&sites, &g, true).unwrap_err();
        assert_eq!(
            e.to_string(),
            "No methylation sites found in the GC distribution"
        );
        let e = bin_means_indexed(&[], &g, true).unwrap_err();
        assert!(e.to_string().contains("No methylation sites"));
    }

    #[test]
    fn windows_are_reduced_by_the_enhancer() {
        let g = vec![w(10, 20, 1), w(100, 110, 2), w(200, 210, 3)];
        let enhancer = vec![(0u32, 105i64, 106i64, Strand::Star)];
        let kept = reduce_windows_by_enhancer(&g, &enhancer, true);
        assert_eq!(kept.len(), 1);
        assert_eq!(kept[0].gc_bin, 2);
        assert_eq!(kept[0].start, 100);
    }
}
