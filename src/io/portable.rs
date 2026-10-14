//! Readers for the portable tab-separated formats (`docs/AGENT_PLAN.md`
//! section 3).
//!
//! | file | columns |
//! |---|---|
//! | `msites.tsv` | `chr start end strand score coverage` |
//! | `tfbs.tsv[.gz]` | `chr start end strand` — original, not resized, original order |
//! | `gc_windows.tsv[.gz]` | `chr start end strand gc_bin` |
//! | `gcfreq.tsv` | no header; one row per bin, `L` columns |
//! | `enhancer.tsv` | `chr start end strand` |
//! | `motifs.tsv` | `motif tfbs_path gcfreq_path`, paths relative to the manifest |
//!
//! Every file has a header line except `gcfreq.tsv`. Coordinates are copied
//! through verbatim: no BED conversion, no half-open ranges.

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use crate::error::{Error, Result};
use crate::io::{count_fields, parse_f64, parse_i64, read_lines};
use crate::model::{ChromTable, GcFreq, GcWindow, MotifAnnotation, Range, Site, parse_strand};

/// One row of the `motifs.tsv` manifest. Paths are resolved against the
/// manifest's own directory when relative, so an annotation directory is
/// relocatable.
#[derive(Clone, Debug)]
pub struct ManifestRow {
    pub motif: String,
    pub tfbs_path: PathBuf,
    pub gcfreq_path: PathBuf,
}

/// Read the `motifs.tsv` manifest, resolving its relative paths.
pub fn read_motifs(manifest: &Path) -> Result<Vec<ManifestRow>> {
    let mut lines = read_lines(manifest)?;
    let (_, header) = lines
        .next_line()?
        .ok_or_else(|| Error::file(manifest, "motifs.tsv is empty"))?;
    let cols: Vec<&str> = header.split('\t').collect();
    let want = ["motif", "tfbs_path", "gcfreq_path"];
    for (i, w) in want.iter().enumerate() {
        if cols.get(i) != Some(w) {
            return Err(Error::parse(
                manifest,
                1,
                format!("expected column {} to be '{w}'", i + 1),
            ));
        }
    }
    let dir = manifest.parent().unwrap_or_else(|| Path::new("."));
    let mut out = Vec::new();
    while let Some((no, line)) = lines.next_line()? {
        if line.is_empty() {
            continue;
        }
        if count_fields(line) < 3 {
            return Err(Error::parse(
                manifest,
                no,
                "expected 3 columns: motif, tfbs_path, gcfreq_path",
            ));
        }
        let f: Vec<&str> = line.split('\t').collect();
        out.push(ManifestRow {
            motif: f[0].to_string(),
            tfbs_path: resolve(dir, f[1]),
            gcfreq_path: resolve(dir, f[2]),
        });
    }
    Ok(out)
}

fn resolve(dir: &Path, p: &str) -> PathBuf {
    let candidate = Path::new(p);
    if candidate.is_absolute() {
        candidate.to_path_buf()
    } else {
        dir.join(candidate)
    }
}

/// Resolve `<name>` or `<name>.gz` inside `dir`, whichever exists.
///
/// Section 3 spells the annotation files `gc_windows.tsv[.gz]`, so the gzip variant
/// is part of the format, not an optional extra: a fixture directory may hold
/// either.
pub fn resolve_existing(dir: &Path, name: &str) -> Result<PathBuf> {
    let plain = dir.join(name);
    if plain.exists() {
        return Ok(plain);
    }
    let gz = dir.join(format!("{name}.gz"));
    if gz.exists() {
        return Ok(gz);
    }
    Err(Error::Missing(format!(
        "{} or {}",
        plain.display(),
        gz.display()
    )))
}

/// Read `msites.tsv` into a methylome.
pub fn read_msites(path: &Path, chroms: &mut ChromTable) -> Result<Vec<Site>> {
    let mut lines = read_lines(path)?;
    let (_, header) = lines
        .next_line()?
        .ok_or_else(|| Error::file(path, "msites.tsv is empty"))?;
    expect_header(
        path,
        1,
        header,
        &["chr", "start", "end", "strand", "score", "coverage"],
    )?;
    let mut out = Vec::with_capacity(1024);
    while let Some((no, line)) = lines.next_line()? {
        if line.is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 6 {
            return Err(Error::parse(
                path,
                no,
                format!("expected 6 columns, found {}", f.len()),
            ));
        }
        out.push(Site {
            chr: chroms.intern(f[0]),
            start: parse_i64(f[1], path, no, "start")?,
            end: parse_i64(f[2], path, no, "end")?,
            strand: parse_strand(f[3], path, no)?,
            score: parse_f64(f[4], path, no, "score")?,
            coverage: parse_f64(f[5], path, no, "coverage")?,
        });
    }
    Ok(out)
}

/// Read a TFBS file: original widths, original row order.
pub fn read_ranges(path: &Path, chroms: &mut ChromTable) -> Result<Vec<Range>> {
    let mut lines = read_lines(path)?;
    let (_, header) = lines
        .next_line()?
        .ok_or_else(|| Error::file(path, "file is empty"))?;
    expect_header(path, 1, header, &["chr", "start", "end", "strand"])?;
    let mut out = Vec::with_capacity(1024);
    while let Some((no, line)) = lines.next_line()? {
        if line.is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 4 {
            return Err(Error::parse(
                path,
                no,
                format!("expected 4 columns, found {}", f.len()),
            ));
        }
        out.push(Range {
            chr: chroms.intern(f[0]),
            start: parse_i64(f[1], path, no, "start")?,
            end: parse_i64(f[2], path, no, "end")?,
            strand: parse_strand(f[3], path, no)?,
        });
    }
    Ok(out)
}

/// Read a GC window file, including the `gc_bin` column.
pub fn read_gc_windows(path: &Path, chroms: &mut ChromTable) -> Result<Vec<GcWindow>> {
    let mut lines = read_lines(path)?;
    let (_, header) = lines
        .next_line()?
        .ok_or_else(|| Error::file(path, "gc_windows.tsv is empty"))?;
    expect_header(
        path,
        1,
        header,
        &["chr", "start", "end", "strand", "gc_bin"],
    )?;
    let mut out = Vec::with_capacity(1024);
    while let Some((no, line)) = lines.next_line()? {
        if line.is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 5 {
            return Err(Error::parse(
                path,
                no,
                format!("expected 5 columns, found {}", f.len()),
            ));
        }
        let gc_bin = parse_i64(f[4], path, no, "gc_bin")?;
        if !(1..=255).contains(&gc_bin) {
            return Err(Error::parse(
                path,
                no,
                format!("gc_bin must be 1..=255, found {gc_bin}"),
            ));
        }
        out.push(GcWindow {
            chr: chroms.intern(f[0]),
            start: parse_i64(f[1], path, no, "start")?,
            end: parse_i64(f[2], path, no, "end")?,
            strand: parse_strand(f[3], path, no)?,
            gc_bin: gc_bin as u8,
        });
    }
    Ok(out)
}

/// Read `gcfreq.tsv`: no header, one row per GC bin, all rows the same length.
pub fn read_gcfreq(path: &Path) -> Result<GcFreq> {
    let mut lines = read_lines(path)?;
    let mut rows: Vec<Vec<f64>> = Vec::new();
    let mut cols = 0usize;
    while let Some((no, line)) = lines.next_line()? {
        if line.is_empty() {
            continue;
        }
        let mut row = Vec::new();
        for field in line.split('\t') {
            row.push(parse_f64(field, path, no, "gcfreq")?);
        }
        if rows.is_empty() {
            cols = row.len();
        } else if row.len() != cols {
            return Err(Error::parse(
                path,
                no,
                format!("row has {} columns, first row has {cols}", row.len()),
            ));
        }
        rows.push(row);
    }
    if rows.is_empty() {
        return Err(Error::file(path, "gcfreq.tsv is empty"));
    }
    let n = rows.len();
    GcFreq::new(n, cols, rows.into_iter().flatten().collect())
}

fn expect_header(path: &Path, line: usize, header: &str, want: &[&str]) -> Result<()> {
    let got: Vec<&str> = header.split('\t').collect();
    if got.len() < want.len() {
        return Err(Error::parse(
            path,
            line,
            format!("expected header {}, found '{header}'", want.join("\t")),
        ));
    }
    for (i, w) in want.iter().enumerate() {
        if got[i] != *w {
            return Err(Error::parse(
                path,
                line,
                format!("expected column {} to be '{w}', found '{}'", i + 1, got[i]),
            ));
        }
    }
    Ok(())
}

/// Read a whole annotation directory: `gc_windows.tsv[.gz]` plus `motifs.tsv`.
///
/// Manifest row order is motif order, and it is preserved here because
/// `AGENT_PLAN.md` section 2.8 makes it the order of the output.
pub fn read_annotation(
    dir: &Path,
    enhancer: Option<&Path>,
    ignore_strand: bool,
    chroms: &mut ChromTable,
) -> Result<crate::model::Annotation> {
    let manifest = resolve_existing(dir, "motifs.tsv")?;
    let windows = resolve_existing(dir, "gc_windows.tsv")?;
    let gc_windows = read_gc_windows(&windows, chroms)?;

    // Two manifest rows may name the same TFBS file (the CLI integration test
    // leans on that), so read each file once and share it.
    let mut cache: HashMap<PathBuf, (Arc<Vec<Range>>, Arc<GcFreq>)> = HashMap::new();
    let mut motifs = Vec::new();
    for row in read_motifs(&manifest)? {
        if !cache.contains_key(&row.tfbs_path) {
            let tfbs = Arc::new(read_ranges(&row.tfbs_path, chroms)?);
            let gcfreq = Arc::new(read_gcfreq(&row.gcfreq_path)?);
            cache.insert(row.tfbs_path.clone(), (tfbs, gcfreq));
        }
        let (tfbs, gcfreq) = &cache[&row.tfbs_path];
        motifs.push(MotifAnnotation {
            name: row.motif,
            tfbs: Arc::clone(tfbs),
            gcfreq: Arc::clone(gcfreq),
        });
    }
    let enhancer = match enhancer {
        Some(p) => Some(read_ranges(p, chroms)?),
        None => None,
    };
    Ok(crate::model::Annotation {
        gc_windows,
        motifs,
        enhancer,
        ignore_strand,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn write(dir: &Path, name: &str, body: &str) -> PathBuf {
        let p = dir.join(name);
        let mut f = std::fs::File::create(&p).unwrap();
        f.write_all(body.as_bytes()).unwrap();
        p
    }

    fn tmpdir() -> PathBuf {
        let d = std::env::temp_dir().join(format!(
            "methyltfr-portable-{}-{:?}",
            std::process::id(),
            std::thread::current().id()
        ));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        d
    }

    #[test]
    fn msites_round_trip() {
        let d = tmpdir();
        let p = write(
            &d,
            "msites.tsv",
            "chr\tstart\tend\tstrand\tscore\tcoverage\nchr1\t10\t11\t+\t0.25\t4\n",
        );
        let mut chroms = ChromTable::new();
        let sites = read_msites(&p, &mut chroms).unwrap();
        assert_eq!(sites.len(), 1);
        assert_eq!(sites[0].width(), 2, "BED width is preserved verbatim");
        assert_eq!(sites[0].score, 0.25);
        assert_eq!(sites[0].coverage, 4.0);
        std::fs::remove_dir_all(&d).ok();
    }

    #[test]
    fn gcfreq_requires_uniform_rows() {
        let d = tmpdir();
        let good = write(&d, "g.tsv", "0.1\t0.2\n0.3\t0.4\n");
        assert_eq!(read_gcfreq(&good).unwrap().cols, 2);
        let bad = write(&d, "b.tsv", "0.1\t0.2\n0.3\n");
        let e = read_gcfreq(&bad).unwrap_err();
        assert!(e.to_string().contains("first row has 2"));
        std::fs::remove_dir_all(&d).ok();
    }

    #[test]
    fn manifest_paths_resolve_against_its_own_directory() {
        let d = tmpdir();
        write(&d, "tfbs.tsv", "chr\tstart\tend\tstrand\nchr1\t1\t5\t+\n");
        write(&d, "gcfreq.tsv", "1\n1\n");
        let m = write(
            &d,
            "motifs.tsv",
            "motif\ttfbs_path\tgcfreq_path\nBATF\ttfbs.tsv\tgcfreq.tsv\n",
        );
        let rows = read_motifs(&m).unwrap();
        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].tfbs_path, d.join("tfbs.tsv"));
        std::fs::remove_dir_all(&d).ok();
    }

    #[test]
    fn missing_column_is_reported_with_the_line() {
        let d = tmpdir();
        let p = write(
            &d,
            "msites.tsv",
            "chr\tstart\tend\tstrand\tscore\tcoverage\nchr1\t1\t2\t+\tnope\t3\n",
        );
        let mut chroms = ChromTable::new();
        let e = read_msites(&p, &mut chroms).unwrap_err();
        assert_eq!(e.line(), Some(2));
        assert!(e.to_string().contains("'nope' is not a number"));
        std::fs::remove_dir_all(&d).ok();
    }

    #[test]
    fn bad_strand_names_the_line() {
        let d = tmpdir();
        let p = write(&d, "tfbs.tsv", "chr\tstart\tend\tstrand\nchr1\t1\t5\t.\n");
        let mut chroms = ChromTable::new();
        let e = read_ranges(&p, &mut chroms).unwrap_err();
        assert_eq!(e.line(), Some(2));
        assert!(e.to_string().contains("strand values must be in"));
        std::fs::remove_dir_all(&d).ok();
    }
}
