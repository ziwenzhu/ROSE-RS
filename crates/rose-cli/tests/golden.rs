//! End-to-end checks against the synthetic fixture in `tests/golden/`.
//!
//! `expected_rose_region_map.tsv` comes from `tests/golden/rose_reference.py`, an independent
//! re-implementation of ROSE 0.1's stitching and signal rules.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::process::Command;

fn golden_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden")
}

fn scratch_dir(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("rose-rs-test-{}-{name}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("create scratch dir");
    dir
}

fn rose(args: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_rose"))
        .args(args)
        .current_dir(golden_dir())
        .output()
        .expect("run rose binary")
}

fn assert_ok(out: &std::process::Output) {
    assert!(
        out.status.success(),
        "rose failed: {}\n{}",
        out.status,
        String::from_utf8_lossy(&out.stderr)
    );
}

/// `(chrom, start, end)` -> `(id, signal columns...)`.
fn read_map(path: &Path, cols: &[usize]) -> BTreeMap<(String, u64, u64), (String, Vec<f64>)> {
    let text = std::fs::read_to_string(path).expect("read region map");
    text.lines()
        .skip(1)
        .filter(|l| !l.is_empty())
        .map(|l| {
            let c: Vec<&str> = l.split('\t').collect();
            let key = (
                c[1].to_owned(),
                c[2].parse().expect("start"),
                c[3].parse().expect("end"),
            );
            let vals = cols
                .iter()
                .map(|&i| c[i].parse().expect("signal"))
                .collect();
            (key, (c[0].to_owned(), vals))
        })
        .collect()
}

fn find_file(dir: &Path, suffix: &str) -> PathBuf {
    std::fs::read_dir(dir)
        .expect("list output dir")
        .filter_map(Result::ok)
        .map(|e| e.path())
        .find(|p| p.to_string_lossy().ends_with(suffix))
        .expect("output file present")
}

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-6 * a.abs().max(1.0)
}

#[test]
fn compat_rose_region_map_matches_rose_reference() {
    let out_dir = scratch_dir("compat");
    let out_s = out_dir.to_string_lossy().to_string();
    assert_ok(&rose(&[
        "run",
        "--compat",
        "rose",
        "-i",
        "peaks.gff",
        "-r",
        "chip.bam",
        "-c",
        "input.bam",
        "-a",
        "hg19_refseq.ucsc",
        "-t",
        "2500",
        "-o",
        &out_s,
    ]));
    let expected = read_map(&golden_dir().join("expected_rose_region_map.tsv"), &[4, 5]);
    let got = read_map(&find_file(&out_dir, "_ENHANCER_REGION_MAP.txt"), &[6, 7]);
    assert_eq!(
        expected.keys().collect::<Vec<_>>(),
        got.keys().collect::<Vec<_>>(),
        "stitched region coordinates differ from ROSE"
    );
    for (k, (id, vals)) in &expected {
        let (gid, gvals) = &got[k];
        assert_eq!(id, gid, "region name differs at {k:?}");
        for (e, g) in vals.iter().zip(gvals) {
            assert!(close(*e, *g), "signal at {k:?}: ROSE {e}, rose-rs {g}");
        }
    }
}

#[test]
fn cram_with_reference_matches_bam() {
    let dir = scratch_dir("cram");
    let bam_out = dir.join("bam.txt").to_string_lossy().to_string();
    let cram_out = dir.join("cram.txt").to_string_lossy().to_string();
    assert_ok(&rose(&[
        "signal",
        "-i",
        "peaks.bed",
        "-r",
        "chip.bam",
        "-o",
        &bam_out,
    ]));
    assert_ok(&rose(&[
        "signal",
        "-i",
        "peaks.bed",
        "-r",
        "chip.cram",
        "--reference",
        "ref.fa",
        "-o",
        &cram_out,
    ]));
    let a = read_map(Path::new(&bam_out), &[6]);
    let b = read_map(Path::new(&cram_out), &[6]);
    assert_eq!(a.len(), b.len());
    for (k, (_, va)) in &a {
        assert!(close(va[0], b[k].1[0]), "CRAM signal differs at {k:?}");
    }
}

#[test]
fn cram_without_reference_is_a_clean_error() {
    let dir = scratch_dir("cram-noref");
    let out = dir.join("x.txt").to_string_lossy().to_string();
    let res = rose(&["signal", "-i", "peaks.bed", "-r", "chip.cram", "-o", &out]);
    assert_eq!(
        res.status.code(),
        Some(1),
        "expected an error exit, not a panic"
    );
    assert!(String::from_utf8_lossy(&res.stderr).contains("--reference"));
}

#[test]
fn nested_region_does_not_steal_reads_from_enclosing_region() {
    let dir = scratch_dir("nested");
    let alone = dir.join("alone.bed");
    let nested = dir.join("nested.bed");
    std::fs::write(&alone, "chr1\t10000\t200000\tbig\n").expect("write bed");
    std::fs::write(
        &nested,
        "chr1\t10000\t200000\tbig\nchr1\t15000\t20000\tsmall\n",
    )
    .expect("write bed");
    let a_out = dir.join("a.txt").to_string_lossy().to_string();
    let n_out = dir.join("n.txt").to_string_lossy().to_string();
    for (bed, out) in [(&alone, &a_out), (&nested, &n_out)] {
        assert_ok(&rose(&[
            "signal",
            "-i",
            &bed.to_string_lossy(),
            "-r",
            "chip.bam",
            "--floor",
            "0",
            "-o",
            out,
        ]));
    }
    let a = read_map(Path::new(&a_out), &[6]);
    let n = read_map(Path::new(&n_out), &[6]);
    let key = ("chr1".to_owned(), 10_000, 200_000);
    assert!(a[&key].1[0] > 0.0);
    assert!(
        close(a[&key].1[0], n[&key].1[0]),
        "nested region changed the enclosing signal"
    );
}
