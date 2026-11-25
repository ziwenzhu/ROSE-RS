//! Pure-Rust UCSC BigWig (`0x888FFC26`) reader and writer powered by `flate2` (`zlib-rs`).
//!
//! Supports spec-compliant B+ contig-name trees, R+ spatial zoom/unzoomed index trees, and
//! `zlib`-compressed type-1 (`bedGraph`), type-2 (`variableStep`), and type-3 (`fixedStep`) data
//! sections without requiring external C `libdeflate` toolchains.

use std::collections::BTreeMap;
use std::fs::File;
use std::io::{Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

use flate2::read::ZlibDecoder;
use flate2::write::ZlibEncoder;
use flate2::Compression;
use rose_core::{integrate_bigwig_values, Interval, StitchedRegion};

use crate::atomic::AtomicWriter;
use crate::error::RoseIoError;

const BIGWIG_MAGIC_LE: u32 = 0x888F_FC26;
const BIGWIG_MAGIC_BE: u32 = 0x26FC_8F88;
const BPT_MAGIC_LE: u32 = 0x78CA_8C91;
const CIR_MAGIC_LE: u32 = 0x2468_ACE0;

/// Metadata and spatial index handle for a UCSC BigWig file.
#[derive(Debug, Clone)]
pub struct BigWigFile {
    path: PathBuf,
    big_endian: bool,
    uncompress_buf_size: u32,
    unzoomed_index_offset: u64,
    chrom_ids: BTreeMap<String, (u32, u32)>,
}

impl BigWigFile {
    /// Opens a BigWig file, validates its header, and loads its B+ chromosome dictionary.
    ///
    /// # Errors
    ///
    /// Returns [`RoseIoError`] if the file is unreadable or not a valid BigWig v3/v4 binary.
    pub fn open(path: impl AsRef<Path>) -> Result<Self, RoseIoError> {
        let p = path.as_ref().to_path_buf();
        let mut f = File::open(&p).map_err(|e| RoseIoError::io(&p, e))?;
        let mut hdr = [0_u8; 64];
        f.read_exact(&mut hdr).map_err(|e| RoseIoError::io(&p, e))?;

        let magic_le = u32::from_le_bytes([hdr[0], hdr[1], hdr[2], hdr[3]]);
        let big_endian = if magic_le == BIGWIG_MAGIC_LE {
            false
        } else if magic_le == BIGWIG_MAGIC_BE {
            true
        } else {
            return Err(RoseIoError::BinaryFormat {
                path: p,
                message: format!("unrecognized BigWig magic 0x{magic_le:08x}"),
            });
        };

        let chrom_tree_offset = read_u64(&hdr[8..16], big_endian);
        let unzoomed_index_offset = read_u64(&hdr[24..32], big_endian);
        let uncompress_buf_size = read_u32(&hdr[52..56], big_endian);

        let chrom_ids = read_bplus_chrom_tree(&mut f, &p, chrom_tree_offset, big_endian)?;
        Ok(Self {
            path: p,
            big_endian,
            uncompress_buf_size,
            unzoomed_index_offset,
            chrom_ids,
        })
    }

    /// Returns contig names present in the BigWig B+ tree.
    #[must_use]
    pub fn contig_names(&self) -> Vec<String> {
        self.chrom_ids.keys().cloned().collect()
    }

    /// Queries overlapping `(Interval, f64)` spans within `[query.start(), query.end())` on `chrom`.
    ///
    /// # Errors
    ///
    /// Returns [`RoseIoError`] on malformed R+ index nodes or corrupted zlib blocks.
    pub fn query_intervals(
        &self,
        chrom: &str,
        query: Interval,
    ) -> Result<Vec<(Interval, f64)>, RoseIoError> {
        let Some(&(chrom_id, _)) = self.chrom_ids.get(chrom) else {
            return Ok(Vec::new());
        };
        let mut f = File::open(&self.path).map_err(|e| RoseIoError::io(&self.path, e))?;
        f.seek(SeekFrom::Start(self.unzoomed_index_offset))
            .map_err(|e| RoseIoError::io(&self.path, e))?;

        let mut r_hdr = [0_u8; 48];
        f.read_exact(&mut r_hdr)
            .map_err(|e| RoseIoError::io(&self.path, e))?;
        let r_magic = read_u32(&r_hdr[0..4], self.big_endian);
        if r_magic != CIR_MAGIC_LE {
            return Err(RoseIoError::BinaryFormat {
                path: self.path.clone(),
                message: format!("invalid R+ tree magic 0x{r_magic:08x}"),
            });
        }

        let root_offset = self.unzoomed_index_offset + 48;
        let q_start = u32::try_from(query.start()).unwrap_or(u32::MAX);
        let q_end = u32::try_from(query.end()).unwrap_or(u32::MAX);
        let mut blocks: Vec<(u64, u64)> = Vec::new();
        collect_rtree_blocks(
            &mut f,
            &self.path,
            root_offset,
            self.big_endian,
            chrom_id,
            q_start,
            q_end,
            &mut blocks,
        )?;

        blocks.sort_unstable();
        blocks.dedup();

        let mut raw_buf = Vec::new();
        let mut dec_buf = Vec::new();
        let mut out = Vec::new();

        for (offset, size) in blocks {
            f.seek(SeekFrom::Start(offset))
                .map_err(|e| RoseIoError::io(&self.path, e))?;
            raw_buf.resize(size as usize, 0);
            f.read_exact(&mut raw_buf)
                .map_err(|e| RoseIoError::io(&self.path, e))?;

            let section_bytes: &[u8] = if self.uncompress_buf_size > 0 {
                dec_buf.clear();
                let mut decoder = ZlibDecoder::new(raw_buf.as_slice());
                decoder
                    .read_to_end(&mut dec_buf)
                    .map_err(|e| RoseIoError::BinaryFormat {
                        path: self.path.clone(),
                        message: format!("zlib decompression failed: {e}"),
                    })?;
                dec_buf.as_slice()
            } else {
                raw_buf.as_slice()
            };

            decode_wig_section(
                section_bytes,
                self.big_endian,
                chrom,
                chrom_id,
                query,
                &mut out,
            )?;
        }
        Ok(out)
    }

    /// Computes width-weighted `(density, total_signal)` across `regions` using resolved contigs.
    ///
    /// # Errors
    ///
    /// Propagates [`RoseIoError`] on corrupted BigWig blocks.
    pub fn compute_region_signals(
        &self,
        regions: &[StitchedRegion],
        chrom_map: &std::collections::HashMap<String, String>,
    ) -> Result<Vec<(f64, f64)>, RoseIoError> {
        let mut results = Vec::with_capacity(regions.len());
        for reg in regions {
            let target_chrom = chrom_map
                .get(&reg.chrom)
                .map_or(reg.chrom.as_str(), String::as_str);
            let vals = self.query_intervals(target_chrom, reg.interval)?;
            results.push(integrate_bigwig_values(reg.interval, &vals));
        }
        Ok(results)
    }
}

fn read_bplus_chrom_tree(
    f: &mut File,
    path: &Path,
    offset: u64,
    be: bool,
) -> Result<BTreeMap<String, (u32, u32)>, RoseIoError> {
    f.seek(SeekFrom::Start(offset))
        .map_err(|e| RoseIoError::io(path, e))?;
    let mut hdr = [0_u8; 32];
    f.read_exact(&mut hdr)
        .map_err(|e| RoseIoError::io(path, e))?;
    let magic = read_u32(&hdr[0..4], be);
    if magic != BPT_MAGIC_LE {
        return Err(RoseIoError::BinaryFormat {
            path: path.to_path_buf(),
            message: format!("invalid B+ tree magic 0x{magic:08x}"),
        });
    }
    let key_size = read_u32(&hdr[8..12], be) as usize;
    let val_size = read_u32(&hdr[12..16], be) as usize;
    let mut out = BTreeMap::new();
    read_bplus_node(f, path, offset + 32, be, key_size, val_size, &mut out)?;
    Ok(out)
}

fn read_bplus_node(
    f: &mut File,
    path: &Path,
    offset: u64,
    be: bool,
    key_size: usize,
    val_size: usize,
    out: &mut BTreeMap<String, (u32, u32)>,
) -> Result<(), RoseIoError> {
    f.seek(SeekFrom::Start(offset))
        .map_err(|e| RoseIoError::io(path, e))?;
    let mut node_hdr = [0_u8; 4];
    f.read_exact(&mut node_hdr)
        .map_err(|e| RoseIoError::io(path, e))?;
    let is_leaf = node_hdr[0] != 0;
    let count = read_u16(&node_hdr[2..4], be) as usize;

    if is_leaf {
        let entry_len = key_size + val_size;
        let mut buf = vec![0_u8; count * entry_len];
        f.read_exact(&mut buf)
            .map_err(|e| RoseIoError::io(path, e))?;
        for chunk in buf.chunks_exact(entry_len) {
            let key_bytes = &chunk[..key_size];
            let nul = key_bytes
                .iter()
                .position(|&b| b == 0)
                .unwrap_or(key_bytes.len());
            let name = String::from_utf8_lossy(&key_bytes[..nul]).into_owned();
            if val_size >= 8 {
                let cid = read_u32(&chunk[key_size..(key_size + 4)], be);
                let clen = read_u32(&chunk[(key_size + 4)..(key_size + 8)], be);
                out.insert(name, (cid, clen));
            }
        }
    } else {
        let entry_len = key_size + 8;
        let mut buf = vec![0_u8; count * entry_len];
        f.read_exact(&mut buf)
            .map_err(|e| RoseIoError::io(path, e))?;
        let children: Vec<u64> = buf
            .chunks_exact(entry_len)
            .map(|c| read_u64(&c[key_size..(key_size + 8)], be))
            .collect();
        for child_off in children {
            read_bplus_node(f, path, child_off, be, key_size, val_size, out)?;
        }
    }
    Ok(())
}

#[expect(
    clippy::too_many_arguments,
    reason = "recursive R+ tree traversal carries query bounds"
)]
fn collect_rtree_blocks(
    f: &mut File,
    path: &Path,
    offset: u64,
    be: bool,
    q_chrom: u32,
    q_start: u32,
    q_end: u32,
    out: &mut Vec<(u64, u64)>,
) -> Result<(), RoseIoError> {
    f.seek(SeekFrom::Start(offset))
        .map_err(|e| RoseIoError::io(path, e))?;
    let mut nh = [0_u8; 4];
    f.read_exact(&mut nh)
        .map_err(|e| RoseIoError::io(path, e))?;
    let is_leaf = nh[0] != 0;
    let count = read_u16(&nh[2..4], be) as usize;
    let entry_len = if is_leaf { 32 } else { 24 };
    let mut buf = vec![0_u8; count * entry_len];
    f.read_exact(&mut buf)
        .map_err(|e| RoseIoError::io(path, e))?;

    let overlaps_key = |s_chr: u32, s_base: u32, e_chr: u32, e_base: u32| -> bool {
        if q_chrom < s_chr || q_chrom > e_chr {
            return false;
        }
        if s_chr == q_chrom && q_end <= s_base {
            return false;
        }
        if e_chr == q_chrom && q_start >= e_base {
            return false;
        }
        true
    };

    if is_leaf {
        for c in buf.as_chunks::<32>().0 {
            let s_chr = read_u32(&c[0..4], be);
            let s_base = read_u32(&c[4..8], be);
            let e_chr = read_u32(&c[8..12], be);
            let e_base = read_u32(&c[12..16], be);
            if overlaps_key(s_chr, s_base, e_chr, e_base) {
                let d_off = read_u64(&c[16..24], be);
                let d_sz = read_u64(&c[24..32], be);
                out.push((d_off, d_sz));
            }
        }
    } else {
        let mut children = Vec::new();
        for c in buf.as_chunks::<24>().0 {
            let s_chr = read_u32(&c[0..4], be);
            let s_base = read_u32(&c[4..8], be);
            let e_chr = read_u32(&c[8..12], be);
            let e_base = read_u32(&c[12..16], be);
            if overlaps_key(s_chr, s_base, e_chr, e_base) {
                children.push(read_u64(&c[16..24], be));
            }
        }
        for child_off in children {
            collect_rtree_blocks(f, path, child_off, be, q_chrom, q_start, q_end, out)?;
        }
    }
    Ok(())
}

fn decode_wig_section(
    bytes: &[u8],
    be: bool,
    chrom_name: &str,
    expected_chrom_id: u32,
    query: Interval,
    out: &mut Vec<(Interval, f64)>,
) -> Result<(), RoseIoError> {
    if bytes.len() < 24 {
        return Ok(());
    }
    let sec_chrom = read_u32(&bytes[0..4], be);
    if sec_chrom != expected_chrom_id {
        return Ok(());
    }
    let chrom_start = u64::from(read_u32(&bytes[4..8], be));
    let item_step = u64::from(read_u32(&bytes[12..16], be));
    let item_span = u64::from(read_u32(&bytes[16..20], be));
    let section_type = bytes[20];
    let item_count = read_u16(&bytes[22..24], be) as usize;
    let body = &bytes[24..];

    match section_type {
        1 => {
            // bedGraph: 12 bytes per item (u32 start, u32 end, f32 val)
            for item in body.as_chunks::<12>().0.iter().take(item_count) {
                let s = u64::from(read_u32(&item[0..4], be));
                let e = u64::from(read_u32(&item[4..8], be));
                let v = f64::from(read_f32(&item[8..12], be));
                if s < e {
                    if let Ok(iv) = Interval::try_new(chrom_name, s, e) {
                        if iv.overlaps(query) {
                            out.push((iv, v));
                        }
                    }
                }
            }
        }
        2 => {
            // variableStep: 8 bytes per item (u32 start, f32 val), width = item_span
            let span = item_span.max(1);
            for item in body.as_chunks::<8>().0.iter().take(item_count) {
                let s = u64::from(read_u32(&item[0..4], be));
                let e = s + span;
                let v = f64::from(read_f32(&item[4..8], be));
                if let Ok(iv) = Interval::try_new(chrom_name, s, e) {
                    if iv.overlaps(query) {
                        out.push((iv, v));
                    }
                }
            }
        }
        3 => {
            // fixedStep: 4 bytes per item (f32 val), start = chrom_start + idx * item_step
            let step = item_step.max(1);
            let span = item_span.max(1);
            for (idx, item) in body.as_chunks::<4>().0.iter().take(item_count).enumerate() {
                let s = chrom_start + (idx as u64) * step;
                let e = s + span;
                let v = f64::from(read_f32(item, be));
                if let Ok(iv) = Interval::try_new(chrom_name, s, e) {
                    if iv.overlaps(query) {
                        out.push((iv, v));
                    }
                }
            }
        }
        _ => {}
    }
    Ok(())
}

/// Writes a spec-compliant little-endian UCSC BigWig v4 file from sorted `(chrom, Interval, f32)`
/// bedGraph records (used by integration tests and track export).
///
/// # Errors
///
/// Returns [`RoseIoError`] if writing fails.
pub fn write_bigwig(
    path: impl AsRef<Path>,
    contigs: &[(String, u32)],
    records: &[(String, Interval, f32)],
) -> Result<(), RoseIoError> {
    let p = path.as_ref();
    let chrom_to_id: BTreeMap<&str, u32> = contigs
        .iter()
        .enumerate()
        .map(|(idx, (name, _))| (name.as_str(), idx as u32))
        .collect();

    // Group records into <= 512-item zlib-compressed bedGraph sections per chromosome.
    let mut grouped: BTreeMap<u32, Vec<(u32, u32, f32)>> = BTreeMap::new();
    for (chr, iv, val) in records {
        if let Some(&cid) = chrom_to_id.get(chr.as_str()) {
            grouped.entry(cid).or_default().push((
                u32::try_from(iv.start()).unwrap_or(0),
                u32::try_from(iv.end()).unwrap_or(1),
                *val,
            ));
        }
    }
    for list in grouped.values_mut() {
        list.sort_by_key(|&(s, e, _)| (s, e));
    }

    let mut sections: Vec<(u32, u32, u32, Vec<u8>)> = Vec::new();
    let mut max_uncomp = 0_usize;
    for (cid, items) in &grouped {
        for chunk in items.chunks(512) {
            let s_min = chunk.first().map_or(0, |x| x.0);
            let e_max = chunk.iter().map(|x| x.1).max().unwrap_or(s_min + 1);
            let mut raw = Vec::with_capacity(24 + chunk.len() * 12);
            raw.extend_from_slice(&cid.to_le_bytes());
            raw.extend_from_slice(&s_min.to_le_bytes());
            raw.extend_from_slice(&e_max.to_le_bytes());
            raw.extend_from_slice(&0_u32.to_le_bytes());
            raw.extend_from_slice(&0_u32.to_le_bytes());
            raw.push(1_u8); // bedGraph
            raw.push(0_u8);
            raw.extend_from_slice(&(chunk.len() as u16).to_le_bytes());
            for &(s, e, v) in chunk {
                raw.extend_from_slice(&s.to_le_bytes());
                raw.extend_from_slice(&e.to_le_bytes());
                raw.extend_from_slice(&v.to_le_bytes());
            }
            max_uncomp = max_uncomp.max(raw.len());
            let mut enc = ZlibEncoder::new(Vec::new(), Compression::fast());
            enc.write_all(&raw).map_err(|e| RoseIoError::io(p, e))?;
            let compressed = enc.finish().map_err(|e| RoseIoError::io(p, e))?;
            sections.push((*cid, s_min, e_max, compressed));
        }
    }

    // Build B+ tree bytes
    let key_size = contigs
        .iter()
        .map(|(n, _)| n.len())
        .max()
        .unwrap_or(4)
        .max(4);
    let mut bpt = Vec::new();
    bpt.extend_from_slice(&BPT_MAGIC_LE.to_le_bytes());
    bpt.extend_from_slice(&(contigs.len().max(1) as u32).to_le_bytes());
    bpt.extend_from_slice(&(key_size as u32).to_le_bytes());
    bpt.extend_from_slice(&8_u32.to_le_bytes());
    bpt.extend_from_slice(&(contigs.len() as u64).to_le_bytes());
    bpt.extend_from_slice(&0_u64.to_le_bytes());
    bpt.push(1_u8); // leaf
    bpt.push(0_u8);
    bpt.extend_from_slice(&(contigs.len() as u16).to_le_bytes());
    let mut sorted_contigs: Vec<(usize, &(String, u32))> = contigs.iter().enumerate().collect();
    sorted_contigs.sort_by(|a, b| a.1 .0.cmp(&b.1 .0));
    for (cid, (name, clen)) in sorted_contigs {
        let mut kb = vec![0_u8; key_size];
        kb[..name.len()].copy_from_slice(name.as_bytes());
        bpt.extend_from_slice(&kb);
        bpt.extend_from_slice(&(cid as u32).to_le_bytes());
        bpt.extend_from_slice(&clen.to_le_bytes());
    }

    let chrom_tree_offset = 64_u64;
    let data_offset = chrom_tree_offset + (bpt.len() as u64);
    let mut data_blob = Vec::new();
    data_blob.extend_from_slice(&(sections.len() as u64).to_le_bytes());
    let mut sec_offsets = Vec::with_capacity(sections.len());
    for (cid, s_min, e_max, comp) in &sections {
        let off = data_offset + (data_blob.len() as u64);
        sec_offsets.push((*cid, *s_min, *e_max, off, comp.len() as u64));
        data_blob.extend_from_slice(comp);
    }

    let index_offset = data_offset + (data_blob.len() as u64);
    let mut rtree = Vec::new();
    rtree.extend_from_slice(&CIR_MAGIC_LE.to_le_bytes());
    rtree.extend_from_slice(&256_u32.to_le_bytes());
    rtree.extend_from_slice(&(sections.len() as u64).to_le_bytes());
    let first_chr = sec_offsets.first().map_or(0, |x| x.0);
    let first_s = sec_offsets.first().map_or(0, |x| x.1);
    let last_chr = sec_offsets.last().map_or(0, |x| x.0);
    let last_e = sec_offsets.last().map_or(0, |x| x.2);
    rtree.extend_from_slice(&first_chr.to_le_bytes());
    rtree.extend_from_slice(&first_s.to_le_bytes());
    rtree.extend_from_slice(&last_chr.to_le_bytes());
    rtree.extend_from_slice(&last_e.to_le_bytes());
    rtree.extend_from_slice(&(index_offset + 48).to_le_bytes()); // endFileOffset placeholder
    rtree.extend_from_slice(&1_u32.to_le_bytes());
    rtree.extend_from_slice(&0_u32.to_le_bytes());
    rtree.push(1_u8); // leaf root
    rtree.push(0_u8);
    rtree.extend_from_slice(&(sec_offsets.len() as u16).to_le_bytes());
    for (cid, s_min, e_max, off, sz) in sec_offsets {
        rtree.extend_from_slice(&cid.to_le_bytes());
        rtree.extend_from_slice(&s_min.to_le_bytes());
        rtree.extend_from_slice(&cid.to_le_bytes());
        rtree.extend_from_slice(&e_max.to_le_bytes());
        rtree.extend_from_slice(&off.to_le_bytes());
        rtree.extend_from_slice(&sz.to_le_bytes());
    }

    let mut hdr = [0_u8; 64];
    hdr[0..4].copy_from_slice(&BIGWIG_MAGIC_LE.to_le_bytes());
    hdr[4..6].copy_from_slice(&4_u16.to_le_bytes()); // version 4
    hdr[8..16].copy_from_slice(&chrom_tree_offset.to_le_bytes());
    hdr[16..24].copy_from_slice(&data_offset.to_le_bytes());
    hdr[24..32].copy_from_slice(&index_offset.to_le_bytes());
    hdr[52..56].copy_from_slice(&(max_uncomp.max(32) as u32).to_le_bytes());

    let mut aw = AtomicWriter::create(p)?;
    aw.write_all(&hdr).map_err(|e| RoseIoError::io(p, e))?;
    aw.write_all(&bpt).map_err(|e| RoseIoError::io(p, e))?;
    aw.write_all(&data_blob)
        .map_err(|e| RoseIoError::io(p, e))?;
    aw.write_all(&rtree).map_err(|e| RoseIoError::io(p, e))?;
    aw.finish()
}

#[inline]
fn read_u16(b: &[u8], be: bool) -> u16 {
    let arr = [b[0], b[1]];
    if be {
        u16::from_be_bytes(arr)
    } else {
        u16::from_le_bytes(arr)
    }
}

#[inline]
fn read_u32(b: &[u8], be: bool) -> u32 {
    let arr = [b[0], b[1], b[2], b[3]];
    if be {
        u32::from_be_bytes(arr)
    } else {
        u32::from_le_bytes(arr)
    }
}

#[inline]
fn read_u64(b: &[u8], be: bool) -> u64 {
    let arr = [b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7]];
    if be {
        u64::from_be_bytes(arr)
    } else {
        u64::from_le_bytes(arr)
    }
}

#[inline]
fn read_f32(b: &[u8], be: bool) -> f32 {
    let arr = [b[0], b[1], b[2], b[3]];
    if be {
        f32::from_be_bytes(arr)
    } else {
        f32::from_le_bytes(arr)
    }
}
