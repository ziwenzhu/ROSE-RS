//! Atomic sibling temporary-file writer and SHA-256 input checksum helper.

use std::fs::{self, File};
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use sha2::{Digest, Sha256};

use crate::error::RoseIoError;

static TMP_SEQ: AtomicU64 = AtomicU64::new(1);

/// Buffered writer that stages output into a sibling temporary file in the same directory and
/// atomically renames it to `target_path` on [`AtomicWriter::finish`].
///
/// If dropped before `finish()`, the temporary file is removed on a best-effort basis so a failed
/// run never leaves partial or stale output behind.
pub struct AtomicWriter {
    target_path: PathBuf,
    temp_path: PathBuf,
    writer: Option<BufWriter<File>>,
}

impl AtomicWriter {
    /// Opens a sibling temporary file next to `target_path`, creating parent directories if needed.
    ///
    /// # Errors
    ///
    /// Returns [`RoseIoError::Io`] if parent directory creation or temp-file creation fails.
    pub fn create(target_path: impl AsRef<Path>) -> Result<Self, RoseIoError> {
        let target = target_path.as_ref().to_path_buf();
        let parent = target.parent().unwrap_or_else(|| Path::new("."));
        if !parent.as_os_str().is_empty() {
            fs::create_dir_all(parent).map_err(|e| RoseIoError::io(parent, e))?;
        }
        let base_name = target
            .file_name()
            .and_then(|s| s.to_str())
            .unwrap_or("output");
        let pid = std::process::id();
        let seq = TMP_SEQ.fetch_add(1, Ordering::Relaxed);
        let temp_path = parent.join(format!(".rose_tmp_{pid}_{seq}_{base_name}"));
        let file = File::create(&temp_path).map_err(|e| RoseIoError::io(&temp_path, e))?;
        Ok(Self {
            target_path: target,
            temp_path,
            writer: Some(BufWriter::with_capacity(64 * 1024, file)),
        })
    }

    /// Flushes, syncs, and atomically renames the staged file to `target_path`.
    ///
    /// # Errors
    ///
    /// Returns [`RoseIoError::Io`] if flushing or renaming fails.
    pub fn finish(mut self) -> Result<(), RoseIoError> {
        if let Some(mut w) = self.writer.take() {
            w.flush().map_err(|e| RoseIoError::io(&self.temp_path, e))?;
            let file = w
                .into_inner()
                .map_err(|e| RoseIoError::io(&self.temp_path, e.into_error()))?;
            let _ = file.sync_all();
        }
        // On Windows, `fs::rename` replaces existing regular files on NTFS; remove target first if
        // a locked/read-only edge case occurs only as fallback.
        if let Err(rename_err) = fs::rename(&self.temp_path, &self.target_path) {
            if self.target_path.exists() {
                let _ = fs::remove_file(&self.target_path);
                fs::rename(&self.temp_path, &self.target_path)
                    .map_err(|e| RoseIoError::io(&self.target_path, e))?;
            } else {
                return Err(RoseIoError::io(&self.target_path, rename_err));
            }
        }
        Ok(())
    }
}

impl Write for AtomicWriter {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        self.writer
            .as_mut()
            .expect("AtomicWriter used after finish")
            .write(buf)
    }

    fn flush(&mut self) -> std::io::Result<()> {
        self.writer
            .as_mut()
            .expect("AtomicWriter used after finish")
            .flush()
    }
}

impl Drop for AtomicWriter {
    fn drop(&mut self) {
        if self.writer.take().is_some() {
            let _ = fs::remove_file(&self.temp_path);
        }
    }
}

/// Atomically writes `bytes` to `target_path`.
///
/// # Errors
///
/// Returns [`RoseIoError::Io`] on write or rename failure.
pub fn write_atomic_bytes(target_path: impl AsRef<Path>, bytes: &[u8]) -> Result<(), RoseIoError> {
    let target = target_path.as_ref();
    let mut aw = AtomicWriter::create(target)?;
    aw.write_all(bytes)
        .map_err(|e| RoseIoError::io(target, e))?;
    aw.finish()
}

/// Computes lowercase hex SHA-256 digest of `path` for provenance (`run.json`).
///
/// # Errors
///
/// Returns [`RoseIoError::Io`] if reading `path` fails.
pub fn sha256_file(path: impl AsRef<Path>) -> Result<String, RoseIoError> {
    let p = path.as_ref();
    let file = File::open(p).map_err(|e| RoseIoError::io(p, e))?;
    let mut reader = BufReader::with_capacity(64 * 1024, file);
    let mut hasher = Sha256::new();
    let mut buf = [0_u8; 32 * 1024];
    loop {
        let n = reader.read(&mut buf).map_err(|e| RoseIoError::io(p, e))?;
        if n == 0 {
            break;
        }
        hasher.update(&buf[..n]);
    }
    let digest = hasher.finalize();
    let mut hex = String::with_capacity(64);
    for byte in digest {
        use std::fmt::Write as _;
        let _ = write!(hex, "{byte:02x}");
    }
    Ok(hex)
}
