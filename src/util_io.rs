//! Pure-Rust I/O helpers: gzip decompression, tar extraction, sha256, and HTTP(S)
//! download — the in-process replacements for the `gunzip` / `tar` / `curl` /
//! `wget` / `sha256sum` shell-outs the DB provisioning + decompression paths used
//! to spawn. `mmseqs` remains the ONLY sanctioned external binary; everything in
//! this module runs without any host tool (flate2 / tar / sha2 / ureq crates).

use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::{Component, Path, PathBuf};

// ------------------------------------------------------------------ sha256
/// Stream-hash a file with sha256; returns the lowercase hex digest.
pub fn sha256_file(path: &Path) -> Result<String, String> {
    use sha2::{Digest, Sha256};
    let f = File::open(path).map_err(|e| format!("open {}: {e}", path.display()))?;
    let mut r = BufReader::new(f);
    let mut hasher = Sha256::new();
    let mut buf = [0u8; 1 << 16];
    loop {
        let n = r.read(&mut buf).map_err(|e| format!("read {}: {e}", path.display()))?;
        if n == 0 {
            break;
        }
        hasher.update(&buf[..n]);
    }
    let mut s = String::with_capacity(64);
    for b in hasher.finalize() {
        s.push_str(&format!("{b:02x}"));
    }
    Ok(s)
}

/// Verify a file's sha256 against an expected lowercase-hex digest.
pub fn verify_sha256(path: &Path, expect: &str) -> Result<(), String> {
    let got = sha256_file(path)?;
    if !got.eq_ignore_ascii_case(expect) {
        return Err(format!("{}: sha256 {got} != expected {expect}", path.display()));
    }
    Ok(())
}

// ------------------------------------------------------------------ gzip
/// Decompress a gzip file `src` (`.gz`) to `dst` (in-process, flate2).
pub fn gunzip_file(src: &Path, dst: &Path) -> Result<(), String> {
    let f = File::open(src).map_err(|e| format!("open {}: {e}", src.display()))?;
    let mut dec = flate2::read::MultiGzDecoder::new(BufReader::new(f));
    let out = File::create(dst).map_err(|e| format!("create {}: {e}", dst.display()))?;
    let mut w = BufWriter::new(out);
    std::io::copy(&mut dec, &mut w).map_err(|e| format!("gunzip {}: {e}", src.display()))?;
    w.flush().map_err(|e| format!("flush {}: {e}", dst.display()))?;
    Ok(())
}

/// Decompress gzip bytes fully into memory (for small `.gz` payloads).
pub fn gunzip_bytes(gz: &[u8]) -> Result<Vec<u8>, String> {
    let mut dec = flate2::read::MultiGzDecoder::new(gz);
    let mut out = Vec::new();
    dec.read_to_end(&mut out).map_err(|e| format!("gunzip bytes: {e}"))?;
    Ok(out)
}

// ------------------------------------------------------------------ tar
/// Reject archive members with absolute paths or `..` traversal; return a
/// dest-relative safe path, or `None` if the entry must be skipped.
fn safe_member_path(p: &Path) -> Option<PathBuf> {
    let mut out = PathBuf::new();
    for comp in p.components() {
        match comp {
            Component::Normal(c) => out.push(c),
            Component::CurDir => {}
            // Absolute roots, prefixes, and `..` are all unsafe → reject the entry.
            Component::RootDir | Component::Prefix(_) | Component::ParentDir => return None,
        }
    }
    if out.as_os_str().is_empty() {
        None
    } else {
        Some(out)
    }
}

/// Extract a `.tar.gz` archive into `dest` in-process, guarding against
/// path-traversal (`../`, absolute members) and WITHOUT preserving ownership
/// (equivalent to `tar --no-same-owner`). Returns the number of members written.
pub fn extract_tar_gz(archive: &Path, dest: &Path) -> Result<usize, String> {
    let f = File::open(archive).map_err(|e| format!("open {}: {e}", archive.display()))?;
    let dec = flate2::read::MultiGzDecoder::new(BufReader::new(f));
    extract_tar_reader(dec, dest)
}

/// Extract a plain (uncompressed) `.tar` into `dest` with the same safety guards.
pub fn extract_tar(archive: &Path, dest: &Path) -> Result<usize, String> {
    let f = File::open(archive).map_err(|e| format!("open {}: {e}", archive.display()))?;
    extract_tar_reader(BufReader::new(f), dest)
}

fn extract_tar_reader<R: Read>(reader: R, dest: &Path) -> Result<usize, String> {
    std::fs::create_dir_all(dest).map_err(|e| format!("mkdir {}: {e}", dest.display()))?;
    let mut ar = tar::Archive::new(reader);
    // Do NOT chown extracted files to the archive's uid/gid (no-same-owner).
    ar.set_preserve_permissions(true);
    ar.set_preserve_ownerships(false);
    let mut n = 0usize;
    let entries = ar.entries().map_err(|e| format!("read tar entries: {e}"))?;
    for entry in entries {
        let mut entry = entry.map_err(|e| format!("tar entry: {e}"))?;
        let path = entry.path().map_err(|e| format!("tar entry path: {e}"))?.into_owned();
        let safe = match safe_member_path(&path) {
            Some(p) => p,
            None => {
                eprintln!("[util_io] skipping unsafe tar member: {}", path.display());
                continue;
            }
        };
        let target = dest.join(&safe);
        if let Some(parent) = target.parent() {
            std::fs::create_dir_all(parent).map_err(|e| format!("mkdir {}: {e}", parent.display()))?;
        }
        entry.unpack(&target).map_err(|e| format!("unpack {}: {e}", target.display()))?;
        n += 1;
    }
    Ok(n)
}

// ------------------------------------------------------------------ http(s)
/// Download `url` to `dst` over HTTP(S) with the pure-Rust `ureq` (rustls) client,
/// streaming to disk. Writes to `<dst>.part` then renames (atomic). Follows
/// redirects (GitHub release assets redirect to a CDN).
pub fn http_download(url: &str, dst: &Path) -> Result<(), String> {
    let resp = ureq::get(url)
        .call()
        .map_err(|e| format!("GET {url}: {e}"))?;
    let part = dst.with_extension(format!(
        "{}part",
        dst.extension().and_then(|e| e.to_str()).map(|e| format!("{e}.")).unwrap_or_default()
    ));
    {
        let out = File::create(&part).map_err(|e| format!("create {}: {e}", part.display()))?;
        let mut w = BufWriter::new(out);
        let mut rdr = resp.into_reader();
        std::io::copy(&mut rdr, &mut w).map_err(|e| format!("download {url}: {e}"))?;
        w.flush().map_err(|e| format!("flush {}: {e}", part.display()))?;
    }
    std::fs::rename(&part, dst).map_err(|e| format!("rename {}: {e}", dst.display()))?;
    Ok(())
}

/// Download `url` to `dst`, then verify its sha256 (removing the file on mismatch).
pub fn http_download_verified(url: &str, dst: &Path, sha256: &str) -> Result<(), String> {
    http_download(url, dst)?;
    if let Err(e) = verify_sha256(dst, sha256) {
        let _ = std::fs::remove_file(dst);
        return Err(e);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sha256_known_vector() {
        let p = std::env::temp_dir().join(format!("bactars_sha_{}.bin", std::process::id()));
        std::fs::write(&p, b"abc").unwrap();
        let h = sha256_file(&p).unwrap();
        let _ = std::fs::remove_file(&p);
        assert_eq!(h, "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
    }

    #[test]
    fn gzip_roundtrip() {
        use flate2::{write::GzEncoder, Compression};
        let dir = std::env::temp_dir().join(format!("bactars_gz_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let gz = dir.join("x.gz");
        let out = dir.join("x");
        let mut enc = GzEncoder::new(File::create(&gz).unwrap(), Compression::default());
        enc.write_all(b"hello gzip world").unwrap();
        enc.finish().unwrap();
        gunzip_file(&gz, &out).unwrap();
        assert_eq!(std::fs::read(&out).unwrap(), b"hello gzip world");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn tar_traversal_is_rejected() {
        // Absolute and parent-dir members must be refused by the sanitizer.
        assert!(safe_member_path(Path::new("/etc/passwd")).is_none());
        assert!(safe_member_path(Path::new("../../evil")).is_none());
        assert!(safe_member_path(Path::new("a/b/c.txt")).is_some());
    }
}
