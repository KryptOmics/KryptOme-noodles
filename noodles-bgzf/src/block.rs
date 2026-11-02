//! Fork-extension: BGZF block layout utilities.
//!
//! This module provides the fixed block geometry needed to inspect encoded
//! BGZF blocks without decoding them.
use std::io;

use crate::gz;

// These values describe the fixed BGZF extra field containing BSIZE.
// XLEN (2)
const GZIP_XLEN_SIZE: usize = 2;

// SI1 (1) + SI2 (1) + SLEN (2) + BSIZE (2)
const BGZF_XLEN: usize = 6;

/// The size of a BGZF block header.
pub const HEADER_SIZE: usize = gz::HEADER_SIZE + GZIP_XLEN_SIZE + BGZF_XLEN;

/// The maximum encoded size of a BGZF block.
pub const MAX_SIZE: usize = 1 << 16;

pub(crate) const MIN_SIZE: usize = HEADER_SIZE + gz::TRAILER_SIZE;

/// Returns the encoded BGZF block size declared by the given block header.
///
/// `src` must contain at least [`HEADER_SIZE`] bytes starting at the
/// beginning of a BGZF block.
pub fn encoded_size(src: &[u8]) -> io::Result<usize> {
    let header = src
        .get(..HEADER_SIZE)
        .ok_or_else(|| io::Error::new(io::ErrorKind::UnexpectedEof, "incomplete BGZF header"))?;

    if !is_valid_header(header) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid BGZF header",
        ));
    }

    // SAFETY: `header.len() == HEADER_SIZE >= size_of::<u16>()`.
    let bsize = header.last_chunk().map(|b| u16::from_le_bytes(*b)).unwrap();

    let block_size = usize::from(bsize) + 1;

    if block_size < MIN_SIZE {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid BGZF block size",
        ));
    }

    Ok(block_size)
}

fn is_valid_header(src: &[u8]) -> bool {
    const BGZF_CM: u8 = 0x08;
    const BGZF_FLG: u8 = 0x04;
    const BGZF_XLEN: [u8; 2] = [0x06, 0x00];
    const BGZF_SI: [u8; 2] = *b"BC";
    const BGZF_SLEN: [u8; 2] = [0x02, 0x00];

    src[0..2] == gz::MAGIC_NUMBER
        && src[2] == BGZF_CM
        && src[3] == BGZF_FLG
        && src[10..12] == BGZF_XLEN
        && src[12..14] == BGZF_SI
        && src[14..16] == BGZF_SLEN
}

#[cfg(test)]
mod tests {
    use super::*;

    const HEADER: [u8; HEADER_SIZE] = [
        0x1f, 0x8b, // ID1, ID2
        0x08, // CM = DEFLATE
        0x04, // FLG = FEXTRA
        0x00, 0x00, 0x00, 0x00, // MTIME = 0
        0x00, // XFL = 0
        0xff, // OS = 255 (unknown)
        0x06, 0x00, // XLEN = 6
        b'B', b'C', // SI1, SI2
        0x02, 0x00, // SLEN = 2
        0x1b, 0x00, // BSIZE = 27
    ];

    #[test]
    fn test_encoded_size() -> io::Result<()> {
        assert_eq!(encoded_size(&HEADER)?, 28);

        let mut header = HEADER;

        // BSIZE = 65535, so the total encoded block size is 65536.
        header[16..18].copy_from_slice(&u16::MAX.to_le_bytes());

        assert_eq!(encoded_size(&header)?, MAX_SIZE);

        Ok(())
    }

    #[test]
    fn test_encoded_size_with_incomplete_header() {
        let src = &HEADER[..HEADER_SIZE - 1];

        let err = encoded_size(src).unwrap_err();

        assert_eq!(err.kind(), io::ErrorKind::UnexpectedEof);
    }

    #[test]
    fn test_encoded_size_with_invalid_header() {
        let mut src = HEADER;

        // Break the gzip magic number.
        src[0] = 0x00;

        let err = encoded_size(&src).unwrap_err();

        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
    }

    #[test]
    fn test_encoded_size_with_invalid_block_size() {
        let mut src = HEADER;

        // BSIZE = 0, so the declared total block size is 1 byte, which cannot
        // contain a BGZF header and trailer.
        src[16..18].copy_from_slice(&0u16.to_le_bytes());

        let err = encoded_size(&src).unwrap_err();

        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
    }
}
