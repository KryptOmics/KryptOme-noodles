//! Fork-extension: BAM index I/O.

use std::io::{self, Read};

use noodles_csi as csi;

use super::Index;
use crate::bai;

const BAI_MAGIC_NUMBER: [u8; 4] = *b"BAI\x01";
const BGZF_MAGIC_NUMBER: [u8; 4] = [0x1f, 0x8b, 0x08, 0x04];

/// A BAM index reader.
pub struct Reader<R> {
    inner: R,
}

impl<R> Reader<R> {
    /// Creates a BAM index reader.
    pub fn new(inner: R) -> Self {
        Self { inner }
    }
}

impl<R> Reader<R>
where
    R: Read,
{
    /// Reads a BAM index.
    ///
    /// The index format is detected from its magic number and decoded as
    /// either BAI or CSI.
    pub fn read_index(&mut self) -> io::Result<Index> {
        let mut magic_number = [0; 4];
        self.inner.read_exact(&mut magic_number)?;

        // Replay the magic number so the format-specific reader receives the
        // complete index from the beginning.
        let source = io::Cursor::new(magic_number).chain(&mut self.inner);
        println!("{magic_number:?}");

        match magic_number {
            BAI_MAGIC_NUMBER => {
                let mut reader = bai::io::Reader::new(source);
                let index = reader.read_index()?;

                Ok(Index::Bai(index))
            }
            BGZF_MAGIC_NUMBER => {
                let mut reader = csi::io::Reader::new(source);
                let index = reader.read_index()?;

                Ok(Index::Csi(index))
            }
            _ => Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid BAM index magic number",
            )),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reads_bai_index() -> io::Result<()> {
        let expected = bai::Index::default();

        let mut data = Vec::new();
        let mut writer = bai::io::Writer::new(&mut data);
        writer.write_index(&expected)?;

        let mut reader = Reader::new(data.as_slice());
        let actual = reader.read_index()?;

        assert_eq!(actual, Index::Bai(expected));

        Ok(())
    }

    #[test]
    fn reads_csi_index() -> io::Result<()> {
        let expected = csi::Index::default();

        let mut writer = csi::io::Writer::new(Vec::new());
        writer.write_index(&expected)?;

        let data = writer.into_inner().finish()?;

        let mut reader = Reader::new(data.as_slice());
        let actual = reader.read_index()?;

        assert_eq!(actual, Index::Csi(expected));

        Ok(())
    }

    #[test]
    fn rejects_invalid_magic_number() {
        let data = b"NDLS";

        let mut reader = Reader::new(data.as_slice());
        let err = reader.read_index().unwrap_err();

        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
    }
}
