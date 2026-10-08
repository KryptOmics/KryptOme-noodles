use std::io::{self, Read, Seek};

use super::{
    Reader,
    container::{read_compression_header, read_container_header, read_slice},
};
use crate::{
    StatsRecord,
    container::{CompressionHeader, Header},
};

impl<R> Reader<R>
where
    R: Read + Seek,
{
    /// Fork-extension: Reads a CRAM container header without materializing its payload.
    ///
    /// On success, the reader is positioned at the start of the container
    /// payload. The returned value is the payload length. A value of `0`
    /// denotes the CRAM EOF container.
    pub fn read_container_header(&mut self, header: &mut Header) -> io::Result<usize> {
        read_container_header(&mut self.inner, header)
    }

    /// Fork-extension: Reads and decodes the compression-header block from
    /// the current container payload.
    ///
    /// The reader must be positioned at the start of the container payload.
    pub fn read_compression_header(
        &mut self,
        header: &Header,
        data_len: usize,
    ) -> io::Result<CompressionHeader> {
        let len = header.landmarks.first().copied().unwrap_or(data_len);

        let mut buf = vec![0; len];
        self.inner.read_exact(&mut buf)?;

        let mut src = buf.as_slice();
        read_compression_header(&mut src)
    }

    /// Fork-extension: Reads one selected CRAM slice and decodes a
    /// mate-complete prefix of stats records.
    ///
    /// The reader must be positioned at the beginning of the selected slice.
    pub fn read_selective_slice_until<F>(
        &mut self,
        slice_len: usize,
        compression_header: &CompressionHeader,
        is_past_boundary: F,
    ) -> io::Result<Vec<StatsRecord>>
    where
        F: FnMut(&StatsRecord) -> bool,
    {
        let mut buf = vec![0; slice_len];
        self.inner.read_exact(&mut buf)?;

        let mut src = buf.as_slice();
        let slice = read_slice(&mut src)?;

        let (core_data_src, external_data_srcs) = slice.decode_blocks()?;

        slice.stats_records_until_boundary(
            compression_header,
            &core_data_src,
            &external_data_srcs,
            is_past_boundary,
        )
    }
}
