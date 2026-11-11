//! Fork-extension: BAM index abstraction.
//!
//! A BAM file can be indexed using either a BAM index (BAI) or coordinate-sorted index (CSI).
//! This module provides a common representation for either index while preserving the existing
//! [`noodles_csi::binning_index::BinningIndex`] interface.
//!
pub mod io;

use noodles_csi::{
    self as csi,
    binning_index::{self, BinningIndex},
};

use crate::bai;

/// A BAM index.
///
/// A BAM index can be either a BAM index (BAI) or coordinate-sorted index (CSI).
/// Both variants implement the common binning-index query interface.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Index {
    /// A BAM index.
    Bai(bai::Index),

    /// A coordinate-sorted index.
    Csi(csi::Index),
}

impl BinningIndex for Index {
    fn min_shift(&self) -> u8 {
        match self {
            Self::Bai(index) => index.min_shift(),
            Self::Csi(index) => index.min_shift(),
        }
    }

    fn depth(&self) -> u8 {
        match self {
            Self::Bai(index) => index.depth(),
            Self::Csi(index) => index.depth(),
        }
    }

    fn header(&self) -> Option<&binning_index::index::Header> {
        match self {
            Self::Bai(index) => index.header(),
            Self::Csi(index) => index.header(),
        }
    }

    fn reference_sequences(
        &self,
    ) -> Box<dyn Iterator<Item = &dyn binning_index::ReferenceSequence> + '_> {
        match self {
            Self::Bai(index) => BinningIndex::reference_sequences(index),
            Self::Csi(index) => BinningIndex::reference_sequences(index),
        }
    }

    fn unplaced_unmapped_record_count(&self) -> Option<u64> {
        match self {
            Self::Bai(index) => index.unplaced_unmapped_record_count(),
            Self::Csi(index) => index.unplaced_unmapped_record_count(),
        }
    }

    fn query(
        &self,
        reference_sequence_id: usize,
        interval: noodles_core::region::Interval,
    ) -> std::io::Result<Vec<binning_index::index::reference_sequence::bin::Chunk>> {
        match self {
            Self::Bai(index) => index.query(reference_sequence_id, interval),
            Self::Csi(index) => index.query(reference_sequence_id, interval),
        }
    }

    fn last_first_record_start_position(&self) -> Option<noodles_bgzf::VirtualPosition> {
        match self {
            Self::Bai(index) => index.last_first_record_start_position(),
            Self::Csi(index) => index.last_first_record_start_position(),
        }
    }
}

impl Index {
    /// Returns the underlying BAI index, if this is a BAI index.
    pub fn as_bai(&self) -> Option<&bai::Index> {
        match self {
            Self::Bai(index) => Some(index),
            Self::Csi(_) => None,
        }
    }

    /// Returns the underlying CSI index, if this is a CSI index.
    pub fn as_csi(&self) -> Option<&csi::Index> {
        match self {
            Self::Bai(_) => None,
            Self::Csi(index) => Some(index),
        }
    }
}

impl From<bai::Index> for Index {
    fn from(index: bai::Index) -> Self {
        Self::Bai(index)
    }
}

impl From<csi::Index> for Index {
    fn from(index: csi::Index) -> Self {
        Self::Csi(index)
    }
}

#[cfg(test)]
mod tests {
    use noodles_csi::binning_index::{BinningIndex, BinningIndexExt as _};

    use super::*;

    #[test]
    fn delegates_bai_binning_index() {
        let inner = bai::Index::default();
        let index = Index::from(inner.clone());

        assert_eq!(index.min_shift(), inner.min_shift());
        assert_eq!(index.depth(), inner.depth());
        assert_eq!(index.header(), inner.header());

        assert_eq!(
            index.reference_sequences().count(),
            inner.reference_sequences().len()
        );

        assert_eq!(
            index.unplaced_unmapped_record_count(),
            inner.unplaced_unmapped_record_count()
        );

        assert_eq!(
            index.last_first_record_start_position(),
            inner.last_first_record_start_position()
        );

        let expected = inner.query_u64(0, 8, 13).unwrap_err();
        let actual = index.query_u64(0, 8, 13).unwrap_err();

        assert_eq!(actual.kind(), expected.kind());
        assert_eq!(actual.to_string(), expected.to_string());
    }

    #[test]
    fn delegates_csi_binning_index() {
        let inner = csi::Index::default();
        let index = Index::from(inner.clone());

        assert_eq!(index.min_shift(), inner.min_shift());
        assert_eq!(index.depth(), inner.depth());
        assert_eq!(index.header(), inner.header());

        assert_eq!(
            index.reference_sequences().count(),
            inner.reference_sequences().len()
        );

        assert_eq!(
            index.unplaced_unmapped_record_count(),
            inner.unplaced_unmapped_record_count()
        );

        assert_eq!(
            index.last_first_record_start_position(),
            inner.last_first_record_start_position()
        );

        let expected = inner.query_u64(0, 8, 13).unwrap_err();
        let actual = index.query_u64(0, 8, 13).unwrap_err();

        assert_eq!(actual.kind(), expected.kind());
        assert_eq!(actual.to_string(), expected.to_string());
    }

    #[test]
    fn gets_bai_index() {
        let index = Index::Bai(bai::Index::default());

        assert!(index.as_bai().is_some());
        assert!(index.as_csi().is_none());
    }

    #[test]
    fn gets_csi_index() {
        let index = Index::Csi(csi::Index::default());

        assert!(index.as_bai().is_none());
        assert!(index.as_csi().is_some());
    }
}
