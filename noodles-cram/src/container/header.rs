use super::ReferenceSequenceContext;

/// A CRAM container header.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Header {
    pub(crate) reference_sequence_context: ReferenceSequenceContext,
    pub(crate) record_count: usize,
    pub(crate) record_counter: u64,
    pub(crate) base_count: u64,
    pub(crate) block_count: usize,
    pub(crate) landmarks: Vec<usize>,
}

impl Header {
    /// Returns the reference sequence context.
    pub fn reference_sequence_context(&self) -> ReferenceSequenceContext {
        self.reference_sequence_context
    }

    /// Returns the number of records.
    pub fn record_count(&self) -> usize {
        self.record_count
    }

    /// Returns the record counter.
    pub fn record_counter(&self) -> u64 {
        self.record_counter
    }

    /// Returns the number of bases.
    pub fn base_count(&self) -> u64 {
        self.base_count
    }

    /// Returns the number of blocks.
    pub fn block_count(&self) -> usize {
        self.block_count
    }

    /// Returns the slice landmarks.
    pub fn landmarks(&self) -> &[usize] {
        &self.landmarks
    }
}
