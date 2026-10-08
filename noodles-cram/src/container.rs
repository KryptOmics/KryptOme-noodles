//! CRAM container and fields.

pub(crate) mod block;
pub mod block_content_encoder_map;
pub mod compression_header;
mod header;
mod reference_sequence_context;
pub(crate) mod slice;

pub use self::header::Header;
pub(crate) use self::reference_sequence_context::ReferenceSequenceContext;
pub use self::{
    block_content_encoder_map::BlockContentEncoderMap, compression_header::CompressionHeader,
};
