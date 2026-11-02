//! Fork-extension: validates BAM decoding from index-selected BGZF byte ranges.

use std::{
    io::{self, Cursor},
    num::NonZero,
    ops::Range,
};

use noodles_bam::{self as bam, bai};
use noodles_bgzf as bgzf;
use noodles_core::{Position, region::Interval};
use noodles_csi::{
    BinningIndex,
    binning_index::{Indexer, index::reference_sequence::bin::Chunk},
};
use noodles_sam::{
    self as sam,
    alignment::{
        RecordBuf,
        io::Write as _,
        record::{
            Flags,
            cigar::{Op, op::Kind},
        },
        record_buf::{Cigar, Sequence},
    },
    header::record::value::{Map, map::ReferenceSequence},
};

fn write_bam(header: &sam::Header, records: &[RecordBuf]) -> io::Result<Vec<u8>> {
    let mut writer = bam::io::Writer::new(Vec::new());

    writer.write_header(header)?;

    for record in records {
        writer.write_alignment_record(header, record)?;
    }

    writer.try_finish()?;

    Ok(writer.into_inner().into_inner())
}

fn build_index(src: &[u8]) -> io::Result<bai::Index> {
    let mut reader = bam::io::Reader::new(src);
    let header = reader.read_header()?;

    let mut indexer = Indexer::default();
    let mut chunk_start = reader.get_ref().virtual_position();

    let mut record = bam::Record::default();

    while reader.read_record(&mut record)? != 0 {
        let chunk_end = reader.get_ref().virtual_position();
        let chunk = Chunk::new(chunk_start, chunk_end);

        let alignment_context = match (
            record.reference_sequence_id().transpose()?,
            record.binning_interval().transpose()?,
        ) {
            (Some(reference_sequence_id), Some(interval)) => {
                let start = interval
                    .start()
                    .expect("BAM binning intervals must be bounded");

                let end = interval
                    .end()
                    .expect("BAM binning intervals must be bounded");

                let is_mapped = !record.flags().is_unmapped();

                Some((reference_sequence_id, start, end, is_mapped))
            }
            _ => None,
        };

        indexer.add_record(alignment_context, chunk)?;

        chunk_start = chunk_end;
    }

    Ok(indexer.build(header.reference_sequences().len()))
}

fn physical_range(src: &[u8], chunk: Chunk) -> io::Result<Range<usize>> {
    let start = usize::try_from(chunk.start().compressed())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    let end = chunk.end();

    let end_compressed = usize::try_from(end.compressed())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    let physical_end = if end.uncompressed() == 0 {
        end_compressed
    } else {
        let header_end = end_compressed
            .checked_add(bgzf::block::HEADER_SIZE)
            .ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "BGZF header range overflow")
            })?;

        let header = src.get(end_compressed..header_end).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "missing terminal BGZF block header",
            )
        })?;

        let block_size = bgzf::block::encoded_size(header)?;

        end_compressed.checked_add(block_size).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "BGZF block range overflow")
        })?
    };

    Ok(start..physical_end)
}

fn rebase_virtual_position(
    position: bgzf::VirtualPosition,
    base: u64,
) -> io::Result<bgzf::VirtualPosition> {
    let compressed = position.compressed().checked_sub(base).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            "virtual position precedes selected BGZF range",
        )
    })?;

    bgzf::VirtualPosition::try_from((compressed, position.uncompressed()))
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn rebase_chunk(chunk: Chunk, base: u64) -> io::Result<Chunk> {
    Ok(Chunk::new(
        rebase_virtual_position(chunk.start(), base)?,
        rebase_virtual_position(chunk.end(), base)?,
    ))
}

#[test]
fn test_decode_records_from_selected_bgzf_range() -> Result<(), Box<dyn std::error::Error>> {
    const RECORD_COUNT: usize = 256;
    const READ_LENGTH: usize = 1024;

    let header = sam::Header::builder()
        .add_reference_sequence(
            "sq0",
            Map::<ReferenceSequence>::new(NonZero::new(1_000_000).unwrap()),
        )
        .build();

    let cigar = Cigar::from(vec![Op::new(Kind::Match, READ_LENGTH)]);

    let sequence = Sequence::from(vec![b'A'; READ_LENGTH]);

    let records = (0..RECORD_COUNT)
        .map(|i| {
            let start = 1 + i * 2_000;

            RecordBuf::builder()
                .set_flags(Flags::empty())
                .set_reference_sequence_id(0)
                .set_alignment_start(Position::try_from(start).unwrap())
                .set_cigar(cigar.clone())
                .set_sequence(sequence.clone())
                .build()
        })
        .collect::<Vec<_>>();

    let src = write_bam(&header, &records)?;
    let index = build_index(&src)?;

    // Pick a record far enough into the file that its index chunk should not
    // require the BAM prefix/header blocks.
    let selected_record_index = 200;
    let expected_start = Position::try_from(1 + selected_record_index * 2_000)?;

    let interval = Interval::from(expected_start..=expected_start);

    let chunks = index.query(0, interval)?;

    assert!(!chunks.is_empty());

    // This is the important range-selection assertion: at least the first
    // selected chunk begins after compressed offset zero.
    assert_ne!(chunks[0].start().compressed(), 0);

    let mut actual_starts = Vec::new();

    for chunk in chunks {
        let range = physical_range(&src, chunk)?;

        let base = u64::try_from(range.start)?;
        let selected = &src[range];

        let rebased_chunk = rebase_chunk(chunk, base)?;

        let cursor = Cursor::new(selected);
        let mut bgzf_reader = bgzf::io::Reader::new(cursor);

        let query = noodles_csi::io::Query::new(&mut bgzf_reader, vec![rebased_chunk]);

        let mut reader = bam::io::Reader::from(query);
        let mut record = bam::Record::default();

        while reader.read_record(&mut record)? != 0 {
            let Some(reference_sequence_id) = record.reference_sequence_id().transpose()? else {
                continue;
            };

            if reference_sequence_id != 0 {
                continue;
            }

            let Some(record_interval) = record.binning_interval().transpose()? else {
                continue;
            };

            if interval.intersects(record_interval) {
                actual_starts.push(
                    record
                        .alignment_start()
                        .transpose()?
                        .expect("selected mapped record must have a position"),
                );
            }
        }
    }

    assert_eq!(actual_starts, [expected_start]);

    Ok(())
}
