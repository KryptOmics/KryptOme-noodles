#![allow(missing_docs)]
use std::{
    fs::File,
    io::{self, SeekFrom},
    path::PathBuf,
};

use noodles_core::Position;
use noodles_cram::{self as cram, StatsRecord, container::Header, io::reader::Container};
use noodles_fasta::{self as fasta, repository::adapters::IndexedReader};

#[test]
/// Tests that direct stats decoding matches the native full-record path.
///
/// Each decoded `StatsRecord` is compared against one derived from the
/// corresponding fully decoded CRAM record.
fn stats_records_match_full_records() -> Result<(), Box<dyn std::error::Error>> {
    let data_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data");

    let cram_path = data_dir.join("NA12878.chr22.cram");
    let reference_path = data_dir.join("hg38.chr22.fasta.gz");

    let reference_reader =
        fasta::io::indexed_reader::Builder::default().build_from_path(reference_path)?;

    let repository = fasta::Repository::new(IndexedReader::new(reference_reader));

    let mut reader = File::open(cram_path).map(cram::io::Reader::new)?;

    let header = reader.read_header()?;

    let mut container = Container::default();
    let mut container_index = 0;

    while reader.read_container(&mut container)? != 0 {
        let compression_header = container.compression_header()?;

        for (slice_index, result) in container.slices().enumerate() {
            let slice = result?;

            let (core_data_src, external_data_srcs) = slice.decode_blocks()?;

            //
            // Native path:
            //
            // Slice::records()
            //     -> cram::Record
            //     -> StatsRecord::try_from_full_record()
            //
            let expected = slice
                .records(
                    repository.clone(),
                    &header,
                    &compression_header,
                    &core_data_src,
                    &external_data_srcs,
                )?
                .iter()
                .map(|record| StatsRecord::try_from_full_record(&header, record))
                .collect::<io::Result<Vec<_>>>()?;

            // New direct stats path:
            //
            // Slice::stats_records()
            //     -> StatsRecord
            //
            let actual =
                slice.stats_records(&compression_header, &core_data_src, &external_data_srcs)?;

            assert_eq!(
                actual.len(),
                expected.len(),
                "record count mismatch at container {container_index}, \
                 slice {slice_index}: actual={}, expected={}",
                actual.len(),
                expected.len(),
            );

            for (record_index, (actual_record, expected_record)) in
                actual.iter().zip(&expected).enumerate()
            {
                assert_stats_record_equivalence(
                    actual_record,
                    expected_record,
                    container_index,
                    slice_index,
                    record_index,
                );
            }
        }

        container_index += 1;
    }

    Ok(())
}

fn assert_stats_record_equivalence(
    actual: &StatsRecord,
    expected: &StatsRecord,
    container_index: usize,
    slice_index: usize,
    record_index: usize,
) {
    let context =
        format!("container {container_index}, slice {slice_index}, record {record_index}");

    assert_eq!(
        actual.bam_flags(),
        expected.bam_flags(),
        "{context}: bam_flags"
    );

    assert_eq!(
        actual.reference_id(),
        expected.reference_id(),
        "{context}: reference_id"
    );

    assert_eq!(
        actual.alignment_start(),
        expected.alignment_start(),
        "{context}: alignment_start"
    );

    // Unlike the native CRAM Record implementation, StatsRecord treats an
    // unmapped record as having no reference alignment span, even when the
    // record is placed on a reference for sorting.
    if actual.bam_flags().is_unmapped() {
        assert_eq!(actual.alignment_span(), None, "{context}: alignment_span");
    } else {
        assert_eq!(
            actual.alignment_span(),
            expected.alignment_span(),
            "{context}: alignment_span"
        );
    }

    assert_eq!(
        actual.mapping_quality(),
        expected.mapping_quality(),
        "{context}: mapping_quality"
    );

    assert_eq!(
        actual.template_length(),
        expected.template_length(),
        "{context}: template_length"
    );

    assert_eq!(
        actual.mate_reference_id(),
        expected.mate_reference_id(),
        "{context}: mate_reference_id"
    );

    assert_eq!(
        actual.mate_alignment_start(),
        expected.mate_alignment_start(),
        "{context}: mate_alignment_start"
    );

    assert_eq!(
        actual.read_length(),
        expected.read_length(),
        "{context}: read_length"
    );

    assert_eq!(
        actual.feature_summary(),
        expected.feature_summary(),
        "{context}: feature_summary"
    );
}

#[test]
/// Tests that `stats_records_until` is identical to full stats decoding
/// when the stopping predicate never signals the end of the query range.
fn stats_records_until_without_boundary_matches_full_records() -> io::Result<()> {
    let data_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data");

    let cram_path = data_dir.join("NA12878.chr22.cram");

    let mut reader = File::open(cram_path).map(cram::io::Reader::new)?;

    let _ = reader.read_header()?;

    let mut container = Container::default();

    while reader.read_container(&mut container)? != 0 {
        let compression_header = container.compression_header()?;

        for result in container.slices() {
            let slice = result?;

            let (core_data_src, external_data_srcs) = slice.decode_blocks()?;

            // Decode all records using the existing path.
            let full_records =
                slice.stats_records(&compression_header, &core_data_src, &external_data_srcs)?;

            // Never signal that the query boundary has been passed, so the selective
            // path must decode the complete slice.
            let select_records = slice.stats_records_until_boundary(
                &compression_header,
                &core_data_src,
                &external_data_srcs,
                |_| false, // no boundary signal
            )?;

            assert_eq!(select_records, full_records);
        }
    }

    Ok(())
}

fn is_past_boundary(record: &StatsRecord, reference_id: usize, alignment_start: Position) -> bool {
    let Some(record_reference_id) = record.reference_id() else {
        return false;
    };

    if record_reference_id > reference_id {
        return true;
    }

    if record_reference_id < reference_id {
        return false;
    }

    record
        .alignment_start()
        .is_some_and(|start| start > alignment_start)
}

fn find_mate_closed_cut(records: &[StatsRecord]) -> Option<usize> {
    if records.len() < 2 {
        return None;
    }

    let mut required_decode_through = 0;

    for i in 0..records.len() - 1 {
        if let Some(mate_distance) = records[i].mate_distance() {
            let mate_index = i.checked_add(mate_distance)?.checked_add(1)?;

            required_decode_through = required_decode_through.max(mate_index);
        }

        // A prefix ending at i is not mate-complete yet.
        if i < required_decode_through {
            continue;
        }

        let current = &records[i];
        let next = &records[i + 1];

        let current_reference_id = current.reference_id()?;
        let next_reference_id = next.reference_id()?;

        let current_start = current.alignment_start()?;
        let next_start = next.alignment_start()?;

        // We need an actual genomic boundary between these two records so
        // `is_past_boundary` can distinguish the next record.
        if current_reference_id < next_reference_id
            || (current_reference_id == next_reference_id && current_start < next_start)
        {
            return Some(i);
        }
    }

    None
}

#[test]
/// Tests that selective stats decoding stops exactly at a mate-complete
/// genomic boundary and returns the same prefix as full-slice decoding.
///
/// full:
/// [R0 R1 R2 R3 R4 R5 ...]
///
///            boundary
///               ↓
/// selected:
/// [R0 R1 R2 R3]
///
/// and R0..R3 are already mate-complete
fn stats_records_until_matches_full_prefix() -> io::Result<()> {
    let data_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let cram_path = data_dir.join("NA12878.chr22.cram");

    let mut reader = File::open(cram_path).map(cram::io::Reader::new)?;

    let _ = reader.read_header()?;

    let mut container = Container::default();
    let mut tested_early_stop = false;

    'containers: while reader.read_container(&mut container)? != 0 {
        let compression_header = container.compression_header()?;

        for result in container.slices() {
            let slice = result?;

            let (core_data_src, external_data_srcs) = slice.decode_blocks()?;

            let full_records =
                slice.stats_records(&compression_header, &core_data_src, &external_data_srcs)?;

            let Some(cut) = find_mate_closed_cut(&full_records) else {
                continue;
            };

            let boundary_record = &full_records[cut];

            let reference_id = boundary_record
                .reference_id()
                .expect("boundary record must have a reference ID");

            let alignment_start = boundary_record
                .alignment_start()
                .expect("boundary record must have an alignment start");

            let selected_records = slice.stats_records_until_boundary(
                &compression_header,
                &core_data_src,
                &external_data_srcs,
                |record| is_past_boundary(record, reference_id, alignment_start),
            )?;

            assert_eq!(selected_records.len(), cut + 1);
            assert_eq!(selected_records, full_records[..=cut]);

            tested_early_stop = true;
            break 'containers;
        }
    }

    assert!(
        tested_early_stop,
        "fixture contains no suitable mate-closed coordinate boundary"
    );

    Ok(())
}

fn mate_closure_end(records: &[StatsRecord], cut: usize) -> Option<usize> {
    let mut required_decode_through = cut;
    let mut record_index = 0;

    while record_index <= required_decode_through {
        let record = records.get(record_index)?;

        if let Some(mate_distance) = record.mate_distance() {
            let mate_index = record_index.checked_add(mate_distance)?.checked_add(1)?;

            if mate_index >= records.len() {
                return None;
            }

            required_decode_through = required_decode_through.max(mate_index);
        }

        record_index += 1;
    }

    Some(required_decode_through)
}

fn find_cut_with_forward_mate_dependency(records: &[StatsRecord]) -> Option<(usize, usize)> {
    if records.len() < 2 {
        return None;
    }

    for cut in 0..records.len() - 1 {
        let current = &records[cut];
        let next = &records[cut + 1];

        let (
            Some(current_reference_id),
            Some(next_reference_id),
            Some(current_start),
            Some(next_start),
        ) = (
            current.reference_id(),
            next.reference_id(),
            current.alignment_start(),
            next.alignment_start(),
        )
        else {
            // This particular boundary cannot be expressed genomically.
            // Keep searching the slice.
            continue;
        };

        let has_coordinate_boundary = current_reference_id < next_reference_id
            || (current_reference_id == next_reference_id && current_start < next_start);

        if !has_coordinate_boundary {
            continue;
        }

        let Some(closure_end) = mate_closure_end(records, cut) else {
            // This candidate cannot form a valid mate-complete prefix.
            // Do not give up on the rest of the slice.
            continue;
        };

        if closure_end > cut {
            return Some((cut, closure_end));
        }
    }

    None
}

#[test]
/// Tests that bounded stats decoding can resolve downstream mates beyond
/// the genomic boundary without returning those downstream records.
///
/// The returned prefix must match full-slice decoding through the boundary,
/// including mate fields and template lengths resolved using later records.
fn stats_records_until_boundary_resolves_mates_beyond_boundary() -> io::Result<()> {
    let data_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let cram_path = data_dir.join("NA12878.chr22.cram");

    let mut reader = File::open(cram_path).map(cram::io::Reader::new)?;

    let _ = reader.read_header()?;

    let mut container = Container::default();
    let mut found_forward_mate = false;
    let mut tested_mate_extension = false;

    'containers: while reader.read_container(&mut container)? != 0 {
        let compression_header = container.compression_header()?;

        for result in container.slices() {
            let slice = result?;

            let (core_data_src, external_data_srcs) = slice.decode_blocks()?;

            let full_records =
                slice.stats_records(&compression_header, &core_data_src, &external_data_srcs)?;

            if full_records
                .iter()
                .any(|record| record.mate_distance().is_some())
            {
                found_forward_mate = true;
            }

            let Some((cut, closure_end)) = find_cut_with_forward_mate_dependency(&full_records)
            else {
                continue;
            };

            assert!(
                closure_end > cut,
                "test candidate must require decoding beyond the query boundary"
            );

            let boundary_record = &full_records[cut];

            let reference_id = boundary_record
                .reference_id()
                .expect("boundary record must have a reference ID");

            let alignment_start = boundary_record
                .alignment_start()
                .expect("boundary record must have an alignment start");

            let selected_records = slice.stats_records_until_boundary(
                &compression_header,
                &core_data_src,
                &external_data_srcs,
                |record| is_past_boundary(record, reference_id, alignment_start),
            )?;

            // The nominal genomic boundary is `cut`, but mate resolution
            // requires the decoder to continue through `closure_end`.
            // assert_eq!(selected_records.len(), closure_end + 1);
            assert_eq!(selected_records.len(), cut + 1);

            // Every returned record, including mate/TLEN fields populated by
            // resolve_stats_mates, must match ordinary full-slice decoding.
            // assert_eq!(selected_records, full_records[..=closure_end],);
            assert_eq!(selected_records, full_records[..=cut],);

            tested_mate_extension = true;
            break 'containers;
        }
    }

    assert!(
        found_forward_mate,
        "fixture must contain CRAM downstream mate relationships (NF)"
    );
    assert!(
        tested_mate_extension,
        "fixture contains no coordinate boundary with a forward mate dependency"
    );

    Ok(())
}

fn find_detached_cut(records: &[StatsRecord]) -> Option<usize> {
    if records.len() < 3 {
        return None;
    }

    // Avoid choosing a trivial boundary at either edge of the slice.
    for cut in 1..records.len() - 1 {
        let current = &records[cut];
        let next = &records[cut + 1];

        let (
            Some(current_reference_id),
            Some(next_reference_id),
            Some(current_start),
            Some(next_start),
        ) = (
            current.reference_id(),
            next.reference_id(),
            current.alignment_start(),
            next.alignment_start(),
        )
        else {
            continue;
        };

        let has_coordinate_boundary = current_reference_id < next_reference_id
            || (current_reference_id == next_reference_id && current_start < next_start);

        if has_coordinate_boundary {
            return Some(cut);
        }
    }

    None
}

#[test]
/// Tests that selectively decoding CRAM records with detached mate metadata
/// can stop exactly at the genomic boundary.
///
/// Because mate reference, position, and template information are stored
/// explicitly, no downstream mate-chain extension is required.
fn stats_records_until_stops_at_boundary_for_detached_mates() -> io::Result<()> {
    let data_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let cram_path = data_dir.join("NA18740.chr22.cram");

    let mut reader = File::open(cram_path).map(cram::io::Reader::new)?;

    let _ = reader.read_header()?;

    let mut container = Container::default();

    let mut found_paired_records = false;
    let mut found_explicit_mates = false;
    let mut tested_detached_boundary = false;

    'containers: while reader.read_container(&mut container)? != 0 {
        let compression_header = container.compression_header()?;

        for result in container.slices() {
            let slice = result?;

            let (core_data_src, external_data_srcs) = slice.decode_blocks()?;

            let full_records =
                slice.stats_records(&compression_header, &core_data_src, &external_data_srcs)?;

            if full_records.is_empty() {
                continue;
            }

            // This fixture is specifically used to exercise detached mate
            // encoding rather than CRAM downstream-mate (NF) relationships.
            if full_records
                .iter()
                .any(|record| record.mate_distance().is_some())
            {
                continue;
            }

            if full_records
                .iter()
                .any(|record| record.bam_flags().is_segmented())
            {
                found_paired_records = true;
            }

            if full_records.iter().any(|record| {
                record.mate_reference_id().is_some() && record.mate_alignment_start().is_some()
            }) {
                found_explicit_mates = true;
            }

            let Some(cut) = find_detached_cut(&full_records) else {
                continue;
            };

            let boundary_record = &full_records[cut];

            let reference_id = boundary_record
                .reference_id()
                .expect("boundary record must have a reference ID");

            let alignment_start = boundary_record
                .alignment_start()
                .expect("boundary record must have an alignment start");

            let selected_records = slice.stats_records_until_boundary(
                &compression_header,
                &core_data_src,
                &external_data_srcs,
                |record| is_past_boundary(record, reference_id, alignment_start),
            )?;

            // Detached mates already carry their mate metadata explicitly,
            // so the decoder should stop exactly at the genomic boundary.
            assert_eq!(selected_records.len(), cut + 1);
            assert_eq!(selected_records, full_records[..=cut]);

            tested_detached_boundary = true;
            break 'containers;
        }
    }

    assert!(found_paired_records, "fixture must contain paired records");

    assert!(
        found_explicit_mates,
        "fixture must contain explicit detached mate information"
    );

    assert!(
        tested_detached_boundary,
        "fixture contains no usable detached-mate coordinate boundary"
    );

    Ok(())
}

#[test]
/// Tests that selective slice reading produces the same stats records as
/// the existing full-container path.
///
/// Each slice is read once from a fully materialized container and once by
/// seeking directly to its landmark and reading only that slice's byte range.
fn selective_stats_records_match_full_container_path() -> io::Result<()> {
    let data_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let cram_path = data_dir.join("NA12878.chr22.cram");

    // Existing full-container path.
    let mut full_reader = File::open(&cram_path).map(cram::io::Reader::new)?;

    // New selective path.
    let mut selective_reader = File::open(cram_path).map(cram::io::Reader::new)?;

    let _ = full_reader.read_header()?;
    let _ = selective_reader.read_header()?;

    let mut container = Container::default();
    let mut container_header = Header::default();

    let mut container_index = 0;

    loop {
        //
        // Existing path:
        //
        // read_container()
        //     -> parses header
        //     -> materializes the entire container payload
        //
        let full_data_len = full_reader.read_container(&mut container)?;

        //
        // Selective path:
        //
        // read_container_header()
        //     -> parses only the container header
        //     -> leaves the payload unread
        //
        let selective_data_len = selective_reader.read_container_header(&mut container_header)?;

        assert_eq!(
            selective_data_len, full_data_len,
            "container payload length mismatch at container {container_index}",
        );

        if full_data_len == 0 {
            break;
        }

        // read_container_header leaves the selective reader exactly at the
        // start of the container content: compression_header.
        let payload_start = selective_reader.seek(SeekFrom::Current(0))?;

        let expected_compression_header = container.compression_header()?;

        let actual_compression_header =
            selective_reader.read_compression_header(&container_header, selective_data_len)?;

        // The compression header occupies payload bytes [0, first_landmark).
        // After reading it, the selective reader must therefore be positioned
        // at the first slice, or at the end of the payload when no slices exist.
        let expected_position = payload_start
            + u64::try_from(
                container_header
                    .landmarks
                    .first()
                    .copied()
                    .unwrap_or(selective_data_len),
            )
            .expect("container payload offset must fit in u64");

        // This test `selective::read_compression_header` moves the reader
        // to the expected position.
        assert_eq!(
            selective_reader.seek(SeekFrom::Current(0))?,
            expected_position,
            "unexpected reader position after compression header \
             at container {container_index}",
        );

        for (slice_index, &landmark) in container_header.landmarks.iter().enumerate() {
            let slice_end = container_header
                .landmarks
                .get(slice_index + 1)
                .copied()
                .unwrap_or(selective_data_len);

            let slice_len = slice_end
                .checked_sub(landmark)
                .expect("slice landmark range must be valid");

            //
            // Existing stats path:
            //
            // The entire container is already in memory, so retrieve the
            // slice from Container::src using its landmark.
            //
            let expected_slice = container.read_slice_at_landmark(
                u64::try_from(landmark).expect("slice landmark must fit in u64"),
            )?;

            let (core_data_src, external_data_srcs) = expected_slice.decode_blocks()?;

            let expected = expected_slice.stats_records(
                &expected_compression_header,
                &core_data_src,
                &external_data_srcs,
            )?;

            //
            // Selective path:
            //
            // Seek directly from the payload start to this slice and read
            // only its physical byte range.
            //
            let slice_offset = payload_start
                .checked_add(u64::try_from(landmark).expect("slice landmark must fit in u64"))
                .expect("slice offset must not overflow");

            selective_reader.seek(SeekFrom::Start(slice_offset))?;

            let actual = selective_reader.read_selective_slice_until(
                slice_len,
                &actual_compression_header,
                |_| false,
            )?;

            assert_eq!(
                actual, expected,
                "selective slice decoding mismatch at container \
                 {container_index}, slice {slice_index}",
            );
        }

        // Unlike read_container(), the selective reader is controlled by
        // explicit seeks. Put it at the next container header before the
        // next iteration.
        let next_container_offset = payload_start
            .checked_add(
                u64::try_from(selective_data_len)
                    .expect("container payload length must fit in u64"),
            )
            .expect("next container offset must not overflow");

        selective_reader.seek(SeekFrom::Start(next_container_offset))?;

        container_index += 1;
    }

    Ok(())
}
