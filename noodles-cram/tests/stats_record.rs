#![allow(missing_docs)]
use std::{fs::File, io, path::PathBuf};

use noodles_cram::{self as cram, StatsRecord, io::reader::Container};
use noodles_fasta::{self as fasta, repository::adapters::IndexedReader};

#[test]
fn stats_records_match_full_records() -> Result<(), Box<dyn std::error::Error>> {
    let data_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data");

    let cram_path = data_dir.join("NA18740.chr22.cram");
    let reference_path = data_dir.join("chr22.fa.gz");

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
