//! Queries a CRAM file for a given region without reference.
//!
use std::{
    env,
    error::Error,
    ffi::OsString,
    fs::File,
    io::{BufWriter, Write},
    path::{Path, PathBuf},
};

use noodles_core::Region;
use noodles_cram::{self as cram, crai};

fn main() -> Result<(), Box<dyn Error>> {
    let mut args = env::args_os();
    let program = args.next().unwrap_or_else(|| OsString::from("query_stats"));

    let cram_path = args.next().map(PathBuf::from).ok_or_else(|| {
        format!(
            "usage: {} <input.cram> <output.tsv> <region>",
            Path::new(&program).display()
        )
    })?;

    let output_path = args.next().map(PathBuf::from).ok_or_else(|| {
        format!(
            "usage: {} <input.cram> <output.tsv> <region>",
            Path::new(&program).display()
        )
    })?;

    let region: Region = args
        .next()
        .ok_or_else(|| {
            format!(
                "usage: {} <input.cram> <output.tsv> <region>",
                Path::new(&program).display()
            )
        })?
        .into_string()
        .map_err(|_| "region must be valid UTF-8")?
        .parse()?;

    if args.next().is_some() {
        return Err(format!(
            "usage: {} <input.cram> <output.tsv> <region>",
            Path::new(&program).display()
        )
        .into());
    }

    let crai_path = crai_path(&cram_path);

    let mut reader = File::open(&cram_path).map(cram::io::Reader::new)?;
    let header = reader.read_header()?;

    let reference_lengths = header
        .reference_sequences()
        .iter()
        .map(|(_, reference_sequence)| reference_sequence.length())
        .collect::<Vec<_>>();

    let index = crai::fs::read(crai_path)?;

    let query = reader.query_stats(&header, &index, &region)?;

    let output = File::create(output_path)?;
    let mut writer = BufWriter::new(output);

    writeln!(
        writer,
        "reference_id\tread_length\tflags\talignment_start\talignment_end\tmapq\ttemplate_length\tfeature_summary"
    )?;

    for result in query {
        let record = result?;

        let Some(reference_id) = record.reference_id() else {
            continue;
        };

        let Some(reference_length) = reference_lengths.get(reference_id) else {
            continue;
        };

        let alignment_start = record
            .alignment_start()
            .map(|position| position.get().to_string())
            .unwrap_or_else(|| ".".into());

        let alignment_end = record
            .alignment_end(reference_length.get())
            .map(|position| position.get().to_string())
            .unwrap_or_else(|| ".".into());

        let mapping_quality = record
            .mapping_quality()
            .map(|mapping_quality| u8::from(mapping_quality).to_string())
            .unwrap_or_else(|| ".".into());

        writeln!(
            writer,
            "{}\t{}\t{:?}\t{}\t{}\t{}\t{}\t{:?}",
            reference_id,
            record.read_length(),
            record.bam_flags(),
            alignment_start,
            alignment_end,
            mapping_quality,
            record.template_length(),
            record.feature_summary(),
        )?;
    }

    writer.flush()?;

    Ok(())
}

fn crai_path(cram_path: &Path) -> PathBuf {
    let mut path = cram_path.as_os_str().to_owned();
    path.push(".crai");
    PathBuf::from(path)
}
