#![allow(missing_docs)]
use std::{
    fs::File,
    io::{self, Cursor, Read},
    path::{Path, PathBuf},
};

use noodles_cram as cram;
use noodles_sam as sam;

fn read_file_header_with_bounded_container(path: &Path) -> io::Result<sam::Header> {
    let mut file = File::open(path)?;

    // Advance to the start of the CRAM header container using the native
    // file-definition reader rather than duplicating its encoded length here.
    {
        let mut reader = cram::io::Reader::new(&mut file);
        reader.read_file_definition()?;
    }

    let prefix_len = cram::io::reader::CONTAINER_HEADER_PREFIX_MAX_LEN;

    let mut container_bytes = vec![0; prefix_len];
    file.read_exact(&mut container_bytes)?;

    let byte_len_upper_bound = {
        let mut reader = cram::io::Reader::new(Cursor::new(&container_bytes));

        reader.read_container_byte_len_upper_bound()?
    };

    if byte_len_upper_bound > container_bytes.len() {
        let already_read = container_bytes.len();

        container_bytes.resize(byte_len_upper_bound, 0);

        file.read_exact(&mut container_bytes[already_read..])?;
    } else {
        container_bytes.truncate(byte_len_upper_bound);
    }

    // The buffer begins exactly at the CRAM header container. The native
    // reader performs the authoritative container-header validation and
    // FILE_HEADER/SAM-header decoding.
    let mut reader = cram::io::Reader::new(Cursor::new(container_bytes));

    reader.read_file_header()
}

fn read_file_header_normally(path: &Path) -> io::Result<sam::Header> {
    let mut reader = cram::io::Reader::new(File::open(path)?);

    reader.read_file_definition()?;
    reader.read_file_header()
}

#[test]
fn read_file_header_from_container_byte_len_upper_bound() -> io::Result<()> {
    let data_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data");

    for name in ["NA12878.chr22.cram", "NA18740.chr22.cram"] {
        let path = data_dir.join(name);

        let expected = read_file_header_normally(&path)?;
        let actual = read_file_header_with_bounded_container(&path)?;

        assert_eq!(actual, expected, "fixture: {name}");
    }

    Ok(())
}
