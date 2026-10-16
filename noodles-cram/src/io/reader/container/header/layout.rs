use std::io::{self, Read, Seek};

use crate::io::reader::num::{read_i32_le, read_itf8, read_itf8_as, read_ltf8};

const ITF8_MAX_LEN: usize = 5;
const LTF8_MAX_LEN: usize = 9;
const CRC32_LEN: usize = size_of::<u32>();

const PREFIX_ITF8_FIELD_COUNT: usize = 6;
const PREFIX_LTF8_FIELD_COUNT: usize = 2;

/// The maximum number of bytes needed from the start of a CRAM container
/// header to read through the landmark count.
pub const CONTAINER_HEADER_PREFIX_MAX_LEN: usize = size_of::<i32>()
    + PREFIX_ITF8_FIELD_COUNT * ITF8_MAX_LEN
    + PREFIX_LTF8_FIELD_COUNT * LTF8_MAX_LEN;

pub(super) fn read_byte_len_upper_bound<R>(reader: &mut R) -> io::Result<usize>
where
    R: Read + Seek,
{
    let start = reader.stream_position()?;

    // CRAM defines this as the exact byte length of the container payload,
    // excluding the encoded container header itself.
    let payload_len = read_i32_le(reader).and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    // Reference sequence context.
    let _ = read_itf8(reader)?;
    let _ = read_itf8(reader)?;
    let _ = read_itf8(reader)?;

    // Record count.
    let _ = read_itf8(reader)?;

    // Record counter and base count.
    let _ = read_ltf8(reader)?;
    let _ = read_ltf8(reader)?;

    // Block count.
    let _ = read_itf8(reader)?;

    // The landmark count determines the maximum remaining encoded length of
    // the container header.
    let landmark_count: usize = read_itf8_as(reader)?;

    let end = reader.stream_position()?;

    let prefix_len = end.checked_sub(start).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid container header reader position",
        )
    })?;

    let prefix_len =
        usize::try_from(prefix_len).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    let landmarks_max_len = landmark_count.checked_mul(ITF8_MAX_LEN).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            "container header landmark length overflow",
        )
    })?;

    let header_len_upper_bound = prefix_len
        .checked_add(landmarks_max_len)
        .and_then(|len| len.checked_add(CRC32_LEN))
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "container header length upper bound overflow",
            )
        })?;

    header_len_upper_bound
        .checked_add(payload_len)
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "container byte length upper bound overflow",
            )
        })
}

#[cfg(test)]
mod tests {
    use std::io::{self, Cursor};

    use super::*;

    #[test]
    fn test_read_byte_len_upper_bound() -> io::Result<()> {
        const START: u64 = 7;
        const PAYLOAD_LEN: i32 = 100;

        let mut src = vec![0; START as usize];

        src.extend_from_slice(&PAYLOAD_LEN.to_le_bytes());

        // All variable-width scalar fields use their one-byte encoding.
        src.extend_from_slice(&[
            0, // reference sequence ID
            0, // alignment start
            0, // alignment span
            0, // record count
            0, // record counter
            0, // base count
            0, // block count
            2, // landmark count
        ]);

        let mut reader = Cursor::new(src);
        reader.set_position(START);

        let actual = read_byte_len_upper_bound(&mut reader)?;

        // Encoded prefix:
        //   4-byte payload length
        //   + 8 one-byte variable-width fields
        // = 12 bytes.
        //
        // Remaining container-header upper bound:
        //   2 landmarks * 5 bytes
        //   + 4-byte CRC
        // = 14 bytes.
        //
        // Complete container byte-length upper bound:
        //   12 + 14 + 100 = 126 bytes.
        assert_eq!(actual, 126);

        // The inspection consumes only through the landmark count.
        assert_eq!(reader.position(), START + 12);

        assert_eq!(CONTAINER_HEADER_PREFIX_MAX_LEN, 52);

        Ok(())
    }
}
