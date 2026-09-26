use std::io;

use noodles_sam::alignment::record::data::field::{Tag, Type};

use crate::{
    container::compression_header::preservation_map::{TagSets, tag_sets::Key},
    io::reader::collections::read_array,
};

pub(super) fn read_tag_sets(src: &mut &[u8]) -> io::Result<TagSets> {
    let mut buf = read_array(src)?;
    read_tag_sets_inner(&mut buf)
}

fn read_tag_sets_inner(src: &mut &[u8]) -> io::Result<TagSets> {
    const NUL: u8 = 0x00;

    let mut sets = Vec::new();

    while let Some(i) = src.iter().position(|&b| b == NUL) {
        let (buf, rest) = src.split_at(i);

        *src = &rest[1..];

        let (chunks, []) = buf.as_chunks() else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid tag set value",
            ));
        };

        let line = chunks
            .iter()
            .map(|[b0, b1, b2]| {
                let tag = Tag::new(*b0, *b1);
                let ty = decode_type(*b2)?;
                Ok(Key::new(tag, ty))
            })
            .collect::<io::Result<_>>()?;

        sets.push(line);
    }

    Ok(sets)
}

fn decode_type(n: u8) -> io::Result<Type> {
    match n {
        b'A' => Ok(Type::Character),
        b'c' => Ok(Type::Int8),
        b'C' => Ok(Type::UInt8),
        b's' => Ok(Type::Int16),
        b'S' => Ok(Type::UInt16),
        b'i' => Ok(Type::Int32),
        b'I' => Ok(Type::UInt32),
        b'f' => Ok(Type::Float),
        b'Z' => Ok(Type::String),
        b'H' => Ok(Type::Hex),
        b'B' => Ok(Type::Array),
        _ => Err(io::Error::from(io::ErrorKind::InvalidData)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_decode_type() -> io::Result<()> {
        assert_eq!(decode_type(b'A')?, Type::Character);
        assert_eq!(decode_type(b'c')?, Type::Int8);
        assert_eq!(decode_type(b'C')?, Type::UInt8);
        assert_eq!(decode_type(b's')?, Type::Int16);
        assert_eq!(decode_type(b'S')?, Type::UInt16);
        assert_eq!(decode_type(b'i')?, Type::Int32);
        assert_eq!(decode_type(b'I')?, Type::UInt32);
        assert_eq!(decode_type(b'f')?, Type::Float);
        assert_eq!(decode_type(b'Z')?, Type::String);
        assert_eq!(decode_type(b'H')?, Type::Hex);
        assert_eq!(decode_type(b'B')?, Type::Array);

        assert!(matches!(
            decode_type(b'n'),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
