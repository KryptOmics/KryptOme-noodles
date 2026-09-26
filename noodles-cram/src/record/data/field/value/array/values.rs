use std::{borrow::Cow, io, marker::PhantomData};

use noodles_sam as sam;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Values<'c, N> {
    src: Cow<'c, [u8]>,
    len: usize,
    _marker: PhantomData<N>,
}

impl<'c, N> Values<'c, N> {
    pub(crate) fn new(src: Cow<'c, [u8]>, len: usize) -> Self {
        Self {
            src,
            len,
            _marker: PhantomData,
        }
    }
}

const OFFSET: usize = 5;

impl<'c> sam::alignment::record::data::field::value::array::Values<'c, i8> for &'c Values<'c, i8> {
    fn len(&self) -> usize {
        self.len
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<i8>> + '_> {
        Box::new(self.src[OFFSET..].iter().copied().map(|n| n as i8).map(Ok))
    }
}

impl<'c> sam::alignment::record::data::field::value::array::Values<'c, u8> for &'c Values<'c, u8> {
    fn len(&self) -> usize {
        self.len
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<u8>> + '_> {
        Box::new(self.src[OFFSET..].iter().copied().map(Ok))
    }
}

impl<'c> sam::alignment::record::data::field::value::array::Values<'c, i16>
    for &'c Values<'c, i16>
{
    fn len(&self) -> usize {
        self.len
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<i16>> + '_> {
        let (chunks, []) = self.src[OFFSET..].as_chunks() else {
            panic!();
        };

        Box::new(
            chunks
                .iter()
                .map(|chunk| i16::from_le_bytes(*chunk))
                .map(Ok),
        )
    }
}

impl<'c> sam::alignment::record::data::field::value::array::Values<'c, u16>
    for &'c Values<'c, u16>
{
    fn len(&self) -> usize {
        self.len
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<u16>> + '_> {
        let (chunks, []) = self.src[OFFSET..].as_chunks() else {
            panic!();
        };

        Box::new(
            chunks
                .iter()
                .map(|chunk| u16::from_le_bytes(*chunk))
                .map(Ok),
        )
    }
}

impl<'c> sam::alignment::record::data::field::value::array::Values<'c, i32>
    for &'c Values<'c, i32>
{
    fn len(&self) -> usize {
        self.len
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<i32>> + '_> {
        let (chunks, []) = self.src[OFFSET..].as_chunks() else {
            panic!();
        };

        Box::new(
            chunks
                .iter()
                .map(|chunk| i32::from_le_bytes(*chunk))
                .map(Ok),
        )
    }
}

impl<'c> sam::alignment::record::data::field::value::array::Values<'c, u32>
    for &'c Values<'c, u32>
{
    fn len(&self) -> usize {
        self.len
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<u32>> + '_> {
        let (chunks, []) = self.src[OFFSET..].as_chunks() else {
            panic!();
        };

        Box::new(
            chunks
                .iter()
                .map(|chunk| u32::from_le_bytes(*chunk))
                .map(Ok),
        )
    }
}

impl<'c> sam::alignment::record::data::field::value::array::Values<'c, f32>
    for &'c Values<'c, f32>
{
    fn len(&self) -> usize {
        self.len
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<f32>> + '_> {
        let (chunks, []) = self.src[OFFSET..].as_chunks() else {
            panic!();
        };

        Box::new(
            chunks
                .iter()
                .map(|chunk| f32::from_le_bytes(*chunk))
                .map(Ok),
        )
    }
}
