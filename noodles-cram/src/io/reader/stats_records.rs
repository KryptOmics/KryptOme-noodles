use std::{
    io::{self, Read},
    vec,
};

use super::{Container, Reader};
use crate::stats_record::StatsRecord;

/// An iterator over stats records of a CRAM reader.
///
/// This is created by calling [`Reader::stats_records`].
pub struct StatsRecords<'r, R>
where
    R: Read,
{
    reader: &'r mut Reader<R>,
    container: Container,
    records: vec::IntoIter<StatsRecord>,
}

impl<'r, R> StatsRecords<'r, R>
where
    R: Read,
{
    pub(crate) fn new(reader: &'r mut Reader<R>) -> Self {
        Self {
            reader,
            container: Container::default(),
            records: Vec::new().into_iter(),
        }
    }

    fn read_container_records(&mut self) -> io::Result<bool> {
        if self.reader.read_container(&mut self.container)? == 0 {
            return Ok(true);
        }

        let compression_header = self.container.compression_header()?;

        self.records = self
            .container
            .slices()
            .map(|result| {
                let slice = result?;

                let (core_data_src, external_data_srcs) = slice.decode_blocks()?;

                slice.stats_records(&compression_header, &core_data_src, &external_data_srcs)
            })
            .collect::<io::Result<Vec<_>>>()?
            .into_iter()
            .flatten()
            .collect::<Vec<_>>()
            .into_iter();

        Ok(false)
    }
}

impl<R> Iterator for StatsRecords<'_, R>
where
    R: Read,
{
    type Item = io::Result<StatsRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            match self.records.next() {
                Some(record) => return Some(Ok(record)),
                None => match self.read_container_records() {
                    Ok(true) => return None,
                    Ok(false) => {}
                    Err(e) => return Some(Err(e)),
                },
            }
        }
    }
}
