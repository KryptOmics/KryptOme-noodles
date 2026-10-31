// Stats queries use the CRAI to prune candidate containers. Each selected
// container is read in full, and records are filtered against the requested
// interval after decoding.

use std::{
    collections::BTreeSet,
    io::{self, Read, Seek, SeekFrom},
    vec,
};

use noodles_core::{position::Position, region::Interval};

use super::{Container, Reader};
use crate::{
    crai,
    stats_record::{StatsRecord, usize_to_u64},
};

/// This is created by calling [`Reader::query_stats`].
pub struct StatsQuery<'r, R>
where
    R: Read + Seek,
{
    reader: &'r mut Reader<R>,
    container_offsets: vec::IntoIter<u64>,
    interval: Interval,
    records: vec::IntoIter<StatsRecord>,
}

impl<'r, R> StatsQuery<'r, R>
where
    R: Read + Seek,
{
    pub(super) fn new(
        reader: &'r mut Reader<R>,
        index: &crai::Index,
        reference_sequence_id: usize,
        interval: Interval,
    ) -> Self {
        let container_offsets =
            prune_duplicate_container_targets(index, reference_sequence_id, interval);

        Self {
            reader,
            container_offsets: container_offsets.into_iter(),
            interval,
            records: Vec::new().into_iter(),
        }
    }

    fn read_next_container(&mut self) -> Option<io::Result<()>> {
        let container_offset = self.container_offsets.next()?;

        if let Err(e) = self.reader.seek(SeekFrom::Start(container_offset)) {
            return Some(Err(e));
        }

        let mut container = Container::default();

        match self.reader.read_container(&mut container) {
            Ok(0) => return None,
            Ok(_) => {}
            Err(e) => return Some(Err(e)),
        }

        let compression_header = match container.compression_header() {
            Ok(compression_header) => compression_header,
            Err(e) => return Some(Err(e)),
        };

        let records = container
            .slices()
            .map(|result| {
                let slice = result?;

                let (core_data_src, external_data_srcs) = slice.decode_blocks()?;

                slice.stats_records(&compression_header, &core_data_src, &external_data_srcs)
            })
            .collect::<io::Result<Vec<_>>>();

        let records = match records {
            Ok(records) => records,
            Err(e) => return Some(Err(e)),
        };

        self.records = records
            .into_iter()
            .flatten()
            .collect::<Vec<_>>()
            .into_iter();

        Some(Ok(()))
    }
}

// Drives the stats query by loading each pruned candidate container on demand,
// decoding all of its slices, and filtering records against the requested
// interval.
impl<R> Iterator for StatsQuery<'_, R>
where
    R: Read + Seek,
{
    type Item = io::Result<StatsRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            match self.records.next() {
                Some(record) => match intersects(&record, self.interval) {
                    Ok(true) => return Some(Ok(record)),
                    Ok(false) => {}
                    Err(e) => return Some(Err(e)),
                },
                None => match self.read_next_container() {
                    Some(Ok(())) => {}
                    Some(Err(e)) => return Some(Err(e)),
                    None => return None,
                },
            }
        }
    }
}

// Applies the final record-level interval filter for the current stats query.
pub(crate) fn intersects(record: &StatsRecord, region_interval: Interval) -> io::Result<bool> {
    let Some(start) = record.alignment_start() else {
        return Ok(false);
    };

    let Some(end) = record.alignment_end()? else {
        return Ok(false);
    };

    if let Some(region_start) = region_interval.start() {
        let region_start = usize_to_u64(region_start.get());

        if end.get() < region_start {
            return Ok(false);
        }
    }

    if let Some(region_end) = region_interval.end() {
        let region_end = usize_to_u64(region_end.get());

        if start.get() > region_end {
            return Ok(false);
        }
    }

    Ok(true)
}

// Returns the unique container offsets whose CRAI records intersect the
// requested reference sequence and interval.
pub(crate) fn prune_duplicate_container_targets(
    index: &crai::Index,
    reference_sequence_id: usize,
    interval: Interval,
) -> Vec<u64> {
    index
        .iter()
        .filter(|record| {
            record.reference_sequence_id() == Some(reference_sequence_id)
                && crai_record_intersects(record, interval)
        })
        .map(crai::Record::offset)
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect()
}

// Returns whether a CRAI slice record overlaps the requested interval.
fn crai_record_intersects(record: &crai::Record, interval: Interval) -> bool {
    let Some(start) = record.alignment_start() else {
        return false;
    };

    let Some(end) = usize::from(start)
        .checked_add(record.alignment_span())
        .and_then(|n| n.checked_sub(1))
        .and_then(Position::new)
    else {
        return false;
    };

    let slice_interval = (start..=end).into();

    interval.intersects(slice_interval)
}

#[cfg(test)]
mod tests {
    use std::num::NonZeroU64;

    use noodles_core::{Position, region::Interval};

    use super::*;

    fn stats_record(start: u64, end: u64) -> StatsRecord {
        StatsRecord {
            alignment_start: NonZeroU64::new(start),
            alignment_span: Some(end - start + 1),
            ..Default::default()
        }
    }

    fn position(n: usize) -> Position {
        Position::try_from(n).unwrap()
    }

    #[test]
    fn test_intersects_bounded_interval() -> io::Result<()> {
        let record = stats_record(10, 20);

        assert!(intersects(
            &record,
            Interval::from(position(10)..=position(20)),
        )?);

        // Inclusive boundary.
        assert!(intersects(
            &record,
            Interval::from(position(20)..=position(30)),
        )?);

        assert!(!intersects(
            &record,
            Interval::from(position(21)..=position(30)),
        )?);

        assert!(!intersects(
            &record,
            Interval::from(position(1)..=position(9)),
        )?);

        Ok(())
    }

    #[test]
    fn test_intersects_unbounded_start() -> io::Result<()> {
        let record = stats_record(10, 20);

        assert!(intersects(&record, Interval::from(..=position(10)),)?);

        assert!(!intersects(&record, Interval::from(..=position(9)),)?);

        Ok(())
    }

    #[test]
    fn test_intersects_unbounded_end() -> io::Result<()> {
        let record = stats_record(10, 20);

        assert!(intersects(&record, Interval::from(position(20)..),)?);

        assert!(!intersects(&record, Interval::from(position(21)..),)?);

        Ok(())
    }

    #[test]
    fn test_intersects_unbounded_interval() -> io::Result<()> {
        let record = stats_record(10, 20);

        assert!(intersects(&record, Interval::from(..))?);

        Ok(())
    }

    #[test]
    fn test_intersects_record_without_alignment_start() -> io::Result<()> {
        let record = StatsRecord::default();

        assert!(!intersects(
            &record,
            Interval::from(position(1)..=position(10)),
        )?);

        Ok(())
    }

    #[test]
    fn test_intersects_zero_alignment_span() -> io::Result<()> {
        let record = StatsRecord {
            alignment_start: NonZeroU64::new(10),
            alignment_span: Some(0),
            ..Default::default()
        };

        assert!(intersects(
            &record,
            Interval::from(position(10)..=position(10)),
        )?);

        assert!(!intersects(
            &record,
            Interval::from(position(11)..=position(20)),
        )?);

        Ok(())
    }
}
