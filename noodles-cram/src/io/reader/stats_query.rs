// TODO: Move target selection and slice planning out of the CRAM reader once
// the direct stats decoding path is validated.

use std::{
    collections::BTreeMap,
    io::{self, Read, Seek, SeekFrom},
    vec,
};

use noodles_core::{position::Position, region::Interval};

use super::{Container, Reader};
use crate::{crai, stats_record::StatsRecord};

/// An iterator over records that intersect a given region.
///
/// This is created by calling [`Reader::query`].
pub struct StatsQuery<'r, R>
where
    R: Read + Seek,
{
    reader: &'r mut Reader<R>,
    target_candidates: vec::IntoIter<TargetCandidate>,
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
        let target_candidates = get_target_candidates(index, reference_sequence_id, interval);

        Self {
            reader,
            target_candidates: target_candidates.into_iter(),
            interval,
            records: Vec::new().into_iter(),
        }
    }

    fn read_next_container(&mut self) -> Option<io::Result<()>> {
        let target = self.target_candidates.next()?;
        if let Err(e) = self.reader.seek(SeekFrom::Start(target.container_offset)) {
            return Some(Err(e));
        }

        let mut container = Container::default();

        match self.reader.read_container(&mut container) {
            Ok(0) => return None,
            Ok(_) => {}
            Err(e) => return Some(Err(e)),
        };

        let compression_header = match container.compression_header() {
            Ok(compression_header) => compression_header,
            Err(e) => return Some(Err(e)),
        };

        let records = target
            .landmarks
            .iter()
            .map(|landmark| {
                let slice = container.read_slice_at_landmark(*landmark)?;

                let (core_data_src, external_data_srcs) = slice.decode_blocks()?;
                // Our refactored stats_records(...) peer to records(...)
                // We build StatsRecord directly
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

impl<R> Iterator for StatsQuery<'_, R>
where
    R: Read + Seek,
{
    type Item = io::Result<StatsRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            match self.records.next() {
                Some(record) => {
                    if intersects(&record, self.interval) {
                        return Some(Ok(record));
                    }
                }
                None => match self.read_next_container() {
                    Some(Ok(())) => {}
                    Some(Err(e)) => return Some(Err(e)),
                    None => return None,
                },
            }
        }
    }
}

pub(crate) fn intersects(record: &StatsRecord, region_interval: Interval) -> bool {
    let (Some(start), Some(span)) = (record.alignment_start, record.alignment_span) else {
        return false;
    };
    let Some(end) = usize::from(start)
        .checked_add(span)
        .and_then(|n| n.checked_sub(1))
        .and_then(Position::new)
    else {
        return false;
    };
    let alignment_interval = (start..=end).into();
    region_interval.intersects(alignment_interval)
}

#[derive(Eq, PartialEq)]
pub(crate) struct TargetCandidate {
    pub(crate) container_offset: u64,
    pub(crate) landmarks: Vec<u64>,
}

pub(crate) fn get_target_candidates(
    index: &crai::Index,
    reference_seq_id: usize,
    interval: Interval,
) -> Vec<TargetCandidate> {
    let mut candidates: BTreeMap<u64, Vec<u64>> = BTreeMap::new();

    for record in index.iter().filter(|record| {
        record.reference_sequence_id() == Some(reference_seq_id)
            && crai_record_intersects(record, interval)
    }) {
        candidates
            .entry(record.offset())
            .or_default()
            .push(record.landmark());
    }

    candidates
        .into_iter()
        .map(|(container_offset, mut landmarks)| {
            landmarks.sort_unstable();
            landmarks.dedup();

            TargetCandidate {
                container_offset,
                landmarks,
            }
        })
        .collect()
}

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
