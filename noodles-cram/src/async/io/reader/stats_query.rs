use std::{io::SeekFrom, vec};

use futures::{Stream, stream};
use noodles_core::region::Interval;
use tokio::io::{self, AsyncRead, AsyncSeek};

use super::Reader;
use crate::{StatsRecord, crai, io::reader::Container};

use crate::io::reader::stats_query::{TargetCandidate, get_target_candidates, intersects};

struct Context<'r, R> {
    reader: &'r mut Reader<R>,
    target_candidates: vec::IntoIter<TargetCandidate>,
    interval: Interval,
    records: vec::IntoIter<StatsRecord>,
}

pub(super) fn query_stats<'r, R>(
    reader: &'r mut Reader<R>,
    index: &crai::Index,
    reference_sequence_id: usize,
    interval: Interval,
) -> impl Stream<Item = io::Result<StatsRecord>> + 'r
where
    R: AsyncRead + AsyncSeek + Unpin,
{
    let target_candidates = get_target_candidates(index, reference_sequence_id, interval);

    let ctx = Context {
        reader,
        target_candidates: target_candidates.into_iter(),
        interval,
        records: Vec::new().into_iter(),
    };

    Box::pin(stream::try_unfold(ctx, |mut ctx| async {
        loop {
            match ctx.records.next() {
                Some(record) => {
                    if intersects(&record, ctx.interval) {
                        return Ok(Some((record, ctx)));
                    }
                }
                None => match read_next_container(&mut ctx).await {
                    Some(Ok(())) => {}
                    Some(Err(e)) => return Err(e),
                    None => return Ok(None),
                },
            }
        }
    }))
}

async fn read_next_container<R>(ctx: &mut Context<'_, R>) -> Option<io::Result<()>>
where
    R: AsyncRead + AsyncSeek + Unpin,
{
    let target = ctx.target_candidates.next()?;

    if let Err(e) = ctx
        .reader
        .seek(SeekFrom::Start(target.container_offset))
        .await
    {
        return Some(Err(e));
    }

    let mut container = Container::default();

    match ctx.reader.read_container(&mut container).await {
        Ok(0) => return None,
        Ok(_) => {}
        Err(e) => return Some(Err(e)),
    }

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

            // Shared sync/async decoding path.
            slice.stats_records(&compression_header, &core_data_src, &external_data_srcs)
        })
        .collect::<io::Result<Vec<_>>>();

    let records = match records {
        Ok(records) => records,
        Err(e) => return Some(Err(e)),
    };

    ctx.records = records
        .into_iter()
        .flatten()
        .collect::<Vec<_>>()
        .into_iter();

    Some(Ok(()))
}
