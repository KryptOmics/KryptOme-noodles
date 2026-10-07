use std::{io::SeekFrom, vec};

use futures::{Stream, stream};
use noodles_core::region::Interval;
use tokio::io::{self, AsyncRead, AsyncSeek};

use super::Reader;
use crate::{StatsRecord, crai, io::reader::Container};

use crate::io::reader::stats_query::{intersects, prune_duplicate_container_targets};

struct Context<'r, R> {
    reader: &'r mut Reader<R>,
    container_targets: vec::IntoIter<u64>,
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
    let container_targets =
        prune_duplicate_container_targets(index, reference_sequence_id, interval);

    let ctx = Context {
        reader,
        container_targets: container_targets.into_iter(),
        interval,
        records: Vec::new().into_iter(),
    };

    Box::pin(stream::try_unfold(ctx, |mut ctx| async {
        loop {
            match ctx.records.next() {
                Some(record) => {
                    if intersects(&record, ctx.interval)? {
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
    let container_target = ctx.container_targets.next()?;

    if let Err(e) = ctx.reader.seek(SeekFrom::Start(container_target)).await {
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

    ctx.records = records
        .into_iter()
        .flatten()
        .collect::<Vec<_>>()
        .into_iter();

    Some(Ok(()))
}
