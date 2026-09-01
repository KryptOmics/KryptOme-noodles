use std::{
    env,
    error::Error,
    ffi::OsString,
    future::Future,
    io::{self, SeekFrom},
    path::Path,
    pin::Pin,
    sync::Arc,
    task::{Context, Poll, ready},
};

use futures::TryStreamExt;
use noodles_core::Region;
use noodles_cram::{self as cram, crai};
use object_store::{ObjectStore, ObjectStoreExt, path::Path as ObjectPath};
use tokio::{
    fs::File,
    io::{AsyncRead, AsyncSeek, AsyncWriteExt, BufWriter, ReadBuf},
};
use url::Url;

// The adapter is only here to demonstrate that the async CRAM stats query can
// operate against a remote object store. A production remote reader would
// likely want a more deliberate caching/range-planning strategy.
//
// In particular, querying CRAM containers causes seeks to CRAI-provided byte
// offsets. For a remote object, those "seeks" are translated into ranged
// object-store requests rather than operating on a local file descriptor.
const READ_CHUNK_SIZE: u64 = 1024 * 1024;

type FetchFuture = Pin<Box<dyn Future<Output = io::Result<(u64, Vec<u8>)>> + Send>>;

enum ReadState {
    Idle,
    Fetching(FetchFuture),
}

/// A minimal buffered `AsyncRead + AsyncSeek` adapter over an object store.
///
/// `AsyncSeek` only changes the logical byte position. The next `AsyncRead`
/// fetches a chunk starting at that position using `ObjectStore::get_range`.
///
/// The read-ahead buffer is important: without it, small CRAM parser reads
/// could cause a separate remote request for every few bytes.
struct ObjectStoreReader {
    store: Arc<dyn ObjectStore>,
    path: ObjectPath,
    len: u64,
    pos: u64,

    buffer: Vec<u8>,
    buffer_start: u64,

    state: ReadState,
    pending_seek: Option<u64>,
}

impl ObjectStoreReader {
    fn new(store: Arc<dyn ObjectStore>, path: ObjectPath, len: u64) -> Self {
        Self {
            store,
            path,
            len,
            pos: 0,
            buffer: Vec::new(),
            buffer_start: 0,
            state: ReadState::Idle,
            pending_seek: None,
        }
    }

    fn clear_buffer(&mut self) {
        self.buffer.clear();
        self.buffer_start = 0;

        // Dropping an in-flight future effectively cancels the pending range
        // read when a seek moves us elsewhere.
        self.state = ReadState::Idle;
    }
}

impl AsyncRead for ObjectStoreReader {
    fn poll_read(
        mut self: Pin<&mut Self>,
        cx: &mut Context<'_>,
        dst: &mut ReadBuf<'_>,
    ) -> Poll<io::Result<()>> {
        loop {
            let this = self.as_mut().get_mut();

            let buffer_end = this.buffer_start + this.buffer.len() as u64;

            // If the current logical position is inside the cached range,
            // serve directly from it.
            if this.pos >= this.buffer_start && this.pos < buffer_end {
                let offset = (this.pos - this.buffer_start) as usize;
                let src = &this.buffer[offset..];

                let n = src.len().min(dst.remaining());

                dst.put_slice(&src[..n]);
                this.pos += n as u64;

                return Poll::Ready(Ok(()));
            }

            if this.pos >= this.len {
                return Poll::Ready(Ok(()));
            }

            if matches!(this.state, ReadState::Idle) {
                let start = this.pos;
                let end = start.saturating_add(READ_CHUNK_SIZE).min(this.len);

                eprintln!("fetch {start}..{end}");

                let store = Arc::clone(&this.store);
                let path = this.path.clone();

                let future = async move {
                    store
                        .get_range(&path, start..end)
                        .await
                        .map(|bytes| (start, bytes.to_vec()))
                        .map_err(io::Error::other)
                };

                this.state = ReadState::Fetching(Box::pin(future));
            }

            let result = match &mut this.state {
                ReadState::Fetching(future) => ready!(future.as_mut().poll(cx)),
                ReadState::Idle => unreachable!(),
            };

            this.state = ReadState::Idle;

            let (start, data) = result?;

            if data.is_empty() {
                return Poll::Ready(Ok(()));
            }

            this.buffer_start = start;
            this.buffer = data;
        }
    }
}
impl AsyncSeek for ObjectStoreReader {
    fn start_seek(mut self: Pin<&mut Self>, position: SeekFrom) -> io::Result<()> {
        if self.pending_seek.is_some() {
            return Err(io::Error::new(
                io::ErrorKind::Other,
                "seek already in progress",
            ));
        }

        let pos = match position {
            SeekFrom::Start(pos) => i128::from(pos),

            SeekFrom::Current(offset) => i128::from(self.pos) + i128::from(offset),

            SeekFrom::End(offset) => i128::from(self.len) + i128::from(offset),
        };

        if pos < 0 || pos > i128::from(u64::MAX) {
            return Err(io::Error::new(io::ErrorKind::InvalidInput, "invalid seek"));
        }

        self.pending_seek = Some(pos as u64);

        Ok(())
    }

    fn poll_complete(mut self: Pin<&mut Self>, _cx: &mut Context<'_>) -> Poll<io::Result<u64>> {
        if let Some(pos) = self.pending_seek.take() {
            self.pos = pos;

            // Any in-flight request was based on the old cursor.
            self.state = ReadState::Idle;
        }

        Poll::Ready(Ok(self.pos))
    }
}

#[tokio::main]
async fn main() -> Result<(), Box<dyn Error>> {
    let mut args = env::args_os();

    let program = args
        .next()
        .unwrap_or_else(|| OsString::from("cram_async_query_stats"));

    let input = args.next().ok_or_else(|| usage(&program))?;

    let output_path = args.next().ok_or_else(|| usage(&program))?;

    let region: Region = args
        .next()
        .ok_or_else(|| usage(&program))?
        .into_string()
        .map_err(|_| "region must be valid UTF-8")?
        .parse()?;

    if args.next().is_some() {
        return Err(usage(&program).into());
    }

    let input = input
        .into_string()
        .map_err(|_| "input URL must be valid UTF-8")?;

    let url = Url::parse(&input)?;

    // The 1000 Genomes bucket is public, so requests do not need AWS
    // credentials or request signing.
    let options = [
        ("aws_region", "us-east-1"),
        ("aws_skip_signature", "true"),
        ("timeout", "3600s"),
    ];

    let (store, path) = object_store::parse_url_opts(&url, options)?;

    let metadata = store.head(&path).await?;

    eprintln!("remote object: {} ({} bytes)", input, metadata.size);

    // Wrap the object store in the minimal seekable reader expected by the
    // current async noodles-cram query implementation.
    let remote_reader = ObjectStoreReader::new(store.into(), path, metadata.size);

    let mut reader = cram::r#async::io::Reader::new(remote_reader);

    // Header decoding is sequential and does not require the reference.
    let header = reader.read_header().await?;

    let reference_lengths = header
        .reference_sequences()
        .iter()
        .map(|(_, reference_sequence)| reference_sequence.length())
        .collect::<Vec<_>>();

    // Infer the index URL in the same way as the local example:
    //
    //     sample.cram
    //     sample.cram.crai
    //
    // CRAI files are small, so unlike the CRAM itself, there is no need for
    // random access. Download the entire index and parse it in memory.
    let crai_url = if let Some(prefix) = input.strip_suffix(".cram") {
        Url::parse(&format!("{prefix}.crai"))?
    } else {
        Url::parse(&format!("{input}.crai"))?
    };

    let (crai_store, crai_path) = object_store::parse_url_opts(&crai_url, options)?;

    let crai_bytes = crai_store.get(&crai_path).await?.bytes().await?;

    let index = {
        let mut reader = crai::io::Reader::new(&crai_bytes[..]);
        reader.read_index()?
    };

    let query = reader.query_stats(&header, &index, &region)?;

    let output = File::create(output_path).await?;
    let mut writer = BufWriter::new(output);

    writer
        .write_all(
            b"reference_id\tread_length\tflags\talignment_start\talignment_end\tmapq\ttemplate_length\tfeature_summary\n",
        )
        .await?;

    tokio::pin!(query);

    while let Some(record) = query.try_next().await? {
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

        let line = format!(
            "{}\t{}\t{:?}\t{}\t{}\t{}\t{}\t{:?}\n",
            reference_id,
            record.read_length(),
            record.bam_flags(),
            alignment_start,
            alignment_end,
            mapping_quality,
            record.template_length(),
            record.feature_summary(),
        );

        writer.write_all(line.as_bytes()).await?;
    }

    writer.flush().await?;

    Ok(())
}

fn usage(program: &OsString) -> String {
    format!(
        "usage: {} <input.cram-url> <output.tsv> <region>",
        Path::new(program).display()
    )
}
