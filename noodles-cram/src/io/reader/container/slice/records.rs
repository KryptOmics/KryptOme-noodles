mod data;
mod external_data_readers;

pub use external_data_readers::ExternalDataReaders;

use std::{borrow::Cow, error, fmt, io};

use noodles_core::Position;
use noodles_sam as sam;

use crate::{
    Record, StatsRecord,
    container::{
        CompressionHeader, ReferenceSequenceContext, block,
        compression_header::{data_series_encodings::DataSeries, preservation_map::tag_sets},
    },
    io::BitReader,
    record::{Feature, Flags, MateFlags, feature},
    stats_record::StatsFeature,
};

#[allow(clippy::enum_variant_names)]
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ReadRecordError {
    MissingDataSeriesEncoding(DataSeries),
    MissingTagEncoding(tag_sets::Key),
}

impl error::Error for ReadRecordError {}

impl fmt::Display for ReadRecordError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingDataSeriesEncoding(data_series) => {
                write!(f, "missing data series encoding: {data_series:?}")
            }
            Self::MissingTagEncoding(key) => write!(f, "missing tag encoding: {key:?}"),
        }
    }
}

pub struct Records<'c, 'ch: 'c> {
    compression_header: &'ch CompressionHeader,
    core_data_reader: BitReader<'c>,
    external_data_readers: ExternalDataReaders<'c>,
    reference_sequence_context: ReferenceSequenceContext,
    id: u64,
    prev_alignment_start: Option<Position>,
}

impl<'c, 'ch: 'c> Records<'c, 'ch> {
    pub fn new(
        compression_header: &'ch CompressionHeader,
        core_data_reader: BitReader<'c>,
        external_data_readers: ExternalDataReaders<'c>,
        reference_sequence_context: ReferenceSequenceContext,
        initial_id: u64,
    ) -> Self {
        let initial_alignment_start = match reference_sequence_context {
            ReferenceSequenceContext::Some(context) => Some(context.alignment_start()),
            _ => None,
        };

        Self {
            compression_header,
            core_data_reader,
            external_data_readers,
            reference_sequence_context,
            id: initial_id,
            prev_alignment_start: initial_alignment_start,
        }
    }

    // Consume the same CRAM data-series stream
    // Populate only the fields StatsRecord needs
    pub fn read_stats_record(&mut self, record: &mut StatsRecord) -> io::Result<()> {
        record.bam_flags = self.read_bam_flags()?;
        let cram_flags = self.read_cram_flags()?;

        self.read_stats_positions(record)?;
        self.skip_names()?;
        self.read_stats_mate(record, cram_flags)?;
        self.skip_data()?;

        if record.bam_flags.is_unmapped() {
            self.read_stats_unmapped_read(record, cram_flags)?;
        } else {
            self.read_stats_mapped_read(record, cram_flags)?;
        }

        self.id += 1;
        self.prev_alignment_start = record.alignment_start;

        Ok(())
    }

    pub fn read_record(&mut self, record: &mut Record<'c>) -> io::Result<()> {
        record.id = self.id;

        record.bam_flags = self.read_bam_flags()?;
        record.cram_flags = self.read_cram_flags()?;

        self.read_positions(record)?;
        self.read_names(record)?;
        self.read_mate(record)?;
        self.read_data(record)?;

        if record.bam_flags.is_unmapped() {
            self.read_unmapped_read(record)?;
        } else {
            self.read_mapped_read(record)?;
        }

        self.id += 1;
        self.prev_alignment_start = record.alignment_start;

        Ok(())
    }

    fn read_bam_flags(&mut self) -> io::Result<sam::alignment::record::Flags> {
        self.compression_header
            .data_series_encodings()
            .bam_flags()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::BamFlags))?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|n| {
                u16::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
            .map(sam::alignment::record::Flags::from)
    }

    fn read_cram_flags(&mut self) -> io::Result<Flags> {
        self.compression_header
            .data_series_encodings()
            .cram_flags()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::CramFlags))?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|n| {
                u8::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
            .map(Flags::from)
    }

    fn read_stats_positions(&mut self, record: &mut StatsRecord) -> io::Result<()> {
        record.reference_id = match self.reference_sequence_context {
            ReferenceSequenceContext::Some(context) => Some(context.reference_sequence_id()),
            ReferenceSequenceContext::None => None,
            ReferenceSequenceContext::Many => self.read_reference_sequence_id()?,
        };

        record.read_length = self.read_read_length()?;
        record.alignment_start = self.read_alignment_start()?;

        // StatsRecord does not use read group, but this data series still
        // belongs to the current CRAM stream and must be consumed
        let _ = self.read_read_group_id()?;

        Ok(())
    }

    fn read_positions(&mut self, record: &mut Record) -> io::Result<()> {
        record.reference_sequence_id = match self.reference_sequence_context {
            ReferenceSequenceContext::Some(context) => Some(context.reference_sequence_id()),
            ReferenceSequenceContext::None => None,
            ReferenceSequenceContext::Many => self.read_reference_sequence_id()?,
        };

        let read_length = self.read_read_length()?;
        record.read_length = read_length;

        record.alignment_start = self.read_alignment_start()?;
        record.read_group_id = self.read_read_group_id()?;

        Ok(())
    }

    fn read_reference_sequence_id(&mut self) -> io::Result<Option<usize>> {
        const UNMAPPED: i32 = -1;

        self.compression_header
            .data_series_encodings()
            .reference_sequence_ids()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::ReferenceSequenceIds))?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|n| match n {
                UNMAPPED => Ok(None),
                _ => usize::try_from(n)
                    .map(Some)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
            })
    }

    fn read_read_length(&mut self) -> io::Result<usize> {
        self.compression_header
            .data_series_encodings()
            .read_lengths()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::ReadLengths))?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|n| {
                usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
    }

    fn read_alignment_start(&mut self) -> io::Result<Option<Position>> {
        let alignment_starts_are_deltas = self
            .compression_header
            .preservation_map()
            .alignment_starts_are_deltas();

        let alignment_start_or_delta = self
            .compression_header
            .data_series_encodings()
            .alignment_starts()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::AlignmentStarts))?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)?;

        let alignment_start = if alignment_starts_are_deltas {
            let prev_alignment_start = i32::try_from(
                self.prev_alignment_start
                    .map(usize::from)
                    .unwrap_or_default(),
            )
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

            prev_alignment_start + alignment_start_or_delta
        } else {
            alignment_start_or_delta
        };

        usize::try_from(alignment_start)
            .map(Position::new)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    }

    fn read_read_group_id(&mut self) -> io::Result<Option<usize>> {
        // § 10.2 "CRAM positional data" (2021-10-15): "-1 for no group".
        const MISSING: i32 = -1;

        self.compression_header
            .data_series_encodings()
            .read_group_ids()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::ReadGroupIds))?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|n| match n {
                MISSING => Ok(None),
                _ => usize::try_from(n)
                    .map(Some)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
            })
    }

    fn skip_names(&mut self) -> io::Result<()> {
        if self
            .compression_header
            .preservation_map()
            .records_have_names()
        {
            let _ = self.read_name()?;
        }

        Ok(())
    }

    fn read_names(&mut self, record: &mut Record<'c>) -> io::Result<()> {
        let preservation_map = self.compression_header.preservation_map();

        // Missing read names are generated when resolving mates.
        if preservation_map.records_have_names() {
            record.name = self.read_name()?;
        }

        Ok(())
    }

    fn read_name(&mut self) -> io::Result<Option<Cow<'c, [u8]>>> {
        const MISSING: &[u8] = b"*\x00";

        self.compression_header
            .data_series_encodings()
            .names()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::Names))?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .map(|buf| match &buf[..] {
                MISSING => None,
                _ => Some(buf),
            })
    }

    fn read_stats_mate(&mut self, record: &mut StatsRecord, cram_flags: Flags) -> io::Result<()> {
        if cram_flags.is_detached() {
            let mate_flags = self.read_mate_flags()?;

            if mate_flags.is_on_negative_strand() {
                record.bam_flags |= sam::alignment::record::Flags::MATE_REVERSE_COMPLEMENTED;
            }

            if mate_flags.is_unmapped() {
                record.bam_flags |= sam::alignment::record::Flags::MATE_UNMAPPED;
            }

            if !self
                .compression_header
                .preservation_map()
                .records_have_names()
            {
                let _ = self.read_name()?;
            }

            record.mate_reference_id = self.read_mate_reference_sequence_id()?;
            record.mate_alignment_start = self.read_mate_alignment_start()?;
            record.template_length = self.read_template_length()?;
        } else if cram_flags.mate_is_downstream() {
            record.mate_distance = self.read_mate_distance().map(Some)?;
        }
        Ok(())
    }

    fn read_mate(&mut self, record: &mut Record<'c>) -> io::Result<()> {
        if record.cram_flags.is_detached() {
            record.mate_flags = self.read_mate_flags()?;

            if record.mate_flags.is_on_negative_strand() {
                record.bam_flags |= sam::alignment::record::Flags::MATE_REVERSE_COMPLEMENTED;
            }

            if record.mate_flags.is_unmapped() {
                record.bam_flags |= sam::alignment::record::Flags::MATE_UNMAPPED;
            }

            let preservation_map = self.compression_header.preservation_map();

            if !preservation_map.records_have_names() {
                record.name = self.read_name()?;
            }

            record.mate_reference_sequence_id = self.read_mate_reference_sequence_id()?;

            record.mate_alignment_start = self.read_mate_alignment_start()?;
            record.template_length = self.read_template_length()?;
        } else if record.cram_flags.mate_is_downstream() {
            record.mate_distance = self.read_mate_distance().map(Some)?;
        }

        Ok(())
    }

    fn read_mate_flags(&mut self) -> io::Result<MateFlags> {
        self.compression_header
            .data_series_encodings()
            .mate_flags()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::MateFlags))?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|n| {
                u8::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
            .map(MateFlags::from)
    }

    fn read_mate_reference_sequence_id(&mut self) -> io::Result<Option<usize>> {
        const UNMAPPED: i32 = -1;

        self.compression_header
            .data_series_encodings()
            .mate_reference_sequence_ids()
            .ok_or_else(|| {
                missing_data_series_encoding_error(DataSeries::MateReferenceSequenceIds)
            })?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|id| match id {
                UNMAPPED => Ok(None),
                _ => usize::try_from(id)
                    .map(Some)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
            })
    }

    fn read_mate_alignment_start(&mut self) -> io::Result<Option<Position>> {
        self.compression_header
            .data_series_encodings()
            .mate_alignment_starts()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::MateAlignmentStarts))?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|n| {
                usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
            .map(Position::new)
    }

    fn read_template_length(&mut self) -> io::Result<i32> {
        self.compression_header
            .data_series_encodings()
            .template_lengths()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::TemplateLengths))?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
    }

    fn read_mate_distance(&mut self) -> io::Result<usize> {
        self.compression_header
            .data_series_encodings()
            .mate_distances()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::MateDistances))?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|n| {
                usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
    }

    // Peer method to read_data. We do not operate on the returned data but
    // only advance the data-series cursor
    fn skip_data(&mut self) -> io::Result<()> {
        let tag_set_id = self.read_tag_set_id()?;

        let tag_set = self
            .compression_header
            .preservation_map()
            .tag_sets()
            .get(tag_set_id)
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "missing tag set"))?;

        for &key in tag_set {
            let id = block::ContentId::from(key);

            let _ = self
                .compression_header
                .tag_encodings()
                .get(&id)
                .ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        ReadRecordError::MissingTagEncoding(key),
                    )
                })?
                // decode() method advance the cursor
                // Not calling self::data::read_value(...) is intentional
                // Malformed tag payload value that would fail read_value will
                // no longer be validated by our stats path here.
                .decode(&mut self.core_data_reader, &mut self.external_data_readers)?;
        }

        Ok(())
    }

    fn read_data(&mut self, record: &mut Record<'c>) -> io::Result<()> {
        let tag_set_id = self.read_tag_set_id()?;

        let tag_set = self
            .compression_header
            .preservation_map()
            .tag_sets()
            .get(tag_set_id)
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "missing tag set"))?;

        record.data.reserve(tag_set.len());

        for &key in tag_set {
            let id = block::ContentId::from(key);

            let value = self
                .compression_header
                .tag_encodings()
                .get(&id)
                .ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        ReadRecordError::MissingTagEncoding(key),
                    )
                })?
                .decode(&mut self.core_data_reader, &mut self.external_data_readers)
                .and_then(|src| self::data::read_value(src, key.ty()))?;

            record.data.push((key.tag(), value));
        }

        Ok(())
    }

    fn read_tag_set_id(&mut self) -> io::Result<usize> {
        self.compression_header
            .data_series_encodings()
            .tag_set_ids()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::TagSetIds))?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|n| {
                usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
    }

    fn read_stats_mapped_read(
        &mut self,
        record: &mut StatsRecord,
        cram_flags: Flags,
    ) -> io::Result<()> {
        let feature_count = self.read_feature_count()?;

        let mut alignment_span = record.read_length;
        let mut prev_position = 0;

        for _ in 0..feature_count {
            let feature = self.read_stats_feature(prev_position, &mut alignment_span)?;

            prev_position = usize::from(feature.position());

            if !matches!(feature, StatsFeature::Ignore(_)) {
                record.features.push(feature);
            }
        }

        record.alignment_span = Some(alignment_span);
        record.mapping_quality = self.read_mapping_quality()?;

        if cram_flags.quality_scores_are_stored_as_array() {
            let _ = self.read_quality_scores(record.read_length)?;
        }

        Ok(())
    }

    fn read_mapped_read(&mut self, record: &mut Record<'c>) -> io::Result<()> {
        let feature_count = self.read_feature_count()?;

        let mut prev_position = 0;

        for _ in 0..feature_count {
            let feature = self.read_feature(prev_position)?;
            prev_position = usize::from(feature.position());
            record.features.push(feature);
        }

        record.mapping_quality = self.read_mapping_quality()?;

        if record.cram_flags.quality_scores_are_stored_as_array() {
            record.quality_scores = self.read_quality_scores(record.read_length)?;
        }

        Ok(())
    }

    fn read_feature_count(&mut self) -> io::Result<usize> {
        self.compression_header
            .data_series_encodings()
            .feature_counts()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::FeatureCounts))?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|n| {
                usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
    }

    // Decodes one CRAM feature and updates the reference alignment span when the
    // feature changes how many reference bases the read spans.
    fn read_stats_feature(
        &mut self,
        prev_position: usize,
        alignment_span: &mut usize,
    ) -> io::Result<StatsFeature> {
        use feature::Code;

        let code = self.read_feature_code()?;

        let delta = self.read_feature_position_delta()?;
        let position = prev_position
            .checked_add(delta)
            .and_then(Position::new)
            .ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "invalid feature position")
            })?;

        match code {
            Code::Bases => {
                let _ = self.read_stretches_of_bases()?;
                Ok(StatsFeature::Ignore(position))
            }
            Code::Scores => {
                let _ = self.read_stretches_of_quality_scores()?;
                Ok(StatsFeature::Ignore(position))
            }
            Code::ReadBase => {
                let _ = self.read_base()?;
                let _ = self.read_quality_score()?;
                Ok(StatsFeature::Ignore(position))
            }
            Code::Substitution => {
                let code = self.read_base_substitution_code()?;
                Ok(StatsFeature::Substitution { position, code })
            }
            Code::Insertion => {
                let bases = self.read_insertion_bases()?;
                let len = bases.len();
                // Update alignment span
                checked_sub_span(alignment_span, len)?;
                Ok(StatsFeature::Insertion { position, len })
            }
            Code::Deletion => {
                let len = self.read_deletion_length()?;
                // Update alignment span
                checked_add_span(alignment_span, len)?;
                Ok(StatsFeature::Deletion { position, len })
            }
            Code::InsertBase => {
                let _ = self.read_base()?;
                // Update alignment span
                checked_sub_span(alignment_span, 1)?;
                Ok(StatsFeature::Ignore(position))
            }
            Code::QualityScore => {
                let _ = self.read_quality_score()?;
                Ok(StatsFeature::Ignore(position))
            }
            Code::ReferenceSkip => {
                let len = self.read_reference_skip_length()?;
                checked_add_span(alignment_span, len)?;
                Ok(StatsFeature::ReferenceSkip { position, len })
            }
            Code::SoftClip => {
                let bases = self.read_soft_clip_bases()?;
                let len = bases.len();
                checked_sub_span(alignment_span, len)?;
                Ok(StatsFeature::SoftClip { position, len })
            }
            Code::Padding => {
                let len = self.read_padding_length()?;
                Ok(StatsFeature::Padding { position, len })
            }
            Code::HardClip => {
                let len = self.read_hard_clip_length()?;
                Ok(StatsFeature::HardClip { position, len })
            }
        }
    }

    fn read_feature(&mut self, prev_position: usize) -> io::Result<Feature<'c>> {
        use feature::Code;

        let code = self.read_feature_code()?;

        let delta = self.read_feature_position_delta()?;
        let position = Position::try_from(prev_position + delta)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        match code {
            Code::Bases => {
                let bases = self.read_stretches_of_bases()?;
                Ok(Feature::Bases { position, bases })
            }
            Code::Scores => {
                let quality_scores = self.read_stretches_of_quality_scores()?;

                Ok(Feature::Scores {
                    position,
                    quality_scores,
                })
            }
            Code::ReadBase => {
                let base = self.read_base()?;
                let quality_score = self.read_quality_score()?;

                Ok(Feature::ReadBase {
                    position,
                    base,
                    quality_score,
                })
            }
            Code::Substitution => {
                let code = self.read_base_substitution_code()?;
                Ok(Feature::Substitution { position, code })
            }
            Code::Insertion => {
                let bases = self.read_insertion_bases()?;
                Ok(Feature::Insertion { position, bases })
            }
            Code::Deletion => {
                let len = self.read_deletion_length()?;
                Ok(Feature::Deletion { position, len })
            }
            Code::InsertBase => {
                let base = self.read_base()?;
                Ok(Feature::InsertBase { position, base })
            }
            Code::QualityScore => {
                let quality_score = self.read_quality_score()?;

                Ok(Feature::QualityScore {
                    position,
                    quality_score,
                })
            }
            Code::ReferenceSkip => {
                let len = self.read_reference_skip_length()?;
                Ok(Feature::ReferenceSkip { position, len })
            }
            Code::SoftClip => {
                let bases = self.read_soft_clip_bases()?;
                Ok(Feature::SoftClip { position, bases })
            }
            Code::Padding => {
                let len = self.read_padding_length()?;
                Ok(Feature::Padding { position, len })
            }
            Code::HardClip => {
                let len = self.read_hard_clip_length()?;
                Ok(Feature::HardClip { position, len })
            }
        }
    }

    fn read_feature_code(&mut self) -> io::Result<feature::Code> {
        self.compression_header
            .data_series_encodings()
            .feature_codes()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::FeatureCodes))?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|id| {
                feature::Code::try_from(id)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
    }

    fn read_feature_position_delta(&mut self) -> io::Result<usize> {
        self.compression_header
            .data_series_encodings()
            .feature_position_deltas()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::FeaturePositionDeltas))?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|n| {
                usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
    }

    fn read_stretches_of_bases(&mut self) -> io::Result<Cow<'c, [u8]>> {
        self.compression_header
            .data_series_encodings()
            .stretches_of_bases()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::StretchesOfBases))?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
    }

    fn read_stretches_of_quality_scores(&mut self) -> io::Result<Cow<'c, [u8]>> {
        self.compression_header
            .data_series_encodings()
            .stretches_of_quality_scores()
            .ok_or_else(|| {
                missing_data_series_encoding_error(DataSeries::StretchesOfQualityScores)
            })?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
    }

    fn read_base(&mut self) -> io::Result<u8> {
        self.compression_header
            .data_series_encodings()
            .bases()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::Bases))?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
    }

    fn read_quality_score(&mut self) -> io::Result<u8> {
        self.compression_header
            .data_series_encodings()
            .quality_scores()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::QualityScores))?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
    }

    fn read_base_substitution_code(&mut self) -> io::Result<u8> {
        self.compression_header
            .data_series_encodings()
            .base_substitution_codes()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::BaseSubstitutionCodes))?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
    }

    fn read_insertion_bases(&mut self) -> io::Result<Cow<'c, [u8]>> {
        self.compression_header
            .data_series_encodings()
            .insertion_bases()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::InsertionBases))?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
    }

    fn read_deletion_length(&mut self) -> io::Result<usize> {
        self.compression_header
            .data_series_encodings()
            .deletion_lengths()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::DeletionLengths))?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|n| {
                usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
    }

    fn read_reference_skip_length(&mut self) -> io::Result<usize> {
        self.compression_header
            .data_series_encodings()
            .reference_skip_lengths()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::ReferenceSkipLengths))?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|n| {
                usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
    }

    fn read_soft_clip_bases(&mut self) -> io::Result<Cow<'c, [u8]>> {
        self.compression_header
            .data_series_encodings()
            .soft_clip_bases()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::SoftClipBases))?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
    }

    fn read_padding_length(&mut self) -> io::Result<usize> {
        self.compression_header
            .data_series_encodings()
            .padding_lengths()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::PaddingLengths))?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|n| {
                usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
    }

    fn read_hard_clip_length(&mut self) -> io::Result<usize> {
        self.compression_header
            .data_series_encodings()
            .hard_clip_lengths()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::HardClipLengths))?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|n| {
                usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
    }

    fn read_mapping_quality(
        &mut self,
    ) -> io::Result<Option<sam::alignment::record::MappingQuality>> {
        self.compression_header
            .data_series_encodings()
            .mapping_qualities()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::MappingQualities))?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|n| {
                u8::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
            .map(sam::alignment::record::MappingQuality::new)
    }

    fn read_stats_unmapped_read(
        &mut self,
        record: &StatsRecord,
        cram_flags: Flags,
    ) -> io::Result<()> {
        // Here, we still materialize the Cow<[u8]> using read_sequence
        // and read_quality_scores. Later lets see if we can skip.
        let _ = self.read_sequence(record.read_length)?;

        if cram_flags.quality_scores_are_stored_as_array() {
            let _ = self.read_quality_scores(record.read_length)?;
        }
        Ok(())
    }

    fn read_unmapped_read(&mut self, record: &mut Record<'c>) -> io::Result<()> {
        record.sequence = self.read_sequence(record.read_length)?;

        if record.cram_flags.quality_scores_are_stored_as_array() {
            record.quality_scores = self.read_quality_scores(record.read_length)?;
        }

        Ok(())
    }

    fn read_sequence(&mut self, read_length: usize) -> io::Result<Cow<'c, [u8]>> {
        let encoding = self
            .compression_header
            .data_series_encodings()
            .bases()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::Bases))?;

        encoding.get().decode_take(
            &mut self.core_data_reader,
            &mut self.external_data_readers,
            read_length,
        )
    }

    fn read_quality_scores(&mut self, read_length: usize) -> io::Result<Cow<'c, [u8]>> {
        const MISSING: u8 = 0xff;

        let encoding = self
            .compression_header
            .data_series_encodings()
            .quality_scores()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::QualityScores))?;

        let src = encoding.get().decode_take(
            &mut self.core_data_reader,
            &mut self.external_data_readers,
            read_length,
        )?;

        if src.iter().all(|&n| n == MISSING) {
            Ok(Cow::from(&[]))
        } else {
            Ok(src)
        }
    }
}

fn missing_data_series_encoding_error(data_series: DataSeries) -> io::Error {
    io::Error::new(
        io::ErrorKind::InvalidData,
        ReadRecordError::MissingDataSeriesEncoding(data_series),
    )
}

// helper for adding to alignment span
fn checked_add_span(span: &mut usize, n: usize) -> io::Result<()> {
    *span = span
        .checked_add(n)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "alignment span overflow"))?;

    Ok(())
}

// helper for subtracting to alignment span
fn checked_sub_span(span: &mut usize, n: usize) -> io::Result<()> {
    *span = span
        .checked_sub(n)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "alignment span underflow"))?;

    Ok(())
}
