// FIXME: fix all of these when mature enough
#![allow(dead_code, unused, missing_docs)]

use std::io;

use crate::Record;
use crate::record::Feature;
use noodles_core::Position;
use noodles_sam::{self as sam, alignment::Record as _};

#[derive(Debug, PartialEq, Eq)]
pub struct StatsRecord {
    pub(crate) bam_flags: sam::alignment::record::Flags,
    pub(crate) reference_id: Option<usize>,
    pub(crate) alignment_start: Option<Position>,
    pub(crate) alignment_span: Option<usize>,
    pub(crate) mapping_quality: Option<sam::alignment::record::MappingQuality>,
    pub(crate) template_length: i32,
    pub(crate) mate_distance: Option<usize>,
    pub(crate) mate_reference_id: Option<usize>,
    pub(crate) mate_alignment_start: Option<Position>,
    pub(crate) read_length: usize,
    pub(crate) features: Vec<StatsFeature>,
}

impl Default for StatsRecord {
    fn default() -> Self {
        Self {
            bam_flags: sam::alignment::record::Flags::UNMAPPED,
            reference_id: None,
            alignment_start: None,
            alignment_span: None,
            mapping_quality: None,
            template_length: 0,
            mate_distance: None,
            mate_reference_id: None,
            mate_alignment_start: None,
            read_length: 0,
            features: Vec::new(),
        }
    }
}

impl StatsRecord {
    // Still reconstruct from the navtive cram::Record<'_>. No need once
    // we have read_stats_record() in which we build StatsRecord directly.
    pub fn try_from_full_record(header: &sam::Header, record: &Record<'_>) -> io::Result<Self> {
        let stats_features = record
            .features
            .iter()
            .filter_map(|feature| match feature {
                Feature::Substitution { position, code } => Some(StatsFeature::Substitution {
                    position: *position,
                    code: *code,
                }),
                Feature::Insertion { position, bases } => Some(StatsFeature::Insertion {
                    position: *position,
                    len: bases.len(),
                }),
                Feature::Deletion { position, len } => Some(StatsFeature::Deletion {
                    position: *position,
                    len: *len,
                }),
                Feature::SoftClip { position, bases } => Some(StatsFeature::SoftClip {
                    position: *position,
                    len: bases.len(),
                }),
                Feature::HardClip { position, len } => Some(StatsFeature::HardClip {
                    position: *position,
                    len: *len,
                }),
                Feature::ReferenceSkip { position, len } => Some(StatsFeature::ReferenceSkip {
                    position: *position,
                    len: *len,
                }),
                Feature::Padding { position, len } => Some(StatsFeature::Padding {
                    position: *position,
                    len: *len,
                }),
                _ => None,
            })
            .collect::<Vec<_>>();

        Ok(Self {
            bam_flags: record.flags()?,
            reference_id: record.reference_sequence_id(header).transpose()?,
            alignment_start: record.alignment_start().transpose()?,
            alignment_span: record.alignment_span().transpose()?,
            mapping_quality: record.mapping_quality().transpose()?,
            template_length: record.template_length()?,
            mate_distance: record.mate_distance,
            mate_reference_id: record.mate_reference_sequence_id(header).transpose()?,
            mate_alignment_start: record.mate_alignment_start().transpose()?,
            read_length: record.read_length,
            features: stats_features,
        })
    }

    pub fn bam_flags(&self) -> sam::alignment::record::Flags {
        self.bam_flags
    }

    pub fn reference_id(&self) -> Option<usize> {
        self.reference_id
    }

    pub fn alignment_start(&self) -> Option<Position> {
        self.alignment_start
    }

    pub fn alignment_span(&self) -> Option<usize> {
        self.alignment_span
    }

    pub(crate) fn raw_alignment_end(&self) -> Option<Position> {
        let start = self.alignment_start?;
        let span = self.alignment_span?;

        let end = span
            .checked_sub(1)
            .and_then(|n| start.get().checked_add(n))?;
        Position::new(end)
    }

    pub fn alignment_end(&self, reference_length: usize) -> Option<Position> {
        let start = self.alignment_start()?.get();
        let span = self.alignment_span()?;

        // A reference length of 0 cannot produce a valid 1-based Position.
        // Position is NonZero.
        Position::new(reference_length)?;

        if start > reference_length {
            return None;
        }

        let end = self.raw_alignment_end()?;

        Position::new(end.get().min(reference_length))
    }

    pub fn mapping_quality(&self) -> Option<sam::alignment::record::MappingQuality> {
        self.mapping_quality
    }

    pub fn template_length(&self) -> i32 {
        self.template_length
    }

    pub fn mate_distance(&self) -> Option<usize> {
        self.mate_distance
    }

    pub fn mate_reference_id(&self) -> Option<usize> {
        self.mate_reference_id
    }

    pub fn mate_alignment_start(&self) -> Option<Position> {
        self.mate_alignment_start
    }

    pub fn read_length(&self) -> usize {
        self.read_length
    }

    pub fn feature_summary(&self) -> Option<FeatureSummary> {
        if self.features.is_empty() {
            return None;
        }
        let mut feature_summary = FeatureSummary::default();
        let read_length = self.read_length();
        for feature in &self.features {
            match feature {
                StatsFeature::Substitution { position, .. } => {
                    feature_summary.n_substitutions += 1;
                    feature_summary.ecnt += 1;
                    feature_summary.count_near_end_event(feature.near_ends(read_length));
                }
                StatsFeature::Insertion { position, len } => {
                    feature_summary.n_insertions += 1;
                    feature_summary.ecnt += 1;
                    feature_summary.inserted_length += len;
                    feature_summary.count_near_end_event(feature.near_ends(read_length));
                }
                StatsFeature::Deletion { position, len } => {
                    feature_summary.n_deletions += 1;
                    feature_summary.ecnt += 1;
                    feature_summary.deleted_length += len;
                    feature_summary.count_near_end_event(feature.near_ends(read_length));
                }
                StatsFeature::SoftClip { len, .. } => {
                    feature_summary.soft_clipped_length += len;
                }
                StatsFeature::HardClip { len, .. } => {
                    feature_summary.hard_clipped_length += len;
                }
                _ => {}
            }
        }

        Some(feature_summary)
    }
}

/// A compact, reference-free summary of a CRAM read feature.
///
/// The `position` is the CRAM feature position, not an absolute reference
/// coordinate. To convert this into a reference coordinate, the caller must
/// walk the feature stream relative to the alignment start.
#[derive(Debug, PartialEq, Eq)]
// cram::record::feature:Feature
pub(crate) enum StatsFeature {
    Substitution { position: Position, code: u8 },

    Insertion { position: Position, len: usize },

    Deletion { position: Position, len: usize },

    SoftClip { position: Position, len: usize },

    HardClip { position: Position, len: usize },

    ReferenceSkip { position: Position, len: usize },

    Padding { position: Position, len: usize },

    Ignore(Position),
}

const ENDS_OF_READ: usize = 10;
impl StatsFeature {
    pub(crate) fn position(&self) -> Position {
        match self {
            Self::Substitution { position, .. }
            | Self::Insertion { position, .. }
            | Self::Deletion { position, .. }
            | Self::SoftClip { position, .. }
            | Self::HardClip { position, .. }
            | Self::ReferenceSkip { position, .. }
            | Self::Padding { position, .. }
            | Self::Ignore(position) => *position,
        }
    }

    fn near_ends(&self, read_length: usize) -> Option<bool> {
        // If the read is too short, "near either end" stops being a meaningful
        // metric because the two end windows overlap or consume the whole read.
        // To compute fraction of near-read-end events, we should not use the
        // total read count as the denominator
        if read_length <= 2 * ENDS_OF_READ {
            return None;
        }

        // Note position is 1-based here
        let near_or_not = |pos: &Position| {
            let pos = pos.get();
            pos <= ENDS_OF_READ || pos > read_length - ENDS_OF_READ
        };

        match self {
            Self::Substitution { position, .. }
            | Self::Insertion { position, .. }
            | Self::Deletion { position, .. } => Some(near_or_not(position)),
            _ => Some(false),
        }
    }
}

#[derive(Debug, Default, PartialEq, Eq)]
pub struct FeatureSummary {
    ecnt: usize,
    n_substitutions: usize,
    n_insertions: usize,
    n_deletions: usize,
    n_near_ends: usize,
    inserted_length: usize,
    deleted_length: usize,
    soft_clipped_length: usize,
    hard_clipped_length: usize,
}

impl FeatureSummary {
    pub fn event_count(&self) -> usize {
        self.ecnt
    }

    pub fn substitution_count(&self) -> usize {
        self.n_substitutions
    }

    pub fn insertion_count(&self) -> usize {
        self.n_insertions
    }

    pub fn deletion_count(&self) -> usize {
        self.n_deletions
    }

    pub fn near_end_count(&self) -> usize {
        self.n_near_ends
    }

    pub fn inserted_length(&self) -> usize {
        self.inserted_length
    }

    pub fn deleted_length(&self) -> usize {
        self.deleted_length
    }

    pub fn soft_clipped_len(&self) -> usize {
        self.soft_clipped_length
    }

    pub fn hard_clipped_length(&self) -> usize {
        self.hard_clipped_length
    }

    fn count_near_end_event(&mut self, near_end: Option<bool>) {
        if let Some(true) = near_end {
            self.n_near_ends += 1;
        }
    }
}
