use std::io;

use crate::Record;
use crate::record::Feature;
use noodles_core::Position;
use noodles_sam::{self as sam, alignment::Record as _};

/// A CRAM record without reference.
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
    // Used by the integration tests to compare reference-backed CRAM decoding
    // against the direct stats-record decoding path.
    #[doc(hidden)]
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

    /// Returns the BAM flags.
    pub fn bam_flags(&self) -> sam::alignment::record::Flags {
        self.bam_flags
    }

    /// Returns the reference sequence ID.
    pub fn reference_id(&self) -> Option<usize> {
        self.reference_id
    }

    /// Returns the alignment start.
    pub fn alignment_start(&self) -> Option<Position> {
        self.alignment_start
    }

    /// Return the alignment span.
    pub fn alignment_span(&self) -> Option<usize> {
        self.alignment_span
    }

    /// Return the alignment end.
    pub fn alignment_end(&self) -> io::Result<Option<Position>> {
        let Some(start) = self.alignment_start else {
            return Ok(None);
        };

        let Some(span) = self.alignment_span else {
            return Ok(Some(start));
        };

        if span == 0 {
            return Ok(Some(start));
        }

        start.checked_add(span - 1).map(Some).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "record alignment end overflow")
        })
    }

    /// Returns mapping quality.
    pub fn mapping_quality(&self) -> Option<sam::alignment::record::MappingQuality> {
        self.mapping_quality
    }

    /// Returns the template length.
    pub fn template_length(&self) -> i32 {
        self.template_length
    }

    /// Returns the mate distance.
    pub fn mate_distance(&self) -> Option<usize> {
        self.mate_distance
    }

    /// Returns the reference ID the mate mapped to.
    pub fn mate_reference_id(&self) -> Option<usize> {
        self.mate_reference_id
    }

    /// Returns the mate's alignment start.
    pub fn mate_alignment_start(&self) -> Option<Position> {
        self.mate_alignment_start
    }

    /// Returns the read length.
    pub fn read_length(&self) -> usize {
        self.read_length
    }

    /// Return the feature summary.
    pub fn feature_summary(&self) -> Option<FeatureSummary> {
        if self.features.is_empty() {
            return None;
        }
        let mut feature_summary = FeatureSummary::default();
        let read_length = self.read_length();
        for feature in &self.features {
            match feature {
                StatsFeature::Substitution { .. } => {
                    feature_summary.n_substitutions += 1;
                    feature_summary.ecnt += 1;
                    feature_summary.count_near_end_event(feature.near_ends(read_length));
                }
                StatsFeature::Insertion { len, .. } => {
                    feature_summary.n_insertions += 1;
                    feature_summary.ecnt += 1;
                    feature_summary.inserted_length += len;
                    feature_summary.count_near_end_event(feature.near_ends(read_length));
                }
                StatsFeature::Deletion { len, .. } => {
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

/// Feature summary for a CRAM record.
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
    /// Returns the number of feature events.
    pub fn event_count(&self) -> usize {
        self.ecnt
    }

    /// Returns the number of substitution events.
    pub fn substitution_count(&self) -> usize {
        self.n_substitutions
    }

    /// Returns the number of insertion events.
    pub fn insertion_count(&self) -> usize {
        self.n_insertions
    }

    /// Returns the number of deletion events.
    pub fn deletion_count(&self) -> usize {
        self.n_deletions
    }

    /// Returns the number of feature event near the end of read.
    pub fn near_end_count(&self) -> usize {
        self.n_near_ends
    }

    /// Returns the total inserted length.
    pub fn inserted_length(&self) -> usize {
        self.inserted_length
    }

    /// Returns the total deleted length.
    pub fn deleted_length(&self) -> usize {
        self.deleted_length
    }

    /// Returns the total soft-clipped length.
    pub fn soft_clipped_length(&self) -> usize {
        self.soft_clipped_length
    }

    /// Returns the total hard-clipped length.
    pub fn hard_clipped_length(&self) -> usize {
        self.hard_clipped_length
    }

    fn count_near_end_event(&mut self, near_end: Option<bool>) {
        if let Some(true) = near_end {
            self.n_near_ends += 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_stats_feature_position() {
        let position = Position::try_from(42).unwrap();

        let features = [
            StatsFeature::Substitution { position, code: 1 },
            StatsFeature::Insertion { position, len: 2 },
            StatsFeature::Deletion { position, len: 3 },
            StatsFeature::SoftClip { position, len: 4 },
            StatsFeature::HardClip { position, len: 5 },
            StatsFeature::ReferenceSkip { position, len: 6 },
            StatsFeature::Padding { position, len: 7 },
            StatsFeature::Ignore(position),
        ];

        for feature in features {
            assert_eq!(feature.position(), position);
        }
    }

    #[test]
    fn test_stats_feature_near_ends() {
        let read_length = 101;

        let near_start = StatsFeature::Substitution {
            position: Position::try_from(5).unwrap(),
            code: 1,
        };

        let middle = StatsFeature::Insertion {
            position: Position::try_from(50).unwrap(),
            len: 2,
        };

        let near_end = StatsFeature::Deletion {
            position: Position::try_from(95).unwrap(),
            len: 3,
        };

        assert_eq!(near_start.near_ends(read_length), Some(true));
        assert_eq!(middle.near_ends(read_length), Some(false));
        assert_eq!(near_end.near_ends(read_length), Some(true));
    }

    #[test]
    fn test_stats_feature_near_ends_for_short_read() {
        let feature = StatsFeature::Substitution {
            position: Position::try_from(5).unwrap(),
            code: 1,
        };

        assert_eq!(feature.near_ends(2 * ENDS_OF_READ), None);
        assert_eq!(feature.near_ends(ENDS_OF_READ), None);
    }

    #[test]
    fn test_feature_summary() {
        let record = StatsRecord {
            read_length: 101,
            features: vec![
                StatsFeature::Substitution {
                    position: Position::try_from(5).unwrap(),
                    code: 1,
                },
                StatsFeature::Substitution {
                    position: Position::try_from(50).unwrap(),
                    code: 2,
                },
                StatsFeature::Insertion {
                    position: Position::try_from(95).unwrap(),
                    len: 3,
                },
                StatsFeature::Deletion {
                    position: Position::try_from(40).unwrap(),
                    len: 4,
                },
                StatsFeature::SoftClip {
                    position: Position::try_from(1).unwrap(),
                    len: 6,
                },
                StatsFeature::HardClip {
                    position: Position::try_from(101).unwrap(),
                    len: 7,
                },
                StatsFeature::ReferenceSkip {
                    position: Position::try_from(60).unwrap(),
                    len: 20,
                },
                StatsFeature::Padding {
                    position: Position::try_from(70).unwrap(),
                    len: 2,
                },
                StatsFeature::Ignore(Position::try_from(5).unwrap()),
            ],
            ..Default::default()
        };

        let summary = record.feature_summary().unwrap();

        assert_eq!(summary.event_count(), 4);
        assert_eq!(summary.substitution_count(), 2);
        assert_eq!(summary.insertion_count(), 1);
        assert_eq!(summary.deletion_count(), 1);

        assert_eq!(summary.near_end_count(), 2);

        assert_eq!(summary.inserted_length(), 3);
        assert_eq!(summary.deleted_length(), 4);
        assert_eq!(summary.soft_clipped_length(), 6);
        assert_eq!(summary.hard_clipped_length(), 7);
    }

    #[test]
    fn test_feature_summary_is_none_without_features() {
        let record = StatsRecord::default();

        assert_eq!(record.feature_summary(), None);
    }
}
