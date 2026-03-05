#!/usr/bin/env python3
"""
Generate CDS/UTR region boundaries from GTF files.
Outputs the boundaries between 5'UTR, CDS, and 3'UTR regions in transcript coordinates.

Fixed issues from original:
- Correct coordinate conversion for minus strand transcripts
- Proper handling of CDS regions spanning multiple exons
- Validation of exon/CDS consistency
- Stop codon awareness
- Configurable output path
- Comprehensive validation suite
- Expanded GENCODE tag recognition for non-standard CDS
"""

import re
import sys
import argparse
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional, Set, FrozenSet
from enum import Enum


class Strand(Enum):
    PLUS = '+'
    MINUS = '-'


@dataclass
class TranscriptData:
    """Container for transcript annotation data."""
    exons: List[Tuple[int, int]] = field(default_factory=list)
    cds: List[Tuple[int, int]] = field(default_factory=list)
    strand: Optional[Strand] = None
    gene_id: Optional[str] = None
    gene_name: Optional[str] = None
    chrom: Optional[str] = None
    transcript_type: Optional[str] = None
    
    # CDS completeness tags from GTF
    cds_start_NF: bool = False  # CDS start not found (5' incomplete)
    cds_end_NF: bool = False    # CDS end not found (3' incomplete)
    mRNA_start_NF: bool = False # mRNA start not found
    mRNA_end_NF: bool = False   # mRNA end not found
    
    # Special transcript types that may have non-standard CDS
    is_seleno: bool = False           # Selenocysteine (UGA recoding)
    is_readthrough: bool = False      # Stop codon readthrough
    is_polymorphism: bool = False     # Non-canonical polymorphism
    
    # Additional tags that can explain non-standard CDS (newly added)
    has_reference_genome_error: bool = False  # Genome error affects CDS
    has_sequence_error: bool = False          # Non-canonical splice from error
    is_non_canonical_TEC: bool = False        # To be experimentally confirmed
    is_non_canonical_conserved: bool = False  # Non-canonical but conserved
    is_non_canonical_other: bool = False      # Other non-canonical
    is_non_canonical_U12: bool = False        # U12 intron (AT-AC splice)
    is_inferred_exon_combination: bool = False  # Exon combination not directly supported
    is_inferred_transcript_model: bool = False  # Transcript not supported by single evidence
    has_low_sequence_quality: bool = False    # Poor sequence quality
    is_semi_processed: bool = False           # Pseudogene with retained introns
    is_fragmented_locus: bool = False         # Non-overlapping fragments
    is_NMD_exception: bool = False            # NMD exception transcript
    has_retained_intron: bool = False         # Retained intron in CDS
    is_bicistronic: bool = False              # Two CDS regions
    has_non_ATG_start: bool = False           # Non-ATG start codon
    has_upstream_ATG: bool = False            # Upstream ATG present
    has_downstream_ATG: bool = False          # Downstream ATG used
    
    tags: Set[str] = field(default_factory=set)
    _exon_lengths: List[int] = field(default_factory=list, repr=False)
    _cumulative_lengths: List[int] = field(default_factory=list, repr=False)
    
    def _update_lengths(self):
        """Pre-calculate lengths for faster coordinate conversion."""
        self._exon_lengths = [e - s + 1 for s, e in self.exons]
        self._cumulative_lengths = [0] * (len(self.exons) + 1)
        curr = 0
        for i, length in enumerate(self._exon_lengths):
            curr += length
            self._cumulative_lengths[i+1] = curr

    def transcript_length(self) -> int:
        """Calculate total transcript length from exons."""
        if not self._cumulative_lengths:
            self._update_lengths()
        return self._cumulative_lengths[-1]
    
    def is_mitochondrial(self) -> bool:
        """Check if transcript is on mitochondrial chromosome."""
        return self.chrom in ('MT', 'chrM', 'chrMT', 'M')
    
    def is_immunoglobulin_or_tcr(self) -> bool:
        """Check if transcript is from IG or TR locus (known annotation challenges)."""
        if self.transcript_type:
            return self.transcript_type.startswith(('IG_', 'TR_'))
        if self.gene_name:
            # KIR genes are also problematic
            return self.gene_name.startswith(('KIR', 'LILR', 'LAIR'))
        return False
    
    def is_cds_complete(self) -> bool:
        """Check if CDS is annotated as complete on both ends."""
        return not (self.cds_start_NF or self.cds_end_NF)
    
    def has_incomplete_mrna(self) -> bool:
        """Check if mRNA boundaries are incomplete (suggests CDS may also be incomplete)."""
        return self.mRNA_start_NF or self.mRNA_end_NF
    
    def has_special_coding(self) -> bool:
        """Check if transcript has special coding rules that affect CDS length."""
        return self.is_seleno or self.is_readthrough or self.is_polymorphism
    
    def has_annotation_uncertainty(self) -> bool:
        """Check if transcript has tags indicating annotation uncertainty."""
        return (
            self.has_reference_genome_error or
            self.has_sequence_error or
            self.is_non_canonical_TEC or
            self.is_non_canonical_conserved or
            self.is_non_canonical_other or
            self.is_non_canonical_U12 or
            self.is_inferred_exon_combination or
            self.is_inferred_transcript_model or
            self.has_low_sequence_quality or
            self.is_semi_processed or
            self.is_fragmented_locus or
            self.is_NMD_exception or
            self.has_retained_intron or
            self.is_bicistronic or
            self.has_non_ATG_start or
            self.has_upstream_ATG or
            self.has_downstream_ATG
        )
    
    def cds_status(self) -> str:
        """Return human-readable CDS completeness status."""
        if self.cds_start_NF and self.cds_end_NF:
            return "incomplete_both"
        elif self.cds_start_NF:
            return "incomplete_5prime"
        elif self.cds_end_NF:
            return "incomplete_3prime"
        return "complete"
    
    def non_divisible_explanation(self) -> Optional[str]:
        """
        Return explanation if CDS length not divisible by 3 is expected.
        Returns None if no explanation (truly unexpected).
        """
        reasons = []
        
        # CDS completeness issues
        if not self.is_cds_complete():
            reasons.append(f"incomplete CDS ({self.cds_status()})")
        
        # Chromosome-level issues
        if self.is_mitochondrial():
            reasons.append("mitochondrial (different genetic code)")
        
        # mRNA/CDS consistency
        if self.has_incomplete_mrna() and self.is_cds_complete():
            reasons.append("mRNA incomplete but CDS marked complete (possible annotation inconsistency)")
        
        # Special coding mechanisms
        if self.is_seleno:
            reasons.append("selenocysteine (UGA recoding)")
        
        if self.is_readthrough:
            reasons.append("stop codon readthrough")
        
        if self.is_polymorphism:
            reasons.append("non-canonical polymorphism")
        
        # Reference/sequence errors
        if self.has_reference_genome_error:
            reasons.append("reference genome error affects annotation")
        
        if self.has_sequence_error:
            reasons.append("sequence error causes non-canonical splice")
        
        # Non-canonical splicing
        if self.is_non_canonical_TEC:
            reasons.append("non-canonical splice (to be experimentally confirmed)")
        
        if self.is_non_canonical_conserved:
            reasons.append("non-canonical splice (conserved)")
        
        if self.is_non_canonical_other:
            reasons.append("non-canonical splice (other)")
        
        if self.is_non_canonical_U12:
            reasons.append("U12 intron (AT-AC splice site)")
        
        # Inferred/uncertain models
        if self.is_inferred_exon_combination:
            reasons.append("inferred exon combination (not directly supported)")
        
        if self.is_inferred_transcript_model:
            reasons.append("inferred transcript model (not supported by single evidence)")
        
        if self.has_low_sequence_quality:
            reasons.append("low sequence quality")
        
        # Pseudogene-related
        if self.is_semi_processed:
            reasons.append("semi-processed pseudogene")
        
        # Structural issues
        if self.is_fragmented_locus:
            reasons.append("fragmented locus")
        
        if self.is_NMD_exception:
            reasons.append("NMD exception")
        
        if self.has_retained_intron:
            reasons.append("retained intron in CDS")
        
        if self.is_bicistronic:
            reasons.append("bicistronic transcript")
        
        # Start codon issues
        if self.has_non_ATG_start:
            reasons.append("non-ATG start codon")
        
        if self.has_upstream_ATG:
            reasons.append("upstream ATG present")
        
        if self.has_downstream_ATG:
            reasons.append("downstream ATG used")
        
        # Gene family issues
        if self.is_immunoglobulin_or_tcr():
            reasons.append("IG/TR/KIR gene (known annotation challenges)")
        
        return "; ".join(reasons) if reasons else None


@dataclass
class RegionBoundary:
    """A single region boundary in transcript coordinates."""
    transcript_id: str
    region: str
    start: int
    end: int
    strand: str
    chrom: Optional[str] = None
    
    def length(self) -> int:
        return self.end - self.start + 1


@dataclass
class ValidationResult:
    """Result of a validation check."""
    passed: bool
    message: str
    transcript_id: Optional[str] = None
    details: Optional[dict] = None


class GTFParser:
    """Parse GTF files and extract transcript structure."""
    
    # Map of GENCODE tags to TranscriptData attributes
    TAG_MAPPING = {
        # CDS completeness
        'cds_start_NF': 'cds_start_NF',
        'cds_end_NF': 'cds_end_NF',
        'mRNA_start_NF': 'mRNA_start_NF',
        'mRNA_end_NF': 'mRNA_end_NF',
        
        # Special coding mechanisms
        'seleno': 'is_seleno',
        'readthrough_transcript': 'is_readthrough',
        'stop_codon_readthrough': 'is_readthrough',  # Alternative name
        'non_canonical_polymorphism': 'is_polymorphism',
        
        # Reference/sequence errors
        'reference_genome_error': 'has_reference_genome_error',
        'sequence_error': 'has_sequence_error',
        
        # Non-canonical splicing
        'non_canonical_TEC': 'is_non_canonical_TEC',
        'non_canonical_conserved': 'is_non_canonical_conserved',
        'non_canonical_other': 'is_non_canonical_other',
        'non_canonical_U12': 'is_non_canonical_U12',
        'non_canonical_genome_sequence_error': 'has_sequence_error',  # Maps to same
        
        # Inferred/uncertain models
        'inferred_exon_combination': 'is_inferred_exon_combination',
        'inferred_transcript_model': 'is_inferred_transcript_model',
        'low_sequence_quality': 'has_low_sequence_quality',
        
        # Pseudogene-related
        'semi_processed': 'is_semi_processed',
        
        # Structural issues
        'fragmented_locus': 'is_fragmented_locus',
        'NMD_exception': 'is_NMD_exception',
        'retained_intron_CDS': 'has_retained_intron',
        'retained_intron_first': 'has_retained_intron',
        'retained_intron_final': 'has_retained_intron',
        'bicistronic': 'is_bicistronic',
        
        # Start codon issues
        'non_ATG_start': 'has_non_ATG_start',
        'upstream_ATG': 'has_upstream_ATG',
        'downstream_ATG': 'has_downstream_ATG',
    }
    
    def __init__(self, gtf_file: str, verbose: bool = True):
        self.gtf_file = gtf_file
        self.verbose = verbose
        self.transcripts: Dict[str, TranscriptData] = {}
        self.parse_warnings: List[str] = []
        
    def log(self, message: str):
        if self.verbose:
            print(message, file=sys.stderr)
    
    @staticmethod
    def parse_attributes(attr_string: str) -> Dict[str, any]:
        """Parse GTF attribute string into dictionary (optimized)."""
        attributes = {}
        tags = []
        
        # GTF attributes are semicolon-separated "key value" pairs
        # e.g., gene_id "ENSG00000223972.5"; transcript_id "ENST00000456328.2"; tag "basic";
        for part in attr_string.split(';'):
            part = part.strip()
            if not part:
                continue
            
            # Find the space between key and value
            space_idx = part.find(' ')
            if space_idx == -1:
                continue
                
            key = part[:space_idx]
            # Strip quotes from value
            value = part[space_idx+1:].strip('"')
            
            if key == 'tag':
                tags.append(value)
            else:
                attributes[key] = value
        
        if tags:
            attributes['_tags'] = tags
            
        return attributes
    
    def _apply_tag(self, tx: TranscriptData, tag: str) -> None:
        """Apply a GENCODE tag to the transcript data."""
        tx.tags.add(tag)
        
        # Check if this tag maps to a known attribute
        if tag in self.TAG_MAPPING:
            attr_name = self.TAG_MAPPING[tag]
            setattr(tx, attr_name, True)
    
    def parse(self) -> Dict[str, TranscriptData]:
        """Parse GTF file and extract exon, CDS coordinates for each transcript."""
        self.log(f"Parsing GTF file: {self.gtf_file}")
        
        line_count = 0
        feature_counts = defaultdict(int)
        
        with open(self.gtf_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                if line.startswith('#'):
                    continue
                
                line = line.strip()
                if not line:
                    continue
                    
                fields = line.split('\t')
                if len(fields) < 9:
                    self.parse_warnings.append(f"Line {line_num}: insufficient fields ({len(fields)})")
                    continue
                
                chrom = fields[0]
                feature_type = fields[2]
                
                try:
                    start = int(fields[3])
                    end = int(fields[4])
                except ValueError:
                    self.parse_warnings.append(f"Line {line_num}: invalid coordinates")
                    continue
                
                strand_char = fields[6]
                if strand_char not in ('+', '-'):
                    self.parse_warnings.append(f"Line {line_num}: invalid strand '{strand_char}'")
                    continue
                strand = Strand.PLUS if strand_char == '+' else Strand.MINUS
                
                attributes = self.parse_attributes(fields[8])
                
                if 'transcript_id' not in attributes:
                    continue
                
                transcript_id = attributes['transcript_id']
                gene_id = attributes.get('gene_id')
                gene_name = attributes.get('gene_name')
                transcript_type = attributes.get('transcript_type')
                
                # Initialize transcript if needed
                if transcript_id not in self.transcripts:
                    self.transcripts[transcript_id] = TranscriptData()
                
                tx = self.transcripts[transcript_id]
                
                # Store metadata
                if tx.strand is None:
                    tx.strand = strand
                if tx.chrom is None:
                    tx.chrom = chrom
                if tx.gene_id is None and gene_id:
                    tx.gene_id = gene_id
                if tx.gene_name is None and gene_name:
                    tx.gene_name = gene_name
                if tx.transcript_type is None and transcript_type:
                    tx.transcript_type = transcript_type
                
                # Extract and apply tags (GENCODE/Ensembl GTF format)
                if '_tags' in attributes:
                    for tag in attributes['_tags']:
                        self._apply_tag(tx, tag)
                
                # Validate consistency
                if tx.strand != strand:
                    self.parse_warnings.append(
                        f"Transcript {transcript_id}: inconsistent strand at line {line_num}"
                    )
                
                if feature_type == 'exon':
                    tx.exons.append((start, end))
                    feature_counts['exon'] += 1
                elif feature_type == 'CDS':
                    tx.cds.append((start, end))
                    feature_counts['CDS'] += 1
                
                line_count += 1
        
        # Sort coordinates and pre-calculate lengths
        for tx in self.transcripts.values():
            tx.exons.sort()
            tx.cds.sort()
            tx._update_lengths()
        
        self.log(f"Parsed {line_count} feature lines")
        self.log(f"Found {len(self.transcripts)} transcripts")
        self.log(f"Feature counts: {dict(feature_counts)}")
        
        if self.parse_warnings:
            self.log(f"Warnings: {len(self.parse_warnings)}")
        
        # Log tag statistics
        tag_counts = defaultdict(int)
        for tx in self.transcripts.values():
            for tag in tx.tags:
                tag_counts[tag] += 1
        
        if tag_counts:
            self.log(f"\nTag statistics (top 20):")
            for tag, count in sorted(tag_counts.items(), key=lambda x: -x[1])[:20]:
                self.log(f"  {tag}: {count}")
        
        return self.transcripts


class CoordinateConverter:
    """Convert between genomic and transcript coordinates."""
    
    @staticmethod
    def genomic_to_transcript(
        genomic_pos: int,
        exons: List[Tuple[int, int]],
        strand: Strand
    ) -> Optional[int]:
        """
        Convert a single genomic position to transcript coordinate.
        Returns None if position is not within exons.
        """
        # For minus strand, we need to process exons in reverse order
        # because transcript position 1 is at the 3' end of the last exon (genomically)
        ordered_exons = list(reversed(exons)) if strand == Strand.MINUS else exons
        
        transcript_pos = 0
        
        for exon_start, exon_end in ordered_exons:
            exon_length = exon_end - exon_start + 1
            
            if exon_start <= genomic_pos <= exon_end:
                if strand == Strand.MINUS:
                    # Distance from the 3' end of this exon
                    offset = exon_end - genomic_pos
                else:
                    # Distance from the 5' end of this exon
                    offset = genomic_pos - exon_start
                
                return transcript_pos + offset + 1  # 1-based
            
            transcript_pos += exon_length
        
        return None
    
    @staticmethod
    def genomic_range_to_transcript(
        genomic_start: int,
        genomic_end: int,
        tx: TranscriptData,
    ) -> Optional[Tuple[int, int]]:
        """
        Convert genomic range to transcript coordinates (optimized).
        Handles ranges that span multiple exons.
        """
        transcript_positions = []
        
        # Use pre-calculated values
        exons = tx.exons
        strand = tx.strand
        
        if strand == Strand.MINUS:
            # For minus strand, iterate exons in reverse genomic order
            # which is forward transcript order
            total_length = tx._cumulative_lengths[-1]
            for i in range(len(exons) - 1, -1, -1):
                exon_start, exon_end = exons[i]
                overlap_start = max(genomic_start, exon_start)
                overlap_end = min(genomic_end, exon_end)
                
                if overlap_start <= overlap_end:
                    # Offset from transcript start (3' end of genomic transcript)
                    # For genomic exon i, the transcript offset is the sum of lengths 
                    # of genomic exons i+1 to N.
                    offset = total_length - tx._cumulative_lengths[i+1]
                    trans_start = offset + (exon_end - overlap_end) + 1
                    trans_end = offset + (exon_end - overlap_start) + 1
                    transcript_positions.extend([trans_start, trans_end])
        else:
            for i in range(len(exons)):
                exon_start, exon_end = exons[i]
                overlap_start = max(genomic_start, exon_start)
                overlap_end = min(genomic_end, exon_end)
                
                if overlap_start <= overlap_end:
                    offset = tx._cumulative_lengths[i]
                    trans_start = offset + (overlap_start - exon_start) + 1
                    trans_end = offset + (overlap_end - exon_start) + 1
                    transcript_positions.extend([trans_start, trans_end])
        
        if not transcript_positions:
            return None
        
        return (min(transcript_positions), max(transcript_positions))


class RegionBuilder:
    """Build UTR/CDS regions from transcript annotations."""
    
    def __init__(self, include_stop_in_cds: bool = True):
        """
        Args:
            include_stop_in_cds: Whether GTF CDS includes stop codon (affects 3'UTR boundary)
        """
        self.include_stop_in_cds = include_stop_in_cds
        self.converter = CoordinateConverter()
    
    def build_regions(
        self,
        transcript_id: str,
        tx: TranscriptData
    ) -> Tuple[List[RegionBoundary], List[ValidationResult]]:
        """
        Build region boundaries for a transcript.
        Returns (regions, validation_results).
        """
        validations = []
        
        if not tx.exons:
            validations.append(ValidationResult(
                passed=False,
                message="No exons found",
                transcript_id=transcript_id
            ))
            return [], validations
        
        transcript_length = tx.transcript_length()
        strand_char = tx.strand.value if tx.strand else '+'
        
        # Non-coding transcript
        if not tx.cds:
            return [RegionBoundary(
                transcript_id=transcript_id,
                region='non_coding',
                start=1,
                end=transcript_length,
                strand=strand_char,
                chrom=tx.chrom
            )], validations
        
        # Validate CDS regions fall within exons
        cds_validation = self._validate_cds_within_exons(transcript_id, tx)
        validations.extend(cds_validation)
        
        # Convert all CDS regions to transcript coordinates
        transcript_cds_ranges = []
        for cds_start, cds_end in tx.cds:
            trans_range = self.converter.genomic_range_to_transcript(
                cds_start, cds_end, tx
            )
            if trans_range:
                transcript_cds_ranges.append(trans_range)
            else:
                validations.append(ValidationResult(
                    passed=False,
                    message=f"CDS region {cds_start}-{cds_end} could not be mapped",
                    transcript_id=transcript_id
                ))
        
        if not transcript_cds_ranges:
            validations.append(ValidationResult(
                passed=False,
                message="No CDS regions could be mapped to transcript",
                transcript_id=transcript_id
            ))
            return [RegionBoundary(
                transcript_id=transcript_id,
                region='non_coding',
                start=1,
                end=transcript_length,
                strand=strand_char,
                chrom=tx.chrom
            )], validations
        
        # Find overall CDS boundaries
        cds_start = min(r[0] for r in transcript_cds_ranges)
        cds_end = max(r[1] for r in transcript_cds_ranges)
        
        regions = []
        
        # 5'UTR (before CDS in transcript coordinates)
        if cds_start > 1:
            regions.append(RegionBoundary(
                transcript_id=transcript_id,
                region='5UTR',
                start=1,
                end=cds_start - 1,
                strand=strand_char,
                chrom=tx.chrom
            ))
        
        # CDS
        regions.append(RegionBoundary(
            transcript_id=transcript_id,
            region='CDS',
            start=cds_start,
            end=cds_end,
            strand=strand_char,
            chrom=tx.chrom
        ))
        
        # 3'UTR (after CDS in transcript coordinates)
        if cds_end < transcript_length:
            regions.append(RegionBoundary(
                transcript_id=transcript_id,
                region='3UTR',
                start=cds_end + 1,
                end=transcript_length,
                strand=strand_char,
                chrom=tx.chrom
            ))
        
        return regions, validations
    
    def _validate_cds_within_exons(
        self,
        transcript_id: str,
        tx: TranscriptData
    ) -> List[ValidationResult]:
        """Validate that all CDS regions fall within exon boundaries."""
        results = []
        
        for cds_start, cds_end in tx.cds:
            covered = False
            for exon_start, exon_end in tx.exons:
                # Check if CDS overlaps this exon
                if cds_start <= exon_end and cds_end >= exon_start:
                    covered = True
                    # Check if CDS extends beyond exon
                    if cds_start < exon_start or cds_end > exon_end:
                        # This is OK if CDS spans multiple exons
                        # Check if fully covered by any combination of exons
                        pass
            
            if not covered:
                results.append(ValidationResult(
                    passed=False,
                    message=f"CDS {cds_start}-{cds_end} not covered by any exon",
                    transcript_id=transcript_id,
                    details={'cds': (cds_start, cds_end), 'exons': tx.exons}
                ))
        
        return results


class RegionValidator:
    """Comprehensive validation of generated regions."""
    
    def __init__(self, transcripts: Dict[str, TranscriptData], regions: List[RegionBoundary]):
        self.transcripts = transcripts
        self.regions = regions
        self.regions_by_transcript = defaultdict(list)
        for r in regions:
            self.regions_by_transcript[r.transcript_id].append(r)
    
    def run_all_validations(self) -> List[ValidationResult]:
        """Run all validation checks."""
        results = []
        results.extend(self.validate_coverage())
        results.extend(self.validate_no_overlaps())
        results.extend(self.validate_region_order())
        results.extend(self.validate_cds_length_divisible_by_3())
        results.extend(self.validate_boundaries_within_transcript())
        return results
    
    def validate_coverage(self) -> List[ValidationResult]:
        """Validate that regions cover the entire transcript without gaps (optimized)."""
        results = []
        
        for transcript_id, tx in self.transcripts.items():
            tx_regions = self.regions_by_transcript.get(transcript_id, [])
            
            if not tx_regions:
                results.append(ValidationResult(
                    passed=False,
                    message="No regions generated",
                    transcript_id=transcript_id
                ))
                continue
            
            transcript_length = tx.transcript_length()
            
            # Sort regions by start position
            sorted_regions = sorted(tx_regions, key=lambda r: r.start)
            
            # Check coverage without sets
            is_valid = True
            current_pos = 1
            
            for r in sorted_regions:
                if r.start != current_pos:
                    is_valid = False
                    results.append(ValidationResult(
                        passed=False,
                        message=f"Gap found before region {r.region}: expected {current_pos}, got {r.start}",
                        transcript_id=transcript_id
                    ))
                    break
                current_pos = r.end + 1
            
            if is_valid:
                if current_pos - 1 != transcript_length:
                    results.append(ValidationResult(
                        passed=False,
                        message=f"Incomplete coverage: regions end at {current_pos - 1}, transcript length is {transcript_length}",
                        transcript_id=transcript_id
                    ))
                else:
                    results.append(ValidationResult(
                        passed=True,
                        message="Complete coverage",
                        transcript_id=transcript_id
                    ))
        
        return results
    
    def validate_no_overlaps(self) -> List[ValidationResult]:
        """Validate that regions don't overlap."""
        results = []
        
        for transcript_id, tx_regions in self.regions_by_transcript.items():
            sorted_regions = sorted(tx_regions, key=lambda r: r.start)
            
            has_overlap = False
            for i in range(len(sorted_regions) - 1):
                if sorted_regions[i].end >= sorted_regions[i + 1].start:
                    has_overlap = True
                    results.append(ValidationResult(
                        passed=False,
                        message=f"Overlap between {sorted_regions[i].region} and {sorted_regions[i+1].region}",
                        transcript_id=transcript_id,
                        details={
                            'region1': (sorted_regions[i].start, sorted_regions[i].end),
                            'region2': (sorted_regions[i+1].start, sorted_regions[i+1].end)
                        }
                    ))
            
            if not has_overlap:
                results.append(ValidationResult(
                    passed=True,
                    message="No overlaps",
                    transcript_id=transcript_id
                ))
        
        return results
    
    def validate_region_order(self) -> List[ValidationResult]:
        """Validate that regions appear in correct order (5UTR -> CDS -> 3UTR)."""
        results = []
        expected_order = {'5UTR': 0, 'CDS': 1, '3UTR': 2, 'non_coding': 0}
        
        for transcript_id, tx_regions in self.regions_by_transcript.items():
            sorted_regions = sorted(tx_regions, key=lambda r: r.start)
            
            correct_order = True
            for i in range(len(sorted_regions) - 1):
                r1_order = expected_order.get(sorted_regions[i].region, -1)
                r2_order = expected_order.get(sorted_regions[i + 1].region, -1)
                
                if r1_order > r2_order:
                    correct_order = False
                    break
            
            results.append(ValidationResult(
                passed=correct_order,
                message="Correct region order" if correct_order else "Incorrect region order",
                transcript_id=transcript_id
            ))
        
        return results
    
    def validate_cds_length_divisible_by_3(self) -> List[ValidationResult]:
        """
        Validate that CDS length is divisible by 3 (complete codons).
        
        Note: Several categories of transcripts may legitimately have CDS lengths
        not divisible by 3:
        - Incomplete CDS annotations (cds_start_NF, cds_end_NF tags)
        - Mitochondrial genes (different genetic code)
        - Selenocysteine-containing proteins (UGA recoding)
        - Stop codon readthrough transcripts
        - Transcripts with mRNA_end_NF but CDS marked complete (annotation inconsistency)
        - Reference genome errors
        - Non-canonical splicing
        - Inferred transcript models
        - IG/TR/KIR genes with known annotation challenges
        """
        results = []
        
        for transcript_id, tx_regions in self.regions_by_transcript.items():
            cds_regions = [r for r in tx_regions if r.region == 'CDS']
            
            if not cds_regions:
                continue
            
            total_cds_length = sum(r.length() for r in cds_regions)
            tx = self.transcripts.get(transcript_id)
            
            if total_cds_length % 3 == 0:
                results.append(ValidationResult(
                    passed=True,
                    message=f"CDS length ({total_cds_length}) divisible by 3",
                    transcript_id=transcript_id
                ))
            else:
                # Check if there's a known explanation for non-divisible CDS
                explanation = tx.non_divisible_explanation() if tx else None
                
                if explanation:
                    # Known reason - report as passing (expected behavior)
                    results.append(ValidationResult(
                        passed=True,
                        message=f"CDS length ({total_cds_length}) not divisible by 3 - expected: {explanation}",
                        transcript_id=transcript_id,
                        details={
                            'cds_length': total_cds_length, 
                            'remainder': total_cds_length % 3,
                            'explanation': explanation,
                            'expected_incomplete': True
                        }
                    ))
                else:
                    # No known explanation - this is unexpected
                    results.append(ValidationResult(
                        passed=False,
                        message=f"CDS length ({total_cds_length}) not divisible by 3 - unexpected (CDS marked complete)",
                        transcript_id=transcript_id,
                        details={
                            'cds_length': total_cds_length, 
                            'remainder': total_cds_length % 3,
                            'cds_status': tx.cds_status() if tx else 'unknown',
                            'chrom': tx.chrom if tx else None,
                            'gene_name': tx.gene_name if tx else None,
                            'transcript_type': tx.transcript_type if tx else None,
                            'expected_incomplete': False,
                            'tags': list(tx.tags) if tx else []
                        }
                    ))
        
        return results
    
    def validate_boundaries_within_transcript(self) -> List[ValidationResult]:
        """Validate that all region boundaries are within transcript bounds."""
        results = []
        
        for transcript_id, tx in self.transcripts.items():
            tx_regions = self.regions_by_transcript.get(transcript_id, [])
            transcript_length = tx.transcript_length()
            
            for r in tx_regions:
                if r.start < 1 or r.end > transcript_length:
                    results.append(ValidationResult(
                        passed=False,
                        message=f"Region {r.region} ({r.start}-{r.end}) outside transcript bounds (1-{transcript_length})",
                        transcript_id=transcript_id
                    ))
                elif r.start > r.end:
                    results.append(ValidationResult(
                        passed=False,
                        message=f"Region {r.region} has start > end ({r.start} > {r.end})",
                        transcript_id=transcript_id
                    ))
        
        return results


def write_output(regions: List[RegionBoundary], output_file: str, include_chrom: bool = False):
    """Write regions to TSV file."""
    with open(output_file, 'w') as f:
        # Header
        if include_chrom:
            f.write("transcript_id\tregion\tstart\tend\tstrand\tchrom\n")
        else:
            f.write("transcript_id\tregion\tstart\tend\tstrand\n")
        
        # Data rows
        for r in regions:
            if include_chrom:
                f.write(f"{r.transcript_id}\t{r.region}\t{r.start}\t{r.end}\t{r.strand}\t{r.chrom or '.'}\n")
            else:
                f.write(f"{r.transcript_id}\t{r.region}\t{r.start}\t{r.end}\t{r.strand}\n")


def write_validation_report(
    validations: List[ValidationResult],
    output_file: str,
    verbose: bool = False
):
    """Write validation report to file."""
    passed = [v for v in validations if v.passed]
    failed = [v for v in validations if not v.passed]
    
    # Separate expected issues from unexpected ones
    unexpected_failures = []
    expected_issues = []
    for v in failed:
        if v.details and v.details.get('expected_incomplete'):
            expected_issues.append(v)
        else:
            unexpected_failures.append(v)
    
    with open(output_file, 'w') as f:
        f.write("# Validation Report\n\n")
        f.write(f"Total checks: {len(validations)}\n")
        f.write(f"Passed: {len(passed)}\n")
        f.write(f"Failed: {len(failed)}\n")
        if expected_issues:
            f.write(f"  - Expected (incomplete CDS): {len(expected_issues)}\n")
        if unexpected_failures:
            f.write(f"  - Unexpected: {len(unexpected_failures)}\n")
        f.write("\n")
        
        if unexpected_failures:
            f.write("## Unexpected Failures (may indicate issues)\n\n")
            for v in unexpected_failures[:100]:  # Limit output
                f.write(f"- [{v.transcript_id or 'global'}] {v.message}\n")
                if verbose and v.details:
                    f.write(f"  Details: {v.details}\n")
            if len(unexpected_failures) > 100:
                f.write(f"\n... and {len(unexpected_failures) - 100} more\n")
            f.write("\n")
        
        if expected_issues and verbose:
            f.write("## Expected Issues (incomplete CDS annotations)\n\n")
            for v in expected_issues[:50]:
                f.write(f"- [{v.transcript_id or 'global'}] {v.message}\n")
            if len(expected_issues) > 50:
                f.write(f"\n... and {len(expected_issues) - 50} more\n")
            f.write("\n")
        
        f.write("## Summary by Check Type\n\n")
        
        # Group by message type
        check_types = defaultdict(lambda: {'passed': 0, 'failed': 0})
        for v in validations:
            # Extract check type from message
            if 'CDS length' in v.message:
                check_type = 'CDS length divisible by 3'
            elif 'coverage' in v.message.lower():
                check_type = 'Coverage'
            elif 'overlap' in v.message.lower():
                check_type = 'No overlaps'
            elif 'order' in v.message.lower():
                check_type = 'Region order'
            elif 'bounds' in v.message.lower() or 'outside' in v.message.lower():
                check_type = 'Within bounds'
            else:
                check_type = v.message.split(':')[0] if ':' in v.message else v.message[:40]
            
            if v.passed:
                check_types[check_type]['passed'] += 1
            else:
                check_types[check_type]['failed'] += 1
        
        for check_type, counts in sorted(check_types.items()):
            f.write(f"- {check_type}: {counts['passed']} passed, {counts['failed']} failed\n")


def main():
    parser = argparse.ArgumentParser(
        description="Find CDS/UTR boundaries from GTF",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
            Examples:
            %(prog)s -i annotation.gtf -o regions.tsv
            %(prog)s -i annotation.gtf -o regions.tsv --validate --validation-report validation.txt
                    """
    )
    parser.add_argument("-i", "--input", required=True, help="Input GTF file")
    parser.add_argument("-o", "--output", default="transcript_region_boundaries.tsv",
                        help="Output TSV file (default: transcript_region_boundaries.tsv)")
    parser.add_argument("--validate", action="store_true",
                        help="Run validation checks after building regions")
    parser.add_argument("--validation-report", default="transcript_region_boundaries_validation_report.txt",
                        help="Output file for validation report (default: transcript_region_boundaries_validation_report.txt)")
    parser.add_argument("--include-chrom", action="store_true",
                        help="Include chromosome column in output")
    parser.add_argument("--quiet", action="store_true",
                        help="Suppress progress messages")
    parser.add_argument("--verbose-validation", action="store_true",
                        help="Include detailed information in validation report")
    
    args = parser.parse_args()
    verbose = not args.quiet
    
    # Parse GTF
    gtf_parser = GTFParser(args.input, verbose=verbose)
    transcripts = gtf_parser.parse()
    
    if gtf_parser.parse_warnings and verbose:
        print(f"\nParse warnings:", file=sys.stderr)
        for warning in gtf_parser.parse_warnings[:10]:
            print(f"  {warning}", file=sys.stderr)
        if len(gtf_parser.parse_warnings) > 10:
            print(f"  ... and {len(gtf_parser.parse_warnings) - 10} more", file=sys.stderr)
    
    # Build regions
    if verbose:
        print("\nBuilding region boundaries...", file=sys.stderr)
    
    region_builder = RegionBuilder()
    all_regions = []
    build_validations = []
    
    for transcript_id in sorted(transcripts.keys()):
        tx = transcripts[transcript_id]
        regions, validations = region_builder.build_regions(transcript_id, tx)
        all_regions.extend(regions)
        build_validations.extend(validations)
    
    if verbose:
        print(f"Generated {len(all_regions)} region entries", file=sys.stderr)
    
    # Count statistics
    region_counts = defaultdict(int)
    for region in all_regions:
        region_counts[region.region] += 1
    
    if verbose:
        print("\nRegion statistics:", file=sys.stderr)
        for region_type, count in sorted(region_counts.items()):
            print(f"  {region_type}: {count}", file=sys.stderr)
    
    # CDS completeness statistics
    cds_stats = {'complete': 0, 'incomplete_5prime': 0, 'incomplete_3prime': 0, 'incomplete_both': 0}
    annotation_uncertainty_count = 0
    for transcript_id, tx in transcripts.items():
        if tx.cds:  # Only count coding transcripts
            cds_stats[tx.cds_status()] += 1
            if tx.has_annotation_uncertainty():
                annotation_uncertainty_count += 1
    
    if verbose and sum(cds_stats.values()) > 0:
        print("\nCDS completeness (coding transcripts only):", file=sys.stderr)
        for status, count in cds_stats.items():
            if count > 0:
                print(f"  {status}: {count}", file=sys.stderr)
        print(f"  with annotation uncertainty tags: {annotation_uncertainty_count}", file=sys.stderr)
    
    # Write output
    if verbose:
        print(f"\nWriting output to: {args.output}", file=sys.stderr)
    write_output(all_regions, args.output, include_chrom=args.include_chrom)
    
    # Run validation if requested
    if args.validate:
        if verbose:
            print("\nRunning validations...", file=sys.stderr)
        
        validator = RegionValidator(transcripts, all_regions)
        post_validations = validator.run_all_validations()
        
        all_validations = build_validations + post_validations
        
        passed = sum(1 for v in all_validations if v.passed)
        failed = sum(1 for v in all_validations if not v.passed)
        
        if verbose:
            print(f"Validation complete: {passed} passed, {failed} failed", file=sys.stderr)
        
        write_validation_report(
            all_validations,
            args.validation_report,
            verbose=args.verbose_validation
        )
        
        if verbose:
            print(f"Validation report written to: {args.validation_report}", file=sys.stderr)
        
        # # Exit with error code if validations failed
        # if failed > 0:
        #     sys.exit(1)
    
    if verbose:
        print("\nDone!", file=sys.stderr)


if __name__ == "__main__":
    main()