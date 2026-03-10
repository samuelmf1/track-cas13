import unittest
import sys
import os
import pandas as pd

# Add src to path to allow import
sys.path.append(os.path.join(os.path.dirname(__file__), '../src'))

from add_cds_annotations import get_region

class TestAddCDSAnnotations(unittest.TestCase):
    def setUp(self):
        # Mock Data
        self.cds_dict = {
            'T1_Full': {'start': 100, 'end': 200},
            'T2_Trimmed': {'start': 100, 'end': 200},
            'T3_NoCDS': {'start': 100, 'end': 200}, # Only in CDS dict but used for testing
        }
        
        self.summary_info = {
            'T1_Full': {'biotype': 'protein_coding', 'is_representative': True},
            'T2_Trimmed': {'biotype': 'protein_coding', 'is_representative': False},
            'T3_lncRNA': {'biotype': 'lncRNA', 'is_representative': True},
            'T4_Other': {'biotype': 'nonsense_mediated_decay', 'is_representative': True},
        }

    def test_unknown_transcript(self):
        row = {'Symbol': 'Unknown', 'Contig_Idx': 100}
        self.assertEqual(get_region(row, self.summary_info, self.cds_dict), "Unknown_Transcript")

    def test_lncrna(self):
        row = {'Symbol': 'T3_lncRNA', 'Contig_Idx': 100}
        self.assertEqual(get_region(row, self.summary_info, self.cds_dict), "lncRNA")

    def test_other_biotype(self):
        row = {'Symbol': 'T4_Other', 'Contig_Idx': 100}
        self.assertEqual(get_region(row, self.summary_info, self.cds_dict), "nonsense_mediated_decay")

    def test_no_cds_found(self):
        # T1_Full has CDS, but let's mock a case where it doesn't
        summary_only = {'T_NoCDS': {'biotype': 'protein_coding', 'is_representative': True}}
        row = {'Symbol': 'T_NoCDS', 'Contig_Idx': 100}
        self.assertEqual(get_region(row, summary_only, self.cds_dict), "No_CDS_Found")

    # --- Representative / Full Sequence Tests ---
    # CDS: [100, 200)

    def test_rep_cds_middle(self):
        # Guide at 150. Length 23 -> [150, 173). Fully inside [100, 200).
        row = {'Symbol': 'T1_Full', 'Contig_Idx': 150}
        self.assertEqual(get_region(row, self.summary_info, self.cds_dict), "CDS")

    def test_rep_5utr(self):
        # Guide at 50. [50, 73). End < Start (100).
        row = {'Symbol': 'T1_Full', 'Contig_Idx': 50}
        self.assertEqual(get_region(row, self.summary_info, self.cds_dict), "5'UTR")

    def test_rep_3utr(self):
        # Guide at 250. [250, 273). Start >= End (200).
        row = {'Symbol': 'T1_Full', 'Contig_Idx': 250}
        self.assertEqual(get_region(row, self.summary_info, self.cds_dict), "3'UTR")
    
    def test_rep_cds_overlap_start_utr(self):
        # Overlap 1bp or 2bp -> UTR. Overlap >= 3bp -> CDS.
        # CDS Start 100.
        
        # Guide End at 101. Start = 101 - 23 = 78. [78, 101). Overlap [100, 101) = 1bp.
        row = {'Symbol': 'T1_Full', 'Contig_Idx': 78}
        self.assertEqual(get_region(row, self.summary_info, self.cds_dict), "5'UTR")

        # Guide End at 102. Start = 79. [79, 102). Overlap 2bp.
        row = {'Symbol': 'T1_Full', 'Contig_Idx': 79}
        self.assertEqual(get_region(row, self.summary_info, self.cds_dict), "5'UTR")

        # Guide End at 103. Start = 80. [80, 103). Overlap 3bp.
        row = {'Symbol': 'T1_Full', 'Contig_Idx': 80}
        self.assertEqual(get_region(row, self.summary_info, self.cds_dict), "CDS")

    def test_rep_cds_overlap_end_utr(self):
        # CDS End 200.
        # Guide Start 197. [197, 220). Overlap [197, 200) = 3bp.
        row = {'Symbol': 'T1_Full', 'Contig_Idx': 197}
        self.assertEqual(get_region(row, self.summary_info, self.cds_dict), "CDS")

        # Guide Start 198. [198, 221). Overlap [198, 200) = 2bp.
        # Should be UTR. Start 198 >= Start 100 -> 3'UTR.
        row = {'Symbol': 'T1_Full', 'Contig_Idx': 198}
        self.assertEqual(get_region(row, self.summary_info, self.cds_dict), "3'UTR")

    # --- Non-Representative / Trimmed Tests ---
    # CDS: [100, 200)
    # Trimmed: FASTA[0] = Genomic[100 - 3] = Genomic[97].
    # Mapping: Genomic = Guide + 97.

    def test_trimmed_cds_start(self):
        # Guide at 0. Genomic = 97.
        # Overlap with [100, 200):
        # Guide [97, 120). Overlap [100, 120) = 20bp.
        # Should be CDS.
        row = {'Symbol': 'T2_Trimmed', 'Contig_Idx': 0}
        self.assertEqual(get_region(row, self.summary_info, self.cds_dict), "CDS")

    def test_trimmed_cds_middle(self):
        # Guide at 53. Genomic = 53 + 97 = 150.
        # Fully inside.
        row = {'Symbol': 'T2_Trimmed', 'Contig_Idx': 53}
        self.assertEqual(get_region(row, self.summary_info, self.cds_dict), "CDS")

    def test_trimmed_3utr(self):
        # Guide at 200. Genomic = 200 + 97 = 297.
        # > 200. 3'UTR.
        row = {'Symbol': 'T2_Trimmed', 'Contig_Idx': 200}
        self.assertEqual(get_region(row, self.summary_info, self.cds_dict), "3'UTR")

if __name__ == '__main__':
    unittest.main()
