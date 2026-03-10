import unittest
import sys
import os
import pandas as pd
import tempfile
import shutil
from io import StringIO

# Add src to path to allow import
sys.path.append(os.path.join(os.path.dirname(__file__), '../src'))

from preproc_reduce_fasta import (
    is_valid_dna,
    parse_header,
    has_no_low_conf_tag,
    has_basic_or_ccds_tag,
    parse_gtf_attributes,
    load_gtf_mapping,
    main
)

class TestPreprocReduceFasta(unittest.TestCase):
    def setUp(self):
        self.test_dir = tempfile.mkdtemp()
        
    def tearDown(self):
        shutil.rmtree(self.test_dir)

    def test_is_valid_dna(self):
        self.assertTrue(is_valid_dna("ACTGN"))
        self.assertTrue(is_valid_dna("actgn"))
        self.assertTrue(is_valid_dna("A"))
        self.assertFalse(is_valid_dna("ACTGR"))
        self.assertFalse(is_valid_dna("ACTG1"))
        self.assertTrue(is_valid_dna(""))

    def test_parse_header(self):
        header = ">ENST00000456328.2|ENSG00000223972.5|OTTHUMG00000000961.2|OTTHUMT00000002844.2|DDX11L1-202|DDX11L1|1653|processed_transcript|"
        expected = ">ENSG00000223972.5|ENST00000456328.2|DDX11L1|processed_transcript"
        self.assertEqual(parse_header(header), expected)

    def test_gtf_tag_parsing(self):
        self.assertTrue(has_no_low_conf_tag("basic,CCDS"))
        self.assertFalse(has_no_low_conf_tag("low_conf,basic"))
        self.assertTrue(has_no_low_conf_tag(""))
        
        self.assertTrue(has_basic_or_ccds_tag("basic"))
        self.assertTrue(has_basic_or_ccds_tag("CCDS"))
        self.assertTrue(has_basic_or_ccds_tag("basic,something_else"))
        self.assertFalse(has_basic_or_ccds_tag("other_tag"))
        self.assertFalse(has_basic_or_ccds_tag(""))

    def test_parse_gtf_attributes(self):
        attr = 'gene_id "ENSG1"; transcript_id "ENST1"; gene_name "GENE1"; transcript_type "protein_coding"; tag "basic"; tag "CCDS";'
        parsed = parse_gtf_attributes(attr)
        self.assertEqual(parsed['gene_id'], "ENSG1")
        self.assertEqual(parsed['transcript_id'], "ENST1")
        self.assertEqual(parsed['tag'], "basic,CCDS")

    def test_load_gtf_mapping(self):
        gtf_content = (
            'chr1\tHAVANA\ttranscript\t1\t100\t.\t+\t.\tgene_id "ENSG1"; transcript_id "ENST1"; gene_name "G1"; transcript_type "protein_coding"; tag "Ensembl_canonical"; tag "basic";\n'
            'chr1\tHAVANA\ttranscript\t200\t300\t.\t+\t.\tgene_id "ENSG1"; transcript_id "ENST2"; gene_name "G1"; transcript_type "protein_coding"; tag "basic";\n'
        )
        gtf_path = os.path.join(self.test_dir, "test.gtf")
        with open(gtf_path, "w") as f:
            f.write(gtf_content)
        
        mapping = load_gtf_mapping(gtf_path)
        self.assertEqual(len(mapping), 2)
        # ENST1 should be representative
        rep = mapping[mapping['transcript_id'] == 'ENST1']
        self.assertTrue(rep.iloc[0]['is_representative'])
        
        non_rep = mapping[mapping['transcript_id'] == 'ENST2']
        self.assertFalse(non_rep.iloc[0]['is_representative'])

    def test_main_with_longest_utr(self):
        # 1. Setup mock GTF
        gtf_content = (
            'chr1\tHAVANA\ttranscript\t1\t1000\t.\t+\t.\tgene_id "ENSG1"; transcript_id "ENST_REP"; gene_name "G1"; transcript_type "protein_coding"; tag "Ensembl_canonical"; tag "basic";\n'
            'chr1\tHAVANA\ttranscript\t1\t1000\t.\t+\t.\tgene_id "ENSG1"; transcript_id "ENST_LONG_UTR"; gene_name "G1"; transcript_type "protein_coding"; tag "basic";\n'
            'chr1\tHAVANA\ttranscript\t1\t1000\t.\t+\t.\tgene_id "ENSG1"; transcript_id "ENST_SHORT"; gene_name "G1"; transcript_type "protein_coding"; tag "basic";\n'
        )
        gtf_path = os.path.join(self.test_dir, "test.gtf")
        with open(gtf_path, "w") as f:
            f.write(gtf_content)

        # 2. Setup mock CDS boundaries
        # ENST_REP: CDS 100-900 (UTR 199)
        # ENST_LONG_UTR: CDS 400-600 (UTR 799) -> Should be marked is_longest_utr
        # ENST_SHORT: CDS 50-950 (UTR 99)
        cds_content = (
            "transcript_id\tregion\tstart\tend\tstrand\n"
            "ENST_REP\t5UTR\t1\t99\t+\n"
            "ENST_REP\tCDS\t100\t900\t+\n"
            "ENST_REP\t3UTR\t901\t1000\t+\n"
            "ENST_LONG_UTR\t5UTR\t1\t399\t+\n"
            "ENST_LONG_UTR\tCDS\t400\t600\t+\n"
            "ENST_LONG_UTR\t3UTR\t601\t1000\t+\n"
            "ENST_SHORT\t5UTR\t1\t49\t+\n"
            "ENST_SHORT\tCDS\t50\t950\t+\n"
            "ENST_SHORT\t3UTR\t951\t1000\t+\n"
        )
        cds_path = os.path.join(self.test_dir, "cds.tsv")
        with open(cds_path, "w") as f:
            f.write(cds_content)

        # 3. Setup mock FASTA
        # Sequence must be same length as transcript (1000)
        seq = "A" * 1000
        fasta_content = (
            f">ENST_REP|ENSG1||||G1||protein_coding\n{seq}\n"
            f">ENST_LONG_UTR|ENSG1||||G1||protein_coding\n{seq}\n"
            f">ENST_SHORT|ENSG1||||G1||protein_coding\n{seq}\n"
        )
        fasta_path = os.path.join(self.test_dir, "test.fa")
        with open(fasta_path, "w") as f:
            f.write(fasta_content)

        output_fasta = os.path.join(self.test_dir, "output.fasta")
        summary_path = os.path.join(self.test_dir, "summary.tsv")

        # 4. Run main
        main(
            input_path=fasta_path,
            output_path=output_fasta,
            cds_boundaries_path=cds_path,
            summary_path=summary_path,
            gtf_path=gtf_path,
            min_length=10
        )

        # 5. Verify results
        summary_df = pd.read_table(summary_path)
        
        # Check flags
        rep_row = summary_df[summary_df['transcript_id'] == 'ENST_REP'].iloc[0]
        long_row = summary_df[summary_df['transcript_id'] == 'ENST_LONG_UTR'].iloc[0]
        short_row = summary_df[summary_df['transcript_id'] == 'ENST_SHORT'].iloc[0]
        
        self.assertTrue(rep_row['is_representative'])
        self.assertFalse(rep_row['is_longest_utr'])
        self.assertEqual(rep_row['cds_length'], 1000) # Preserved

        self.assertFalse(long_row['is_representative'])
        self.assertTrue(long_row['is_longest_utr'])
        self.assertEqual(long_row['cds_length'], 1000) # Preserved

        self.assertFalse(short_row['is_representative'])
        self.assertFalse(short_row['is_longest_utr'])
        # Short should be trimmed: 950 - (50-3) + 1 = 904 if using start-3:end logic in script 
        # (Start is 50, End is 950 -> full_seq[50-3:950] which is indices 47 to 950 -> length 903)
        # Wait, the script does: cds_seq = full_seq[start-3:end]
        # start=50, end=950 -> full_seq[47:950] -> length 903.
        self.assertEqual(short_row['cds_length'], 903) 

        # Verify output FASTA has 3 sequences
        with open(output_fasta, 'r') as f:
            lines = f.readlines()
            # 3 headers + 3 sequences = 6 lines
            self.assertEqual(len(lines), 6)

if __name__ == '__main__':
    unittest.main()
