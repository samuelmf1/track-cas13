import pytest
import pandas as pd
import sys
import os
import tempfile
import re

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))
from reduce_utr_library import SelectionEngine, SelectionConfig, GuideSelector, load_utr_data

class TestUTRGuideSelection:
    @pytest.fixture
    def mock_utr_data(self):
        # tx_id -> {'five_prime_UTR': len, 'three_prime_UTR': len, 'generic': len}
        return {
            'ENST1': {'five_prime_UTR': 100, 'three_prime_UTR': 200, 'generic': 0},
            'ENST2': {'five_prime_UTR': 160, 'three_prime_UTR': 200, 'generic': 0}, # 60bp diff in 5'
            'ENST3': {'five_prime_UTR': 150, 'three_prime_UTR': 210, 'generic': 0}, # 50bp diff in 5' from ENST1
            'ENST4': {'five_prime_UTR': 100, 'three_prime_UTR': 300, 'generic': 0}, # 100bp diff in 3' from ENST1
            'ENST_NC': {'five_prime_UTR': 0, 'three_prime_UTR': 0, 'generic': 500} # Non-coding
        }

    def test_disjoint_selection_with_delta(self, mock_utr_data):
        # Create a mock dataframe
        data = {
            'gene_id': ['G1']*5,
            'transcript_id_group': ['ENST1', 'ENST2', 'ENST3', 'ENST4', 'ENST_NC'],
            'tiger_score': [0.9, 0.8, 0.7, 0.6, 0.5],
            'position': [10, 20, 30, 40, 50],
            'biotype': ['protein_coding']*5
        }
        df = pd.DataFrame(data)
        
        config = SelectionConfig(target_n=1)
        selector = GuideSelector(target_n=1, utr_data=mock_utr_data)
        res_df, summary = selector.process_gene('G1', df, config)
        
        processed_groups = res_df['targeted_transcript'].unique()
        # The new logic groups into ensembles
        assert 'Short_UTR_Ensemble' in processed_groups
        assert 'Long_UTR_Ensemble' in processed_groups

    def test_gtf_parsing(self):
        # Create a mock GTF file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.gtf', delete=False) as f:
            f.write('chr1\t.\tCDS\t100\t201\t.\t+\t.\ttranscript_id "T1";\n')
            f.write('chr1\t.\tUTR\t50\t80\t.\t+\t.\ttranscript_id "T1";\n') # 5' UTR, len 31
            f.write('chr1\t.\tUTR\t220\t250\t.\t+\t.\ttranscript_id "T1";\n') # 3' UTR, len 31
            
            f.write('chr1\t.\tCDS\t100\t200\t.\t-\t.\ttranscript_id "T2";\n')
            f.write('chr1\t.\tUTR\t50\t80\t.\t-\t.\ttranscript_id "T2";\n') # 3' UTR, len 31
            f.write('chr1\t.\tUTR\t220\t250\t.\t-\t.\ttranscript_id "T2";\n') # 5' UTR, len 31
            
            f.write('chr1\t.\tUTR\t1000\t1100\t.\t+\t.\ttranscript_id "T_NC";\n') # Generic UTR, len 101
            gtf_path = f.name

        try:
            utr_data = load_utr_data(gtf_path)
            assert utr_data['T1']['five_prime_UTR'] == 31
            assert utr_data['T1']['three_prime_UTR'] == 31
            assert utr_data['T2']['five_prime_UTR'] == 31
            assert utr_data['T2']['three_prime_UTR'] == 31
            assert utr_data['T_NC']['generic'] == 101
        finally:
            os.remove(gtf_path)

@pytest.mark.parametrize("strand, start, end, cds_min, cds_max, expected", [
    ('+', 50, 80, 100, 200, '5'),
    ('+', 220, 250, 100, 200, '3'),
    ('-', 50, 80, 100, 200, '3'),
    ('-', 220, 250, 100, 200, '5'),
])
def test_utr_inference_logic(strand, start, end, cds_min, cds_max, expected):
    # This matches the logic in the script
    if strand == '+':
        if end <= cds_min: res = '5'
        elif start >= cds_max: res = '3'
        else: res = None
    else:
        if start >= cds_max: res = '5'
        elif end <= cds_min: res = '3'
        else: res = None
    assert res == expected
