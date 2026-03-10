import pytest
import pandas as pd
import numpy as np
import sys
import os

# Add src to path
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))

# Mock matplotlib and seaborn before they are imported by preproc_expression
from unittest.mock import MagicMock
sys.modules['matplotlib'] = MagicMock()
sys.modules['matplotlib.pyplot'] = MagicMock()
sys.modules['matplotlib.ticker'] = MagicMock()
sys.modules['seaborn'] = MagicMock()

from preproc_expression import expression
from reduce_collapse import collapse_guides

def test_effective_tpm_logic():
    # 1. Create Mock Data
    # Gene A has 2 transcripts: T1 (75% of mass) and T2 (25% of mass)
    # Total Gene TPM in Cell Line 1 = 100
    # Total Gene TPM in Cell Line 2 = 100
    # In DepMap format (log2(TPM+1)): log2(101) = 6.658
    
    # We'll set Linear TPMs directly in our logic test
    # T1 = 75, T2 = 25
    # log2(76) = 6.248, log2(26) = 4.700
    
    mapping_df = pd.DataFrame({
        'transcript_id': ['T1.1', 'T2.1'],
        'gene_id': ['GENE_A', 'GENE_A'],
        'gene_name': ['A', 'A'],
        'biotype': ['protein_coding', 'protein_coding'],
        'transcript_support_level': ['1', '1'],
        'tag': ['basic', 'basic']
    })
    
    # Simulate 10 cell lines, all with same expression (top 10% will be these values)
    omic_df = pd.DataFrame({
        'ModelID': [f'M{i}' for i in range(10)],
        'T1.1': [6.248] * 10,
        'T2.1': [4.700] * 10
    })
    
    # 2. Run Expression Preprocessing
    expr_results = expression(
        mapping_df, omic_df, 
        biotype='protein_coding', 
        tpm_threshold=1.0, 
        n_cell_lines=1
    )
    
    # Gene Operational TPM should be log2(75 + 25 + 1) = log2(101) = 6.658
    # sum_expression_norm for T1 = 75/(75+25) = 0.75
    # sum_expression_norm for T2 = 25/(75+25) = 0.25
    
    t1_res = expr_results[expr_results['transcript_id'] == 'T1.1'].iloc[0]
    t2_res = expr_results[expr_results['transcript_id'] == 'T2.1'].iloc[0]
    
    assert np.isclose(t1_res['gene_log_median_expr'], 6.658, atol=1e-2)
    assert np.isclose(t1_res['sum_expression_norm'], 0.75, atol=1e-2)
    assert np.isclose(t2_res['sum_expression_norm'], 0.25, atol=1e-2)
    
    # Effective TPM = Fraction * GeneExpression
    # T1 Effective TPM = 0.75 * 6.658 = 4.9935
    # T2 Effective TPM = 0.25 * 6.658 = 1.6645
    assert np.isclose(t1_res['effective_tpm'], 4.9935, atol=1e-2)
    assert np.isclose(t2_res['effective_tpm'], 1.6645, atol=1e-2)
    
    # 3. Test Collapse (Additivity)
    # A guide G1 targeting BOTH T1 and T2 should have guide_expression = 4.9935 + 1.6645 = 6.658
    # A guide G2 targeting ONLY T1 should have guide_expression = 4.9935
    
    expanded_df = pd.DataFrame({
        'Guide Sequence': ['G1', 'G1', 'G2'],
        'Target Sequence': ['TSEG1', 'TSEG1', 'TSEG2'],
        'Symbol': ['T1.1', 'T2.1', 'T1.1'],
        'Gene': ['GENE_A', 'GENE_A', 'GENE_A'],
        'effective_tpm': [t1_res['effective_tpm'], t2_res['effective_tpm'], t1_res['effective_tpm']],
        'is_expressed': [True, True, True],
        'pce': [100.0, 100.0, 100.0],
        'n_expressed': [10, 10, 10],
        'exon_id': ['E1', 'E2', 'E1'],
        'region': ['CDS', 'CDS', 'CDS'],
        'tags': ['tag1', 'tag2', 'tag1'],
        'Guide Score': [0.9, 0.9, 0.9],
        'Title': ['A', 'A', 'A'],
        'Contig_Idx': [100, 100, 200]
    })
    
    collapsed = collapse_guides(expanded_df)
    
    g1_coll = collapsed[collapsed['Guide Sequence'] == 'G1'].iloc[0]
    g2_coll = collapsed[collapsed['Guide Sequence'] == 'G2'].iloc[0]
    
    # Guide G1 hits both, should equal gene expression
    assert np.isclose(g1_coll['guide_expression'], 6.658, atol=1e-2)
    # Guide G2 hits only T1
    assert np.isclose(g2_coll['guide_expression'], 4.9935, atol=1e-2)
    
    # Gene Max Expression should be 6.658 (G1)
    # guide_expression_norm for G1 = 1.0
    # guide_expression_norm for G2 = 4.9935 / 6.658 = 0.75
    assert np.isclose(g1_coll['guide_expression_norm'], 1.0, atol=1e-2)
    assert np.isclose(g2_coll['guide_expression_norm'], 0.75, atol=1e-2)

if __name__ == "__main__":
    pytest.main([__file__])
