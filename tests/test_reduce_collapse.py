import pytest
import pandas as pd
import numpy as np
import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))
from reduce_collapse import collapse_guides

@pytest.fixture
def expanded_df():
    data = {
        'Guide Sequence': ['G1', 'G1', 'G2'],
        'Target Sequence': ['TSEQ1', 'TSEQ1', 'TSEQ2'],
        'Guide Score': [0.9, 0.9, 0.8],
        'Title': ['PC', 'PC', 'PC'],
        'Gene': ['G1', 'G1', 'G1'],
        'Symbol': ['T1', 'T2', 'T1'],
        'Contig_Idx': [100, 200, 100],
        'exon_id': ['E1', 'E2', 'E1'],
        'region': ['CDS', 'CDS', 'CDS'],
        'tags': ['tag1', 'tag2', 'tag1'],
        'median_expression_norm': [0.7, 0.3, 0.7]
    }
    return pd.DataFrame(data)

def test_collapse_guides_aggregation(expanded_df):
    coll_df = collapse_guides(expanded_df)
    
    # Check shape
    assert len(coll_df) == 2
    
    # Check G1 aggregation
    g1 = coll_df[coll_df['Guide Sequence'] == 'G1'].iloc[0]
    assert g1['sum_median_expression_norm'] == pytest.approx(1.0)
    assert g1['ntargeted_tx'] == 2
    assert 'T1' in g1['Symbol'] and 'T2' in g1['Symbol']
    
    # Check G2 aggregation
    g2 = coll_df[coll_df['Guide Sequence'] == 'G2'].iloc[0]
    assert g2['sum_median_expression_norm'] == pytest.approx(0.7)
    assert g2['ntargeted_tx'] == 1

def test_collapse_guides_no_expression(expanded_df):
    df_no_expr = expanded_df.drop(columns=['median_expression_norm'])
    coll_df = collapse_guides(df_no_expr)
    
    assert 'sum_median_expression_norm' in coll_df.columns
    assert (coll_df['sum_median_expression_norm'] == 0.0).all()
