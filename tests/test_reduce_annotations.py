import pytest
import pandas as pd
import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))
from reduce_expression_annotations import add_expression_info

@pytest.fixture
def guides_df():
    data = {
        'Guide Sequence': ['G1', 'G1', 'G2'],
        'Symbol': ['T1', 'T2', 'T3']
    }
    return pd.DataFrame(data)

@pytest.fixture
def expr_df():
    data = {
        'transcript_id': ['T1', 'T2'],
        'median_expression_norm': [0.7, 0.3]
    }
    return pd.DataFrame(data)

def test_add_expression_info(guides_df, expr_df):
    result = add_expression_info(guides_df, expr_df)
    
    # T1 -> 0.7
    # T2 -> 0.3
    # T3 -> 0.0 (default)
    
    assert result.loc[0, 'median_expression_norm'] == 0.7
    assert result.loc[1, 'median_expression_norm'] == 0.3
    assert result.loc[2, 'median_expression_norm'] == 0.0
