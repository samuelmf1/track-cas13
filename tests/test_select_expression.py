import pytest
import pandas as pd
import numpy as np
import io
from pathlib import Path
import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))

from select_expression import (
    extract_stable_id, 
    load_transcript_summary, 
    expression, 
    extract_gencode_version
)

# --- Fixtures for Mock Data ---

@pytest.fixture
def mock_mapping_df():
    """Simulates a transcript summary file (e.g., GENCODE v49)."""
    data = {
        'transcript_id': ['ENST000001.1', 'ENST000002.5', 'ENST000003.1_PAR_Y'],
        'gene_id': ['ENSG000001', 'ENSG000002', 'ENSG000003'],
        'gene_name': ['GENE1', 'GENE2', 'GENE3'],
        'biotype': ['protein_coding', 'protein_coding', 'lncRNA'],
        'transcript_support_level': ['1', '1', '1'],
        'tag': ['basic', 'basic', 'basic'],
        'canonical_ranking': [1, 1, 1]
    }
    return pd.DataFrame(data)

@pytest.fixture
def mock_omic_df():
    """Simulates expression data (e.g., DepMap/CCLE parquet)."""
    # TPM values are log2(TPM + 1)
    # ENST000001: Expressed (log2(1+1)=1)
    # ENST000002: Not expressed (0)
    data = {
        'ModelID': ['M1', 'M2', 'M3'],
        'ENST000001.5': [1.0, 1.2, 0.9], 
        'ENST000002.5': [0.0, 0.0, 0.0],
        'ENST000003.1_PAR_Y': [0.5, 0.6, 0.4]
    }
    return pd.DataFrame(data)

# --- Tests for ID Extraction & Regex ---

@pytest.mark.parametrize("input_id, expected", [
    ("ENST00000452144.6", "ENST00000452144"),
    ("ENST00000452144.6_PAR_Y", "ENST00000452144_PAR_Y"),
    ("ENST000001", "ENST000001"),
])
def test_extract_stable_id(input_id, expected):
    assert extract_stable_id(input_id) == expected

@pytest.mark.parametrize("filename, expected", [
    ("gencode.v49.metadata.tsv", "v49"),
    ("expression_v38.parquet", "v38"),
    ("no_version_file.txt", "v38"), # Default case
])
def test_extract_gencode_version(filename, expected):
    assert extract_gencode_version(filename) == expected

# --- Tests for Data Loading ---

def test_load_transcript_summary_missing_cols(tmp_path):
    # Create a TSV missing 'biotype'
    d = tmp_path / "bad_summary.tsv"
    d.write_text("gene_id\ttranscript_id\nENSG1\tENST1")
    
    with pytest.raises(ValueError, match="Missing required columns"):
        load_transcript_summary(str(d))

# --- Tests for Expression Logic ---

def test_expression_filtering(mock_mapping_df, mock_omic_df):
    """Verify thresholding and version mismatch flagging."""
    # Run for protein_coding, threshold 0.5, requires 2 cell lines
    result = expression(
        mapping_df=mock_mapping_df,
        omic_df=mock_omic_df,
        biotype='protein_coding',
        tpm_threshold=0.5,
        n_cell_lines=2
    )
    
    # ENST000001 should be is_expressed=True
    # (Values are ~1.0, which is > log2(0.5 + 1) ~ 0.58)
    enst1 = result[result['transcript_id'] == 'ENST000001.1'].iloc[0]
    assert enst1['is_expressed'] == True
    
    # Check version mismatch flag
    # Omic had ENST...01.5, Mapping had ENST...01.1
    assert enst1['version_mismatch'] == True

def test_expression_duplicate_handling(caplog):
    """Ensure duplicates in mapping/expression are handled gracefully."""
    # Add the columns the function explicitly tries to subset
    mapping_dups = pd.DataFrame({
        'transcript_id': ['ENST1.1', 'ENST1.2'],
        'gene_id': ['G1', 'G1'],
        'biotype': ['protein_coding', 'protein_coding'],
        'gene_name': ['G1', 'G1'],
        'transcript_support_level': ['1', '1'],
        'tag': ['basic', 'basic'],
        'canonical_ranking': [1, 1]
    })
    omic = pd.DataFrame({
        'ModelID': ['M1'],
        'ENST1.1': [2.0]
    })

    # Set logging level to capture the warning
    import logging
    caplog.set_level(logging.WARNING)
    
    result = expression(mapping_dups, omic, 'protein_coding', 1.0, 1)
    
    assert "Found 1 stable IDs with multiple versions" in caplog.text
    assert len(result) == 1

# --- Integration / Main Workflow Test ---

def test_main_cli_flow(mocker, tmp_path, mock_mapping_df, mock_omic_df):
    """Tests the full main() execution using mocks for file I/O."""
    
    # Create dummy files
    expr_path = tmp_path / "expr.parquet"
    mock_omic_df.to_parquet(expr_path)
    
    model_path = tmp_path / "model.csv"
    pd.DataFrame({'ModelID': ['M1', 'M2', 'M3'], 'OncotreeLineage': ['A', 'B', 'C']}).to_csv(model_path)
    
    summary_path = tmp_path / "gencode.v49.summary.tsv"
    mock_mapping_df.to_csv(summary_path, sep='\t', index=False)
    
    # Mock CLI arguments
    mocker.patch('argparse.ArgumentParser.parse_args', return_value=type('Args', (), {
        'expr': str(expr_path),
        'model': str(model_path),
        'transcript_summary': str(summary_path),
        'output_dir': str(tmp_path),
        'tpm_pc': 1.0,
        'tpm_lnc': 0.5,
        'n_cell_lines': 1,
        'abundance_threshold': 0.1
    }))
    
    from select_expression import main
    main()
    
    # Check if output file exists
    output_file = tmp_path / "gencode.v49.expr.tsv"
    assert output_file.exists()
    
    output_df = pd.read_csv(output_file, sep='\t')
    assert 'is_abundant' in output_df.columns
    assert 'version_mismatch' in output_df.columns

import pandas as pd
import numpy as np
from select_expression import relexpr, expression

@pytest.fixture
def logic_test_data():
    """
    Creates a scenario where:
    - Gene A has 2 transcripts: one highly expressed (TPM 30), one low (TPM 10).
    - Gene B has 1 transcript: expressed at TPM 20.
    """
    # Note: Input omic data is log2(TPM + 1)
    # log2(31) ≈ 4.95, log2(11) ≈ 3.46, log2(21) ≈ 4.39
    omic_df = pd.DataFrame({
        'ModelID': ['CL1', 'CL2'],
        'ENST_A_HIGH.1': [4.95, 4.95], 
        'ENST_A_LOW.1':  [3.46, 3.46],
        'ENST_B_ONLY.1': [4.39, 4.39]
    })
    
    mapping_df = pd.DataFrame({
        'transcript_id': ['ENST_A_HIGH.1', 'ENST_A_LOW.1', 'ENST_B_ONLY.1'],
        'gene_id': ['GENE_A', 'GENE_A', 'GENE_B'],
        'gene_name': ['A', 'A', 'B'],
        'biotype': ['protein_coding'] * 3,
        'transcript_support_level': ['1'] * 3,
        'tag': ['basic'] * 3,
        'canonical_ranking': [1, 2, 1]
    })
    return mapping_df, omic_df

def test_mathematical_sanity_linear(logic_test_data):
    mapping_df, omic_df = logic_test_data
    expr_df = expression(mapping_df, omic_df, 'protein_coding', 1.0, 1)
    final_df = relexpr(expr_df, omic_df, n_cell_lines=1, abundance_threshold=0.25)
    
    # Linear TPMs were 30 (High) and 10 (Low). Total = 40.
    # Expected Norm High = 30 / 40 = 0.75
    # Expected Norm Low = 10 / 40 = 0.25
    high_norm = final_df.loc[final_df['transcript_id'] == 'ENST_A_HIGH.1', 'median_expression_norm'].item()
    low_norm = final_df.loc[final_df['transcript_id'] == 'ENST_A_LOW.1', 'median_expression_norm'].item()
    
    assert np.isclose(high_norm, 0.75, atol=1e-2)
    assert np.isclose(low_norm, 0.25, atol=1e-2)

def test_mathematical_sanity(logic_test_data):
    """Verify ranking and relative abundance values."""
    mapping_df, omic_df = logic_test_data
    
    # 1. Test basic expression filtering & ranking
    expr_df = expression(
        mapping_df, omic_df, 
        biotype='protein_coding', 
        tpm_threshold=1.0, 
        n_cell_lines=1
    )
    
    # Check Ranking: HIGH should be 1, LOW should be 2 for GENE_A
    gene_a_res = expr_df[expr_df['gene_id'] == 'GENE_A'].set_index('transcript_id')
    assert gene_a_res.loc['ENST_A_HIGH.1', 'expression_rank'] == 1
    assert gene_a_res.loc['ENST_A_LOW.1', 'expression_rank'] == 2
    
    # 2. Test Relative Abundance (relexpr)
    # threshold 0.25 (meaning 25% of gene total)
    final_df = relexpr(expr_df, omic_df, n_cell_lines=1, abundance_threshold=0.25)
    
    # Verification of GENE_A proportions
    # Total median expression = 4.95 + 3.46 = 8.41
    # Expected Norm High = 4.95 / 8.41 ≈ 0.588
    # Expected Norm Low = 3.46 / 8.41 ≈ 0.411
    high_norm = final_df.loc[final_df['transcript_id'] == 'ENST_A_HIGH.1', 'median_expression_norm'].values[0]
    low_norm = final_df.loc[final_df['transcript_id'] == 'ENST_A_LOW.1', 'median_expression_norm'].values[0]
    
    assert np.isclose(high_norm + low_norm, 1.0)
    assert high_norm > low_norm
    
    # Check Abundance Flags (m10 = 10% threshold)
    # Both are > 25%, so both should be True for m25
    # .item() is best for getting a single scalar out of a single-element Series
    assert final_df.loc[final_df['transcript_id'] == 'ENST_A_LOW.1', 'is_m25_abundant'].item() is True

def test_non_expressed_behavior(logic_test_data):
    """Ensure transcripts below threshold don't get ranked but stay in output."""
    mapping_df, omic_df = logic_test_data
    
    # Set threshold very high so everything is 'not expressed'
    expr_df = expression(
        mapping_df, omic_df, 
        biotype='protein_coding', 
        tpm_threshold=100.0, 
        n_cell_lines=1
    )
    
    # All ranks should be 0, is_expressed should be False
    assert (expr_df['expression_rank'] == 0).all()
    assert (expr_df['is_expressed'] == False).all()
    
    # Run through relexpr
    final_df = relexpr(expr_df, omic_df, n_cell_lines=1, abundance_threshold=0.1)
    
    # Abundance boolean flags should be False for non-expressed items
    assert (final_df['is_m10_abundant'] == False).all()


if __name__ == "__main__":
    import pytest
    import sys
    # This ensures that running the script directly invokes pytest
    sys.exit(pytest.main([__file__]))