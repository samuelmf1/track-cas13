import pytest
import pandas as pd
import io
import os
import sys

# Ensure the src directory is in the path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))

from reduce_annotate import main as run_annotate

@pytest.fixture
def mock_inputs(tmp_path):
    """Creates dummy files for testing the annotation logic."""
    # 1. Cleared Guides
    guides_csv = tmp_path / "test_guides.csv"
    pd.DataFrame({
        'Contig_Idx': [10, 100],
        'Target Sequence': ['ATGC', 'GGGG'], 
        'Symbol': ['T1', 'T1']
    }).to_csv(guides_csv, index=False)

    # 2. Exons
    exons_tsv = tmp_path / "exons.tsv"
    pd.DataFrame({
        'transcript_id': ['T1'],
        'exon_id': ['E1'],
        'exon_sequence': ['AAAAGGGGCCCC']
    }).to_csv(exons_tsv, sep='\t', index=False)

    # 3. Gencode Summary
    summary_tsv = tmp_path / "summary.tsv"
    pd.DataFrame({
        'transcript_id': ['T1'],
        'biotype': ['protein_coding'],
        'is_representative': [True],
        'tag': ['basic;appris_p1']
    }).to_csv(summary_tsv, sep='\t', index=False)

    # 4. CDS Boundaries
    cds_tsv = tmp_path / "cds.tsv"
    pd.DataFrame({
        'transcript_id': ['T1'],
        'region': ['CDS'],
        'start': [50],
        'end': [150]
    }).to_csv(cds_tsv, sep='\t', index=False)

    # 5. Expression
    expr_tsv = tmp_path / "expr.tsv"
    pd.DataFrame({
        'transcript_id': ['T1'],
        'sum_expression_norm': [10.5],
        'is_expressed': [True]
    }).to_csv(expr_tsv, sep='\t', index=False)

    return {
        'csv': str(guides_csv),
        'exons': str(exons_tsv),
        'cds': str(cds_tsv),
        'summary': str(summary_tsv),
        'expr': str(expr_tsv),
        'out_prefix': str(tmp_path / "test_output")
    }

def test_annotation_logic(mock_inputs, monkeypatch):
    """Verifies that all three output stages are produced and accurate."""
    args = [
        "reduce_annotate.py",
        "--csv", mock_inputs['csv'],
        "--exons", mock_inputs['exons'],
        "--cds", mock_inputs['cds'],
        "--summary", mock_inputs['summary'],
        "--expr", mock_inputs['expr'],
        "--out_prefix", mock_inputs['out_prefix']
    ]
    monkeypatch.setattr("sys.argv", args)
    run_annotate()

    final_file = f"{mock_inputs['out_prefix']}.annots.csv"
    assert os.path.exists(final_file)
    df = pd.read_csv(final_file)
    
    assert df.iloc[0]['exon_id'] == 'E1'
    assert df.iloc[0]['region'] == 'CDS'
    assert df.iloc[0]['sum_expression_norm'] == 10.5
    assert "basic|appris_p1" in df.iloc[0]['tags']

def test_inadmissible_filtering(mock_inputs, monkeypatch):
    """Checks that guides with no exon match are filtered out of the a2/annots files."""
    args = [
        "reduce_annotate.py", "--csv", mock_inputs['csv'], "--exons", mock_inputs['exons'],
        "--cds", mock_inputs['cds'], "--summary", mock_inputs['summary'],
        "--expr", mock_inputs['expr'], "--out_prefix", mock_inputs['out_prefix']
    ]
    monkeypatch.setattr("sys.argv", args)
    run_annotate()

    a1 = pd.read_csv(f"{mock_inputs['out_prefix']}.a1.csv")
    annots = pd.read_csv(f"{mock_inputs['out_prefix']}.annots.csv")

    assert len(a1) == 2
    assert len(annots) == 1

def test_lncrna_and_representative_logic(mock_inputs, monkeypatch):
    """Verifies biotype override for lncRNA."""
    guides_df = pd.DataFrame({
        'Contig_Idx': [10, 10], 
        'Target Sequence': ['GGGG', 'GGGG'],
        'Symbol': ['T_LNC', 'T_NONREP']
    })
    guides_df.to_csv(mock_inputs['csv'], index=False)

    pd.DataFrame({
        'transcript_id': ['T_LNC', 'T_NONREP'],
        'biotype': ['lncRNA', 'protein_coding'],
        'is_representative': [False, False],
        'tag': ['basic', 'basic']
    }).to_csv(mock_inputs['summary'], sep='\t', index=False)

    pd.DataFrame({
        'transcript_id': ['T_LNC', 'T_NONREP'],
        'exon_id': ['E_LNC', 'E_NON'],
        'exon_sequence': ['AAAAGGGGCCCC', 'AAAAGGGGCCCC']
    }).to_csv(mock_inputs['exons'], sep='\t', index=False)

    monkeypatch.setattr("sys.argv", [
        "reduce_annotate.py", "--csv", mock_inputs['csv'], "--exons", mock_inputs['exons'],
        "--cds", mock_inputs['cds'], "--summary", mock_inputs['summary'],
        "--expr", mock_inputs['expr'], "--out_prefix", mock_inputs['out_prefix']
    ])
    
    run_annotate()
    df = pd.read_csv(f"{mock_inputs['out_prefix']}.annots.csv")

    assert df[df['Symbol'] == 'T_LNC']['region'].iloc[0] == 'lncRNA'
    assert df[df['Symbol'] == 'T_NONREP']['region'].iloc[0] != 'lncRNA'

def test_utr_boundary_overlap(mock_inputs, monkeypatch):
    """Validates the 3nt overlap rule with representative +3 offset."""
    # T1 is representative (+3 offset). CDS is [50, 150)
    # 1. Idx 24 -> Eff 27. End 50. Overlap 0. -> 5'UTR
    # 2. Idx 26 -> Eff 29. End 52. Overlap 2 (50, 51). -> 5'UTR (< 3)
    # 3. Idx 27 -> Eff 30. End 53. Overlap 3 (50, 51, 52). -> CDS (>= 3)
    
    guides_df = pd.DataFrame({
        'Contig_Idx': [24, 26, 27], 
        'Target Sequence': ['GGGG', 'GGGG', 'GGGG'],
        'Symbol': ['T1', 'T1', 'T1']
    })
    guides_df.to_csv(mock_inputs['csv'], index=False)

    monkeypatch.setattr("sys.argv", [
        "reduce_annotate.py", "--csv", mock_inputs['csv'], "--exons", mock_inputs['exons'],
        "--cds", mock_inputs['cds'], "--summary", mock_inputs['summary'],
        "--expr", mock_inputs['expr'], "--out_prefix", mock_inputs['out_prefix']
    ])
    
    run_annotate()
    df = pd.read_csv(f"{mock_inputs['out_prefix']}.annots.csv").sort_values('Contig_Idx')

    assert df.iloc[0]['region'] == "5'UTR"
    assert df.iloc[1]['region'] == "5'UTR"
    assert df.iloc[2]['region'] == "CDS"

def test_missing_metadata_handling(mock_inputs, monkeypatch):
    """Ensures resilience when transcript IDs are missing from metadata."""
    guides_df = pd.DataFrame({
        'Contig_Idx': [10], 
        'Target Sequence': ['GGGG'],
        'Symbol': ['T_MISSING']
    })
    guides_df.to_csv(mock_inputs['csv'], index=False)

    pd.DataFrame({
        'transcript_id': ['T_MISSING'],
        'exon_id': ['E_MISSING'],
        'exon_sequence': ['GGGGGGGGGGGG']
    }).to_csv(mock_inputs['exons'], sep='\t', index=False)

    monkeypatch.setattr("sys.argv", [
        "reduce_annotate.py", "--csv", mock_inputs['csv'], "--exons", mock_inputs['exons'],
        "--cds", mock_inputs['cds'], "--summary", mock_inputs['summary'],
        "--expr", mock_inputs['expr'], "--out_prefix", mock_inputs['out_prefix']
    ])
    
    run_annotate()
    df = pd.read_csv(f"{mock_inputs['out_prefix']}.annots.csv")

    assert df.iloc[0]['region'] == "Unknown"
    assert df.iloc[0]['sum_expression_norm'] == -1