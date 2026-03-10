
import os
import shutil
import pytest
from Bio import SeqIO
from src.ontarget_chunk_fasta import split_fasta

# Sample FASTA content
SAMPLE_FASTA = """
>1|ENSG000001|ENST000001|GeneA|lncRNA
ACGT
>2|ENSG000001|ENST000002|GeneA|lncRNA
CGTA
>3|ENSG000001|ENST000003|GeneA|lncRNA
GTAC
>4|ENSG000002|ENST000004|GeneB|lncRNA
TACG
>5|ENSG000002|ENST000005|GeneB|lncRNA
ACGT
>6|ENSG000003|ENST000006|GeneC|lncRNA
AAAA
>7|ENSG000004|ENST000007|GeneD|lncRNA
GGGG
>8|ENSG000004|ENST000008|GeneD|lncRNA
CCCC
>9|ENSG000005|ENST000009|GeneE|lncRNA
TTTT
>10|ENSG000005|ENST000010|GeneE|lncRNA
AAAA
"""

@pytest.fixture
def run_dir(tmp_path):
    """Setup a temporary directory for the test run"""
    d = tmp_path / "test_run"
    d.mkdir()
    return d

def test_split_fasta_grouping(run_dir, monkeypatch):
    """Test that split_fasta groups same gene_ids into the same chunk"""
    
    # Run in the temporary directory so output files are created there
    monkeypatch.chdir(run_dir)
    
    input_fasta = run_dir / "test.fasta"
    input_fasta.write_text(SAMPLE_FASTA.strip())
    
    # We have 5 genes (A, B, C, D, E) and 10 records.
    # GeneA: 3 records
    # GeneB: 2 records
    # GeneC: 1 record
    # GeneD: 2 records
    # GeneE: 2 records
    
    # Try splitting into 3 chunks.
    # Total records = 10. Target size = ceil(10/3) = 4.
    
    # Expected behavior:
    # Chunk 1: GeneA (3) + maybe GeneB(2)? No, 3+2=5 > 4. 
    # But strictly speaking, my greedy logic was: 
    # if current >= target: split.
    # So:
    # 1. GeneA added. Current size = 3. 3 < 4.
    # 2. GeneB added. Current size = 3+2 = 5.
    # 3. Next iteration: 5 >= 4, creating new chunk for NEXT items.
    # Wait, let's look at logic again.
    # for gene_id in ordered_gene_ids:
    #    check if current >= target. If so, flush current and start new.
    #    Add gene to current.
    
    # Trace:
    # Target = 4.
    # 1. GeneA (3). Current=0. Valid. Add A. Current=3.
    # 2. GeneB (2). Current=3. 3 < 4. Valid. Add B. Current=5.
    # 3. GeneC (1). Current=5. 5 >= 4. FLUSH Chunk 1 (A, B). Start new Chunk 2. Add C. Current=1.
    # 4. GeneD (2). Current=1. 1 < 4. Add D. Current=3.
    # 5. GeneE (2). Current=3. 3 < 4. Add E. Current=5.
    # 6. End loop. Flush Chunk 2 (C, D, E).
    
    # So we expect 2 chunks?
    # Chunk 1: A, B (5 records)
    # Chunk 2: C, D, E (5 records)
    
    # But wait, n_chunks=3.
    # Loop Trace again:
    # ...
    # 3. GeneC (1). Current=5. 5>=4 and len(chunks) < 2. FLUSH 1. Chunk 1 = [A, B]. Chunks=[[A,B]]. Current=[] -> Add C. Current=[C]. Size=1.
    # 4. GeneD (2). Current=1 < 4. Add D. Size=3.
    # 5. GeneE (2). Current=3 < 4. Add E. Size=5.
    # End. Flush Current. Chunks=[[A,B], [C,D,E]].
    # Only 2 chunks created even though 3 requested. This is acceptable for this logic as long as genes are grouped.
    
    # Let's try n_chunks=4. Target = ceil(10/4) = 3.
    # 1. GeneA (3). Add A. Size=3.
    # 2. GeneB (2). Size=3 >= 3. FLUSH 1 [A]. Chunks=[[A]]. New Current. Add B. Size=2.
    # 3. GeneC (1). Size=2 < 3. Add C. Size=3.
    # 4. GeneD (2). Size=3 >= 3. FLUSH 2 [B, C]. Chunks=[[A], [B,C]]. New Current. Add D. Size=2.
    # 5. GeneE (2). Size=2 < 3. Add E. Size=4.
    # End. Flush [D, E]. Chunks=[[A], [B,C], [D,E]].
    # 3 chunks created.
    
    # Let's run with n_chunks=4.
    split_fasta(str(input_fasta), start_seq_id=1, n_chunks=4)
    
    # Verify outputs
    chunk_files = sorted(list(run_dir.glob("test.chunk_*.fasta")))
    
    # Check that we have created some chunks
    assert len(chunk_files) > 0
    
    seen_genes = set()
    
    for cf in chunk_files:
        records = list(SeqIO.parse(cf, "fasta"))
        chunk_genes = set()
        for r in records:
            parts = r.id.split("|")
            gene_id = parts[1]
            chunk_genes.add(gene_id)
            
            # Check if this gene was seen in a previous chunk
            assert gene_id not in seen_genes, f"Gene {gene_id} split across chunks!"
        
        # Add to seen
        seen_genes.update(chunk_genes)

    # Also verify total records
    total_output_records = 0
    for cf in chunk_files:
        total_output_records += len(list(SeqIO.parse(cf, "fasta")))
    assert total_output_records == 10
