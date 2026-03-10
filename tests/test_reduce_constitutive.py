import pytest
import pandas as pd
import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))
from reduce_constitutive import SelectionEngine, SelectionConfig, GuideSelector

class TestGuideSelection:

    @pytest.fixture
    def selector(self):
        return GuideSelector(target_n=2)

    @pytest.fixture
    def config(self):
        return SelectionConfig(target_n=2)

    # --- RANK WINDOW TESTS ---

    def test_rank_window_preference_logic(self, selector, config):
        """Verify that a higher expression group wins if it is within 4 ranks."""
        df = pd.DataFrame({
            'position': [100, 115, 300, 315], # Both fail 23bp (R0), hit 12bp (R1)
            'tiger_score': [0.9, 0.9, 0.9, 0.9],
            'transcript_id_group': ['Low_Exp', 'Low_Exp', 'High_Exp', 'High_Exp'],
            'region': ['CDS']*4,
            'ttm_gene_norm': [0.1, 0.1, 1.0, 1.0],
            'expression_content': [0.1, 0.1, 0.1, 0.1],
            'sum_log_median_expr_norm': [0.1, 0.1, 1.0, 1.0],
            'Guide Sequence': ['S1', 'S2', 'S3', 'S4'],
            'biotype': ['protein_coding']*4
        })
        res_df, summary = selector.process_gene('G1', df, config)
        assert summary['transcript_id_group_targeted'] == 'High_Exp'

    def test_rank_window_cutoff(self, selector, config):
        """Verify that a high expression group is IGNORED if it is >4 ranks away."""
        df = pd.DataFrame({
            'position': [100, 200, 300, 400],
            'tiger_score': [0.9, 0.9, 0.3, 0.3], # Q1 (Rank 0-11) vs Q3 (Rank 24+)
            'transcript_id_group': ['Q1_Low', 'Q1_Low', 'Q3_High', 'Q3_High'],
            'region': ['CDS']*4,
            'ttm_gene_norm': [0.1, 0.1, 1.0, 1.0],
            'expression_content': [0.1, 0.1, 0.1, 0.1],
            'sum_log_median_expr_norm': [0.1, 0.1, 1.0, 1.0],
            'Guide Sequence': ['A1', 'A2', 'B1', 'B2'],
            'biotype': ['protein_coding']*4
        })
        res_df, summary = selector.process_gene('G2', df, config)
        assert summary['transcript_id_group_targeted'] == 'Q1_Low'

    # --- REGION TESTS ---

    def test_utr_limit_enforced(self, selector, config):
        """Verify that a 'UTR' strategy respects the 2-guide limit for UTR regions."""
        df = pd.DataFrame({
            'position': [100, 200, 300, 400],
            'tiger_score': [0.9, 0.85, 0.8, 0.75],
            'transcript_id_group': ['T1']*4,
            'region': ['3UTR', '5UTR', '3UTR', 'CDS'],
            'ttm_gene_norm': [1.0]*4,
            'expression_content': [0.5]*4,
            'sum_log_median_expr_norm': [1.0]*4,
            'Guide Sequence': ['U1', 'U2', 'U3', 'C1']
        })
        # target_n=3. 
        # The selector iterates UTR 0..5.
        # When it hits max_utr=2, it should pick U1, U2, and then fail to pick U3 (limit 2), forcing it to pick C1 instead (assuming C1 score is ok? C1 score is 0.75, which is < 0.8 min score for Q1/Q2 strategies...)
        # Wait, the strategies loop through score cutoffs [0.8, 0.6] then [0.4, 0.2].
        # Guides: U1(0.9), U2(0.85), U3(0.8), C1(0.75).
        # Strategy Q1 (0.8): Hits U1, U2, U3. 
        # If max_utr=0: Selects only non-UTR? No non-UTR >= 0.8. Selects []. Fails.
        # If max_utr=1: Selects U1. Fails (target_n=3).
        # If max_utr=2: Selects U1, U2. Fails (target_n=3).
        # If max_utr=3: Selects U1, U2, U3. Success.
        
        # But we want to test that the LIMIT works.
        # Let's verify that with max_utr=2 passed manually to 'select_guides', it respects it.
        
        # NOTE: The full integration 'process_gene' will now FIND the solution with 3 UTRs if needed.
        # So we should test SelectionEngine.select_guides directly to verify the limiter logic.
        
        config.target_n = 3
        res_df, stats = SelectionEngine.select_guides(df, config, 23, 0.6, max_utr=2)

        # With min_score=0.6, Candidates: U1, U2, U3, C1.
        # max_utr=2. 
        # Pick U1 (0.9, UTR) -> OK. utr=1
        # Pick U2 (0.85, UTR) -> OK. utr=2
        # Pick U3 (0.8, UTR) -> Skip (utr >= max).
        # Pick C1 (0.75, CDS) -> OK.
        # Result: U1, U2, C1.
        
        assert len(res_df) == 3
        regions = res_df['region'].tolist()
        utr_count = sum(1 for r in regions if 'UTR' in r)
        assert utr_count == 2
        assert 'C1' in res_df['Guide Sequence'].values

    # --- ENGINE DISTANCE TESTS ---

    def test_engine_respects_distance(self):
        """Verify SelectionEngine rejects guides too close together."""
        df = pd.DataFrame({
            'position': [100, 110, 200], # 100 and 110 are 10bp apart
            'tiger_score': [0.9, 0.85, 0.8],
            'position_col': [100, 110, 200]
        })
        config = SelectionConfig(target_n=2, pos_col='position_col')
        # Min dist 23: Should pick 100 and 200 (skipping 110)
        res_df, stats = SelectionEngine.select_guides(df, config, 23, 0.0, max_utr=5)
        assert list(res_df['position_col']) == [100, 200]

    # --- TIE BREAKER TESTS ---

    def test_transcript_breadth_tie_breaker(self, selector, config):
        """If expression is identical, the group with more transcripts wins."""
        df = pd.DataFrame({
            'position': [100, 200, 300, 400],
            'tiger_score': [0.9]*4,
            'transcript_id_group': ['T1', 'T1', 'T2|T3', 'T2|T3'],
            'region': ['CDS']*4,
            'ttm_gene_norm': [1.0]*4,
            'expression_content': [0.5]*4,
            'sum_log_median_expr_norm': [1.0]*4,
            'Guide Sequence': ['S1', 'S2', 'S3', 'S4'],
            'biotype': ['protein_coding']*4
        })
        res_df, summary = selector.process_gene('G_TIE', df, config)
        assert summary['transcript_id_group_targeted'] == 'T2|T3'
        assert summary['n_targeted_tx'] == 2

    def test_summary_includes_tags(self, selector, config):
        """Verify that tags are included in the summary info."""
        df = pd.DataFrame({
            'position': [100, 200],
            'tiger_score': [0.9, 0.9],
            'transcript_id_group': ['T1', 'T1'],
            'region': ['CDS', 'CDS'],
            'ttm_gene_norm': [1.0, 1.0],
            'expression_content': [0.5, 0.5],
            'sum_log_median_expr_norm': [1.0, 1.0],
            'Guide Sequence': ['S1', 'S2'],
            'tags': ['tag1|tag2', 'tag1|tag2']
        })
        res_df, summary = selector.process_gene('G_TAGS', df, config)
        assert summary['tags'] == 'tag1|tag2'

    def test_summary_tags_none_if_missing(self, selector, config):
        """Verify that tags are "None" if missing from input."""
        df = pd.DataFrame({
            'position': [100, 200],
            'tiger_score': [0.9, 0.9],
            'transcript_id_group': ['T1', 'T1'],
            'region': ['CDS', 'CDS'],
            'ttm_gene_norm': [1.0, 1.0],
            'expression_content': [0.5, 0.5],
            'sum_log_median_expr_norm': [1.0, 1.0],
            'Guide Sequence': ['S1', 'S2']
            # 'tags' column missing
        })
        res_df, summary = selector.process_gene('G_NO_TAGS', df, config)
        assert summary['tags'] == 'None'

    def test_failure_summary_columns(self, selector, config):
        """Verify that failures return a complete summary_info dictionary."""
        df = pd.DataFrame({
            'position': [100],
            'tiger_score': [0.1], # Low score, won't be selected if target_n=2 or min_score higher
            'transcript_id_group': ['T1'],
            'region': ['CDS'],
            'ttm_gene_norm': [1.0],
            'expression_content': [0.5],
            'sum_log_median_expr_norm': [1.0],
            'Guide Sequence': ['S1'],
            'biotype': ['protein_coding'],
            'tags': ['None']
        })
        # target_n=2 will cause failure because only 1 guide available
        res_df, summary = selector.process_gene('G_FAIL', df, config)
        assert res_df.empty
        assert summary['gene_id'] == 'G_FAIL'
        assert summary['biotype'] == 'protein_coding'
        assert summary['result_class'] == '999_failed'
        assert summary['n_targeted_tx'] == 0
        assert summary['transcript_ids_not_targeted'] == 'T1'
        assert summary['tags'] == 'None'

    # --- CANONICAL FALLBACK & CONDITION PRIORITY TESTS ---

    def test_canonical_fallback_fragmentation_fix(self, selector, config):
        """Verify that guides targetting CAN are aggregated even if they also target LUTR."""
        # target_n=2
        # T1 = Canonical, T2 = LUTR
        df = pd.DataFrame({
            'position': [100, 200], 
            'tiger_score': [0.9, 0.9],
            'transcript_id_group': ['T1|T2', 'T1'], # Fragmented: T1|T2 vs T1
            'region': ['CDS', 'CDS'],
            'sum_log_median_expr_norm': [1.0, 1.0],
            'expression_content': [0.1, 0.1], # Low expression
            'ttm_gene_norm': [0.5, 0.5],
            'Guide Sequence': ['S1', 'S2'],
            'tags': ['Ensembl_canonical', 'Ensembl_canonical']
        })
        # Without fallback, T1|T2 has 1 guide, T1 has 1 guide. Both < target_n=2. Fails.
        # With fallback, f"FALLBACK_T1" will have 2 guides.
        res_df, summary = selector.process_gene('G_FALLBACK', df, config, canonical_tx='T1')
        assert not res_df.empty
        assert summary['transcript_id_group_targeted'] == 'Canonical_Fallback'
        assert len(res_df) == 2

    def test_expressed_non_canonical_preferred(self, selector, config):
        """Verify that high-expression non-canonical (L2) beats low-expression non-canonical (L0)
        when there is no canonical group available."""
        # target_n=2
        # T1 = Non-Canonical (Low Exp), T2 = Non-Canonical (High Exp)
        # canonical_tx='T_MISSING' so neither group overlaps canonical
        df = pd.DataFrame({
            'position': [100, 150, 200, 250], 
            'tiger_score': [0.9, 0.9, 0.9, 0.9],
            'transcript_id_group': ['T1', 'T1', 'T2', 'T2'],
            'region': ['CDS']*4,
            'sum_log_median_expr_norm': [0.1, 0.1, 1.0, 1.0],
            'expression_content': [0.1, 0.1, 0.5, 0.5], # T2 group max expression = 0.5 (>0.4)
            'ttm_gene_norm': [0.2, 0.2, 1.0, 1.0], # T2 has much better priority
            'Guide Sequence': ['S1', 'S2', 'S3', 'S4'],
            'tags': ['None', 'None', 'None', 'None']
        })
        res_df, summary = selector.process_gene('G_EXPR', df, config, canonical_tx='T_MISSING')
        assert summary['transcript_id_group_targeted'] == 'T2'

    def test_low_expression_canonical_priority(self, selector, config):
        """Verify that if expression is low (<0.4), canonical wins even if lower priority."""
        # target_n=2
        # T1 = Canonical (Lower Priority), T2 = Non-Canonical (Higher Priority)
        df = pd.DataFrame({
            'position': [100, 150, 200, 250], 
            'tiger_score': [0.9, 0.9, 0.9, 0.9],
            'transcript_id_group': ['T1', 'T1', 'T2', 'T2'],
            'region': ['CDS']*4,
            'sum_log_median_expr_norm': [0.1, 0.1, 0.2, 0.2],
            'expression_content': [0.1, 0.1, 0.3, 0.3], # Max expression 0.3 (<0.4)
            'ttm_gene_norm': [0.5, 0.5, 1.0, 1.0], # T2 has higher priority
            'Guide Sequence': ['S1', 'S2', 'S3', 'S4'],
            'tags': ['Ensembl_canonical', 'Ensembl_canonical', 'None', 'None']
        })
        res_df, summary = selector.process_gene('G_LOW', df, config, canonical_tx='T1')
        assert summary['transcript_id_group_targeted'] == 'T1'

    def test_worthiness_hierarchy(self, selector):
        """Verify the Worthiness Level hierarchy (0-3)."""
        config = SelectionConfig(target_n=2, lookahead=4)
        
        # Scenario:
        # We need a dataframe with 3 potential groups:
        # 1. Real Canonical (Level 3) -> should win if present
        # 2. Real High-Expression Non-Canonical (Level 2)
        # 3. Canonical Fallback (Level 1)
        # 4. Real Low-Expression Non-Canonical (Level 0)
        
        df = pd.DataFrame({
            'position': [100, 110, 300, 310, 500, 510, 700, 710], 
            'tiger_score': [0.9]*8,
            'transcript_id_group': ['T_CAN', 'T_CAN', 'T_HIGH_NON', 'T_HIGH_NON', 'T_OTHER', 'T_OTHER', 'T_LOW_NON', 'T_LOW_NON'], 
            'region': ['CDS']*8,
            'sum_log_median_expr_norm': [1.0]*8,
            'expression_content': [0.1, 0.1, 0.5, 0.5, 0.1, 0.1, 0.1, 0.1], # T_HIGH_NON is 0.5, others low
            'ttm_gene_norm': [1.0]*8,
            'Guide Sequence': ['A1', 'A2', 'B1', 'B2', 'C1', 'C2', 'D1', 'D2'],
            'biotype': ['protein_coding']*8,
            'tags': ['Ensembl_canonical', 'Ensembl_canonical', 'None', 'None', 'None', 'None', 'None', 'None']
        })
        
        # Test 1: Real Canonical (L3) vs Real Low-Exp Non-Canonical (L0) -> Real Canonical wins
        res_df, summary = selector.process_gene('G1', df[df['Guide Sequence'].isin(['A1', 'A2', 'C1', 'C2'])], config, canonical_tx='T_CAN')
        assert summary['transcript_id_group_targeted'] == 'T_CAN'
        
        # Test 2: Real High-Exp Non-Canonical (L2) vs Fallback (L1) -> Real High-Exp wins
        # Fragments are split across two DIFFERENT groups so no single real group
        # can reach target_n=2, but the FALLBACK aggregation (str.contains T_CAN) can.
        subset2 = df[df['Guide Sequence'].isin(['B1', 'B2'])].copy()
        frag = pd.DataFrame({
            'position': [500, 510], 'tiger_score': [0.9]*2,
            'transcript_id_group': ['T_CAN|T_FRAG', 'T_CAN'],  # Two different groups, 1 guide each
            'region': ['CDS']*2, 'sum_log_median_expr_norm': [1.0]*2,
            'expression_content': [0.1]*2, 'ttm_gene_norm': [1.0]*2,
            'Guide Sequence': ['E1', 'E2'], 'biotype': ['protein_coding']*2,
            'tags': ['None', 'None']
        })
        res_df, summary = selector.process_gene('G2', pd.concat([subset2, frag]), config, canonical_tx='T_CAN')
        assert summary['transcript_id_group_targeted'] == 'T_HIGH_NON'
        
        # Test 3: Fallback (L1) vs Real Low-Exp Non-Canonical (L0) -> Fallback wins
        subset3 = df[df['Guide Sequence'].isin(['D1', 'D2'])].copy()
        res_df, summary = selector.process_gene('G3', pd.concat([frag, subset3]), config, canonical_tx='T_CAN')
        assert summary['transcript_id_group_targeted'] == 'Canonical_Fallback'

    # --- REGION FALLBACK TESTS ---

    def test_intron_only_guides_succeed_via_region_fallback(self, selector, config):
        """Verify that genes with only intron-region guides can still succeed."""
        df = pd.DataFrame({
            'position': [100, 200],
            'tiger_score': [0.9, 0.9],
            'transcript_id_group': ['T1', 'T1'],
            'region': ['intron', 'intron'],  # No CDS/UTR/lncRNA
            'ttm_gene_norm': [1.0, 1.0],
            'expression_content': [0.5, 0.5],
            'sum_log_median_expr_norm': [1.0, 1.0],
            'Guide Sequence': ['S1', 'S2'],
            'biotype': ['protein_coding', 'protein_coding'],
            'tags': ['Ensembl_canonical', 'Ensembl_canonical']
        })
        res_df, summary = selector.process_gene('G_INTRON', df, config)
        assert not res_df.empty, "Gene with only intron guides should succeed via region fallback"
        assert len(res_df) == 2

    def test_nan_region_guides_succeed_via_region_fallback(self, selector, config):
        """Verify that guides with NaN region can still be selected via fallback."""
        df = pd.DataFrame({
            'position': [100, 200],
            'tiger_score': [0.9, 0.9],
            'transcript_id_group': ['T1', 'T1'],
            'region': [float('nan'), float('nan')],
            'ttm_gene_norm': [1.0, 1.0],
            'expression_content': [0.5, 0.5],
            'sum_log_median_expr_norm': [1.0, 1.0],
            'Guide Sequence': ['S1', 'S2'],
            'biotype': ['protein_coding', 'protein_coding'],
            'tags': ['Ensembl_canonical', 'Ensembl_canonical']
        })
        res_df, summary = selector.process_gene('G_NAN_REGION', df, config)
        assert not res_df.empty, "Gene with NaN region should succeed via region fallback"

    # --- FALLBACK TAG FIX TESTS ---

    def test_fallback_group_gets_canonical_tags(self, selector, config):
        """Verify FALLBACK group winner reports canonical-like tags, not inherited non-canonical tags."""
        df = pd.DataFrame({
            'position': [100, 200],
            'tiger_score': [0.9, 0.9],
            'transcript_id_group': ['T1|T2', 'T1'],  # Fragmented canonical
            'region': ['CDS', 'CDS'],
            'sum_log_median_expr_norm': [1.0, 1.0],
            'expression_content': [0.1, 0.1],
            'ttm_gene_norm': [0.5, 0.5],
            'Guide Sequence': ['S1', 'S2'],
            'tags': ['non_canonical_tag', 'non_canonical_tag']  # No canonical tags
        })
        res_df, summary = selector.process_gene('G_FBTAG', df, config, canonical_tx='T1')
        assert not res_df.empty
        # The fallback group should win AND have canonical-like tags
        assert 'canonical' in summary['tags'].lower(), \
            f"FALLBACK winner should have canonical tags, got: {summary['tags']}"

    # --- OVERLAPS_CANONICAL TESTS ---

    def test_overlaps_canonical_in_summary(self, selector, config):
        """Verify overlaps_canonical is True when canonical transcript is targeted."""
        df = pd.DataFrame({
            'position': [100, 200],
            'tiger_score': [0.9, 0.9],
            'transcript_id_group': ['T1', 'T1'],
            'region': ['CDS', 'CDS'],
            'ttm_gene_norm': [1.0, 1.0],
            'expression_content': [0.5, 0.5],
            'sum_log_median_expr_norm': [1.0, 1.0],
            'Guide Sequence': ['S1', 'S2'],
            'tags': ['Ensembl_canonical', 'Ensembl_canonical']
        })
        _, summary = selector.process_gene('G_OC', df, config, canonical_tx='T1')
        assert summary['overlaps_canonical'] == True

    def test_overlaps_canonical_false_when_not_targeting(self, selector, config):
        """Verify overlaps_canonical is False when winner doesn't target canonical."""
        df = pd.DataFrame({
            'position': [100, 200],
            'tiger_score': [0.9, 0.9],
            'transcript_id_group': ['T2', 'T2'],
            'region': ['CDS', 'CDS'],
            'ttm_gene_norm': [1.0, 1.0],
            'expression_content': [0.5, 0.5],
            'sum_log_median_expr_norm': [1.0, 1.0],
            'Guide Sequence': ['S1', 'S2'],
            'tags': ['None', 'None']
        })
        _, summary = selector.process_gene('G_OC2', df, config, canonical_tx='T1')
        assert summary['overlaps_canonical'] == False

    def test_failed_gene_has_overlaps_canonical_false(self, selector, config):
        """Verify failed genes have overlaps_canonical=False."""
        df = pd.DataFrame({
            'position': [100],
            'tiger_score': [0.1],
            'transcript_id_group': ['T1'],
            'region': ['CDS'],
            'ttm_gene_norm': [1.0],
            'expression_content': [0.5],
            'sum_log_median_expr_norm': [1.0],
            'Guide Sequence': ['S1'],
            'tags': ['Ensembl_canonical']
        })
        _, summary = selector.process_gene('G_FAIL2', df, config)
        assert summary['overlaps_canonical'] == False

    # --- DIST=0 OVERLAP RECOVERY TEST ---

    def test_dist_zero_recovers_clustered_guides(self, selector, config):
        """Verify that heavily clustered guides can succeed with dist=0 fallback."""
        # All guides within 2bp of each other — fails dist=3 and dist=6 and dist=23
        df = pd.DataFrame({
            'position': [100, 101],
            'tiger_score': [0.9, 0.85],
            'transcript_id_group': ['T1', 'T1'],
            'region': ['CDS', 'CDS'],
            'ttm_gene_norm': [1.0, 1.0],
            'expression_content': [0.5, 0.5],
            'sum_log_median_expr_norm': [1.0, 1.0],
            'Guide Sequence': ['S1', 'S2'],
            'tags': ['Ensembl_canonical', 'Ensembl_canonical']
        })
        res_df, summary = selector.process_gene('G_CLUSTER', df, config)
        assert not res_df.empty, "Clustered guides should succeed via dist=0 fallback"
        assert len(res_df) == 2

if __name__ == "__main__":
    pytest.main([__file__])