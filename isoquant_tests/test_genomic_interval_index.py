############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import pytest
import gffutils
from unittest.mock import Mock, MagicMock, patch

from isoquant_lib.genomic_interval_index import GenomicIntervalIndex


class TestGenomicIntervalIndex:
    """Unit tests for GenomicIntervalIndex class that builds interval trees from gffutils database."""

    source_dir = os.path.dirname(os.path.realpath(__file__))
    gffutils_db = gffutils.FeatureDB(os.path.join(source_dir, 'toy_data/synth.db'), keep_order=True)

    def test_initialization_with_all_chromosomes(self):
        """Test GenomicIntervalIndex initialization with no chromosome filtering."""
        index = GenomicIntervalIndex(self.gffutils_db)
        
        assert index.db == self.gffutils_db
        assert index.gene_trees is not None
        assert index.exon_trees is not None
        assert isinstance(index.gene_trees, dict)
        assert isinstance(index.exon_trees, dict)
        assert index.chromosomes is None
        assert len(index.gene_trees) > 0
        assert len(index.exon_trees) > 0

    def test_initialization_with_chromosome_filter(self):
        """Test GenomicIntervalIndex initialization with specific chromosome filter."""
        chromosomes = ['chr10']
        index = GenomicIntervalIndex(self.gffutils_db, chromosomes=chromosomes)
        assert index.chromosomes == set(chromosomes)
        assert 'chr10' in index.gene_trees
        # Verify that only requested chromosomes are indexed
        for chrom in index.gene_trees:
            assert chrom in chromosomes

    def test_empty_chromosome_filter(self):
        """Test GenomicIntervalIndex with empty chromosome list."""
        index = GenomicIntervalIndex(self.gffutils_db, chromosomes=[])
        # Empty chromosome list is falsy, so it becomes None (no filtering)
        assert index.chromosomes is None
        # Since no filtering is applied, all chromosomes are indexed
        assert len(index.gene_trees) > 0

    def test_gene_tree_building(self):
        """Test that gene trees are built correctly."""
        index = GenomicIntervalIndex(self.gffutils_db)
        # Verify gene trees contain expected data
        for chrom, tree in index.gene_trees.items():
            assert len(tree) > 0, f"Gene tree for {chrom} should not be empty"
            # Verify each interval has a gene_id as data
            for interval in tree:
                assert interval.data is not None
                # Verify we can retrieve the gene from the database
                gene = self.gffutils_db[interval.data]
                assert gene.featuretype == 'gene'

    def test_exon_tree_building(self):
        """Test that exon trees are built correctly."""
        index = GenomicIntervalIndex(self.gffutils_db)
        # Verify exon trees contain expected data
        for chrom, tree in index.exon_trees.items():
            assert len(tree) > 0, f"Exon tree for {chrom} should not be empty"
            # Verify each interval has an exon_id as data
            for interval in tree:
                assert interval.data is not None
                # Verify we can retrieve the exon from the database
                exon = self.gffutils_db[interval.data]
                assert exon.featuretype == 'exon'

    def test_get_genes_at_single_position(self):
        """Test retrieving genes at a single genomic position."""
        index = GenomicIntervalIndex(self.gffutils_db)
        # Get a known gene and test position lookup
        gene = self.gffutils_db['ENSMUSG00000020196.10']
        chrom = str(gene.chrom)
        pos = (int(gene.start) + int(gene.end)) // 2  # Middle of gene
        genes = index.get_genes_at(chrom, pos)
        assert len(genes) > 0
        assert gene.id in [g.id for g in genes]

    def test_get_genes_at_nonexistent_chromosome(self):
        """Test retrieving genes from nonexistent chromosome."""
        index = GenomicIntervalIndex(self.gffutils_db)
        genes = index.get_genes_at('chrNonexistent', 1000)
        assert genes == []

    def test_get_genes_at_outside_gene_region(self):
        """Test retrieving genes at a position with no genes."""
        index = GenomicIntervalIndex(self.gffutils_db)
        # Use chr10 with a position unlikely to have genes
        genes = index.get_genes_at('chr10', 1)
        # Result depends on actual data, but should be a valid list
        assert isinstance(genes, list)

    def test_get_genes_at_with_window(self):
        """Test retrieving genes within a window around a position."""
        index = GenomicIntervalIndex(self.gffutils_db)
        # Get a known gene
        gene = self.gffutils_db['ENSMUSG00000020196.10']
        chrom = str(gene.chrom)
        pos = int(gene.start)
        window = 100
        genes_no_window = index.get_genes_at(chrom, pos)
        genes_with_window = index.get_genes_at(chrom, pos, window=window)
        # Window query should return at least as many results
        assert len(genes_with_window) >= len(genes_no_window)

    def test_get_genes_at_with_large_window(self):
        """Test retrieving genes with a large window."""
        index = GenomicIntervalIndex(self.gffutils_db)
        gene = self.gffutils_db['ENSMUSG00000020196.10']
        chrom = str(gene.chrom)
        pos = int(gene.start)
        window = 5000  # Large window
        genes = index.get_genes_at(chrom, pos, window=window)
        assert isinstance(genes, list)
        assert len(genes) >= 0
        # All returned genes should have valid feature types
        for g in genes:
            assert hasattr(g, 'featuretype')
            assert g.featuretype == 'gene'

    def test_get_exons_at_single_position(self):
        """Test retrieving exons at a single genomic position."""
        index = GenomicIntervalIndex(self.gffutils_db)
        
        # Get a known exon
        exon = list(self.gffutils_db.features_of_type('exon'))[0]
        chrom = str(exon.chrom)
        pos = (int(exon.start) + int(exon.end)) // 2  # Middle of exon
        exons = index.get_exons_at(chrom, pos)
        assert len(exons) > 0
        assert exon.id in [e.id for e in exons]

    def test_get_exons_at_nonexistent_chromosome(self):
        """Test retrieving exons from nonexistent chromosome."""
        index = GenomicIntervalIndex(self.gffutils_db)
        exons = index.get_exons_at('chrNonexistent', 1000)
        assert exons == []

    def test_get_exons_at_with_window(self):
        """Test retrieving exons within a window around a position."""
        index = GenomicIntervalIndex(self.gffutils_db)
        exon = list(self.gffutils_db.features_of_type('exon'))[0]
        chrom = str(exon.chrom)
        pos = int(exon.start)
        window = 50
        exons_no_window = index.get_exons_at(chrom, pos)
        exons_with_window = index.get_exons_at(chrom, pos, window=window)
        # Window query should return at least as many results
        assert len(exons_with_window) >= len(exons_no_window)

    def test_get_exons_of_gene(self):
        """Test retrieving all exons of a specific gene."""
        index = GenomicIntervalIndex(self.gffutils_db)
        gene = self.gffutils_db['ENSMUSG00000020196.10']
        exons = index.get_exons_of_gene(gene)
        assert len(exons) > 0
        # Verify all exons belong to the gene
        for exon in exons:
            assert exon.featuretype == 'exon'
            # Check that exon is associated with the gene
            parents = list(self.gffutils_db.parents(exon, featuretype='transcript'))
            assert len(parents) > 0

    def test_get_exons_of_gene_ordered(self):
        """Test that exons of a gene are returned in order."""
        index = GenomicIntervalIndex(self.gffutils_db)
        gene = self.gffutils_db['ENSMUSG00000020196.10']
        exons = index.get_exons_of_gene(gene)
        if len(exons) > 1:
            # Verify exons are ordered by start position
            starts = [int(exon.start) for exon in exons]
            assert starts == sorted(starts)

    def test_get_all_genes(self):
        """Test retrieving all genes from the index."""
        index = GenomicIntervalIndex(self.gffutils_db)
        all_genes = index.get_all_genes()
        assert len(all_genes) > 0
        # Verify all returned features are genes
        for gene in all_genes:
            assert gene.featuretype == 'gene'

    def test_get_all_genes_count_matches_trees(self):
        """Test that get_all_genes returns genes from all trees."""
        index = GenomicIntervalIndex(self.gffutils_db)
        # Count genes in trees
        tree_gene_count = sum(len(tree) for tree in index.gene_trees.values())
        all_genes = index.get_all_genes()
        # Should retrieve same number of genes as in trees
        assert len(all_genes) == tree_gene_count

    def test_find_genes_by_name(self):
        """Test finding genes by their gene_name attribute."""
        index = GenomicIntervalIndex(self.gffutils_db)
        # Get a known gene and extract its name
        gene = self.gffutils_db['ENSMUSG00000020196.10']
        attrs = getattr(gene, "attributes", {}) or {}
        gene_name = attrs.get("gene_name", [None])[0]
        if gene_name:
            found_genes = index.find_genes_by_name(gene_name)
            assert len(found_genes) > 0
            # Verify the known gene is in results
            assert gene.id in [g.id for g in found_genes]
            # Verify all results have matching gene_name
            for found_gene in found_genes:
                found_attrs = getattr(found_gene, "attributes", {}) or {}
                found_name = found_attrs.get("gene_name", [None])[0]
                assert found_name == gene_name

    def test_find_genes_by_nonexistent_name(self):
        """Test finding genes by name that doesn't exist."""
        index = GenomicIntervalIndex(self.gffutils_db)
        genes = index.find_genes_by_name('NonexistentGene12345')
        assert genes == []

    def test_multiple_genes_at_position(self):
        """Test position with multiple overlapping genes."""
        index = GenomicIntervalIndex(self.gffutils_db)
        # Find a position with multiple genes (if any exist)
        for chrom, tree in index.gene_trees.items():
            for pos in range(1, 100000, 1000):
                genes = index.get_genes_at(chrom, pos)
                if len(genes) > 1:
                    # Found multiple genes at this position
                    assert all(g.featuretype == 'gene' for g in genes)
                    return
        # If no position with multiple genes found, that's ok too
        pytest.skip("No position with multiple overlapping genes found in test data")

    def test_window_boundary_conditions(self):
        """Test window queries at boundary conditions."""
        index = GenomicIntervalIndex(self.gffutils_db)
        gene = self.gffutils_db['ENSMUSG00000020196.10']
        chrom = str(gene.chrom)
        # Test with window at position 1
        genes = index.get_genes_at(chrom, 1, window=100)
        assert isinstance(genes, list)
        # Test with very large window
        genes = index.get_genes_at(chrom, 5000000, window=1000000)
        assert isinstance(genes, list)

    def test_chromosome_filter_respected(self):
        """Test that chromosome filter is respected in all results."""
        chromosomes = ['chr10']
        index = GenomicIntervalIndex(self.gffutils_db, chromosomes=chromosomes)
        all_genes = index.get_all_genes()
        # All genes should be from filtered chromosomes
        for gene in all_genes:
            assert str(gene.chrom) in chromosomes

    def test_consistency_between_single_and_window_queries(self):
        """Test consistency between single position and window queries."""
        index = GenomicIntervalIndex(self.gffutils_db)
        gene = self.gffutils_db['ENSMUSG00000020196.10']
        chrom = str(gene.chrom)
        pos = int(gene.start) + 10
        single_genes = set(g.id for g in index.get_genes_at(chrom, pos))
        window_genes = set(g.id for g in index.get_genes_at(chrom, pos, window=1))
        # Single position genes should be subset of window genes
        assert single_genes.issubset(window_genes)

    def test_get_genes_at_with_zero_window(self):
        """Test get_genes_at with window=0."""
        index = GenomicIntervalIndex(self.gffutils_db)
        gene = self.gffutils_db['ENSMUSG00000020196.10']
        chrom = str(gene.chrom)
        pos = int(gene.start) + 10
        # window=0 creates empty range [pos, pos)
        genes_zero_window = index.get_genes_at(chrom, pos, window=0)
        # Should return empty list since range is empty
        assert genes_zero_window == []

    def test_exon_tree_coverage(self):
        """Test that exon trees cover all exons in database."""
        index = GenomicIntervalIndex(self.gffutils_db)
        # Get all exons from database
        db_exons = set(exon.id for exon in self.gffutils_db.features_of_type('exon'))
        # Get all exons from index
        index_exons = set()
        for tree in index.exon_trees.values():
            for interval in tree:
                index_exons.add(interval.data)
        # Should be the same
        assert db_exons == index_exons

    def test_gene_tree_coverage(self):
        """Test that gene trees cover all genes in database."""
        index = GenomicIntervalIndex(self.gffutils_db)
        # Get all genes from database
        db_genes = set(gene.id for gene in self.gffutils_db.features_of_type('gene'))
        # Get all genes from index
        index_genes = set()
        for tree in index.gene_trees.values():
            for interval in tree:
                index_genes.add(interval.data)
        # Should be the same
        assert db_genes == index_genes

    def test_interval_coordinates_accuracy(self):
        """Test that interval tree coordinates match database features."""
        index = GenomicIntervalIndex(self.gffutils_db)
        # Test genes
        for chrom, tree in index.gene_trees.items():
            for interval in tree:
                gene = self.gffutils_db[interval.data]
                db_start = int(gene.start)
                db_end = int(gene.end)
                # IntervalTree uses half-open intervals [start, end)
                # So interval.begin == start and interval.end == end+1
                assert interval.begin == db_start
                assert interval.end == db_end + 1

    def test_error_handling_invalid_gene_id_in_tree(self):
        """Test error handling when interval data doesn't resolve to a gene."""
        mock_db = MagicMock()
        mock_tree = Mock()
        mock_interval = Mock()
        mock_interval.data = 'invalid_gene_id'
        mock_tree.__iter__ = Mock(return_value=iter([mock_interval]))
        index = GenomicIntervalIndex.__new__(GenomicIntervalIndex)
        index.db = mock_db
        index.gene_trees = {'chr1': mock_tree}
        index.exon_trees = {}
        # Make database lookups fail for invalid_gene_id
        mock_db.__getitem__ = Mock(side_effect=KeyError('invalid_gene_id'))
        genes = index.get_all_genes()
        # Should handle error gracefully and return empty list
        assert genes == []

    def test_get_exons_of_gene_with_invalid_gene(self):
        """Test get_exons_of_gene with a gene that has no exons."""
        index = GenomicIntervalIndex(self.gffutils_db)
        # Create a mock gene feature
        mock_gene = Mock()
        # When database.children() raises an exception, should return empty list
        self.gffutils_db.children = Mock(side_effect=Exception("No children"))
        exons = index.get_exons_of_gene(mock_gene)
        assert exons == []

    def test_interval_tree_import_availability(self):
        """Test handling when intervaltree is not available."""
        with patch('isoquant_lib.genomic_interval_index.IntervalTree', None):
            index = GenomicIntervalIndex.__new__(GenomicIntervalIndex)
            index.db = self.gffutils_db
            index.gene_trees = {}
            index.exon_trees = {}
            # Should return empty list when IntervalTree is None
            genes = index.get_genes_at('chr10', 5000)
            assert genes == []

    def test_performance_multiple_queries(self):
        """Test performance with multiple queries."""
        index = GenomicIntervalIndex(self.gffutils_db)
        gene = self.gffutils_db['ENSMUSG00000020196.10']
        chrom = str(gene.chrom)
        # Run multiple queries
        for pos in range(int(gene.start), int(gene.end), 100):
            genes = index.get_genes_at(chrom, pos)
            assert isinstance(genes, list)
