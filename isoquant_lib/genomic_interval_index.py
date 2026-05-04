import logging
from typing import Optional

try:
    from intervaltree import IntervalTree
except ImportError:
    IntervalTree = None

logger = logging.getLogger('IsoQuant')

class GenomicIntervalIndex:
    """Interval-tree index for mapping genomic coordinates to biological entities.

    Pre-builds interval trees for genes and exons to eliminate repeated database queries.
    """

    def __init__(self, genedb, chromosomes: Optional[list[str]] = None) -> None:
        self.db = genedb
        self.gene_trees: dict = {}
        self.exon_trees: dict = {}
        self.chromosomes = set(chromosomes) if chromosomes else None
        self._build_indices()

    def _build_indices(self) -> None:
        """Build interval trees for all genes and exons in the database."""
        logger.info("Building genomic interval trees (memory-efficient mode)...")
        gene_count = 0
        try:
            for gene in self.db.features_of_type('gene'):
                chrom = str(gene.chrom)
                if self.chromosomes and chrom not in self.chromosomes:
                    continue

                if chrom not in self.gene_trees:
                    self.gene_trees[chrom] = IntervalTree()
                start = int(gene.start)
                end = int(gene.end)
                gene_id = gene.id
                self.gene_trees[chrom][start:end+1] = gene_id
                gene_count += 1
                if gene_count % 10000 == 0:
                    logger.debug(f"  Indexed {gene_count} genes...")
        except Exception as e:
            logger.warning(f"Error building gene index: {e}")
        exon_count = 0
        try:
            for exon in self.db.features_of_type('exon'):
                chrom = str(exon.chrom)
                if self.chromosomes and chrom not in self.chromosomes:
                    continue
                if chrom not in self.exon_trees:
                    self.exon_trees[chrom] = IntervalTree()
                start = int(exon.start)
                end = int(exon.end)
                exon_id = exon.id
                self.exon_trees[chrom][start:end+1] = exon_id
                exon_count += 1
                if exon_count % 50000 == 0:
                    logger.debug(f"  Indexed {exon_count} exons...")
        except Exception as e:
            logger.warning(f"Error building exon index: {e}")
        logger.info(f"Built gene trees for {len(self.gene_trees)} chromosomes with {gene_count} genes")
        logger.info(f"Built exon trees for {len(self.exon_trees)} chromosomes with {exon_count} exons")

    def get_genes_at(self, chrom: str, pos: int, window: Optional[int] = None) -> list:
        """Return all gene features overlapping ``pos`` (or ``[pos-window, pos+window]``) on ``chrom``."""
        if IntervalTree is None or chrom not in self.gene_trees:
            return []
        if window is None:
            intervals = self.gene_trees[chrom][pos]
        else:
            start = max(1, pos - window)
            end = pos + window
            intervals = self.gene_trees[chrom][start:end]
        genes = []
        for iv in intervals:
            try:
                gene_id = iv.data
                gene = self.db[gene_id]
                genes.append(gene)
            except Exception:
                logger.debug(
                    "Failed to resolve gene %r from interval tree at %s:%s",
                    iv.data, chrom, pos,
                    exc_info=True,
                )
        return genes

    def get_exons_at(self, chrom: str, pos: int, window: Optional[int] = None) -> list:
        """Return all exon features overlapping ``pos`` (or ``[pos-window, pos+window]``) on ``chrom``."""
        if IntervalTree is None or chrom not in self.exon_trees:
            return []
        if window is None:
            intervals = self.exon_trees[chrom][pos]
        else:
            start = max(1, pos - window)
            end = pos + window
            intervals = self.exon_trees[chrom][start:end]
        exons = []
        for iv in intervals:
            try:
                exon_id = iv.data
                exon = self.db[exon_id]
                exons.append(exon)
            except Exception:
                logger.debug(
                    "Failed to resolve exon %r from interval tree at %s:%s",
                    iv.data, chrom, pos,
                    exc_info=True,
                )
        return exons

    def get_exons_of_gene(self, gene_feature) -> list:
        """Return all exon children for a specific gene feature, ordered by start position."""
        try:
            return list(self.db.children(gene_feature, featuretype="exon", order_by="start"))
        except Exception:
            return []

    def get_all_genes(self) -> list:
        """Return all gene features known to the index."""
        genes = []
        for chrom, tree in self.gene_trees.items():
            for interval in tree:
                try:
                    gene_id = interval.data
                    gene = self.db[gene_id]
                    genes.append(gene)
                except Exception:
                    pass
        return genes

    def find_genes_by_name(self, gene_name: str) -> list:
        """Return all gene features whose ``gene_name`` attribute matches ``gene_name``."""
        matching_genes = []
        for chrom, tree in self.gene_trees.items():
            for interval in tree:
                try:
                    gene_id = interval.data
                    gene = self.db[gene_id]
                    attrs = getattr(gene, 'attributes', {}) or {}
                    if attrs.get('gene_name', [None])[0] == gene_name:
                        matching_genes.append(gene)
                except Exception:
                    pass
        return matching_genes
