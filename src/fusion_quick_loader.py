import os
import json
import logging
from collections import defaultdict

import gffutils
from .gene_info import GeneInfo

from .intron_graph import VERTEX_polya, VERTEX_polyt, VERTEX_read_end, VERTEX_read_start
from .graph_based_model_construction import GraphBasedModelConstructor
from .id_policy import SimpleIDDistributor

logger = logging.getLogger('IsoQuant')


class _PolyaStub:
    def __init__(self):
        self.external_polya_pos = -1
        self.internal_polya_pos = -1
        self.external_polyt_pos = -1
        self.internal_polyt_pos = -1


class _AssignmentTypeStub:
    # minimal interface used by downstream code
    def is_unique(self): return True
    def is_unassigned(self): return False
    def is_inconsistent(self): return False
    def is_consistent(self): return True


class FusionSyntheticAssignment:
    """
    Minimal assignment object compatible with GraphBasedModelConstructor.
    Fields used by the pipeline:
      - read_id
      - corrected_exons (list of (start,end))
      - corrected_introns (list of (start,end))
      - multimapper (bool)
      - mapping_quality (int)
      - polya_info (object with polyA positions)
      - strand ( '+' | '-' | '.' )
      - read_group (string)
      - assignment_type (object with is_unique / is_consistent etc.)
      - isoform_matches (list)  -- kept empty for synthetic fusion reads
      - gene_info (gffutils DB) -- set later by caller when needed
    """
    def __init__(self, read_id, corrected_exons, corrected_introns, mapping_quality=255):
        self.read_id = read_id
        self.corrected_exons = corrected_exons
        self.corrected_introns = corrected_introns
        self.multimapper = False
        self.mapping_quality = mapping_quality
        self.polya_info = _PolyaStub()
        self.strand = '.'
        self.read_group = 'fusion'
        self.assignment_type = _AssignmentTypeStub()
        self.isoform_matches = []
        self.gene_info = None


def _load_bundle(bundle_dir):
    bundle_dir = os.path.abspath(bundle_dir)
    fusions = []
    chr_record = {}
    multimapped = {}
    try:
        with open(os.path.join(bundle_dir, "fusion_candidates.json"), "r") as fh:
            fusions = json.load(fh)
    except Exception:
        logger.debug("No fusion_candidates.json found or failed to parse in %s" % bundle_dir)

    try:
        with open(os.path.join(bundle_dir, "chr_record.json"), "r") as fh:
            chr_record = json.load(fh)
    except Exception:
        chr_record = {}

    try:
        with open(os.path.join(bundle_dir, "multimapped.json"), "r") as fh:
            multimapped = json.load(fh)
    except Exception:
        multimapped = {}

    return fusions, chr_record, multimapped


def build_synthetic_assignments_from_bundle(bundle_dir):
    # Return (assignments_list, chr_record) built from fusion_candidates.json.
    fusions, chr_record, multimapped = _load_bundle(bundle_dir)
    assignments = []
    for fusion in fusions:
        bp = fusion.get("consensus_bp")
        if not bp:
            continue
        left_chr, left_pos, right_chr, right_pos = bp
        # create small exon windows around breakpoints so graph code forms introns
        left_exon = (max(1, left_pos - 50), left_pos)
        right_exon = (right_pos, right_pos + 50)
        corrected_exons = [left_exon, right_exon]
        # single junction between the two exons (may be cross-chromosome)
        corrected_introns = [(left_exon[1], right_exon[0])]
        for rid in fusion.get("supporting_reads", []):
            a = FusionSyntheticAssignment(read_id=rid,
                                          corrected_exons=corrected_exons,
                                          corrected_introns=corrected_introns,
                                          mapping_quality=255)
            assignments.append(a)
    return assignments, chr_record


def process_fusion_bundle(bundle_dir, gffutils_db_path, params, transcript_counter, gene_counter, id_distributor=None):
    """
    Build synthetic assignments from a fusion bundle and run GraphBasedModelConstructor.process()
    on them. Returns the instantiated GraphBasedModelConstructor on success.

    - bundle_dir: directory produced by FusionDetector.export_assignment_bundle()
    - gffutils_db_path: path to gene DB (.db) used to create GeneInfo (must exist)
    - params: pipeline parameter object (args)
    - transcript_counter / gene_counter: objects expected by GraphBasedModelConstructor (can be dummies)
    - id_distributor: optional SimpleIDDistributor; if None a new one is created.

    This function keeps all code out of isoquant.py; call it from isoquant when you want to
    hand the fusion candidates to IsoQuant's model construction.
    """
    if id_distributor is None:
        id_distributor = SimpleIDDistributor()

    # load gffutils DB
    try:
        db = gffutils.FeatureDB(gffutils_db_path, keep_order=True)
    except Exception as e:
        logger.error("Failed to open gffutils DB at %s: %s" % (gffutils_db_path, str(e)))
        raise

    assignments, chr_record = build_synthetic_assignments_from_bundle(bundle_dir)
    if not assignments:
        logger.info("No synthetic fusion assignments created from bundle %s" % bundle_dir)
        return None

    # Collect gene objects overlapping fusion breakpoints to build GeneInfo
    gene_db_list = []
    seen = set()
    try:
        with open(os.path.join(bundle_dir, "fusion_candidates.json"), "r") as fh:
            fusions = json.load(fh)
    except Exception:
        fusions = []

    for fusion in fusions:
        bp = fusion.get("consensus_bp")
        if not bp:
            continue
        left_chr, left_pos, right_chr, right_pos = bp
        for chr_id, pos in ((left_chr, left_pos), (right_chr, right_pos)):
            try:
                for gene in db.region(seqid=chr_id, start=max(1, pos - 1), end=pos, featuretype="gene"):
                    if gene.id not in seen:
                        gene_db_list.append(gene)
                        seen.add(gene.id)
            except Exception:
                continue

    # Fallback: if no genes found, create GeneInfo from a small region around first breakpoint
    if not gene_db_list and fusions:
        first = fusions[0].get("consensus_bp")
        if first:
            chr_id, pos = first[0], first[1]
            gene_info = GeneInfo.from_region(chr_id, max(1, pos - 1000), pos + 1000, delta=getattr(params, "delta", 0))
        else:
            gene_info = GeneInfo([], None, delta=getattr(params, "delta", 0))
    else:
        gene_info = GeneInfo(gene_db_list, db, delta=getattr(params, "delta", 0))

    # attach gene_info reference to synthetic assignments (some downstream code expects this)
    for a in assignments:
        a.gene_info = gene_info

    # instantiate and run model constructor with GeneInfo (not raw FeatureDB)
    try:
        gmc = GraphBasedModelConstructor(gene_info, chr_record, params, transcript_counter, gene_counter, id_distributor)
        gmc.process(assignments)
        logger.info("GraphBasedModelConstructor processed %d fusion assignments" % len(assignments))
        return gmc
    except Exception as e:
        logger.exception("GraphBasedModelConstructor failed on fusion assignments: %s" % str(e))
        raise