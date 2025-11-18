import re
import pysam
import gffutils
from collections import defaultdict
import mappy as mp
import json
import os
import tempfile
import logging
from src.fusion_quick_loader import process_fusion_bundle

logger = logging.getLogger('IsoQuant')

class FusionDetector:
    def __init__(self, bam_path, gene_db_path, reference_fasta):
        self.bam_path = bam_path
        # keep the DB path for downstream consumers
        self.genedb_path = gene_db_path
        self.db = gffutils.FeatureDB(gene_db_path, keep_order=True)
        self.fusion_candidates = defaultdict(set)
        # store breakpoint counts per fusion key: {(chr1,pos1,chr2,pos2): count}
        self.fusion_breakpoints = defaultdict(lambda: defaultdict(int))
        # optional reference and aligner for soft-clip realignment
        self.reference_fasta = reference_fasta
        self.aligner = mp.Aligner(reference_fasta)
        # collected metadata for downstream filtering / model construction
        # fusion_metadata[fusion_key] = {
        #   "supporting_reads": set(),
        #   "consensus_bp": (chr1,pos1,chr2,pos2),
        #   "left_gene": str, "right_gene": str,
        #   "support": int
        # }
        self.fusion_metadata = {}

    def safe_reference_start(self, read):
        try:
            return read.reference_start + 1 if read.reference_start is not None else None
        except Exception:
            return None

    def compute_aligned_length(self, read):
        try:
            if getattr(read, "cigartuples", None):
                return self.aligned_len_from_cigartuples(read.cigartuples)
            return self.aligned_len_from_cigarstring(getattr(read, "cigarstring", None))
        except Exception:
            return 0

    def parse_sa_fields(self, sa):
        # Parse a single SA tag entry and return (chrom, pos, cigar).
        # format: chrom,pos,strand,CIGAR,mapQ,NM
        fields = sa.split(",")
        chrom = fields[0]
        pos = int(fields[1])
        sa_cigar = fields[3] if len(fields) > 3 else None
        return chrom, pos, sa_cigar
    
    def get_context(self, chrom, pos):
        try:
            genes = list(self.db.region(region=(chrom, pos, pos), featuretype='gene'))
            if genes:
                return genes[0].attributes.get('gene_name', [genes[0].id])[0]
            exons = list(self.db.region(region=(chrom, pos, pos), featuretype='exon'))
            if exons:
                return "intronic"
            return "intergenic"
        except Exception:
            return "unknown"

    def detect_softclip(self, read):
        """
        Detect large soft-clips (>=50 bp) on either end of the read.
        Returns (clip_side, clip_len) or (None, 0) if no significant clip.
        """
        cigartuples = getattr(read, "cigartuples", None)
        clip_side, clip_len = None, 0

        if cigartuples:
            # Leading soft-clip
            if len(cigartuples) > 0 and cigartuples[0][0] == 4 and cigartuples[0][1] >= 50:
                clip_side, clip_len = "left", cigartuples[0][1]
            # Trailing soft-clip
            elif len(cigartuples) > 0 and cigartuples[-1][0] == 4 and cigartuples[-1][1] >= 50:
                clip_side, clip_len = "right", cigartuples[-1][1]
        else:
            # Fallback: parse cigarstring
            cigarstr = getattr(read, "cigarstring", "") or ""
            m1 = re.match(r'^(\d+)S', cigarstr)
            m2 = re.search(r'(\d+)S$', cigarstr)
            if m1 and int(m1.group(1)) >= 50:
                clip_side, clip_len = "left", int(m1.group(1))
            elif m2 and int(m2.group(1)) >= 50:
                clip_side, clip_len = "right", int(m2.group(1))
        return clip_side, clip_len

    def aligned_len_from_cigartuples(self, cigartuples):
        # count read-aligned operations: M (=0), = (=7), X (=8)
        if not cigartuples:
            return 0
        return sum(l for op, l in cigartuples if op in (0, 7, 8))

    def aligned_len_from_cigarstring(self, cigar):
        if not cigar:
            return 0
        total = 0
        for m in re.finditer(r'(\d+)([MIDNSHP=XB])', cigar):
            length = int(m.group(1)); op = m.group(2)
            if op in ('M', '=', 'X'):
                total += length
        return total
    
    def realign_clipped_seq(self, seq):
        # Realign a clipped sequence using mappy aligner. 
        try:
            for hit in self.aligner.map(seq):
                hit_len = (getattr(hit, "r_en", None) - getattr(hit, "r_st", None)) if getattr(hit, "r_en", None) else getattr(hit, "mlen", None) or getattr(hit, "alen", None) or 0
                if hit_len and hit_len >= 50:
                    return hit  # first good hit
            return None
        except Exception:
            return None

    def record_fusion(self, context1, context2, read_name, chrom1, pos1, chrom2, pos2):
        fusion_key = "--".join(sorted([context1, context2]))
        self.fusion_candidates[fusion_key].add(read_name)
        bp = (chrom1, int(pos1), chrom2, int(pos2))
        self.fusion_breakpoints[fusion_key][bp] += 1
        # also maintain supporting_reads set for quick access
        meta = self.fusion_metadata.setdefault(fusion_key, {"supporting_reads": set(), "consensus_bp": None,
                                                            "left_gene": None, "right_gene": None, "support": 0})
        meta["supporting_reads"].add(read_name)
        meta["support"] = len(meta["supporting_reads"])
        # consensus_bp will be calculated later in build_metadata()

    def build_metadata(self, min_support=1):
        # Compute consensus breakpoint and gene names for all fusion keys,
        for fusion_key, reads in self.fusion_candidates.items():
            support = len(reads)
            meta = self.fusion_metadata.setdefault(fusion_key, {"supporting_reads": set(), "consensus_bp": None,
                                                                "left_gene": None, "right_gene": None, "support": 0})
            meta["supporting_reads"].update(reads)
            meta["support"] = len(meta["supporting_reads"])
            if meta["support"] < min_support:
                continue
            bp_counts = self.fusion_breakpoints.get(fusion_key, {})
            if not bp_counts:
                continue
            consensus_bp, _ = max(bp_counts.items(), key=lambda item: item[1])
            meta["consensus_bp"] = consensus_bp
            left_chr, left_pos, right_chr, right_pos = consensus_bp
            meta["left_gene"] = self.get_context(left_chr, left_pos)
            meta["right_gene"] = self.get_context(right_chr, right_pos)

    def get_fusion_metadata(self, min_support=1):
        """
        Returns list of metadata dicts for fusions with support >= min_support.
        Each dict contains keys: fusion_key, left_gene, right_gene, consensus_bp, support, supporting_reads
        """
        self.build_metadata(min_support=min_support)
        results = []
        for fusion_key, meta in self.fusion_metadata.items():
            if meta.get("support", 0) < min_support:
                continue
            results.append({
                "fusion_key": fusion_key,
                "left_gene": meta.get("left_gene"),
                "right_gene": meta.get("right_gene"),
                "consensus_bp": meta.get("consensus_bp"),
                "support": meta.get("support"),
                "supporting_reads": set(meta.get("supporting_reads", set()))
            })
        return results

    def detect_fusions(self):
        bam = pysam.AlignmentFile(self.bam_path, "rb")
        for read in bam:
            if read.has_tag("SA"):
                self._process_supplementary_alignment(read)
            else:
                self._process_softclip(read)
        bam.close()

        # export assignment bundle and run downstream IsoQuant model construction
        try:
            bundle_dir = tempfile.mkdtemp(prefix="fusion_bundle_")
            self.export_assignment_bundle(bundle_dir)

            # minimal dummy counters/params for process_fusion_bundle
            class _DummyCounter:
                def add_read_info_raw(self, *a, **kw): pass
                def add_unassigned(self, *a, **kw): pass
                def add_confirmed_features(self, *a, **kw): pass

            params = type("P", (), {})()  # empty params object; pipeline may override as needed
            transcript_counter = _DummyCounter()
            gene_counter = _DummyCounter()

            try:
                process_fusion_bundle(bundle_dir, self.genedb_path, params, transcript_counter, gene_counter)
            except Exception as e:
                logger.warning("process_fusion_bundle() failed: %s" % str(e))

            # delegate TSV generation to report()
            out_tsv = os.path.join(bundle_dir, "fusion_candidates_postprocessed.tsv")
            self.report(output_path=out_tsv, min_support=3)
            logger.info("Fusion bundle exported to %s; postprocessed TSV written to %s" % (bundle_dir, out_tsv))
        except Exception as e:
            logger.error("Fusion post-processing failed: %s" % str(e))

    def _process_supplementary_alignment(self, read):
        sa_tag = read.get_tag("SA")
        primary_chrom = read.reference_name
        primary_pos = self.safe_reference_start(read)
        primary_context = self.get_context(primary_chrom, primary_pos)
        primary_al_len = self.compute_aligned_length(read)
        # iterate each supplementary alignment entry and handle them independently;
        # record a fusion immediately when a valid cross-context supplementary alignment is found
        for sa in sa_tag.split(";"):
            if not sa:
                continue
            chrom, pos, sa_cigar = self.parse_sa_fields(sa)
            sa_al_len = self.aligned_len_from_cigarstring(sa_cigar)

            # require sufficient aligned length on both sides
            if primary_al_len < 100 or sa_al_len < 100:
                continue

            context = self.get_context(chrom, pos)
            if context != primary_context:
                self.record_fusion(primary_context, context, read.query_name, primary_chrom, primary_pos, chrom, pos)

    def _process_softclip(self, read):
        clip_side, clip_len = self.detect_softclip(read)
        if not clip_side or self.aligner is None:
            return
        seq = read.query_sequence
        if not seq:
            return
        clipped_seq = seq[:clip_len] if clip_side == "left" else seq[-clip_len:]
        best_hit = self.realign_clipped_seq(clipped_seq)
        if not best_hit:
            return
        sa_chrom, sa_pos = getattr(best_hit, "ctg", None), getattr(best_hit, "r_st", None)
        if sa_chrom is None or sa_pos is None:
            return
        sa_pos += 1
        primary_chrom = read.reference_name
        primary_pos = self.safe_reference_start(read)
        primary_al_len = self.compute_aligned_length(read)
        if primary_al_len < 100:
            return
        primary_context = self.get_context(primary_chrom, primary_pos)
        sa_context = self.get_context(sa_chrom, sa_pos)
        if primary_context != sa_context:
            self.record_fusion(primary_context, sa_context, read.query_name, primary_chrom, primary_pos, sa_chrom, sa_pos)

    def report(self, output_path="fusion_candidates.tsv", min_support=3):
        #  write a TSV with one consensus breakpoint per reported fusion.
        # prefer metadata-driven output (consensus breakpoint + supporting reads)
        self.build_metadata(min_support=min_support)
        with open(output_path, "w") as f:
            f.write("LeftGene\tLeftChromosome\tLeftBreakpoint\tRightGene\tRightChromosome\tRightBreakpoint\tSupportingReads\tFusionName\n")
            for fusion_key, meta in sorted(self.fusion_metadata.items(), key=lambda x: -x[1].get("support", 0)):
                support = meta.get("support", 0)
                if support < min_support:
                    continue
                consensus_bp = meta.get("consensus_bp")
                if not consensus_bp:
                    continue
                left_chr, left_pos, right_chr, right_pos = consensus_bp
                left_gene = meta.get("left_gene") or self.get_context(left_chr, left_pos)
                right_gene = meta.get("right_gene") or self.get_context(right_chr, right_pos)
                fusion_name = fusion_key
                f.write(f"{left_gene}\t{left_chr}\t{left_pos}\t{right_gene}\t{right_chr}\t{right_pos}\t{support}\t{fusion_name}\n")

    def export_assignment_bundle(self, out_dir):
        os.makedirs(out_dir, exist_ok=True)
        fusions = []
        # ensure metadata is built
        self.build_metadata(min_support=1)
        for fk, meta in self.fusion_metadata.items():
            fusions.append({
                "fusion_key": fk,
                "consensus_bp": meta.get("consensus_bp"),
                "support": meta.get("support", 0),
                "supporting_reads": list(meta.get("supporting_reads", []))
            })

        with open(os.path.join(out_dir, "fusion_candidates.json"), "w") as fh:
            json.dump(fusions, fh, indent=2)

        # write chr record if reference fasta known
        chr_record = {}
        if getattr(self, "reference_fasta", None):
            try:
                import pysam
                fa = pysam.FastaFile(self.reference_fasta)
                for ctg in fa.references:
                    chr_record[ctg] = fa.get_reference_length(ctg)
                fa.close()
            except Exception:
                chr_record = {}
        with open(os.path.join(out_dir, "chr_record.json"), "w") as fh:
            json.dump(chr_record, fh, indent=2)

        # placeholder multimapped dictionary (fusion detection may produce cross-chromosome multimappers)
        multimapped = {}
        with open(os.path.join(out_dir, "multimapped.json"), "w") as fh:
            json.dump(multimapped, fh, indent=2)

        return out_dir