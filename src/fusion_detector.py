import re
import pysam
import gffutils
from collections import defaultdict
import mappy as mp

class FusionDetector:
    def __init__(self, bam_path, gene_db_path, reference_fasta):
        self.bam_path = bam_path
        self.db = gffutils.FeatureDB(gene_db_path, keep_order=True)
        self.fusion_candidates = defaultdict(set)
        # store breakpoint counts per fusion key: {(chr1,pos1,chr2,pos2): count}
        self.fusion_breakpoints = defaultdict(lambda: defaultdict(int))
        # optional reference and aligner for soft-clip realignment
        self.reference_fasta = reference_fasta
        self.aligner = mp.Aligner(reference_fasta)

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
        """
        Parse a single SA tag entry and return (chrom, pos, cigar).
        SA format: chrom,pos,strand,CIGAR,mapQ,NM
        """
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
        """
        Realign a clipped sequence using mappy aligner.
        Returns the best hit or None if no suitable alignment found.
        """
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

    def detect_fusions(self):
        bam = pysam.AlignmentFile(self.bam_path, "rb")
        for read in bam:
            if read.has_tag("SA"):
                self._process_supplementary_alignment(read)
            else:
                self._process_softclip(read)
        bam.close()

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
        with open(output_path, "w") as f:
            f.write("LeftGene\tLeftChromosome\tLeftBreakpoint\tRightGene\tRightChromosome\tRightBreakpoint\tSupportingReads\tFusionName\n")
            for fusion, reads in self.fusion_candidates.items():
                support = len(reads)
                if support < min_support:
                    continue
                # choose consensus breakpoint = bp with maximum supporting read count
                bp_counts = self.fusion_breakpoints.get(fusion, {})
                if not bp_counts:
                    # fallback: no stored breakpoints -> skip
                    continue
                consensus_bp, _ = max(bp_counts.items(), key=lambda item: item[1])
                left_chrom, left_pos, right_chrom, right_pos = consensus_bp
                # obtain gene names/contexts at consensus positions
                left_gene = self.get_context(left_chrom, left_pos)
                right_gene = self.get_context(right_chrom, right_pos)
                fusion_name = fusion  # keep existing fusion key/name
                f.write(f"{left_gene}\t{left_chrom}\t{left_pos}\t{right_gene}\t{right_chrom}\t{right_pos}\t{support}\t{fusion_name}\n")