import pysam
import gffutils
from collections import defaultdict

class FusionDetector:
    def __init__(self, bam_path, gene_db_path):
        self.bam_path = bam_path
        self.db = gffutils.FeatureDB(gene_db_path, keep_order=True)
        self.fusion_candidates = defaultdict(set)
        # store breakpoint counts per fusion key: {(chr1,pos1,chr2,pos2): count}
        self.fusion_breakpoints = defaultdict(lambda: defaultdict(int))

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

    def detect_fusions(self):
        bam = pysam.AlignmentFile(self.bam_path, "rb")
        for read in bam:
            if read.has_tag("SA"):
                sa_tag = read.get_tag("SA")
                contexts = set()
                # primary alignment coordinates 
                primary_chrom = read.reference_name
                primary_pos = None
                try:
                    primary_pos = read.reference_start + 1 if read.reference_start is not None else None
                except Exception:
                    primary_pos = None

                primary_context = self.get_context(primary_chrom, primary_pos if primary_pos is not None else 0)
                contexts.add(primary_context)

                for sa in sa_tag.split(";"):
                    if not sa:
                        continue
                    fields = sa.split(",")
                    chrom, pos = fields[0], int(fields[1])
                    context = self.get_context(chrom, pos)
                    contexts.add(context)

                    # if this read indicates a cross-context event, record it and the bp
                    if primary_chrom and primary_pos is not None and context and primary_context and (context != primary_context):
                        fusion_key = "--".join(sorted([primary_context, context]))
                        self.fusion_candidates[fusion_key].add(read.query_name)
                        # store/increment breakpoint count (use 1-based positions)
                        bp = (primary_chrom, int(primary_pos), chrom, int(pos))
                        self.fusion_breakpoints[fusion_key][bp] += 1

        bam.close()

    def report(self, output_path="fusion_candidates.tsv", min_support=3):
        """
        Write a TSV with one consensus breakpoint per reported fusion.
        """
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