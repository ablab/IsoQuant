import pysam
import gffutils
from collections import defaultdict

class FusionDetector:
    def __init__(self, bam, gene_db):
        self.bam_path = bam
        self.db = gffutils.FeatureDB(gene_db, keep_order=True)
        # store supporting read names as a set to count unique supporting reads
        self.fusion_candidates = defaultdict(set)

    def get_gene_name(self, chrom, pos):
        try:
            for feature in self.db.region(region=(chrom, pos, pos), featuretype='gene'):
                return feature.attributes.get('gene_name', [feature.id])[0]
        except Exception:
            return None

    def detect_fusions(self):
        bam = pysam.AlignmentFile(self.bam_path, "rb")
        for read in bam:
            if read.has_tag("SA"):  # Supplementary alignment tag
                sa_tag = read.get_tag("SA")
                for sa in sa_tag.split(";"):
                    if not sa:
                        continue
                    fields = sa.split(",")
                    chrom, pos = fields[0], int(fields[1])
                    gene1 = self.get_gene_name(read.reference_name, read.reference_start)
                    gene2 = self.get_gene_name(chrom, pos)
                    if gene1 and gene2 and gene1 != gene2:
                        fusion_key = f"{gene1}--{gene2}"
                        self.fusion_candidates[fusion_key].add(read.query_name)
        bam.close()

    def report(self, output_path="fusion_candidates.tsv"):
        MIN_SUPPORT = 3
        with open(output_path, "w") as f:
            f.write("Fusion\tSupporting_Reads\n")
            for fusion, reads in self.fusion_candidates.items():
                support = len(reads)
                if support >= MIN_SUPPORT:
                    f.write(f"{fusion}\t{support}\n")


