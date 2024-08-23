############################################################################
# Copyright (c) 2024 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################


from .graph_based_model_construction import GraphBasedModelConstructor


class SimpleIDDistributor(object):
    def __init__(self):
        self.value = 0

    def increment(self):
        self.value += 1
        return self.value


class ExcludingIdDistributor(SimpleIDDistributor):
    def __init__(self, genedb, chr_id):
        SimpleIDDistributor.__init__(self)
        self.forbidden_ids = set()
        if not genedb:
            return

        for g in genedb.region(seqid=chr_id, start=1, featuretype="gene"):
            if g.id.startswith(GraphBasedModelConstructor.novel_gene_prefix):
                try:
                    gene_num = int(g.id.split("_")[-1])
                    self.forbidden_ids.add(gene_num)
                except IndexError:
                    pass
                except ValueError:
                    pass

        transcript_num_start_pos = len(GraphBasedModelConstructor.transcript_prefix)
        for t in genedb.region(seqid=chr_id, start=1, featuretype=("transcript", "mRNA")):
            if t.id.startswith(GraphBasedModelConstructor.transcript_prefix):
                try:
                    transcript_num = int(t.id.split(".")[0][transcript_num_start_pos:])
                    self.forbidden_ids.add(transcript_num)
                except IndexError:
                    pass
                except ValueError:
                    pass

    def increment(self):
        self.value += 1
        while self.value in self.forbidden_ids:
            self.value += 1
        return self.value

