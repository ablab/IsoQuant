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


class FeatureIdStorage:
    def __init__(self, id_distributor, genedb = None, chr_id = None, feature = "exon"):
        self.id_distributor = id_distributor
        self.id_dict = {}
        self.feature_name = feature
        if not genedb or not chr_id:
            return

        id_attribute = feature + "_id"
        for f in genedb.region(seqid=chr_id, start=1, featuretype=feature):
            if id_attribute in f.attributes:
                feature_tuple = (chr_id, f.start, f.end, f.strand)
                try:
                    self.id_dict[feature_tuple] = f.attributes[id_attribute][0]
                except IndexError:
                    pass

    def get_id(self, chr_id, feature, strand):
        feature_tuple = (chr_id, feature[0], feature[1], strand)
        if feature_tuple not in self.id_dict:
            feature_id = self.id_distributor.increment()
            self.id_dict[feature_tuple] = chr_id + ".%d" % feature_id
        else:
            feature_id =  self.id_dict[feature_tuple]

        return feature_id
