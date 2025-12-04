import pandas as pd
import numpy as np
from scipy.signal import find_peaks, peak_prominences, peak_widths
from scipy import stats
from sklearn.metrics import confusion_matrix, f1_score, roc_auc_score, precision_recall_curve
from xgboost import XGBClassifier
import xgboost
import pickle
import matplotlib.pyplot as plt
import logging
import os
from collections import defaultdict, OrderedDict
from enum import Enum, unique

from .isoform_assignment import (
    ReadAssignmentType,
    ReadAssignment,
    IsoformMatch
)
from .gene_info import FeatureInfo, GeneInfo
from .read_groups import AbstractReadGrouper
from .convert_grouped_counts import GROUP_COUNT_CUTOFF, convert_to_mtx, convert_to_matrix

logger = logging.getLogger('IsoQuant')

from .long_read_counter import AbstractCounter



class PolyACounter(AbstractCounter):
    def __init__(self, output_prefix, read_groups, ignore_read_groups=False):
        super().__init__(output_prefix, ignore_read_groups=False)
        self.ignore_read_groups = ignore_read_groups
        self.read_groups = read_groups
        self.transcripts = {}    
        self.gene_info = None

        # get read groups

        self.first = True

        self.group_numeric_ids = {}
        self.ordered_groups = sorted(self.read_groups)
        if self.ordered_groups:
            for i, g in enumerate(self.ordered_groups):
                self.group_numeric_ids[g] = i
        
    
    def add_read_info(self, read_assignment: ReadAssignment):
        if self.gene_info is None:
            self.gene_info: GeneInfo = read_assignment.gene_info
        group_id = read_assignment.read_group
        group_id = self.group_numeric_ids[group_id]
        # add a single read_assignment to features dataframe


        if read_assignment.polyA_found == True and read_assignment.assignment_type in [ReadAssignmentType.inconsistent_non_intronic, ReadAssignmentType.unique, ReadAssignmentType.unique_minor_difference, ReadAssignmentType.inconsistent]:
            ok = False
            if read_assignment.strand == '-' and read_assignment.polya_info.external_polyt_pos != -1:
                polya_pos = read_assignment.polya_info.external_polyt_pos
                ok = True
            elif read_assignment.strand == '+' and read_assignment.polya_info.external_polya_pos != -1:
                polya_pos = read_assignment.polya_info.external_polya_pos
                ok = True
            if ok:                
                isoform_match: IsoformMatch = read_assignment.isoform_matches[0]
                if isoform_match.assigned_transcript not in self.transcripts:

                    if read_assignment.strand == '+':
                        annot = read_assignment.gene_info.all_isoforms_exons[isoform_match.assigned_transcript][-1][1] + 1
                    else:
                        annot = read_assignment.gene_info.all_isoforms_exons[isoform_match.assigned_transcript][0][0] - 1
                    self.transcripts[isoform_match.assigned_transcript] = {'chr': read_assignment.chr_id, 'gene_id': isoform_match.assigned_gene, 'data': [], 'annotated': annot}
                    if not self.ignore_read_groups:
                        for i in self.group_numeric_ids.values():
                            self.transcripts[isoform_match.assigned_transcript][i] = []

                
                self.transcripts[isoform_match.assigned_transcript]['data'].append(int(polya_pos))
                if not self.ignore_read_groups:
                    self.transcripts[isoform_match.assigned_transcript][group_id].append(int(polya_pos))
      

    def dump(self):
        self.model = XGBClassifier()
        self.model.load_model('src/model.json')
        self.df = pd.DataFrame({'transcript_id': self.transcripts.keys()})
        self.df['chromosome'] = self.df['transcript_id'].apply(lambda x: self.transcripts[x]['chr'])
        self.df['gene_id'] = self.df['transcript_id'].apply(lambda x: self.transcripts[x]['gene_id'])
        self.df['start'] = self.df['transcript_id'].apply(lambda x: np.min(self.transcripts[x]['data']))
        self.df['max'] = self.df['transcript_id'].apply(lambda x: stats.mode(self.transcripts[x]['data'])[0])-self.df.start
        self.df['mean'] = self.df['transcript_id'].apply(lambda x: np.mean(self.transcripts[x]['data']))-self.df.start
        self.df['histogram'] = self.df['transcript_id'].apply(lambda x: [0]*10 + list(np.histogram(self.transcripts[x]['data'], bins = 1+max(self.transcripts[x]['data']) - min(self.transcripts[x]['data']))[0])+[0]*10)
        self.df['mean_height'] = self.df['histogram'].apply(lambda x: np.mean(x))
        self.df['var'] = self.df['transcript_id'].apply(lambda x: np.var(self.transcripts[x]['data']))
        self.df['range'] = self.df['transcript_id'].apply(lambda x: np.max(self.transcripts[x]['data']) - np.min(self.transcripts[x]['data'])) + 1
        self.df['skew'] = self.df['transcript_id'].apply(lambda x: stats.skew(self.transcripts[x]['data']))
        self.df['peak_count'] = self.df['histogram'].apply(lambda x: len(find_peaks(x, distance=10)[0]))
        self.df['peak_info'] = self.df['histogram'].apply(lambda x: find_peaks(x, distance=10, threshold=(None, None), height=(None, None))[1])
        self.df['peak_location'] = self.df['histogram'].apply(lambda x: [int(j - 10) for j in find_peaks(x, distance=10)[0]])
        self.df['peak_prominence'] = self.df['histogram'].apply(lambda x: peak_prominences(x, find_peaks(x, distance=10)[0])[0])
        self.df['peak_width'] = self.df['histogram'].apply(lambda x: peak_widths(x, find_peaks(x, distance=10)[0], rel_height=0.98)[0])
        self.df['peak_left'] = self.df['histogram'].apply(lambda x: [int(j - 10) for j in peak_widths(x, find_peaks(x, distance=10)[0], rel_height=0.98)[2]])
        self.df['peak_right'] = self.df['histogram'].apply(lambda x: [int(j - 10) for j in peak_widths(x, find_peaks(x, distance=10)[0], rel_height=0.98)[3]])
        self.df['entropy'] = self.df['histogram'].apply(lambda x: stats.entropy(x))
        self.df['annotated'] = self.df['transcript_id'].apply(lambda x: self.transcripts[x]['annotated'])

        self.dfResult = pd.DataFrame()
        if 0 in list(self.df['peak_count']):
            self.dfResult = self.df.loc[self.df['peak_count'] == 0].copy()
            self.dfResult['prediction'] = self.dfResult['max']
            self.dfResult['peak_location'] = self.dfResult['max']
            self.dfResult = self.dfResult.drop('max', axis = 1)
            self.df = self.df.drop(self.df[self.df.peak_count==0].index, axis = 0).reset_index(drop=True)
            self.dfResult['peak_heights'] = self.dfResult.apply(lambda x: x.histogram[x.prediction+10], axis = 1)        


        keys = list(self.df.peak_info[0].keys())
        peaks = self.df.drop('max', axis = 1).copy()
        for i in keys:
            peaks[i] = peaks.peak_info.apply(lambda x: x[i])
        peaks = peaks.drop(['left_thresholds', 'right_thresholds'], axis = 1)
        peaks = peaks.apply(lambda x: self.sort_peaks(x), axis = 1)
        keys.remove('left_thresholds')
        keys.remove('right_thresholds')
        peaks = peaks.explode(keys+['peak_location', 'rank', 'peak_prominence', 'peak_width', 'peak_left', 'peak_right']).drop(['peak_info'], axis = 1).reset_index(drop=True)
        peaks[['mean', 'var', 'range', 'skew', 'peak_count', 'peak_location', 'entropy', 
            'peak_heights', 'peak_width', 'peak_prominence',
            'rank']] = peaks[['mean', 'var', 'range', 'skew',
            'peak_count', 'peak_location', 'entropy', 
            'peak_heights', 'peak_width', 'peak_prominence', 'rank']].astype(float, errors='ignore') 
        
    
        peaks['prediction'] = self.model.predict(peaks.drop(['peak_left', 'peak_right', 'histogram', 
                                                             'annotated', 'rank', 'transcript_id', 'gene_id', 'start', 'chromosome'], axis = 1).astype(float, errors='ignore'))
    
        peaks = peaks[peaks.prediction ==True].reset_index(drop=True)
        peaks['prediction'] = peaks['peak_location']
        self.dfResult = pd.concat([self.dfResult, peaks], axis=0).reset_index(drop=True)
        self.dfResult['counts'] = self.dfResult.apply(lambda x: self.counts(x), axis = 1)

        
        self.dfResult['prediction'] += self.dfResult['start']
        self.dfResult['prediction'] = self.dfResult['prediction'].astype(int)
        self.dfResult['flag'] = self.dfResult.apply(lambda x: self.flag(x), axis = 1)
        self.dfResult['peak_heights'] = self.dfResult['peak_heights'].astype(int)
        if 'peak_info' in self.dfResult:
            self.dfResult = self.dfResult.drop('peak_info', axis = 1)


        if not self.ignore_read_groups:
            self.dfResult['counts_groups'] = self.dfResult.apply(lambda x: self.counts_byGroup(x), axis = 1)
            self.dfResult['counts_byGroup'] = self.dfResult['counts_groups'].apply(lambda x: x[0])
            self.dfResult['peak_heights_byGroup'] = self.dfResult['counts_groups'].apply(lambda x: x[1])
            self.dfResult.drop('counts_groups', inplace=True, axis = 1)
            self.dfResult = self.dfResult[['chromosome', 'transcript_id', 'gene_id', 'prediction', 'counts', 'counts_byGroup', 'flag']]
            self.dfResult['group_id'] = [self.ordered_groups for _ in range(len(self.dfResult))]
            self.dfResult = self.dfResult.explode(['counts_byGroup', 'group_id']).reset_index(drop=True)
        else:
            self.dfResult = self.dfResult[['chromosome', 'transcript_id', 'gene_id', 'prediction', 'counts', 'flag']]

        
        self.transcripts = {}
        
        if self.first:
            self.dfResult.to_csv(self.output_prefix, sep="\t", index=False, mode="w", header=True)
        else:
            self.dfResult.to_csv(self.output_prefix, sep="\t", index=False, mode="a", header=False)
        
        self.first = False
    
    def finalize(self, args=None):
        pass    


    def flag(self, x):
        if abs(x.prediction - x.annotated) > 10:
            return 'Novel'
        return 'Known'



    def counts(self, x):
        low = int(10+x['peak_left'])
        high = int(10+x['peak_right']+1)
        return int(np.array(x.histogram[low:high]).sum())
    

    def counts_byGroup(self, x):
        counts = []
        heights = []
        for i in self.group_numeric_ids.values():
            hist = np.histogram(self.transcripts[x['transcript_id']][i], bins = 1+np.max(self.transcripts[x['transcript_id']]['data']) - np.min(self.transcripts[x['transcript_id']]['data']), range=(np.min(self.transcripts[x['transcript_id']]['data']), np.max(self.transcripts[x['transcript_id']]['data'])))[0]
            low = int(x['peak_left'])
            high = int(x['peak_right']+1)
            if low < 0:
                low = 0
            if high > len(hist):
                high = len(hist)
                
            counts.append(int(np.array(hist[low:high]).sum()))
            heights.append(int(hist[int(x['peak_location'])]))
        return (counts, heights)

    def sort_peaks(self, x):
        if x.peak_count > 1:
            sorted_indices = np.argsort(-x['peak_prominence'])  
            x['peak_location'] = np.array(x['peak_location'])[sorted_indices].tolist()
            x['peak_prominence'] = x['peak_prominence'][sorted_indices]
            x['peak_width'] = x['peak_width'][sorted_indices]
            x['peak_heights'] = x['peak_heights'][sorted_indices]
            x['peak_left'] = np.array(x['peak_left'])[sorted_indices].tolist()
            x['peak_right'] = np.array(x['peak_right'])[sorted_indices].tolist()
            x['rank'] = list(range(1, len(sorted_indices)+1))
        else:
            x['rank'] = 0
        return x


    def result(self, dff):
        dff['count'] = 1

        matrix = pd.pivot_table(dff, 
                            values='count', 
                            index='rank', 
                            columns=['prediction', 'true_peak'], 
                            aggfunc='count', 
                            fill_value=0)

        print(matrix)





