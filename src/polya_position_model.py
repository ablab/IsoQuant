import pandas as pd
import numpy as np
from scipy.signal import find_peaks, peak_prominences, peak_widths
from scipy import stats
from sklearn.metrics import confusion_matrix, f1_score, roc_auc_score, precision_recall_curve
from xgboost import XGBClassifier
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
    def __init__(self, output_prefix, ignore_read_groups=False):
        super().__init__(output_prefix, ignore_read_groups=False)
        self.transcripts = {}    
        self.gene_info = None
        with open('src/model.pkl', 'rb') as file:
            self.model = pickle.load(file)
        
    
    def add_read_info(self, read_assignment: ReadAssignment):
        if self.gene_info is None:
            self.gene_info: GeneInfo = read_assignment.gene_info
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
                    self.transcripts[isoform_match.assigned_transcript] = {'chr': read_assignment.chr_id, 'gene_id': isoform_match.assigned_gene, 'data': []}
                
                self.transcripts[isoform_match.assigned_transcript]['data'].append(int(polya_pos))

      

    def dump(self):
        # compute polyA peaks, filtering etc, save to self.output_file and clean features dataframe 

        self.df = pd.DataFrame({
            'transcript_id': self.transcripts.keys()     
        })

            
        self.df['chromosome'] = self.df['transcript_id'].apply(lambda x: self.transcripts[x]['chr'])
        self.df['gene_id'] = self.df['transcript_id'].apply(lambda x: self.transcripts[x]['gene_id'])
        self.df['start'] = self.df['transcript_id'].apply(lambda x: np.min(self.transcripts[x]['data']))
        self.df['max'] = self.df['transcript_id'].apply(lambda x: stats.mode(self.transcripts[x]['data'])[0])-self.df.start
        self.df['mean'] = self.df['transcript_id'].apply(lambda x: np.mean(self.transcripts[x]['data']))-self.df.start
        self.df['median'] = self.df['transcript_id'].apply(lambda x: np.median(self.transcripts[x]['data']))-self.df.start
        self.df['var'] = self.df['transcript_id'].apply(lambda x: np.var(self.transcripts[x]['data']))
        self.df['range'] = self.df['transcript_id'].apply(lambda x: np.max(self.transcripts[x]['data']) - np.min(self.transcripts[x]['data'])) + 1
        self.df['skew'] = self.df['transcript_id'].apply(lambda x: stats.skew(self.transcripts[x]['data']))
        self.df['peak_count'] = self.df['transcript_id'].apply(lambda x: len(find_peaks(np.histogram(self.transcripts[x]['data'], bins = 1+np.max(self.transcripts[x]['data']) - np.min(self.transcripts[x]['data']))[0], distance=10)[0]))
        self.df['peak_info'] = self.df['transcript_id'].apply(lambda x: find_peaks(np.histogram(self.transcripts[x]['data'], bins = 1+np.max(self.transcripts[x]['data']) - np.min(self.transcripts[x]['data']))[0], distance=10, threshold=(None, None), height=(None, None))[1])
        self.df['peak_location'] = self.df['transcript_id'].apply(lambda x: find_peaks(np.histogram(self.transcripts[x]['data'], bins = 1+np.max(self.transcripts[x]['data']) - np.min(self.transcripts[x]['data']))[0], distance=10)[0])
        self.df['peak_prominence'] = self.df['transcript_id'].apply(lambda x: peak_prominences(np.histogram(self.transcripts[x]['data'], bins = 1+np.max(self.transcripts[x]['data']) - np.min(self.transcripts[x]['data']))[0], find_peaks(np.histogram(self.transcripts[x]['data'], bins = 1+np.max(self.transcripts[x]['data']) - np.min(self.transcripts[x]['data']))[0], distance=10)[0])[0])
        self.df['peak_width'] = self.df['transcript_id'].apply(lambda x: peak_widths(np.histogram(self.transcripts[x]['data'], bins = 1+np.max(self.transcripts[x]['data']) - np.min(self.transcripts[x]['data']))[0], find_peaks(np.histogram(self.transcripts[x]['data'], bins = 1+np.max(self.transcripts[x]['data']) - np.min(self.transcripts[x]['data']))[0], distance=10)[0])[0])
        self.df['entropy'] = self.df['transcript_id'].apply(lambda x: stats.entropy(np.histogram(self.transcripts[x]['data'], bins = 1+np.max(self.transcripts[x]['data']) - np.min(self.transcripts[x]['data']))[0]))



        self.dfResult = self.df.loc[self.df['peak_count'] == 0].copy()
        self.dfResult['prediction'] = self.dfResult['max']
        self.dfResult['peak_location'] = self.dfResult['max']
        self.dfResult = self.dfResult.drop('max', axis = 1)
        self.df = self.df.drop(self.df[self.df.peak_count==0].index, axis = 0).reset_index(drop=True)

        self.dfResult['peak_heights'] = self.dfResult.apply(lambda x: np.histogram(self.transcripts[x['transcript_id']]['data'], bins = 1+np.max(self.transcripts[x['transcript_id']]['data']) - np.min(self.transcripts[x['transcript_id']]['data']))[0][x.prediction], axis = 1)

        


        keys = list(self.df.peak_info[0].keys())


        peaks = self.df.drop('max', axis = 1).copy()
        for i in keys:
            peaks[i] = peaks.peak_info.apply(lambda x: x[i])

        
        peaks = peaks.drop(['left_thresholds', 'right_thresholds'], axis = 1)
        peaks = peaks.apply(lambda x: self.sort_peaks(x), axis = 1)
        keys.remove('left_thresholds')
        keys.remove('right_thresholds')
        peaks = peaks.explode(keys+['peak_location', 'rank', 'peak_prominence', 'peak_width']).drop(['peak_info'], axis = 1).reset_index(drop=True)
        
        peaks[['mean', 'var', 'range', 'skew',
            'peak_count', 'peak_location', 'entropy', 
            'peak_heights', 'peak_width', 'peak_prominence',
            'rank']] = peaks[['mean', 'var', 'range', 'skew',
            'peak_count', 'peak_location', 'entropy', 
            'peak_heights', 'peak_width', 'peak_prominence', 'rank']].astype(float, errors='ignore') 

    
        peaks['prediction'] = self.model.predict(peaks.drop(['rank', 'transcript_id', 'gene_id', 'start', 'chromosome'], axis = 1).astype(float, errors='ignore'))
        peaks = peaks[peaks.prediction ==True].reset_index(drop=True)
        peaks['prediction'] = peaks['peak_location']
        self.dfResult = pd.concat([self.dfResult, peaks], axis=0).reset_index(drop=True)
        self.dfResult['counts'] = self.dfResult.apply(lambda x: self.test(x), axis = 1)
        self.dfResult = self.dfResult.drop(['median', 'mean', 'var', 'range', 'skew', 'peak_prominence', 'entropy', 'peak_count', 'peak_info', 'peak_width', 'peak_location', 'rank'], axis = 1)
    
        self.dfResult['prediction'] += self.dfResult['start']
        self.dfResult = self.dfResult.drop('start', axis = 1)


        # self.gene_info.all_isoforms_exons[transcript_id][-1][1] + 1

        self.transcripts = {}
        self.gene_info = None
        self.dfResult['prediction'] = self.dfResult['prediction'].astype(int)
        self.dfResult['peak_heights'] = self.dfResult['peak_heights'].astype(int)
        self.dfResult.to_csv(
        self.output_prefix,
        sep="\t",
        index=False,
        mode="a",
        header=not os.path.exists(self.output_prefix)
    )
    
    def finalize(self, args=None):
        pass    


    def test(self, x):
        low = int(x['peak_location']-5)
        high = int(x['peak_location']+6)
        data = np.histogram(self.transcripts[x['transcript_id']]['data'], bins = 1+np.max(self.transcripts[x['transcript_id']]['data']) - np.min(self.transcripts[x['transcript_id']]['data']))[0]
        if low < 0:
            low = 0
        if high > len(data):
            high = len(data)
            
        return data[low:high].sum()
    

    def sort_peaks(self, x):
        if x.peak_count > 1:
            sorted_indices = np.argsort(-x['peak_prominence'])  
            x['peak_location'] = x['peak_location'][sorted_indices]
            x['peak_prominence'] = x['peak_prominence'][sorted_indices]
            x['peak_width'] = x['peak_width'][sorted_indices]
            x['peak_heights'] = x['peak_heights'][sorted_indices]
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





