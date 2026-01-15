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
import random
import sys

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



class TerminalCounter(AbstractCounter):
    def __init__(self, output_prefix, read_groups, ignore_read_groups=False):
        super().__init__(output_prefix, ignore_read_groups=False)
        self.ignore_read_groups = ignore_read_groups
        self.read_groups = read_groups
        self.transcripts = {} 
        self.transcripts_additional = {}   
        self.gene_info = None
        self.first = True
        self.group_numeric_ids = {}
        self.ordered_groups = sorted(self.read_groups)
        if self.ordered_groups:
            for i, g in enumerate(self.ordered_groups):
                self.group_numeric_ids[g] = i

    def peak_target(self, x):
        ok = False
        delta = 10
        for i in x['target']:
            ok = ok or ((x['peak_location'] >= i - x['start'] - delta) & (x['peak_location'] <= i - x['start'] + delta))
        return ok

    def flag(self, x):
        if abs(x.prediction - x.annotated) > 10:
            return 'Novel'
        return 'Known'

    def counts(self, x):
        
        low = int(10+x['peak_left'])
        high = int(10+x['peak_right']+1)
        
        return int(np.array(x['histogram'])[low:high].sum())
    
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

    def collect_data(self):
        if self.transcripts == {}:
            self.df = pd.DataFrame({'chromosome':[], 'transcript_id':[], 'gene_id':[], 'prediction':[], 'counts':[], 'counts_byGroup':[], 'flag':[], 'group_id':[]})
            if self.first:
                self.df.to_csv(self.output_prefix, sep="\t", index=False, mode="w", header=True)
            else:
                self.df.to_csv(self.output_prefix, sep="\t", index=False, mode="a", header=False)
            return
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
        self.df['skew'] = self.df['transcript_id'].apply(lambda x: stats.skew(self.transcripts[x]['data']) if (len(self.transcripts[x]['data']) > 3 and np.var(self.transcripts[x]['data']) >= 1e-12) else None)
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
        self.peaks = self.df.drop('max', axis = 1).copy()
        for i in keys:
            self.peaks[i] = self.peaks.peak_info.apply(lambda x: x[i])
        self.peaks = self.peaks.drop(['left_thresholds', 'right_thresholds'], axis = 1)
        self.peaks = self.peaks.apply(lambda x: self.sort_peaks(x), axis = 1)
        keys.remove('left_thresholds')
        keys.remove('right_thresholds')
        self.peaks = self.peaks.explode(keys+['peak_location', 'rank', 'peak_prominence', 'peak_width', 'peak_left', 'peak_right']).drop(['peak_info'], axis = 1).reset_index(drop=True)
             
            
        self.peaks[['mean', 'var', 'range', 'skew', 'peak_count', 'peak_location', 'entropy', 
            'peak_heights', 'peak_width', 'peak_prominence',
            'rank']] = self.peaks[['mean', 'var', 'range', 'skew',
            'peak_count', 'peak_location', 'entropy', 
            'peak_heights', 'peak_width', 'peak_prominence', 'rank']].astype(float, errors='ignore') 
        
        self.peaks['target'] = self.peaks['transcript_id'].apply(lambda x: [self.transcripts[x]['annotated']])

        k = 0
        for id in self.transcripts:
            if "_" in id:
                k += 1
                if self.peaks['transcript_id'].eq(id.split("_")[0]).any():
                    self.peaks[self.peaks['transcript_id']==id.split("_")[0]]['target'] += [self.transcripts[id]['annotated']]

        print(k)

        self.peaks['true_peak'] = self.peaks.apply(lambda x: self.peak_target(x) , axis = 1)

        # self.peaks = self.peaks.drop('target', axis = 1)

        if self.first:
            self.peaks.to_csv("src/model_df.csv", index=False, header=True)
        else:
            self.peaks.to_csv("src/model_df.csv", mode="a", index=False, header=False)
        
        self.peaks = self.peaks.drop('true_peak', axis = 1)
        
    def train_model(self):
        self.model = XGBClassifier()
        self.model.load_model('src/model.json')
        params = self.model.get_xgb_params()
        # {'booster': 'gbtree', 'learning_rate': 0.1711139976081603, 'gamma': 3.213024131696218, 'max_depth': 8, 'min_child_weight': 7.031066410000928, 'reg_lambda': 0.00035253812038082424, 'reg_alpha': 0.010639993402868728, 'subsample': 0.6310944088222122, 'colsample_bytree': 0.7856502544324856, 'n_estimators': 56, 'scale_pos_weight': 0.3971102549981997}
        # {'booster': 'dart', 'learning_rate': 0.09115397231222125, 'gamma': 0.7911797085921515, 'max_depth': 6, 'min_child_weight': 8.890261356320382, 'reg_lambda': 3.4438281320977343e-07, 'reg_alpha': 4.583918034010515e-06, 'subsample': 0.9139332107255594, 'colsample_bytree': 0.9804235433943478, 'n_estimators': 95, 'scale_pos_weight': 0.41416529852133854}
        params = {'booster': 'dart', 'learning_rate': 0.25989024187450366, 'gamma': 4.087502269121285, 'max_depth': 6, 'min_child_weight': 5.722951381412213, 'reg_lambda': 4.84758559371017, 'reg_alpha': 0.03795296107941127, 'subsample': 0.8809954546952385, 'colsample_bytree': 0.7483035972922923, 'n_estimators': 111, 'scale_pos_weight': 0.297061977156673}
        self.model_new = XGBClassifier(**params)
        
        self.model_df = pd.read_csv("src/model_df.csv").drop('target', axis = 1)

        X = self.model_df.drop('true_peak', axis = 1).copy().astype(float, errors='ignore') 
        y = self.model_df['true_peak'].copy()
        
        chrs = list(X.chromosome.unique())
        n = len(chrs)//2
        test_chrs = random.sample(chrs, n)
        train_chrs = [i for i in chrs if i not in test_chrs]

        X_train = X[X["chromosome"].isin(train_chrs)].drop('chromosome', axis = 1).copy().reset_index(drop=True)
        y_train = y[X["chromosome"].isin(train_chrs)].copy().reset_index(drop=True)

        X_test = X[X["chromosome"].isin(test_chrs)].drop('chromosome', axis = 1).copy().reset_index(drop=True)
        y_test = y[X["chromosome"].isin(test_chrs)].copy().reset_index(drop=True)

        self.model_new.fit(X_train.drop(['peak_left', 'peak_right', 'histogram','annotated', 'rank', 'transcript_id', 'gene_id', 'start'], axis = 1).astype(float, errors='ignore'), y_train)
        
        X_test['prediction'] = self.model_new.predict(X_test.drop(['peak_left', 'peak_right', 'histogram', 'annotated', 'rank', 'transcript_id', 'gene_id', 'start'], axis = 1).astype(float, errors='ignore'))

        X_test['true_peak'] = y_test
        self.result(X_test)
        self.model_new.save_model('src/model_new.json')

    def create_df(self):
        if self.transcripts == {}:
            self.df = pd.DataFrame({'chromosome':[], 'transcript_id':[], 'gene_id':[], 'prediction':[], 'counts':[], 'counts_byGroup':[], 'flag':[], 'group_id':[]})
            if self.first:
                self.df.to_csv(self.output_prefix, sep="\t", index=False, mode="w", header=True)
            else:
                self.df.to_csv(self.output_prefix, sep="\t", index=False, mode="a", header=False)
            return
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
        self.df['skew'] = self.df['transcript_id'].apply(lambda x: stats.skew(self.transcripts[x]['data']) if (len(self.transcripts[x]['data']) > 3 and np.var(self.transcripts[x]['data']) >= 1e-12) else None)
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
        self.peaks = self.df.drop('max', axis = 1).copy()


        for i in keys:
            self.peaks[i] = self.peaks.peak_info.apply(lambda x: x[i])
        self.peaks = self.peaks.drop(['left_thresholds', 'right_thresholds'], axis = 1)
        self.peaks = self.peaks.apply(lambda x: self.sort_peaks(x), axis = 1)



        


        keys.remove('left_thresholds')
        keys.remove('right_thresholds')
        self.peaks = self.peaks.explode(keys+['peak_location', 'rank', 'peak_prominence', 'peak_width', 'peak_left', 'peak_right']).drop(['peak_info'], axis = 1).reset_index(drop=True)
             


        self.peaks[['mean', 'var', 'range', 'skew', 'peak_count', 'peak_location', 'entropy', 
            'peak_heights', 'peak_width', 'peak_prominence',
            'rank']] = self.peaks[['mean', 'var', 'range', 'skew',
            'peak_count', 'peak_location', 'entropy', 
            'peak_heights', 'peak_width', 'peak_prominence', 'rank']].astype(float, errors='ignore') 
        
        

        self.peaks['prediction'] = self.model.predict(self.peaks.drop(['peak_left', 'peak_right', 'histogram', 
            'annotated', 'rank', 'transcript_id', 'gene_id', 'start', 'chromosome'], axis = 1).astype(float, errors='ignore'))
        
        self.peaks['prediction'] = self.peaks['prediction'].apply(lambda x: True if x == 1 else False)
        if self.first:
            self.peaks.to_csv("src/peaks_1.csv", sep="\t", index=False, mode="w", header=True)
        else:
            self.peaks.to_csv("src/peaks_1.csv", sep="\t", index=False, mode="a", header=False)


        
        # self.peaks = self.peaks[self.peaks.prediction ==True].reset_index(drop=True) 
        # 
        self.peaks = self.peaks[self.peaks.peak_location > 40].reset_index(drop=True)   
        
        
        
        
        
        self.peaks['prediction'] = self.peaks['peak_location']
        self.dfResult = pd.concat([self.dfResult, self.peaks], axis=0).reset_index(drop=True)

        
        if self.dfResult.empty:
            self.df = pd.DataFrame({'chromosome':[], 'transcript_id':[], 'gene_id':[], 'prediction':[], 'counts':[], 'counts_byGroup':[], 'flag':[], 'group_id':[]})
            if self.first:
                self.df.to_csv(self.output_prefix, sep="\t", index=False, mode="w", header=True)
            else:
                self.df.to_csv(self.output_prefix, sep="\t", index=False, mode="a", header=False)
            return

        print(len(self.dfResult))
        


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

        
        
        if self.first:
            self.dfResult.to_csv(self.output_prefix, sep="\t", index=False, mode="w", header=True)
        else:
            self.dfResult.to_csv(self.output_prefix, sep="\t", index=False, mode="a", header=False)

    def test_model(self):
        self.model = XGBClassifier()
        self.model.load_model('src/model.json')
        params = self.model.get_xgb_params()
        params = {'booster': 'dart', 'learning_rate': 0.09115397231222125, 'gamma': 0.7911797085921515, 'max_depth': 6, 'min_child_weight': 8.890261356320382, 'reg_lambda': 3.4438281320977343e-07, 'reg_alpha': 4.583918034010515e-06, 'subsample': 0.9139332107255594, 'colsample_bytree': 0.9804235433943478, 'n_estimators': 95, 'scale_pos_weight': 0.41416529852133854}
        self.model_new = XGBClassifier(**params)

        self.model_df = pd.read_csv("src/model_df.csv")

        X = self.model_df.drop('true_peak', axis = 1).copy().astype(float, errors='ignore') 
        y = self.model_df['true_peak'].copy()
        
        chrs = list(X.chromosome.unique())
        n = len(chrs)//2
        test_chrs = random.sample(chrs, n)

        X_test = X[X["chromosome"].isin(test_chrs)].drop('chromosome', axis = 1).copy().reset_index(drop=True)
        y_test = y[X["chromosome"].isin(test_chrs)].copy().reset_index(drop=True)

        
        X_test['prediction'] = self.model.predict(X_test.drop(['peak_left', 'peak_right', 'annotated', 'rank', 'transcript_id', 'gene_id', 'start'], axis = 1).astype(float, errors='ignore'))

        X_test['true_peak'] = y_test
        self.result(X_test)
        


class TSSCounter(TerminalCounter):
    def __init__(self, output_prefix, read_groups, ignore_read_groups=False):
        super().__init__(output_prefix, read_groups, ignore_read_groups=False)
        
    def add_read_info(self, read_assignment: ReadAssignment):
        if self.gene_info is None:
            self.gene_info: GeneInfo = read_assignment.gene_info
        group_id = read_assignment.read_group
        group_id = self.group_numeric_ids[group_id]
        # add a single read_assignment to features dataframe

        if read_assignment.assignment_type in [ReadAssignmentType.inconsistent_non_intronic, ReadAssignmentType.unique, ReadAssignmentType.unique_minor_difference, ReadAssignmentType.inconsistent]:
            ok = False
            if read_assignment.strand == '-':
                start_pos = read_assignment.corrected_exons[-1][1]
                ok = True
            elif read_assignment.strand == '+':
                start_pos = read_assignment.corrected_exons[0][0]
                ok = True
            
            if ok:
                isoform_match: IsoformMatch = read_assignment.isoform_matches[0]
                if isoform_match.assigned_transcript not in self.transcripts:

                    if read_assignment.strand == '+':
                        annot_start = read_assignment.gene_info.all_isoforms_exons[isoform_match.assigned_transcript][0][0] - 1
                    else:
                        annot_start = read_assignment.gene_info.all_isoforms_exons[isoform_match.assigned_transcript][-1][1] + 1
                    self.transcripts[isoform_match.assigned_transcript] = {'chr': read_assignment.chr_id, 'gene_id': isoform_match.assigned_gene, 'data': [], 'annotated': annot_start}
                    
                    
                    if not self.ignore_read_groups:
                        for i in self.group_numeric_ids.values():
                            self.transcripts[isoform_match.assigned_transcript][i] = []
                self.transcripts[isoform_match.assigned_transcript]['data'].append(int(start_pos))

                if not self.ignore_read_groups:
                    self.transcripts[isoform_match.assigned_transcript][group_id].append(int(start_pos))
        
    def dump(self):
        self.model = XGBClassifier()
        self.model.load_model('src/model_new.json')

        self.create_df()
        # self.collect_data()

        self.transcripts = {}
        self.first = False

    def finalize(self, args=None):
        # self.train_model()
        pass    


class PolyACounter(TerminalCounter):
    def __init__(self, output_prefix, read_groups, ignore_read_groups=False):
        super().__init__(output_prefix, read_groups, ignore_read_groups=False)

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
                        annot_polya = read_assignment.gene_info.all_isoforms_exons[isoform_match.assigned_transcript][-1][1] + 1
                    else:
                        annot_polya = read_assignment.gene_info.all_isoforms_exons[isoform_match.assigned_transcript][0][0] - 1
                    self.transcripts[isoform_match.assigned_transcript] = {'chr': read_assignment.chr_id, 'gene_id': isoform_match.assigned_gene, 'data': [], 'annotated': annot_polya}
                    
                    if not self.ignore_read_groups:
                        for i in self.group_numeric_ids.values():
                            self.transcripts[isoform_match.assigned_transcript][i] = []
                self.transcripts[isoform_match.assigned_transcript]['data'].append(int(polya_pos))

                if not self.ignore_read_groups:
                    self.transcripts[isoform_match.assigned_transcript][group_id].append(int(polya_pos))
     
    def dump(self):
        self.model = XGBClassifier()
        self.model.load_model('src/model.json')
        # self.model.load_model('src/model_new.json')


        # self.create_df()
        # self.collect_data()

        self.transcripts = {}
        self.first = False

    def finalize(self, args=None):
        # self.test_model()
        # self.train_model()
        pass

# to train a model: 
# in the dump function of the wanted class comment self.create_df() and uncomment self.collect_data() 
# and in the finalize function uncomment self.train_model()