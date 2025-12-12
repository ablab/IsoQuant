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
        self.transcripts_additional = {}   
        self.gene_info = None

        # get read groups

        self.first = True

        self.group_numeric_ids = {}
        self.ordered_groups = sorted(self.read_groups)
        if self.ordered_groups:
            for i, g in enumerate(self.ordered_groups):
                self.group_numeric_ids[g] = i


    def result(self, dff):
        dff['count'] = 1

        matrix = pd.pivot_table(dff, 
                            values='count', 
                            index='rank', 
                            columns=['prediction', 'true_peak'], 
                            aggfunc='count', 
                            fill_value=0)

        return matrix

    def peak_target(self, x):
        ok = False
        delta = 10
        for i in x['target']:
            ok = ok or ((x['peak_location'] >= i - delta) & (x['peak_location'] <= i + delta))
        return ok

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
                start_pos = read_assignment.corrected_exons[-1][1]
                ok = True
            elif read_assignment.strand == '+' and read_assignment.polya_info.external_polya_pos != -1:
                polya_pos = read_assignment.polya_info.external_polya_pos
                start_pos = read_assignment.corrected_exons[0][0]
                ok = True
            if ok:                
                isoform_match: IsoformMatch = read_assignment.isoform_matches[0]
                if isoform_match.assigned_transcript not in self.transcripts:

                    if read_assignment.strand == '+':
                        annot_polya = read_assignment.gene_info.all_isoforms_exons[isoform_match.assigned_transcript][-1][1] + 1
                        annot_start = read_assignment.gene_info.all_isoforms_exons[isoform_match.assigned_transcript][0][0] - 1
                    else:
                        annot_polya = read_assignment.gene_info.all_isoforms_exons[isoform_match.assigned_transcript][0][0] - 1
                        annot_start = read_assignment.gene_info.all_isoforms_exons[isoform_match.assigned_transcript][-1][1] + 1
                    self.transcripts[isoform_match.assigned_transcript] = {'chr': read_assignment.chr_id, 'gene_id': isoform_match.assigned_gene, 'data_polya': [], 'annotated_polya': annot_polya, 'data_start': [], 'annotated_start': annot_start}
                    
                    
                    if not self.ignore_read_groups:
                        for i in self.group_numeric_ids.values():
                            self.transcripts[isoform_match.assigned_transcript][i] = {}
                            self.transcripts[isoform_match.assigned_transcript][i]['data_polya'] = []
                            self.transcripts[isoform_match.assigned_transcript][i]['data_start'] = []
                self.transcripts[isoform_match.assigned_transcript]['data_start'].append(int(start_pos))
                self.transcripts[isoform_match.assigned_transcript]['data_polya'].append(int(polya_pos))

                if not self.ignore_read_groups:
                    self.transcripts[isoform_match.assigned_transcript][group_id]['data_start'].append(int(start_pos))
                    self.transcripts[isoform_match.assigned_transcript][group_id]['data_polya'].append(int(polya_pos))
                
        
        
        elif read_assignment.polyA_found == False and read_assignment.assignment_type in [ReadAssignmentType.inconsistent_non_intronic, ReadAssignmentType.unique, ReadAssignmentType.unique_minor_difference, ReadAssignmentType.inconsistent]:
            isoform_match: IsoformMatch = read_assignment.isoform_matches[0]
            if read_assignment.strand == '-':
                start_pos = read_assignment.corrected_exons[-1][1]
            else:
                start_pos = read_assignment.corrected_exons[0][0]
            
            if isoform_match.assigned_transcript not in self.transcripts_additional:
                if read_assignment.strand == '+':
                    annot_start = read_assignment.gene_info.all_isoforms_exons[isoform_match.assigned_transcript][0][0] - 1
                else:
                    annot_start = read_assignment.gene_info.all_isoforms_exons[isoform_match.assigned_transcript][-1][1] + 1
                self.transcripts_additional[isoform_match.assigned_transcript] = {'chr': read_assignment.chr_id, 'gene_id': isoform_match.assigned_gene, 'data_start': [], 'annotated_start': annot_start}
                if not self.ignore_read_groups:
                    for i in self.group_numeric_ids.values():
                        self.transcripts_additional[isoform_match.assigned_transcript][i] = {}
                        self.transcripts_additional[isoform_match.assigned_transcript][i]['data_start'] = []
            
            self.transcripts_additional[isoform_match.assigned_transcript]['data_start'].append(int(start_pos))
            if not self.ignore_read_groups:
                self.transcripts_additional[isoform_match.assigned_transcript][group_id]['data_start'].append(int(start_pos))
    

    def collect_data(self):

        self.peaks['target'] = self.peaks['transcript_id'].apply(lambda x: [self.transcripts[x]['annotated_start']])

        for id in self.transcripts:
            if "_" in id:
                if self.peaks['transcript_id'].eq(id.split("_")[0]).any():
                    self.peaks[self.peaks['transcript_id']==id.split("_")[0]]['target'] += [self.transcripts[id]['annotated_start']]


        self.peaks['true_peak'] = self.peaks.apply(lambda x: self.peak_target(x) , axis = 1)

        self.peaks = self.peaks.drop('target', axis = 1)

        
        if self.first:
            self.peaks.to_csv("src/model_df.csv", index=False, header=True)
        else:
            self.peaks.to_csv("src/model_df.csv", mode="a", index=False, header=False)

        self.peaks = self.peaks.drop('true_peak', axis = 1)
        




    def train_model(self):
        self.model = XGBClassifier()
        self.model.load_model('src/model.json')
    
        params = self.model.get_xgb_params()
        self.model_start = XGBClassifier(**params)

        self.model_df = pd.read_csv("src/model_df.csv")

        X = self.model_df.copy().astype(float, errors='ignore').drop('true_peak', axis = 1)
        y = self.model_df['true_peak'].copy()
        
        chrs = list(X.chromosome.unique())
        n = len(chrs)//2
        test_chrs = random.sample(chrs, n)
        train_chrs = [i for i in chrs if i not in test_chrs]

        X_train = X[X["chromosome"].isin(train_chrs)].drop('chromosome', axis = 1).copy().reset_index(drop=True)
        y_train = y[X["chromosome"].isin(train_chrs)].copy().reset_index(drop=True)

        X_test = X[X["chromosome"].isin(test_chrs)].drop('chromosome', axis = 1).copy().reset_index(drop=True)
        y_test = y[X["chromosome"].isin(test_chrs)].copy().reset_index(drop=True)


        self.model_start.fit(X_train.drop(['peak_left', 'peak_right', 'histogram', 'annotated', 'rank', 'transcript_id', 'gene_id', 'start'], axis = 1).astype(float, errors='ignore'), y_train)
        
        
        rank = X_test['rank']
        X_test.drop('rank', axis = 1, inplace=True)
        X_test['prediction'] = self.model_start.predict(X_test.drop(['peak_left', 'peak_right', 'histogram', 'annotated', 'transcript_id', 'gene_id', 'start'], axis = 1).astype(float, errors='ignore'))

        X_test['rank'] = rank
        X_test['true_peak'] = y_test
        self.result(X_test)
        self.model_start.save_model('src/model_start.json')



    def create_df(self, data, annotated, file_name):
        self.df = pd.DataFrame({'transcript_id': self.transcripts.keys()})
        self.df['chromosome'] = self.df['transcript_id'].apply(lambda x: self.transcripts[x]['chr'])
        self.df['gene_id'] = self.df['transcript_id'].apply(lambda x: self.transcripts[x]['gene_id'])
        self.df['start'] = self.df['transcript_id'].apply(lambda x: np.min(self.transcripts[x][data]))
        self.df['max'] = self.df['transcript_id'].apply(lambda x: stats.mode(self.transcripts[x][data])[0])-self.df.start
        self.df['mean'] = self.df['transcript_id'].apply(lambda x: np.mean(self.transcripts[x][data]))-self.df.start
        self.df['histogram'] = self.df['transcript_id'].apply(lambda x: [0]*10 + list(np.histogram(self.transcripts[x][data], bins = 1+max(self.transcripts[x][data]) - min(self.transcripts[x][data]))[0])+[0]*10)
        self.df['mean_height'] = self.df['histogram'].apply(lambda x: np.mean(x))
        self.df['var'] = self.df['transcript_id'].apply(lambda x: np.var(self.transcripts[x][data]))
        self.df['range'] = self.df['transcript_id'].apply(lambda x: np.max(self.transcripts[x][data]) - np.min(self.transcripts[x][data])) + 1
        self.df['skew'] = self.df['transcript_id'].apply(lambda x: stats.skew(self.transcripts[x][data]))
        self.df['peak_count'] = self.df['histogram'].apply(lambda x: len(find_peaks(x, distance=10)[0]))
        self.df['peak_info'] = self.df['histogram'].apply(lambda x: find_peaks(x, distance=10, threshold=(None, None), height=(None, None))[1])
        self.df['peak_location'] = self.df['histogram'].apply(lambda x: [int(j - 10) for j in find_peaks(x, distance=10)[0]])
        self.df['peak_prominence'] = self.df['histogram'].apply(lambda x: peak_prominences(x, find_peaks(x, distance=10)[0])[0])
        self.df['peak_width'] = self.df['histogram'].apply(lambda x: peak_widths(x, find_peaks(x, distance=10)[0], rel_height=0.98)[0])
        self.df['peak_left'] = self.df['histogram'].apply(lambda x: [int(j - 10) for j in peak_widths(x, find_peaks(x, distance=10)[0], rel_height=0.98)[2]])
        self.df['peak_right'] = self.df['histogram'].apply(lambda x: [int(j - 10) for j in peak_widths(x, find_peaks(x, distance=10)[0], rel_height=0.98)[3]])
        self.df['entropy'] = self.df['histogram'].apply(lambda x: stats.entropy(x))
        self.df['annotated'] = self.df['transcript_id'].apply(lambda x: self.transcripts[x][annotated])

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
    
        self.peaks = self.peaks[self.peaks.prediction ==True].reset_index(drop=True)
        self.peaks['prediction'] = self.peaks['peak_location']
        self.dfResult = pd.concat([self.dfResult, self.peaks], axis=0).reset_index(drop=True)
        
        self.dfResult['counts'] = self.dfResult.apply(lambda x: self.counts(x), axis = 1)
        self.dfResult['prediction'] += self.dfResult['start']
        self.dfResult['prediction'] = self.dfResult['prediction'].astype(int)
        self.dfResult['flag'] = self.dfResult.apply(lambda x: self.flag(x), axis = 1)
        self.dfResult['peak_heights'] = self.dfResult['peak_heights'].astype(int)
        if 'peak_info' in self.dfResult:
            self.dfResult = self.dfResult.drop('peak_info', axis = 1)


        if not self.ignore_read_groups:
            self.dfResult['counts_groups'] = self.dfResult.apply(lambda x: self.counts_byGroup(x, data), axis = 1)
            self.dfResult['counts_byGroup'] = self.dfResult['counts_groups'].apply(lambda x: x[0])
            self.dfResult['peak_heights_byGroup'] = self.dfResult['counts_groups'].apply(lambda x: x[1])
            self.dfResult.drop('counts_groups', inplace=True, axis = 1)
            self.dfResult = self.dfResult[['chromosome', 'transcript_id', 'gene_id', 'prediction', 'counts', 'counts_byGroup', 'flag']]
            self.dfResult['group_id'] = [self.ordered_groups for _ in range(len(self.dfResult))]
            self.dfResult = self.dfResult.explode(['counts_byGroup', 'group_id']).reset_index(drop=True)
        else:
            self.dfResult = self.dfResult[['chromosome', 'transcript_id', 'gene_id', 'prediction', 'counts', 'flag']]

        
        
        if self.first:
            self.dfResult.to_csv(file_name, sep="\t", index=False, mode="w", header=True)
        else:
            self.dfResult.to_csv(file_name, sep="\t", index=False, mode="a", header=False)
    



    def dump(self):
        # self.model = XGBClassifier()
        # self.model.load_model('src/model.json')

        # self.create_df('data_polya', 'annotated_polya', self.output_prefix)

        # for i in self.transcripts_additional:
        #     if i in self.transcripts:
        #         self.transcripts_additional[i]['data_start'].extend(self.transcripts[i]['data_start'])
        #         if not self.ignore_read_groups:
        #             for j in self.group_numeric_ids.values():
        #                 self.transcripts_additional[i][j]['data_start'].extend(self.transcripts[i][j]['data_start'])



        # self.transcripts = self.transcripts | self.transcripts_additional

        self.model = XGBClassifier()
        self.model.load_model('src/model_start.json')

        self.create_df('data_start', 'annotated_start', self.output_prefix+"_start")
                    
        
        self.transcripts = {}
        self.first = False


        
        
    
    def finalize(self, args=None):
        self.train_model()
        pd.set_option('display.max_columns', None)
        pd.set_option('display.width', None) 
        print(self.model_df)


    def flag(self, x):
        if abs(x.prediction - x.annotated) > 10:
            return 'Novel'
        return 'Known'



    def counts(self, x):
        
        low = int(10+x['peak_left'])
        high = int(10+x['peak_right']+1)
        
        return int(np.array(x.histogram[low:high]).sum())
    

    def counts_byGroup(self, x, data):
        counts = []
        heights = []
        for i in self.group_numeric_ids.values():
            hist = np.histogram(self.transcripts[x['transcript_id']][i][data], bins = 1+np.max(self.transcripts[x['transcript_id']][data]) - np.min(self.transcripts[x['transcript_id']][data]), range=(np.min(self.transcripts[x['transcript_id']][data]), np.max(self.transcripts[x['transcript_id']][data])))[0]
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





