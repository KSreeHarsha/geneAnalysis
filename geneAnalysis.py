#!/usr/bin/env python

"""
@file geneAnalysis.py
@author K Sree Harsha
@brief Gene analysis for various cell types
"""

#import modules
import sys
from pandas import DataFrame
import pandas as pd
import numpy as np
import itertools
import optparse




class GeneAnalysis(object):


    def create_writer(self):
        ''' Creates result file for given dataset file
        '''
        self.resfile='result_'+self.filename
        self.file = open('./'+self.resfile, 'w')
        self.file.write('@author K Sree Harsha\n')
        self.file.write('Gene analysis results\n\n')

    def get_cell_lines(self):
        ''' Returns the cell types in the given data set with time points '''
        cell_lines ={}
        cell_lines['HL60']= [0, .5, 4, 24]
        cell_lines['U937']=[0, .5, 4, 24]
        cell_lines['NB4']=[0, 5.5, 24, 48,72]
        cell_lines['Jurkat']=[0, .5, 4, 24]
        return cell_lines

    def import_file(self,filename, sep='\t', header=0):
        ''' import data from a file and return dataframe'''
       
        #Verify input parameters
        if not isinstance(filename, str):
            raise ValueError('Invalid filename.')
        if not isinstance(sep, str):
            raise ValueError('Invalid sep param.')
        if not isinstance(header, int):
            raise ValueError('Invalid header param.')
        self.filename=filename
        self.dataframe=pd.read_csv(filename, sep=sep, header=header)
        self.create_writer()
    
    def writer(self, string):
        ''' Write results to file'''
        #Verify input parameters
        if not isinstance(string, str):
            raise ValueError('Invalid string param')

        self.file.write('--------------------------\n')
        self.file.write(string)
        self.file.write('\n\n')



    def get_unique_genes(self):
        ''' Computes the number of distinct gene acession numbers in the dataset. '''

        cnt_unique_genes=len(set(self.dataframe['Gene Accession Number']))

        result=('Total number of distinct Gene Accession Numbers '
                   'represented in the data set are {}'.format(cnt_unique_genes))
        self.writer(result)

    def get_time_corr(self):
        ''' Compute highly correlated time points for all cell lines'''

        cell_line_time_corr={}
        cell_lines=self.get_cell_lines()

        for cell in cell_lines:
            cell_keys = ['{}_{}_hrs'.format(cell, time) for time in cell_lines[cell]]

            time_pairs=list(itertools.combinations(cell_keys, 2))
            time_corrs = dict(((t1.split('_')[1], t2.split('_')[1]),
                         abs(self.dataframe[t1].corr(self.dataframe[t2])))
						 for t1, t2 in time_pairs)

            max_corr_timepoints = max(time_corrs, key=time_corrs.get)
            cell_line_time_corr[cell]=max_corr_timepoints

        result=[]

        for cell_line,time_corrs in cell_line_time_corr.items():
            time1,time2= time_corrs
            result.append('Time points {} hr and {} hr are highly correlated for Cell line {} '.format(
                              time1, time2,cell_line))

        self.writer('\n'.join(result))


    def get_cell_corr(self):
        ''' Compute correlation between cell_lines. This function creates single vector of all time
              points for a cell line to compute correlation between  cell lines. First 4 time points for
              all cell lines are taken into account for computing correlation between cells. Note: Fifth
 	      time point (72 hrs) for cell NB4 is not taken
        '''

        cell_lines=self.get_cell_lines()
        cell_line_timepoints={}

        for cell in cell_lines:
            cell_keys = ['{}_{}_hrs'.format(cell, time) for time in cell_lines[cell] if time!=72]
            cell_line_timepoints[cell]=self.dataframe[cell_keys].unstack().values


        cell_pairs = list(itertools.combinations(cell_lines.keys(), 2))

        cell_corrs = dict(((cell1, cell2),
					abs(np.corrcoef(cell_line_timepoints[cell1],cell_line_timepoints[cell2]))[0, 1])
					for cell1, cell2 in cell_pairs)

        cell_line1, cell_line2 = max(cell_corrs, key=cell_corrs.get)
        result=('Most correlated cell lines are {} and {}'.format(cell_line1,cell_line2))
        self.writer(result)

    def done(self):
        print 'Results written to '+self.resfile

if __name__ == '__main__':

    ''' run the program as python geneAnalysis.py 'data_set_HL60_U937_NB4_Jurkat.txt' '''
    # Parse command line arguments and options.
    usage = 'usage: %prog fileName'
    description = 'Gene analysis for cell lines, results are written to result_dataset--.txt'
    p = optparse.OptionParser(usage=usage, description=description)

    (opts, args) = p.parse_args()

    if len(args) < 1:
        p.print_usage()
        sys.exit(1)

    filename=args[0]

    analyzer = GeneAnalysis()
    analyzer.import_file(filename)
    analyzer.get_unique_genes()
    analyzer.get_time_corr()
    analyzer.get_cell_corr()
    analyzer.done()


		

