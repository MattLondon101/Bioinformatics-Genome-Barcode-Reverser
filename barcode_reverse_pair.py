#!/usr/bin/env python

# Required imports
import os
import sys
import subprocess
import time
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import collections
from collections import Counter
from itertools import islice
import csv
import numpy as np
import pandas as pd
import psutil

# Set timer to output program execution time
starttime=time.time()
process = psutil.Process(os.getpid())

# Initiate class that takes user input
class code_pairs():
    def __init__(self,fastq_file,start_position,stop_position):
        self.a=fastq_file
        self.b=start_position
        self.c=stop_position
        
# Convert .fastq / .fq file to list
# O(n**2)
    def fastq_to_list(self):
        try:
            fn = self.a
        except IndexError as ie:
            raise SystemError("Error: Specify file name\n")
        if not os.path.exists(fn):
            raise SystemError("Error: File does not exist\n")
        
        start = int(self.b)
        stop = int(self.c)
        n = 4
        liner = []
        with open(fn, 'r') as fh:
            lines = []
            for line in fh:
                lines.append(line.rstrip())
                if len(lines) == n:
                    newlines = islice(lines, 1, None, 4)
                    liny=list(newlines)
                    liner.append(liny)
                    lines = []
        return (liner)

# Extract barcodes and create reverse complements (rev_comp)
# O(n**2)
    def list_to_barcodes(self, listt=None):
        barcodes = []
        rev_comps = []
        start = int(self.b)
        stop = int(self.c)
        listt = self.fastq_to_list()
        for i in listt:
            for j in i:
                gen=j[start:stop]
                s=Seq(gen,IUPAC.unambiguous_dna)
                barcodes.append(gen)
                rev=s.reverse_complement()
                revstr=str(rev)
                rev_comps.append(revstr)
        return (barcodes, rev_comps)

# Create lists for barcode reads and counts
# O(n)
    def barcode_counts(self, bars=None, revs=None):
        bars, revs = self.list_to_barcodes()
        counted_bars = dict(Counter(bars))
        bk = list(counted_bars.keys())
        bv = list(counted_bars.values())
        return (bk, bv)

# Create lists for rev_comp and counts
# O(n log n)
    def rev_counts(self, bars=None, revs=None):
        bars, revs = self.list_to_barcodes()
        counter = []
        for i in revs:
            j = i, bars.count(i)
            counter.append(j)
        counted_revs = dict(counter)
        rk = list(counted_revs.keys())
        rv = list(counted_revs.values())
        return (rk, rv)

# Create dataframe of barcode and rev_comps lists
# O(n)
    def dataframe_main(self, bdk=None, bdv=None, rdk=None, rdv=None):
        bdk, bdv = self.barcode_counts()
        rdk, rdv = self.rev_counts()
        df_main = pd.DataFrame(list(zip(bdk, bdv, rdk, rdv)), 
                  columns=['Barcode 1', 'Number times Barcode 1 found',
                  'Barcode 2 (reverse complement of Barcode 1)',
                  'Number times Barcode 2 found'])
        return (df_main)

# Create specific dataframe for paired barcodes and rev_comp
# O(n)
    def dataframe_pairs(self, datafr=None):
        datafr = self.dataframe_main()
        no_zero_df = datafr[datafr['Number times Barcode 2 found'] != 0]
        min_pair = no_zero_df.sort_values(by=['Number times Barcode 1 found',
                                             'Number times Barcode 2 found'])
        return (min_pair)

# Create specific dataframe for unpaired barcodes
# O(n)
    def dataframe_unpaired(self, dataf=None):
        dataf = self.dataframe_main()
        zero_df = dataf[dataf['Number times Barcode 2 found'] == 0]
        nopair = zero_df[['Barcode 1', 'Number times Barcode 1 found']]
        min_nopair = nopair.sort_values(by=['Number times Barcode 1 found'])
        return (min_nopair)

# Write dataframe to CSV files
# O(n)        
    def df_to_csv(self, pair=None, nopair=None):
        pair = self.dataframe_pairs()
        nopair = self.dataframe_unpaired()
        pair_csv = pair.to_csv('paired_barcodes_revcomps.csv', index=False)
        nopair_csv = nopair.to_csv('unpaired_barcodes.csv', index=False)
        return (pair_csv, nopair_csv)
        
# User input commands        
a=str(input("Enter path to .fastq file: "))
b=int(input("Enter start position of sequence: "))
c=int(input("Enter stop position of sequence: "))

# Initiatin of class and final module
action=code_pairs(a,b,c)
action.df_to_csv()

# Output for timer
print ("Program execution time = ", time.time()-starttime,"s")
print ("Memory usage of program = ", process.memory_info().rss, "MB")





