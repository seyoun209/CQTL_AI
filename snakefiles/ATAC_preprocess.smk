#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import glob
import shutil
import os
from utils.namer import namer
#from snakemake.io import *
#from snakefiles.utils.namer import namer

## Read in samplesheet
samples = pd.read_csv(config["samplesheet"], sep='\t')
#samples = pd.read_csv("samplesheet_NA_removed.txt",sep='\t')
## Convert all columns to strings
samples = samples.astype(str)

## Concatenate the sequencing directory to Read1 and Read2 for full paths
samples['Read1'] = samples[['Sequencing_Directory', 'Read1']].apply(lambda row: os.path.join(*row), axis=1)
samples['Read2'] = samples[['Sequencing_Directory', 'Read2']].apply(lambda row: os.path.join(*row), axis=1)

## Concatenate columns to identify which groups to run (i.e. Seq_Rep will be run together)
samples['mn'] = samples[config['mergeBy']].agg('_'.join, axis=1)
#mergeBy = ['Proj','Donor','Condition','Time']
#samples['mn'] = samples[mergeBy].agg('_'.join, axis=1)

## Group by mn and extract Read1 & Read2
read1 = samples.groupby('mn')['Read1'].apply(list).to_dict()
read2 = samples.groupby('mn')['Read2'].apply(list).to_dict()

## Set run summary name using helper script
runName = namer(samples, config['mergeBy'])
#runName = namer(samples, samples[mergeBy])
