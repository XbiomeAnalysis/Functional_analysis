import argparse
from pprint import pprint
import msea
from msea import SetLibrary
import argparse
import pandas as pd
import sys
import random

## Define arguments
parser = argparse.ArgumentParser(description='Run Microbe-set Enrichment Analysis (MSEA)')

parser.add_argument('--input',
                    type = str,
                    help = 'Input text file containing all DA genus, one genus per line.')

parser.add_argument('--output',
                    type = str,
                    default = './output.csv',
                    help = 'Output csv table of enrichment analysis result. Default ouput.csv in current directory.')

parser.add_argument('--PerturbationTimes',
                    default = 50,
                    type = int,
                    help = 'Output csv table of enrichment analysis result. Default ouput.csv in current directory.')

args = parser.parse_args()

## Define functions
def Read_input(InputGenusFile):       ## Read in genus list in text format, one genus per line
        GenusList = pd.read_table(InputGenusFile, header=None)
        return GenusList        ## Return a python list of input Genus


def Run_MSEA(GenusList, PerturbationTimes):
    gmt_filepath = \
    'https://bitbucket.org/wangz10/msea/raw/aee6dd184e9bde152b4d7c2f3c7245efc1b80d23/msea/data/human_genes_associated_microbes/set_library.gmt'

    print('\nReading database file\n')
    set_library = SetLibrary.load(gmt_filepath)     ## Read human-gene associated microbes library from url.
    random.seed(1)
    set_library.get_empirical_ranks(n=PerturbationTimes)      ## Randomly sampled (n times), without replacement, a universe of microbes under consideration to remove the bias of Fisher's exact test. n is the number of perturbation.

    microbe_set_input = set(GenusList)
    
    print('\nPerforming MSEA with adjustment\n')
    msea_result_adj = set_library.enrich(microbe_set_input, adjust=True, universe=1000)     ## Perform MSEA with adjustment.
    
    print('Showing the head of result table:')
    print(msea_result_adj.head())
    msea_result_adj.to_csv(args.output, sep = '\t')
    
    print('\nFinished !!\n')


## Check if input is 
if (args.input == None):
    parser.print_usage()
    sys.exit(1)
else:
    pass

## Main code
GenusList = Read_input(args.input)
GenusList = GenusList.iloc[:,0].to_list()

## Run MSEA
Run_MSEA(GenusList, args.PerturbationTimes)
