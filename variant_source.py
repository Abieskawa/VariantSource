import pandas as pd
import os
import argparse as par
import subprocess as sbp
import sys
import dateutil.parser as psr

from Plot_Variant import plot_variant


#create parser
parser = par.ArgumentParser(
    prog = 'Visual the variant source',
    description="The grogram is to produce digest vcf and categorize the variant\
                with heatmap")
# Add an argument to the parser
parser.add_argument("config_files", nargs='+', help="Paths to the configuration files")
# Parse the arguments
args = parser.parse_args()
print("Path to config file: " + args.config_files[0])

with open(args.config_files[0],mode = 'r')as configure:
    arglines = configure.readlines()

#Data preprocessing
class visualize_variant_source(object):
    def __init__(self, argument):
        self.arg = argument
        self.mode = argument[0].split('\t')[0]
        self.vcf_path = argument[1].split('\t')[0]
        self.parent = argument[2].split('\t')[0].split(',')
        self.ind = argument[3].split('\t')[0].split(',')
        self.ref_path = argument[4].split('\t')[0]
        self.range = argument[5].split('\t')[0].split(',')
        self.out = argument[6].split('\t')[0]
        self.skip = argument[7].split('\t')[0].split(',')
        self.gff = argument[8].split('\t')[0]

    
    def run(self):
        plot_variant(self).run()

visualize_variant_source(arglines).run()
