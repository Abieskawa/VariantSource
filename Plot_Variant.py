import pandas as pd
import os
import argparse as par
import subprocess as sbp
import sys
import dateutil.parser as psr
import gzip
import re
from utils import time_stamp
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import matplotlib.patches as mpatches


class plot_variant(object):
    def __init__(self, args):
        self.arg = args
        self.mode = args.mode
        self.vcf_path = args.vcf_path
        self.parent = args.parent
        self.ind = args.ind
        self.ref_path = args.ref_path
        self.range = args.range
        self.out = args.out
        self.skip = args.skip
        self.gff = args.gff

    def readData(self):
        with open(self.vcf_path, 'r') as file:
            data = file.readlines()  # or any other method to read the file
        return data

    def get_field(self, vcf):
        for line in vcf:
            if re.match(r'[^#]', line): 
                fields = line.split()[8].split(':')
                
                try:
                    GT_pos = fields.index('GT')
                except ValueError:
                    sys.stderr.write(('{} No GT field'
                                      ' in your VCF!!\n').format(time_stamp()))
                    sys.exit(1)
                break
                
            if re.match(r'#CHROM',line):
                cols = line.split('\t')
                parent_pos = {}
                target_pos = {}
                
                for element in cols:
                    if 'parent' in element:
                        parent_pos[element.strip()] = cols.index(element)
                    if '#CHROM' in element:
                        chrom_pos = cols.index(element)
                    if 'POS' in element:
                        variant_pos = cols.index(element)
                    if 'INFO' in element:
                        info_pos = cols.index(element)

                for ind in self.ind:
                    target_pos[ind.strip()] = cols.index(ind)
                
                    
        return GT_pos, parent_pos, target_pos, chrom_pos, variant_pos, info_pos
    
    def process_data(self, vcf, GT_pos, parent_pos, target_ind_pos, chrom_pos, variant_pos, info_pos,
                      chr, pos_range_low = -1, pos_range_up = -1, list_of_skips = []):
        if "ANN" in cols[info_pos]:
            self.ann = 1
            data = {}
        else:
            self.ann = 0
            data = []
        feature_id_list = []
        for line in vcf:
            if './.' not in line:
                cols = line.split('\t')
                if re.match(r'[^#]', line): 
                    if cols[parent_pos[self.parent[0]]].split(':')[GT_pos] \
                        in ['1|1', '1/1']:
                        chrom = cols[chrom_pos]
                        pos = int(cols[variant_pos])
                        if chrom == chr:
                            if pos_range_low != -1:
                                if pos < pos_range_low or pos > pos_range_up :
                                    continue
                            
                            if len(list_of_skips) != 0:
                                t = []
                                for skip in list_of_skips:
                                    if pos == skip[1]:
                                        t.append(-1)
                                if len(t) != 0:
                                    continue 
                                        
                                
                            genotypes = []
                            for ind in self.ind:
                                genotypes.append(cols[target_ind_pos[ind]].split(':')[GT_pos]) 
                            
                            
                            if self.ann == 1:
                                details = cols[info_pos].split(';')[-1].split('|')
                                Annotation_Impact = details[2]
                                Gene_ID = details[4]
                                feature_id = details[6]
                                if feature_id not in feature_id_list:
                                    feature_id_list.append(feature_id)
                                    data[feature_id] = []
                                    data[feature_id].append((chrom, pos, genotypes, Annotation_Impact, Gene_ID))
                                else:
                                    data[feature_id].append((chrom, pos, genotypes, Annotation_Impact, Gene_ID))
                            else:
                                data.append((chrom, pos, genotypes))
                         
        return data

    def create_matrix(self, data):
        length = len(data)
        matrix = pd.DataFrame('white', index=self.ind, columns=range(0, length))
        impact_matrix = pd.DataFrame('white', index=self.ind, columns=range(0, length))
        matrix_list = {}
        rename_list = {}
        
        if self.ann == 1:
            for key, value in data.items():
                for i, record in enumerate(value):
                    chrom, pos, genotypes,Annotation_Impact, Gene_ID, feature_id  = record
                    rename_list[i] = pos
                    for j, Annotation_Impact in enumerate(genotypes):
                        if Annotation_Impact in ['HIGH']:
                            color = 1 
                        elif Annotation_Impact in ['MODERATE']:
                            color = 2
                        elif Annotation_Impact in ['LOW']:
                            color = 3
                        elif Annotation_Impact in ['MODIFIER']:
                            color = 4
                        matrix.at[self.ind[j], i] = color

                    for j, genotype in enumerate(genotypes):
                        if genotype in ['1|1', '1/1']:
                            color = 1 
                        elif genotype in ['0|0', '0/0']:
                            color = 2
                        elif genotype in ['1|0', '1/0', '0/1', '0|1']:
                            color = 3
                    matrix.at[self.ind[j], i] = color
                
                self.gene_id = Gene_ID
                self.feature_id = feature_id

                #rename the column to variant name
                impact_matrix.rename(columns=rename_list, inplace=True)
                matrix.rename(columns=rename_list, inplace=True)

            matrix_list[key] = [matrix,impact_matrix]

            return matrix_list

        else:
            for i, record in enumerate(data):
                chrom, pos, genotypes = record
                rename_list[i] = pos
                for j, genotype in enumerate(genotypes):
                    if genotype in ['1|1', '1/1']:
                        color = 1 
                    elif genotype in ['0|0', '0/0']:
                        color = 2
                    elif genotype in ['1|0', '1/0', '0/1', '0|1']:
                        color = 3
                    matrix.at[self.ind[j], i] = color
            #rename the column to variant name
            matrix.rename(columns=rename_list, inplace=True)
            return matrix


    def load_gff3(self, file_path):
        # Load the GFF3 file into a DataFrame
        columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
        gff3_df = pd.read_csv(file_path, sep='\t', comment='#', names=columns, low_memory=False)
        return gff3_df

    def parse_attributes(self, attributes_str):
        attributes = {}
        for attribute in attributes_str.split(';'):
            if '=' in attribute:
                key, value = attribute.strip().split('=', 1)
                attributes[key] = value
        return attributes
 
    def extract_transcripts(self,gff3_df):
        # Filter for transcript features
        transcript_df = gff3_df[gff3_df['type'].isin(['transcript', 'mRNA'])]

        # Parse attributes column
        transcript_df['attributes_dict'] = transcript_df['attributes'].apply(self.parse_attributes)

        # Extract relevant attributes into separate columns
        for attr in ['ID', 'Parent', 'Name']:
            transcript_df[attr] = transcript_df['attributes_dict'].apply(lambda x: x.get(attr))

        return transcript_df
    
    def modify_column_names(self,columns):
        new_columns = []
        for col in columns:
            col_str = str(col)  # Convert to string for checking
            try:
                if 'K' in col_str:
                    col_int = int(col_str.replace('K', '')) * 1000  # Convert '3K' to 3000
                elif 'M' in col_str:
                    col_int = int(col_str.replace('M', '')) * 1000000  # Convert '3M' to 3000000
                else:
                    col_int = int(col_str)  # Convert to integer for comparison

                if col_int >= 1000000:
                    new_columns.append(f'{col_int // 1000000}M')
                elif col_int >= 1000:
                    new_columns.append(f'{col_int // 1000}K')
                else:
                    new_columns.append(str(col_int))
            except ValueError:
                new_columns.append(col_str) 
        return new_columns
    
    def plot_heatmap(self, matrix, chromosome, target_range, impact_matrix):
        if isinstance(matrix, list):
            matrix = matrix[0]
            impact_matrix = matrix[1]

        # Select specific colors from the PiYG colormap
        cmap = plt.get_cmap('viridis')
        color_list = [cmap(1.2), '#1E90FF', cmap(0.7)]  # Adjust the indices to pick desired colors
        custom_cmap = mcolors.ListedColormap(color_list)
        if self.ann == 1:
            cmap_impact = plt.get_cmap('inferno')
            

        matrix_numeric = matrix.map(lambda x: pd.to_numeric(x, errors='coerce'))
        matrix_numeric.columns = self.modify_column_names(matrix_numeric.columns)
        
        if self.ann == 1:
            impact_matrix_numeric = matrix.map(lambda x: pd.to_numeric(x, errors='coerce'))
            impact_matrix_numeric.columns = self.modify_column_names(matrix_numeric.columns)
        
        fig, ax = plt.subplots(figsize=(20, 5))
        #cax = ax.imshow(matrix_numeric.values, aspect='auto', cmap=custom_cmap, vmin=1, vmax=3)
        ax.pcolormesh(matrix_numeric.values, cmap=custom_cmap, vmin=1, vmax=3, edgecolors='none')
    
        # Color text based on impact
        if self.ann == 1:
            for i,impacts in enumerate(impact_matrix_numeric.values):
                color = 'black'  # default color
                if 1 in impacts:
                    color = "red"  # HIGH impact
                elif 2 in impacts:
                    color = cmap_impact(0.8)  # MODERATE impact
                elif 3 in impacts:
                    color = cmap_impact(0.4) #MODERATE impact
                elif 4 in impacts:
                    color = cmap_impact(0)
                ax.text(-0.5, i + 0.5, ha='right', va='center', color=color, fontsize=8, weight='bold')
        # Show only a subset of the x-axis labels
        x_labels = matrix_numeric.columns
        
        
        ## Add gridlines: x
        seen_m_labels = set()
        if self.mode == "chrom_level":
            visible_labels = []
            visible_ticks = []
            for i, label in enumerate(x_labels):
                if 'M' in label:
                    m_value = int(label.replace('M', ''))
                    if m_value % 5 == 0 and m_value not in seen_m_labels:
                        visible_labels.append(label)
                        visible_ticks.append(i)
                        seen_m_labels.add(m_value)
                    else:
                        visible_labels.append('')
                else:
                    visible_labels.append('')

            # Set the x-ticks to only the positions of the visible labels
            ax.set_xticks(visible_ticks)
            ax.set_xticklabels([label for label in visible_labels if label != ''], rotation=90)
            for x in visible_ticks:
                ax.axvline(x, color='black', linewidth=0.5)
        
        elif self.mode == "gene_level":
            plt.xticks(ticks=np.arange(len(matrix.columns)) + 0.5, labels=matrix.columns, rotation=90)
            for x in np.arange(len(matrix.columns)):
                ax.axvline(x, color='black', linewidth=0.5)
        
        else:
            sys.stderr.write(('{} There is no this kind of mode'
                                      ' type error!!\n').format(time_stamp()))
            sys.exit(1)

        
        ## Add gridlines: y
        for y in np.arange(matrix.shape[0]):
            if y == matrix.shape[0]:
                continue
            ax.axhline(y, color='black', linewidth=0.5)
        ax.set_yticks(np.arange(len(matrix_numeric.index))+0.5)
        ax.set_yticklabels(matrix_numeric.index)
        
        #label and title
        ax.set_xlabel('Position')
        ax.set_ylabel('Individual')
        ax.set_title(f'Chromosome {chromosome} Variants')
        
        # Add legend for custom colors
        colors = {'Parent 1 (V9)': color_list[0], 
                  'Parent 2 (ref)': color_list[1], 
                  'Heterozygous': color_list[2]}
        labels = list(colors.keys())
        handles = [plt.Rectangle((0,0),1,1, color=colors[label]) for label in labels]
        ax.legend(handles, labels, loc='upper center', bbox_to_anchor=(1.15, 1), borderaxespad=0.)

        # Add legend for text colors
        variant_legend = plt.legend(handles, labels, loc='upper right')
        high_patch = mpatches.Patch(color='red', label='HIGH impact')
        moderate_patch = mpatches.Patch(color='orange', label='MODERATE impact')
        plt.legend(handles=[high_patch, moderate_patch], loc='lower right')
        plt.gca().add_artist(variant_legend)


        #Adjust the canvas
        plt.tight_layout(rect=[0.1, 0, 0.9, 1])  

        #save figure to output folder
        plt.savefig(self.out + target_range + '.png', bbox_inches='tight')
        plt.close(fig) 

    def run(self):
        vcf = self.readData()
        GT_pos, parent_pos, target_ind_pos, chrom_pos, variant_pos, info_pos = self.get_field(vcf)
        list_of_skips = []
        if self.skip[0] != '':
            for skip_element in self.skip:
                chr_skip = skip_element.split(':')[0]
                pos_skip = int(skip_element.split(':')[1])
                list_of_skips.append([chr_skip, pos_skip])
    
        for element in self.range:
            if ":" in element:
                #limit its range
                chr = element.split(':')[0]
                pos_range_low = int(element.split(':')[1].split('-')[0])
                pos_range_up = int(element.split(':')[1].split('-')[1])
                list_of_skips_this_time = []
                for skip in list_of_skips:
                    if chr == skip[0] and skip[1] >= pos_range_low and \
                    skip[1] <= pos_range_up:
                        list_of_skips_this_time.append([skip[0],skip[1]])

                data = self.process_data(vcf, GT_pos, parent_pos, target_ind_pos, chrom_pos, variant_pos,info_pos,
                                         chr,pos_range_low,pos_range_up,
                                         list_of_skips = list_of_skips_this_time)
                if self.ann == 1:
                    matrix_list = self.create_matrix(data)
                    for feature_id, matrix_complex in matrix_list.items():
                        print("plotting feature_id", feature_id)
                        self.plot_heatmap(matrix_complex, chr, element)
                else:
                    matrix = self.create_matrix(data)
                    self.plot_heatmap(matrix, chr, element)
            else :
                chr = element
                data = self.process_data(vcf, GT_pos, parent_pos, target_ind_pos, chrom_pos, variant_pos,info_pos,
                                         chr)
                
                matrix = self.create_matrix(data)
                self.plot_heatmap(matrix, chr, element)
                
            
            
            

