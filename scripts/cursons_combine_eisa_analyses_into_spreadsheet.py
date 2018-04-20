#!/usr/bin/env python
'''
@author: katherine
'''
from argparse import ArgumentParser
from collections import OrderedDict, Counter
from os import path
import time
import matplotlib.pyplot as pl
import numpy as np
import pandas as pd
import seaborn as sns

def do_eisa_targetscan_gene_ranks_plot(df, ts_score_col, intron_col, exon_col, diff_col):
    '''Produces plot comparing gene rankings from different EISA data types (from EISA paper)'''
    plot_labels = OrderedDict()
    plot_labels[intron_col] = r'$\Delta$intron'
    plot_labels[exon_col] = r'$\Delta$exon'
    plot_labels[diff_col] = r'$\Delta$exon-$\Delta$intron'
    
    eisa_plot_df = pd.DataFrame(data=None, index=df.index)    
    eisa_plot_df['target_in_utr'] = False
    eisa_plot_df.loc[df[ts_score_col] < 0, 'target_in_utr'] = True
    print "Genes with targets in UTR: %r" % Counter(eisa_plot_df['target_in_utr'])
    
    for col in plot_labels.iterkeys():
        sorted_data_for_type = df[col].dropna().sort_values(inplace=False)
        eisa_plot_df['%s_ranks' % col] = sorted_data_for_type.rank(ascending=True) #Most negative have highest ranking
        
    x_range = range(1, 1000 + 1)
    sns.set(font_scale=2.0)
    sns.set_style("white")
    pl.figure(figsize=(6,6))
    for data_type in plot_labels.keys():
               
        data_type_percentages = pd.Series(index=x_range)
        for i in x_range:
            gene_mask = eisa_plot_df["%s_ranks" % data_type] <= i
            data_for_genes = eisa_plot_df.loc[gene_mask, ['target_in_utr']]
            percent_of_genes_with_target = float(sum(data_for_genes.values)) / float(len(data_for_genes)) * 100
                   
            data_type_percentages[i] = percent_of_genes_with_target
               
        pl.plot(x_range, data_type_percentages, label=plot_labels[data_type], linewidth=4)
           
    ax = pl.gca()

    [i.set_linewidth(1.5) for i in ax.spines.itervalues()]
    [i.set_edgecolor('black') for i in ax.spines.itervalues()]
    box = ax.get_position()
    # ax.set_position([box.x0, box.y0 + 0.05 * box.height, box.width * 0.7, box.height])
    xticks = [0, 500, 1000]
    ax.xaxis.set_ticks(xticks)
    ax.set_xticklabels(xticks)#, fontsize=20)
    yticks = [0, 50, 100]
    ax.yaxis.set_ticks(yticks)
    ax.set_yticklabels(yticks)#, fontsize=20)
    pl.xlim(xmin=0) #, xmax=num_genes_on_plot)
    pl.xlabel('Top-ranked genes')
    pl.ylabel('% with seed match')
    ax.legend(loc='upper center', bbox_to_anchor=(0.68, 1.0))
#     pl.text(910, 92, "3'UTR", horizontalalignment='right', verticalalignment='top', fontsize=35) #Title
#     ax.legend(loc='upper center', bbox_to_anchor=(0.67, 0.80))
    ax.set_position([box.x0 + 0.05, box.y0 + 0.04, box.width * 0.98, box.height])
    return ax
    
def do_eisa_hexbin_plot(final_df, intron_col, exon_col, ts_score_col, axis_label):
        '''Produce EISA hexbin of scatterplot
        Input: df with columns: TS scores, mean_intron, mean_exon '''
        sns.set_style("white")
        final_df[ts_score_col] = final_df[ts_score_col].fillna(0)
        pl.figure(figsize=(8,7))
        x = final_df[intron_col]
        y = final_df[exon_col]
        
        pl.hexbin(x, y, C=final_df[ts_score_col], gridsize=24, mincnt=1, cmap='YlOrRd_r', extent=[-6, 6, -6, 6])
        pl.hlines(0, xmin=-6, xmax=6, alpha=0.5)
        pl.vlines(0, ymin=-6, ymax=6, alpha=0.5)
        pl.plot([-5, 5], [-5, 5], '--k', alpha=0.8)
        pl.xlim((-5, 5))
        cb = pl.colorbar()
        cb.set_label('Mean TargetScan score')
        pl.xlabel('%s $\Delta$intron' % axis_label)
        pl.ylabel('%s $\Delta$exon' % axis_label)
        return pl.gca()

def cufflinks_identify_diff_for_2_samples(gene_exp_diff_df, sample_a, sample_b):
    '''Identify the sub-df with expression change data for a pair of samples.
    This is useful because there is no obvious rule regarding which sample will be value_1 or value_2
    '''
    df = gene_exp_diff_df[(gene_exp_diff_df['sample_1'] == sample_a) & (gene_exp_diff_df['sample_2'] == sample_b)]
    
    if len(df) == 0:
        df = gene_exp_diff_df[(gene_exp_diff_df['sample_1'] == sample_b) & (gene_exp_diff_df['sample_2'] == sample_a)]
        if len(df) == 0:
            raise ValueError('sample_1 &/or sample_2 are not in dataframe')
    
    return df.copy()
        
def get_dist_between_confidence_intervals(mean_a, ci_a, mean_b, ci_b):
    '''ci_x is should be stdev or something of the sort, a value to +- to mean to make a ci'''
    if mean_a < mean_b:
        if (mean_a + ci_a) < (mean_b - ci_b):
            return (mean_a + ci_a) - (mean_b - ci_b)
        else:
            return 0
    elif mean_b <= mean_a:
        if (mean_b + ci_b) < (mean_a - ci_a):
            return (mean_a - ci_a) - (mean_b + ci_b) 
        else:
            return 0
    else:
        'Error, whats going on here', mean_a, ci_a, mean_b, ci_b

#Argument parser
def handle_args():
    parser = ArgumentParser(description='Collate EISA analyses into a spreadsheet')
    parser.add_argument("--outdir", required=True, type=str, help="Directory to use for output")
    parser.add_argument("--basedir", required=True, type=str, help="Base directory with the dirs from the pairwise expts in it")
    parser.add_argument("--expt-dir-strs", required=True, type=str, help="Comma-separated list of directory names to use in analysis. e.g. r1,r2,r3")
    parser.add_argument("--project-name", required=True, type=str, help="One_word name for this project.")
    parser.add_argument("--ctrl-name", required=True, type=str, help="Which group was control - for plot label")
    parser.add_argument("--expt-name", required=True, type=str, help="Which group was experimental - for plot label")
    return parser.parse_args()

####Main Script ###
if __name__ == '__main__':
    
    args = handle_args()
    insertion_pos = 6
    outdir = args.outdir
    basedir = args.basedir
    eisa_analyses_dirs = args.expt_dir_strs.split(',')
    print "Analysis dirs: %r" % eisa_analyses_dirs
    project_name = args.project_name
    control_name = args.ctrl_name
    experimental_name = args.expt_name
    
    timestr = time.strftime("%Y%m%d")
    final_df_file = path.join(outdir, '%s_project_EISA_spreadsheet-%s.tsv' % (project_name, timestr))
    final_excel_file = path.join(outdir, '%s_project_EISA_spreadsheet-%s.xlsx' % (project_name, timestr))
    
    #Combine output from multiple analyses.
    final_df = pd.DataFrame()
    eisa_data_files = [path.join(outdir, '%s/%s_EISA_data.tsv' % (d, d)) for d in eisa_analyses_dirs]
    
    for expt_name in eisa_analyses_dirs:
        eisa_data_file = path.join(basedir, expt_name, '%s_EISA_data.tsv' % expt_name)
        one_analysis_df = pd.read_table(eisa_data_file, index_col=0)
        
        delta_intron_data = one_analysis_df['delta-intron']
        delta_intron_data.name = 'delta-intron-%s' % expt_name 
        delta_exon_data = one_analysis_df['delta-exon']
        delta_exon_data.name = 'delta-exon-%s' % expt_name 
        final_df = final_df.join(delta_intron_data, how='outer')
        final_df = final_df.join(delta_exon_data, how='outer')
        
        for col in one_analysis_df.columns:
            if (col.startswith('exon') | col.startswith('intron')):
                if not col.endswith('raw_counts'):
                    if not (col in final_df.columns):                        
                        final_df = final_df.join(one_analysis_df[col].copy(), how='outer')
                        final_df = final_df.join(one_analysis_df["%s_raw_counts" % col].copy(), how='outer') #add column to end
                        
                    else:
                        for gene_name, norm_number in one_analysis_df[col].iteritems():
                            if norm_number >= 0:
                                final_df.set_value(gene_name, col, norm_number) #fill nan holes
                                
                            raw_count_number = one_analysis_df.get_value(gene_name, "%s_raw_counts" % col)
                            if raw_count_number >= 0:
                                final_df.set_value(gene_name, "%s_raw_counts" % col, raw_count_number) #fill nan holes
                                    
    #Messing around to get columns in desired order
    ordered_cols = [] 
    for eisa_analysis_expt in eisa_analyses_dirs:
        ordered_cols.append('delta-exon-%s' % eisa_analysis_expt)
        ordered_cols.append('delta-intron-%s' % eisa_analysis_expt)
    
    exon_cols = [c for c in final_df.columns if (c.startswith('exon') & ~c.endswith('raw_counts'))]
    for exon_col in sorted(exon_cols):
        ordered_cols.append(exon_col)
        
    intron_cols = [c for c in final_df.columns if (c.startswith('intron') & ~c.endswith('raw_counts'))]
    for intron_col in sorted(intron_cols):
        ordered_cols.append(intron_col)
    
    exon_cols = [c for c in final_df.columns if (c.startswith('exon') & ~c.endswith('raw_counts'))]
    for exon_col in sorted(exon_cols):
        ordered_cols.append("%s_raw_counts" % exon_col)
        
    intron_cols = [c for c in final_df.columns if (c.startswith('intron') & ~c.endswith('raw_counts'))]
    for intron_col in sorted(intron_cols):
        ordered_cols.append("%s_raw_counts" % intron_col)
    
    final_df = final_df[ordered_cols]
    
    #Add summary data
    #Calculate means and stdev for exon and intron data
    exon_cols = [c for c in final_df.columns if c.startswith('delta-exon')]
    intron_cols = [c for c in final_df.columns if c.startswith('delta-intron')]
    final_df.insert(0, 'mean_exon', final_df[exon_cols].mean(axis=1))
    final_df.insert(1, 'std_exon', final_df[exon_cols].std(axis=1))
    final_df.insert(2, 'mean_intron', final_df[intron_cols].mean(axis=1))
    final_df.insert(3, 'std_intron', final_df[intron_cols].std(axis=1))
    
    #Calculate minimum distance between samples
    final_df.insert(4, 'ci_gap', np.nan)
    for gene_name in final_df.index:
        gene_data = final_df.loc[gene_name, :]
        distance = get_dist_between_confidence_intervals(gene_data['mean_exon'], gene_data['std_exon'], gene_data['mean_intron'], gene_data['std_intron'])
        final_df.set_value(gene_name, 'ci_gap', distance)
        
    #Calculate diff
    final_df.insert(4, 'diff_e-i', final_df['mean_exon'] - final_df['mean_intron'])
        
    ### Whole experiment EISA Scatterplot ###
    sns.set_style("white")
    pl.figure(figsize=(7,7))
    x = final_df['mean_intron']
    y = final_df['mean_exon']
    sns.jointplot(x=x, y=y, kind='reg')
#     pl.scatter(x, y, alpha=0.07, s=10, c='k', edgecolor='')
    pl.hlines(0, xmin=-7, xmax=7, alpha=0.5)
    pl.vlines(0, ymin=-10, ymax=10, alpha=0.5)
    pl.plot([-5, 5], [-5, 5], '--k', alpha=0.8)
    pl.xlim(xmin=-7, xmax=7)
    pl.ylim(ymin=-10, ymax=10)
    pl.xlabel('%s - %s mean $\Delta$intron' % (experimental_name, control_name))
    pl.ylabel('%s - %s mean $\Delta$exon' % (experimental_name, control_name))
    pl.savefig(path.join(outdir, '%s_project_mean_data_scatterplot.png' % project_name), dpi=150)
    pl.close() 
    
    #Print spreadsheet
    print final_df.head()
    final_df.index.name = 'gene_id' #Set first column header
    final_df.to_csv(final_df_file, sep='\t')
    final_df.to_excel(final_excel_file, engine='xlsxwriter')
    
