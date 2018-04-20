#!/usr/bin/env python
'''
@author: katherine
'''
import HTSeq
from collections import OrderedDict
import configargparse
from os import path
from pyreference.reference import Reference
from pyreference.utils.file_utils import mk_path
from scipy.stats import linregress

import matplotlib.pyplot as pl
import numpy as np
import pandas as pd


class ReferenceArgumentParser(configargparse.ArgumentParser):
    def __init__(self, *args, **kwargs):
        super(ReferenceArgumentParser, self).__init__(*args, **kwargs)
        self.add("--build", env_var="BUILD", help='Use [build] section of config file.')

    def parse_args(self):
        ''' get args from command line, adding 'reference' field set to PyReference instance '''
        args = super(ReferenceArgumentParser, self).parse_args()
        args.reference = Reference(**args.__dict__)
        return args


def get_overlapping_gene_names(reference, is_stranded):
    '''is_stranded: Whether a transcript on the opposite strand should be considered overlapping. If is_stranded=True, only consider genes on same strand as overlapping. 
    Returns: lists of names of genes that overlap another gene or do not overlap another gene (strand-specific or not as specified in input)'''
    all_gene_names = set()
    gene_gas = HTSeq.GenomicArrayOfSets(chroms='auto', stranded=is_stranded)
    for gene_name, gene in reference.genes.iteritems():
        all_gene_names.add(gene_name)
        
        for transcript in gene.transcripts:
            if gene_name == 'ZEB2-AS1':
                if transcript.iv.strand == '-': #ZEB2 AS transcript strand was misannotated in iGenomes UCSC hg19.
                    print 'Switching ZEB2-AS1 transcript strand to + (misannotated!)'
                    transcript.iv.strand = '+'
                else:
                    print "ZEB2-AS1 is on '+' strand, %r" % transcript
            gene_gas[transcript.iv] += gene_name
    
    overlapping_genes = set()
    for gene_name, gene in reference.genes.iteritems():
        for transcript in gene.transcripts:
            
            for _, genes_in_region in gene_gas[transcript.iv].steps():
                if len(genes_in_region) > 1:
                    overlapping_genes.update(genes_in_region)
    
    non_overlapping_genes = all_gene_names - overlapping_genes
    
    return overlapping_genes, non_overlapping_genes


def normalise_counts_for_expt(ordered_sample_grouping, count_file_regex, norm_factors_file, sample_table, PSEUDOCOUNT, genes_in_all_analyses, MINIMUM_COUNT, outdir, expt_name, four_col_input, args):
    #Normalisation is performed using all samples in the analysis.
    #This means it must be redone each time a new sample is included/excluded.
     
    #Calculate normalisation factors for all samples, for intron and exonic counts separately
     
    #Ensure all samples have same number of genes, otherwise averages will be out.
    #Exclude genes where counts could not be obtained for all samples.
    
    exon_df = pd.DataFrame()
    intron_df = pd.DataFrame()
    sample_names = []
    for _, group_sample_names in ordered_sample_grouping.iteritems():
        for sample_name in group_sample_names:
            sample_names.append(sample_name)
            stranded_str = sample_table.get_value(sample_name, 'stranded')
            
            count_df = pd.read_table(count_file_regex % (sample_name, stranded_str), index_col=0)
            
            #If strand-specific, counts have 4 columns so discard the columns relating to the wrong strand and rename
            if four_col_input:
                strand_cols = [c for c in count_df.columns if c.endswith('%s_strand' % args.choose_strand)]
                count_df = count_df.loc[:, strand_cols].copy()
                count_df.columns = [c.split('-')[0] for c in count_df.columns]
            else:
                print "%s input only has 2 columns" % sample_name
            #drop genes which were not included in counting in all samples. i.e. if one sample was stranded and another wasn't 
            count_df = count_df.loc[genes_in_all_analyses, :]
            
            exon_df = exon_df.join(pd.Series(count_df['exon_counts'], name=sample_name), how='outer')
            intron_df = intron_df.join(pd.Series(count_df['intron_counts'], name=sample_name), how='outer')
    
    #Fill missing data with zeros - because we have already removed genes that weren't included in counting, it's safe to conclude absence means 0 reads were present at that locus.
    exon_df = exon_df.fillna(0)
    intron_df = intron_df.fillna(0)
    
    assert len(sample_names) == len(set(sample_names)) #sample names must be unique
    
    #Analyse distributions of the data
    pl.figure()
    for col in intron_df.columns:
        hist_data = intron_df[intron_df[col] > 0][col]
        pl.hist(np.log2(hist_data + 1), bins=50, range=(0, 15), alpha=0.2, label=col)
    pl.legend()
    pl.savefig(path.join(outdir, '%s_intron_raw_counts_histogram.png' % expt_name))
    pl.close()
    
    #This is the total number of reads mapped to each region
    norm_factors_df = pd.DataFrame(index=sample_names, columns=['exon', 'intron'])
    
    #Calculate total number of reads mapped to genes.
    norm_factors_df.loc[sample_names, 'exon'] = exon_df.sum(axis=0)
    norm_factors_df.loc[sample_names, 'intron'] = intron_df.sum(axis=0)
    #Convert to per million reads mapped
    norm_factors_df = norm_factors_df / 1000000.0
    print 'Norm factors: Number of reads mapped (x10^6)'
    print norm_factors_df
    
    norm_factors_df.to_csv(norm_factors_file, sep='\t')
    
    ### Apply expression cutoffs to - exclude genes which don't have good expression in at least 1 sample.
    '''Genes should be expressed above a threshold in at least 1 group.'''
    # Select genes with sufficient counts - >= 24 in at least 1 group for both exonic and intronic counts.
    print "Minimum count required for both exonic and intronic in at least 1 group: %d" % MINIMUM_COUNT
    #Genes where exon passes threshold and intron passes threshold for group 1
    sample_1_passes_threshold_mask = (exon_df[sample_names[0]] >= MINIMUM_COUNT) & (intron_df[sample_names[0]] >= MINIMUM_COUNT)
    sample_2_passes_threshold_mask = (exon_df[sample_names[1]] >= MINIMUM_COUNT) & (intron_df[sample_names[1]] >= MINIMUM_COUNT)
    
    at_least_1_sample_passes_threshold_mask = sample_1_passes_threshold_mask | sample_2_passes_threshold_mask
    print "%s: %d genes pass the threshold" % (sample_names[0], sum(sample_1_passes_threshold_mask))
    print "%s: %d genes pass the threshold" % (sample_names[1], sum(sample_2_passes_threshold_mask))
    
    good_coverage_exon_df = exon_df.loc[at_least_1_sample_passes_threshold_mask, :]
    good_coverage_intron_df = intron_df.loc[at_least_1_sample_passes_threshold_mask, :]
    print "There are %d genes which pass the minimum cutoffs for both exon and intron in at least 1 sample" % sum(at_least_1_sample_passes_threshold_mask)
    
    # Normalise counts for library size - divide each sample by total number of reads (x10^6).
    norm_factors_df = pd.read_table(norm_factors_file, index_col=0)
        
    #Use standardised ratio so that we can compare between experiments
    average_norm_factor = {'exon': 40.0, 'intron': 2.0} # Approximately 20-fold diff between intron and exon counts observed for polyA data.
    
    normalised_data_df = pd.DataFrame()
    norm_pos = 0
    for sample_type, type_df in zip(['exon', 'intron'], [good_coverage_exon_df, good_coverage_intron_df]):
        pl.figure()
        for sample_name in sample_names:
            norm_factor_for_sample_and_type = norm_factors_df.loc[sample_name, sample_type]
            ave_norm_factor_for_type = average_norm_factor[sample_type]
#             ave_norm_factor_for_type = 50.0 #ie express as values per 50 million reads
            normalised_counts_for_type = (type_df[sample_name] / norm_factor_for_sample_and_type) * ave_norm_factor_for_type
            
            #plot normalised data
            pl.hist(np.log2(normalised_counts_for_type + 1), bins=50, range=(0, 15), alpha=0.2, label=sample_name)
            
            # Add pseudocount and calculate log2 expression levels.
            normalised_counts_for_type_plus_pseudocount = normalised_counts_for_type + PSEUDOCOUNT
            normalised_data_df.insert(norm_pos, '%s-%s' % (sample_type, sample_name), np.log2(normalised_counts_for_type_plus_pseudocount))
            count_data_for_type = type_df.loc[normalised_counts_for_type.index, sample_name]
            normalised_data_df.insert(len(normalised_data_df.columns), '%s-%s_raw_counts' % (sample_type, sample_name), count_data_for_type)
            norm_pos += 1
        
        pl.legend()
        pl.savefig(path.join(outdir, '%s_%s_normalised_counts_histogram.png' % (expt_name, sample_type)))
        pl.close()
    
    pl.figure()
    for sample_name in sample_names:
        pl.hist(normalised_data_df['intron-%s' % sample_name], bins=50, range=(0, 15), alpha=0.2, label=sample_name)
    pl.legend()
    pl.savefig(path.join(outdir, '%s_intron_normalised_and_pseudocount_histogram.png' % expt_name))
    pl.close()
    
    return normalised_data_df

def get_sample_groups_from_sample_table(sample_table):
    #Get order of samples from table
    ordered_sample_grouping = OrderedDict()
    for group, group_sample_table in sample_table.groupby('group'):
        ordered_sample_grouping[group] = group_sample_table.index.tolist()
    
    sample_groups = ordered_sample_grouping.keys()
    if len(sample_groups) != 2:
        assert ValueError('There cannot be more than 2 groups!, not %r' % sample_groups)
    
    return ordered_sample_grouping
    
#Argument parser
def handle_args():
    parser = ReferenceArgumentParser(description='Perform EISA analysis on bam files')
    parser.add_argument("--outdir", required=True, type=str, help="Base directory to put output dir in")
    parser.add_argument("--countdir", required=True, type=str, help="Directory with counts that match sample table")
    parser.add_argument("--expt-name", required=True, type=str, help="One_word name for this combination of samples, e.g. r1.")
    parser.add_argument("--control", required=True, type=str, help="Which group in the table should be treated as the control?")
    parser.add_argument("--choose-strand", required=False, type=str, help="If data is stranded, which strand to choose? Options: r1 (if read1 is on transcribed strand) or r2")
    parser.add_argument("--min-count", required=False, type=int, default=24, help="Minimum number of reads for a gene (in introns/exons)")
    parser.add_argument("table", help='config table')
    return parser.parse_args()

####Main Script ###

if __name__ == '__main__':
    
    args = handle_args()
    
    sample_table = pd.read_table(args.table, index_col=None)
    sample_table = sample_table.set_index('name', verify_integrity=True) #Names must be unique
    
    base_dir = args.outdir
    expt_name = args.expt_name
    count_dir = args.countdir
    control_group = args.control
    MINIMUM_COUNT = args.min_count #24 #Paper: 24 - This is (2^5) - 8(pseudocount)
            
    outdir = path.join(base_dir, expt_name) 
    mk_path(outdir)
    print "Outdir: %s" % outdir
    print "Count files are in %s" % count_dir
    
    #Strand specific data has 4 columns
    if args.choose_strand is not None:
        four_col_input = True
        assert ((args.choose_strand == "r1") | (args.choose_strand == "r2"))
        print "Expt is stranded, using strand %s" % args.choose_strand
    else:
        four_col_input = False
        print "Expt is unstranded"
        
    ### Parameters ###
    SIDE_PADDING = 10 #Paper: 10
    PSEUDOCOUNT = 8 #Paper: 8
    overlapping_genes_file = path.join(outdir, 'overlapping_genes-%s.txt')
    non_overlapping_genes_file = path.join(outdir, 'nonoverlapping_genes-%s.txt')
    count_file_regex = path.join(count_dir, '%s-%s.tsv') # format: sample_name, stranded_str
    exonic_gas_pkl_file = path.join(outdir, "exonic_gas-%s.pkl")
    intronic_gas_pkl_file = path.join(outdir, "intronic_gas-%s.pkl")
    norm_factors_file = path.join(outdir, 'normalisation_factors.tsv')
    normalised_count_file = path.join(outdir, 'normalised_counts.tsv') 
    mean_df_file = path.join(outdir, 'mean_data_for_groups.tsv')
    final_eisa_df_file = path.join(outdir, '%s_EISA_data.tsv' % expt_name)
    stranded_to_boolean = {'stranded': True, 'unstranded': False}

    ### Checking setup is OK ###
    print 'Sample table'                        
    print sample_table
    #Check that there are only 2 samples 
    assert len(sample_table) == 2

    ordered_sample_grouping = get_sample_groups_from_sample_table(sample_table) #Keep order in table - not sure whether to do this or change to put control first
    
    sample_groups = ordered_sample_grouping.keys()
    if len(sample_groups) != 2:
        assert ValueError('There cannot be more than 2 groups!, not %r' % sample_groups)
        
    if not control_group in set(sample_groups):
        assert ValueError('The control label given ("%s") is not in the sample table groups: %r' % (control_group, sample_groups))
    
    #Identify control and experimental groups and sample names
    for sample_name, sample_data in sample_table.iterrows():
        if sample_data['group'] == control_group:
            control_sample = sample_name
        else:
            experimental_sample = sample_name
            experimental_group = sample_data['group']
    print 'Control sample: %s, control group: %s' % (control_sample, control_group)
    print 'Experimental sample: %s, experimental group: %s' % (experimental_sample, experimental_group)
    
    #Check count files exist
    assert path.isfile(count_file_regex % (control_sample, sample_table.get_value(control_sample, 'stranded'))), count_file_regex % (control_sample, sample_table.get_value(control_sample, 'stranded'))
    assert path.isfile(count_file_regex % (experimental_sample, sample_table.get_value(experimental_sample, 'stranded')))
        
    ### Main Analysis ###
    pd.set_option("display.large_repr", "info") #Pretty table printing
    
    print "Printing a copy of the sample table used to the output dir as a record\n"
    sample_table.to_csv(path.join(outdir, '%s_sample_table.tsv' % expt_name), sep='\t')
    
    reference = Reference()
    
    genes_in_all_analyses = set(reference.genes.keys())
    #If some samples are counted strand-specifically and others unstranded, there will be some genes that are not counted in one or the other.
    #Only want to work with genes which are counted in ALL analyses.
    #Identify genes counted in all analyses:
    for stranded_str in set(sample_table['stranded']): 
        is_stranded = stranded_to_boolean[stranded_str]
        
        overlapping_genes, non_overlapping_genes = get_overlapping_gene_names(reference, is_stranded)
        print "\nThere are %d overlapping and %d nonoverlapping genes when analysed as %s" % (len(overlapping_genes), len(non_overlapping_genes), stranded_str)
        
        genes_in_all_analyses = genes_in_all_analyses.intersection(non_overlapping_genes)
    
    print "There are %d genes in all analyses" % len(genes_in_all_analyses)
    
    print "\nNormalising count data"
    # Normalise counts for this experiment    
    normalised_data_df = normalise_counts_for_expt(ordered_sample_grouping, count_file_regex, norm_factors_file, sample_table, PSEUDOCOUNT, genes_in_all_analyses, MINIMUM_COUNT, outdir, expt_name, four_col_input, args)
    
    # Calculate change in exon (or intron) = the difference between log2 exonic counts between experimental conditions.    
    print "\nComparing samples: %s - %s\n" % (experimental_sample, control_sample)
    change_in_exon = normalised_data_df['exon-%s' % experimental_sample] - normalised_data_df['exon-%s' % control_sample] 
    change_in_intron = normalised_data_df['intron-%s' % experimental_sample] - normalised_data_df['intron-%s' % control_sample]
    posttranscriptional_reg = change_in_exon - change_in_intron
    
    change_between_groups_df = pd.DataFrame({'delta-intron': change_in_intron,
                  'delta-exon': change_in_exon,
                  'diff_bet_delt-exon_and_delt-intron': posttranscriptional_reg})
    
    final_df = change_between_groups_df.join(normalised_data_df)
    print "Final df:"
    print final_df
    final_df.to_csv(final_eisa_df_file, sep='\t')
    print "Finished printing final df to file: %s" % final_eisa_df_file
    
    ### EISA scatterplot ###
    pl.figure(figsize=(7,7))
    x = final_df['delta-intron']
    y = final_df['delta-exon']
    slope, intercept, r, _, _ = linregress(x,y)
    print "slope: %.2f, intercept: %.2f, r: %.2f" % (slope, intercept, r)
    pl.scatter(x, y, label='R=%.2f' % r, alpha=0.15, s=4, c='k')
    pl.hlines(0, xmin=-5, xmax=5, alpha=0.5)
    pl.vlines(0, ymin=-5, ymax=5, alpha=0.5)
    pl.plot([-5, 5], [-5, 5], '--k', alpha=0.8)
    pl.xlabel('%s - %s ($\Delta$intron)' % (experimental_group, control_group))
    pl.ylabel('%s - %s ($\Delta$exon)' % (experimental_group, control_group))
    pl.legend(loc='upper left', scatterpoints=1)
    pl.title('%s, n=%d' % (expt_name, len(final_df)))
    pl.ylim(ymin=-11, ymax=11)
    pl.xlim(xmin=-8, xmax=8)
    pl.savefig(path.join(outdir, '%s_EISA_scatterplot.png' % expt_name))
    pl.close()
    

    
