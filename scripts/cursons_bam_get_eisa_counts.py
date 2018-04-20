#!/usr/bin/env python
'''
@author: katherine
'''
import HTSeq
from collections import Counter
import configargparse
from os import path
from pyreference import Reference
from pyreference.utils.file_utils import mk_path, mk_path_for_file
from pyreference.utils.genomics_utils import opposite_strand
import sys
import pandas as pd

class ReferenceArgumentParser(configargparse.ArgumentParser):
    def __init__(self, *args, **kwargs):
        super(ReferenceArgumentParser, self).__init__(*args, **kwargs)
        self.add("--build", env_var="BUILD", help='Use [build] section of config file.')
        self.add("--genes_json", help='gtf.json.gz file created by pyreference_gtf_to_json')

    def parse_args(self):
        ''' get args from command line, adding 'reference' field set to PyReference instance '''
        args = super(ReferenceArgumentParser, self).parse_args()
        args.reference = Reference(**args.__dict__)
        return args


##### Utility functions from SACGF libraries #####
def iv_pad(iv, side_padding):
    return HTSeq.GenomicInterval(iv.chrom, iv.start - side_padding, iv.end + side_padding, iv.strand)

def iv_opposite_strand(iv):
    iv = iv.copy()
    iv.strand = opposite_strand(iv.strand)
    return iv

def array_to_file(file_name, array):
    mk_path_for_file(file_name)
    with open(file_name, "w") as f:
        for line in array:
            f.write(line + "\n")

class Reader_format_chrom(object):
    reader_clazz = None
    
    def __init__(self, alignment_filename, want_chr):
        self.reader = self.reader_clazz(alignment_filename)
        self.has_chr = sam_or_bam_has_chrom(self.reader)
        self.want_chr = want_chr
        
        if self.has_chr == self.want_chr: #Use standard reader class if possible
            self.iter_method = self.reader.__iter__
            self.getitem_method = self.reader.__getitem__
        else:
            self.iter_method = self._format_chrom_iter
            self.getitem_method = self._format_chrom_getitem
            
    def _format_chrom_iter(self):
        for aln in self.reader:
            if aln.aligned:
                aln.iv = iv_format_chrom(aln.iv, self.want_chr)
            yield aln
    
    def _format_chrom_getitem(self, iv):
        iv = iv_format_chrom(iv, self.has_chr) # to not get ValueError
        for aln in self.reader[iv]:
            if aln.aligned:
                aln.iv.chrom = format_chrom(aln.iv.chrom, self.want_chr)
            yield aln

    def __iter__(self):
        return self.iter_method()
        
    def __getitem__(self, iv):
        return self.getitem_method(iv)
        
class BAM_Reader_format_chrom(Reader_format_chrom):
    ''' Formats alignments to have/not have 'chr' in chromosome ID
        file name
        Allows random access
        want_chr: Boolean - whether you want "chr" at the beginning of chrom
        return: yields formatted alignment 
    '''
    reader_clazz = HTSeq.BAM_Reader

def sam_or_bam_has_chrom(sam_or_bam):
    ''' return true if sam or bam has chr in header reference names
        warning: this will die if you have BOTH chr and NO chr
    '''
    header_dict = sam_or_bam.get_header_dict()
    
    sq = header_dict['SQ']
    has_chr = 0
    no_chr = 0
    for row in sq:
        chrom = row['SN']
        if chrom.startswith('chr'):
            has_chr += 1
        else:
            no_chr += 1
            
    if has_chr and no_chr:
        raise ValueError("Header has both starting with 'chr' and non-starting with 'chr' reference name!")
    
    return bool(has_chr)

def iv_format_chrom(iv, want_chr):
    iv = iv.copy()
    iv.chrom = format_chrom(iv.chrom, want_chr)
    return iv

def format_chrom(chrom, want_chr):
    ''' Pass in a chromosome (unknown format), return in your format
        @param chrom: chromosome ID (eg 1 or 'chr1')
        @param want_chr: Boolean - whether you want "chr" at the beginning of chrom
        @return: "chr1" or "1" (for want_chr True/False) 
    '''
    if want_chr:   
        if chrom.startswith("chr"):
            return chrom
        else:
            return "chr%s" % chrom
    else:
        return chrom.replace("chr", "")

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
                    #print 'Switching ZEB2-AS1 transcript strand to + (misannotated!)'
                    transcript.iv.strand = '+'
                #else:
                #    print "ZEB2-AS1 is on '+' strand, %r" % transcript
            gene_gas[transcript.iv] += gene_name
    
    overlapping_genes = set()
    for gene_name, gene in reference.genes.iteritems():
        for transcript in gene.transcripts:
            
            for _, genes_in_region in gene_gas[transcript.iv].steps():
                if len(genes_in_region) > 1:
                    overlapping_genes.update(genes_in_region)
    
    non_overlapping_genes = all_gene_names - overlapping_genes
    
    return overlapping_genes, non_overlapping_genes

###### main code #####

def prepare_gas_for_exons_and_introns(is_stranded, reference, overlapping_genes, SIDE_PADDING):
    intronic_gas = HTSeq.GenomicArrayOfSets(chroms='auto', stranded=is_stranded)
    exonic_gas = HTSeq.GenomicArrayOfSets(chroms='auto', stranded=is_stranded)
    
    for gene_name, gene in reference.genes.iteritems():    
        #exclude overlapping genes
        if gene_name in overlapping_genes:
            continue
        
        for transcript in gene.transcripts:
            intron_ivs = transcript.get_intron_ivs()
            
            #Set all introns as intronic
            for intron_iv in intron_ivs:
                intronic_gas[intron_iv] += gene_name
    
    for gene_name, gene in reference.genes.iteritems():
    
        #exclude overlapping genes
        if gene_name in overlapping_genes:
            continue
        
        for transcript in gene.transcripts:
            exon_ivs = [e.iv for e in transcript.exons]
    
            for exon_iv in exon_ivs:
                # For intron counts, clear regions around exon coordinates by 10 bp on both sides
                # to ensure overhanging reads arent counted as intronic.
                #Note: This MUST be done AFTER ALL introns are added to array
                extended_exon_iv = iv_pad(exon_iv, SIDE_PADDING)
                intronic_gas[extended_exon_iv] = set() #Reset
                
                #Set all exons in exonic array
                exonic_gas[exon_iv] += gene_name
    
    return exonic_gas, intronic_gas

def get_exon_and_intron_counts_unstranded(exonic_gas, intronic_gas, bam, want_chr):
    
    bam_reader = BAM_Reader_format_chrom(bam, want_chr=want_chr)
    
    exonic_gene_counts = Counter()
    intronic_gene_counts = Counter()
    
    chromosomes = exonic_gas.chrom_vectors.keys()
    unmatched_chroms = set()
    for i, aln in enumerate(bam_reader):
#         is_exonic = False 
#         is_intronic = False 
        
        if not aln.aligned:
            continue
        
        if not (aln.iv.chrom in chromosomes): #to catch 'chrM'
            unmatched_chroms.add(aln.iv.chrom)
            continue
        
        #If paired end, only use first read in pair, and only uniquely mapping reads
        if (((aln.pe_which == 'not_paired_end') | (aln.pe_which == 'first')) & (aln.aQual >= 50)):
            #Want uniquely mapped reads. These have mapping quality (aln.aQual):
                # 50 for Tophat
                # 255 for STAR
            
            ### Exonic counts ### Count reads that start within any exon of a gene
#             assert len(exonic_gas[aln.iv.start_d_as_pos]) <= 1
            
            for gene_in_exon_region in exonic_gas[aln.iv.start_d_as_pos]:
                exonic_gene_counts[gene_in_exon_region] += 1
#                 is_exonic = True
            
            ### Intronic counts ###
            # Count reads within the gene body that do not overlap any of the annotated exons (intronic).
            overlapping_regions = list(intronic_gas[aln.iv].steps())
            
            if len(overlapping_regions) == 1: #If overlaps and entirely contained within an intron
                genes_in_intron_region = overlapping_regions[0][1]
                
#                 assert len(genes_in_intron_region) <= 1
                
                for gene_in_intron_region in genes_in_intron_region:
                    intronic_gene_counts[gene_in_intron_region] += 1
#                     is_intronic = True
            
#             assert not (is_exonic & is_intronic) #An alignment should not be able to be counted in both exon and intron
        
        if not (i % 1000000):
            print 'Finished %d reads' % (i + 1) 
    print "Finished all. Could not match these chromosome names: %r" % unmatched_chroms
    return exonic_gene_counts, intronic_gene_counts


def get_exon_and_intron_counts_stranded(exonic_gas, intronic_gas, bam, want_chr):
    
    bam_reader = BAM_Reader_format_chrom(bam, want_chr=want_chr)
    
    read1_strand_exonic_gene_counts = Counter()
    read1_strand_intronic_gene_counts = Counter()
    read2_strand_exonic_gene_counts = Counter()
    read2_strand_intronic_gene_counts = Counter()
    
    chromosomes = exonic_gas.chrom_vectors.keys()
    unmatched_chroms = set()
    for i, aln in enumerate(bam_reader):
#         is_exonic = False 
#         is_intronic = False 
        
        if not aln.aligned:
            continue
        
        if not (aln.iv.chrom in chromosomes): #to catch 'chrM'
            unmatched_chroms.add(aln.iv.chrom)
            continue
        
        #If paired end, only use first read in pair, and only uniquely mapping reads
        if (((aln.pe_which == 'not_paired_end') | (aln.pe_which == 'first')) & (aln.aQual >= 50)):
            #Want uniquely mapped reads. These have mapping quality (aln.aQual):
                # 50 for Tophat
                # 255 for STAR
            
            ### Exonic counts ### Count reads that start within any exon of a gene
#             assert len(exonic_gas[aln.iv.start_d_as_pos]) <= 1
            
            aln_opposite_strand_iv = iv_opposite_strand(aln.iv)
            
            #Count as though read1 is on transcribed strand (i.e. allow for different stranded sequencing methods)
            for gene_in_exon_region in exonic_gas[aln.iv.start_d_as_pos]:
                read1_strand_exonic_gene_counts[gene_in_exon_region] += 1
            #Count as though transcribed strand is opposite to read1 strand (= same strand as read2) (i.e. allow for different stranded sequencing methods)
            for gene_in_exon_region in exonic_gas[aln_opposite_strand_iv.start_d_as_pos]:
                read2_strand_exonic_gene_counts[gene_in_exon_region] += 1
            
            ### Intronic counts ###
            # Count reads within the gene body that do not overlap any of the annotated exons (intronic).
            read1_strand_overlapping_regions = list(intronic_gas[aln.iv].steps())
            read2_strand_overlapping_regions = list(intronic_gas[aln_opposite_strand_iv].steps())
            
            if len(read1_strand_overlapping_regions) == 1: #If overlaps and entirely contained within an intron
                read1_strand_genes_in_intron_region = read1_strand_overlapping_regions[0][1]                
#                 assert len(genes_in_intron_region) <= 1
                
                for read1_strand_gene_in_intron_region in read1_strand_genes_in_intron_region:
                    read1_strand_intronic_gene_counts[read1_strand_gene_in_intron_region] += 1
           
            if len(read2_strand_overlapping_regions) == 1: #If overlaps and entirely contained within an intron
                read2_strand_genes_in_intron_region = read2_strand_overlapping_regions[0][1]                
#                 assert len(genes_in_intron_region) <= 1
                
                for read2_strand_gene_in_intron_region in read2_strand_genes_in_intron_region:
                    read2_strand_intronic_gene_counts[read2_strand_gene_in_intron_region] += 1
           
#             assert not (is_exonic & is_intronic) #An alignment should not be able to be counted in both exon and intron
        
        if not (i % 5000000):
            print 'Finished %d reads' % (i + 1) 
    print "Finished all. Could not match these chromosome names: %r" % unmatched_chroms
    return read1_strand_exonic_gene_counts, read2_strand_exonic_gene_counts, read1_strand_intronic_gene_counts, read2_strand_intronic_gene_counts
   
   
#Argument parser
def handle_args():
    parser = ReferenceArgumentParser(description='Produce count file for EISA analysis from bam. Works for Tophat and STAR alignments')
    parser.add_argument("--outdir", required=True, help="Directory to put count file and other stuff into.")
    parser.add_argument("--sample-name", required=True, help="One_word sample name to use for output count file")
    parser.add_argument("--stranded", action='store_true', default=False, help="Should counts be done strand-specifically? Default: unstranded. Set if yes, data is strand-specific.")
    parser.add_argument("bam", help='bam file path')
    return parser.parse_args()

if __name__ == '__main__':
    
    #Parameters
    SIDE_PADDING = 10
    
    args = handle_args()
    outdir = args.outdir 
    bam = args.bam
    sample_name = args.sample_name
    is_stranded = args.stranded
    
    mk_path(outdir)
    print "Outdir: %s" % outdir
    
    boolean_to_stranded_str = {True: 'stranded', False: 'unstranded'}
    stranded_str = boolean_to_stranded_str[is_stranded] 
    print "Working on %s, treating as %s." % (sample_name, stranded_str)
    
    count_file = path.join(outdir, '%s-%s.tsv' % (sample_name, stranded_str))
    print "Count file will be %s" % count_file
    #check if count file exists
    if path.exists(count_file):
        sys.exit("Exiting: %s, count file exists: %s" % (sample_name, count_file))
    
    reference = Reference()
    
    #Identify overlapping genes and prepare genomic array for exonic and intronic regions.
    #To avoid complexities, decided to get the overlapping genes and genomic arrays afresh for every bam (since the time taken to do this is trivial compared to bam processing time) 
    overlapping_genes, non_overlapping_genes = get_overlapping_gene_names(reference, is_stranded)
    print "\nThere are %d overlapping and %d nonoverlapping genes when analysed as %s" % (len(overlapping_genes), len(non_overlapping_genes), stranded_str)
    array_to_file(path.join(outdir, 'overlapping_genes-%s.txt' % stranded_str), overlapping_genes)
    array_to_file(path.join(outdir, 'nonoverlapping_genes-%s.txt' % stranded_str), non_overlapping_genes)
    exonic_gas, intronic_gas = prepare_gas_for_exons_and_introns(is_stranded, reference, overlapping_genes, SIDE_PADDING)
    
    #TODO: Filter gtf - only transcripts that map to a unique position in the genome.
    
    print 'Getting counts for bam: %s' % bam
    if is_stranded:
        read1_strand_exonic_gene_counts, read2_strand_exonic_gene_counts, read1_strand_intronic_gene_counts, read2_strand_intronic_gene_counts = get_exon_and_intron_counts_stranded(exonic_gas, intronic_gas, bam, reference.has_chr)
        count_df = pd.DataFrame({'exon_counts-r1_strand':read1_strand_exonic_gene_counts, 
                                 'intron_counts-r1_strand':read1_strand_intronic_gene_counts,
                                 'exon_counts-r2_strand':read2_strand_exonic_gene_counts, 
                                 'intron_counts-r2_strand':read2_strand_intronic_gene_counts})

    else:        
        exonic_gene_counts, intronic_gene_counts = get_exon_and_intron_counts_unstranded(exonic_gas, intronic_gas, bam, reference.has_chr)
        count_df = pd.DataFrame({'exon_counts':exonic_gene_counts, 'intron_counts':intronic_gene_counts})
    
    count_df = count_df.fillna(0)
    print "Number of reads"
    print count_df.sum(axis=0)
    count_df.to_csv(count_file, sep='\t')
    
    
