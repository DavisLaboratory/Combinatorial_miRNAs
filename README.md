# Combinatorial_miRNAs in breast cancer EMT
This code repository accompanies the manuscript **Post-Transcriptional Control Of EMT Is Coordinated Through Combinatorial Targeting By Multiple microRNAs** which is currently under review at Cell Systems.

A pre-print article of this work is available from:
* https://doi.org/10.1101/138024

For further information please contact:
* Joe Cursons (cursons.j (at) wehi.edu.au) 
   * General enquiries and computational scripts for reproducing results presented within the corresponding manuscript
* Katherine Pillman (kapillman (at) gmail.com) 
   * Alignment and quantification procedures, including implementation of the exon-intron split analysis
* Cameron Bracken (cameron.bracken (at) sa.gov.au)
   * General and scientific enquiries, in particular relating to experimental protocols
* Melissa Davis (davis.m (at) wehi.edu.au)
   * General and scientific enquiries, in particular relating to the computational analysis 

## Overview
This documentation has been structured as following:
1. RNA-seq data
2. RNA-seq alignment and quantification (general)
3. RNA-seq quantification using the exon-intron split analysis 
4. Other data
5. Computational scripts to reproduce the figures and analysis of the associated manuscript

## 1. RNA-seq data
As described in the associated manuscript, this work involves involves a number of RNA-seq data files that are collected from HMLE cells, mesenchymal HMLE (MesHMLE) cells and mesenchymal HMLE cells transfected with a high concentration of miR-200c-3p (MesHMLE+miR200c). These data are hosted upon the European Nucleotide Archive with the following study accession numbers:
* PRJEB8225
* PRJEB25061
* PRJEB25042
* PRJEB25116
   * small RNA-seq for the quantification of miRNAs
   * **NB**: these data are only available for the HMLE and MesHMLE cells, and small RNA-seq was not performed for the miR-200c transfected MesHMLE cells

   
## 2. RNA-seq alignment and quantification (general)


## 3. RNA-seq quantification using the exon-intron split analysis 
For further details on the exon-intron split analysis (EISA) method, please refer to the associated manuscript:
* Gaidatzis D, Burger L, Florescu M, & Stadler MB. (2015). Analysis of intronic and exonic reads in RNA-seq data characterizes transcriptional and post-transcriptional regulation. *Nature Biotechnology*. **33**(7): pp. 722-9.
* [Link via DOI](http://dx.doi.org/10.1038/nbt.3269): 10.1038/nbt.3269

More detailed instructions and scripts to reproduce this analysis are given in 'EISA_instructions.txt', briefly:
* This work requires the 'pyreference' package which can be installed via pip package manager
~~~
sudo pip install pyreference
~~~

* To improve run times and resource allocation, the reference genome gtf file needs to be converted to the json format 
~~~
pyreference_gtf_to_json.py genes.gtf
~~~

* For the specified data files (note that this analysis uses polyA-depleted RNA-seq data to improve coverage of intronic regions) the numbers of intronic and exonic reads can be extracted using the script provided:
~~~
cursons_bam_get_eisa_counts.py --outdir <outdir> --sample-name <sample_name> --stranded <bam_file>
~~~

* Using sample tables provided (please see the config folder), perform the EISA quantification:
~~~
for i in 1 2 3 4 5 6; do mkdir r${i}; cursons_counts_do_eisa_analysis_on_pair_of_samples.py --build hg19 --outdir <output_directory> --countdir <path_to_directory_with_count_files> --expt-name r${i} --control HMLE --choose-strand r2 --min-count 24 <path_to_repo>/config/HMLE_vs_MesHMLE_sample_table-r${i}.tsv > r${i}/run.log ; done
~~~

* Combine the resulting data and calculate mean values:
~~~
cursons_combine_eisa_analyses_into_spreadsheet.py --outdir <output_directory> --basedir <directory_with_output_from_previous_step> --expt-dir-strs r1,r2,r3,r4,r5,r6 --project-name HMLE_vs_MesHMLE --ctrl-name HMLE --expt-name MesHMLE 
~~~

* **NB**: As demonstrated in 'EISA_instructions.txt' this process is repeated separately for the HMLE-to-MesHMLE and the MesHMLE-to-MesHMLE+miR200c comparisons, as EISA calculates the *change* in intronic and exonic reads between condition pairs.