# CHONK_functions
----------------------

## Training Set: 

### Alignments

You can work with the platinum genomes from the 1000Genomes project. This was the training set I used for my analysis.

They will be in this location

`/projects/ps-gleesonlab5/reference_panels/1kgp/platinumgenomes/`

* GRCh37 (hg19)
* 27 samples with SV calls
* PCR-free 
* 250bp reads
* 500bp insert size
* ~50X Coverage

[Ethnic groups](http://www.internationalgenome.org/category/population/):
```
ACB ASW BEB CDX 
CEU CHB CHS CLM 
ESN FIN GBR GIH 
GWD IBS ITU JPT 
KHV LWK MSL MXL 
PEL PJL PUR STU 
TSI YRI
```

### SV calls

Phase 3 SV calls from 2,504 low-coverage genomes ([Sudmant et al. 2015](https://www.nature.com/articles/nature15394))
```
/projects/ps-gleesonlab5/reference_panels/1kgp/platinumgenomes/platinumgenomes.del.bed  
/projects/ps-gleesonlab5/reference_panels/1kgp/platinumgenomes/platinumgenomes.dup.bed

/projects/ps-gleesonlab5/reference_panels/1kgp/platinumgenomes/platinumgenomes.sv.v8.vcf
```

In the BED files: 
* 0-base
* The genotypes will be your training labels
  * 0|0 = 0
  * 0/1 (1/0) = 1
  * 1/1 = 2
* for now, just focus on chr1-22 (autosomes) 

From the README:
```
Here are some basic statistics about the structural variant call set.

samples: 2504
records: 68818

ALUs: 12748
LINE1s: 3048
SVAs: 835
DELs: 40975
DEL_ALUs: 1238
DEL_LINE1s: 56
DEL_SVAs: 9
DEL_HERVs: 1
DUPs: 6025
CNVs: 2929
INVs: 786
INSs: 168

All male genotypes are unknown on chrX for the union deletions
 because the genotype encoding was not picked up correctly by the phasing
 software. Please refer to the unphased union deletion genotypes for chrX
 where 0/1 encodes copy-number 1 and 1/1 encodes copy-number 0. The unphased
 union deletions are available here:
 ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130723_phase3_wg/union_gls/ALL.wgs.phase3_dels_merged_genome_strip.20130502.dels.low_coverage.genotypes.vcf.gz
```

------------

## Features
[Medvedev et al. 2009](https://www.nature.com/articles/nmeth.1374)
[Alkan et al. 2011](https://www.nature.com/articles/nrg2958)

### Depth of Coverage

Papers:

[CNVnator: Abyzov et al. 2011](https://genome.cshlp.org/content/21/6/974.short)
[Duphold: Pedersen et al. 2018](https://www.biorxiv.org/content/early/2018/11/08/465385.full.pdf)

#### Features

The first three features are taken from the duphold program. If you can install it, you can compare your estimates with dupold. 

* doc_fc
  * fold-change for the variant depth relative to the rest of the chromosome the variant was found on
  * chromosome doc:
    * determine the average read length (you can do this once) = *L*
    * count the reads (readname + mate_1 similar to your split-read code) = *N*
    * count the number of base-pairs you scan = *G*
    * chromosome doc = *L* X *N* / *G* 
  * SV doc:
    * same thing, but now you just do the SV region. You don't have to recalculate the average read length  
    * `for Aln in bam.fetch(reference=chrom, start=SV_start, end=SV_end): read_length.append(Aln.reference_length)`
    * SV doc = *L* X *N* / *SV length*
  * doc_fc = SV doc / chromosome doc
  
* doc_gc
  * fold-change for the variant depth relative to bins in the genome with similar GC-content.
  * bins sizes: 100bp, 500bp, 1kbp
  * in the same regions for the null distribution in `doc_fc`, now determine the doc for 500-1000 bins at different GC context
  * for example, chr1:1-500 , 0.62 GC; doc=0.8. then for all bins at 0.62GC, calculate the median (mean and check if different) 
    * this median/mean DoC at specific GC context will be used, and save it to a file 
    * do this for each chromosome
  * doc_gc = SV_doc  / median doc of the bin with the same GC context (so if the SV has 0.42 GC, you use the 0.42 GC bin) 
  
* doc_flank
  * fold-change for the variant depth relative to flanking regions.

* doc_mad
  * fold-change for the variant median absolute deviation (MAD) of depth relative to the MAD for the rest of the chromosome
* doc_mad_gc
  * fold-change for the variant MAD of depth to bins in the genome with similar GC-content
* doc_mad_flank
  * fold-change for the variant MAD of depth to MAD bins in flanking regions.

#### Defintions

* depth of coverage:
  * *number_of_reads* * *average_read_length* / *base-pairs parsed*
  * for example. chr1:1-300 DEL. 20 reads * 150bp / 300bp = 10X

* (median absolute deviation)[https://en.wikipedia.org/wiki/Median_absolute_deviation]

* fold-change
  * take the depth of coverage of the variant and divide it by the chromosome or bins

* GC-content: the fraction of GC nucleotides for a given region, like 0.65 

* bins
  * experiment with bin sizes. read the methods of CNVnator. Recommend 100bp, 500bp, or 1kbp

* flank
  * experiment with flanking sizes: 500bp or 1kbp

#### Tips

To calculate GC content, use the `bedtools nuc` command or use `samtools faidx` and count it yourself (the latter is likely to be faster)

The command below makes 100bp windows. You only have to do this once and then you can make another file with the chrom:start-end and the GC content 
```
bedtools makewindows -w 100 -b /projects/ps-gleesonlab5/reference_panels/1kgp/platinumgenomes/cn2.regions.grch37.masked.bed >cn2.regions.grch37.masked.windowed100bp.bed

bedtools nuc -fi /home/dantakli/ref/human_g1k_v37.fasta -bed cn2.regions.grch37.masked.windowed100bp.bed | cut -f 1-3,5 >cn2.regions.grch37.masked.windowed100bp.gc.bed
```
Make sure you round to 2 decimal places. 



You only need to calculate the chromosome coverage and binned coverage once. Save the data in a `metadata.txt` file that can be read when
extracting features again. To make this step run faster, build the null distribution of coverage in regions that are intolerant to copy
number variants. `/projects/ps-gleesonlab5/reference_panels/1kgp/platinumgenomes/cn2.regions.grch37.masked.bed`

For the null model:
  * chromosome depth of coverage
  * average read length (pysam Aln.reference_length)
  * chromosome MAD of coverage in bins
  * median (or mean see if there is a significant difference) coverage in bins at different GC content
  * MAD of coverage in bins at different GC content

For the null model you should check the coverage only within regions intolerant to copy number change. 
For you should try to have at least 500 bins at each GC content but no more than 1000 to save computation time.

##### Exploring the Data

You can make plots like histograms and look at the distribution of features with respect to genotype

````
doc    sample   geno
1.2    HG000096  0/0
0.5    HG000096  0/1
0.01   HG000096  1/1
````

```
import seaborn as sns
import pandas as pd
df = pd.read_table("...")
sns.violinplot(y="doc",x="geno",data=df) 
```
https://seaborn.pydata.org/api.html#categorical-plots

Look at the doc_gc distribution with respect to SV length and bin size 


