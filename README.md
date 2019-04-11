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



