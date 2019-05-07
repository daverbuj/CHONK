# Somatic Training Set
----------------------

## Objective:

Make a machine learning training set for somatic SVs. 

Overview:

Somatic: these are mutations that are present in some cells in an organism. This is in contrast to germline mutations which are present in all cells.

SV: these mutations change the structure of a chromosome. The basic types are deletions, duplications, insertions, and inversions. SVs are typically >=50bp.

-----------------

## Methodology:

The machine learning classifier will distinguish between true somatic SVs vs. false positive calls. 

SV predictions are extremely noisy and are enriched for false positives. Therefore we will use machine learning to distinguish between signal and noise. 

### Simulation of Somatic SVs

This approach will implement digital spike-ins, which are simply, taking the sequence reads from one individual and mixing them into another individual. 

For example, a somatic deletion in 1/4 of cells can be simply simulated by mixing a heterozygous carrier with a non-carrier.

Let's break this down:

For a non-carrier (no deletion; two normal copies) the coverage of a region is 30X. For a heterozygous carrier (one normal copy) the coverage would be half, so 15X. 

Mixing the two genomes you would get 60X for sequence flanking the deletion, but 45X (30X + 15X) within the deletion. 45x/60x = 3/4 so 1/4 deletion which is somatic. 

### Datasets

* 1000 Genomes Platinum Genomes (27 samples)
  * `/projects/ps-gleesonlab5/reference_panels/1kgp/platinumgenomes/*bam`

* SV gold-standard callset
  * `/home/dantakli/chonk/somatic_training_set/platinumgenomes.del.bed`
  * `/home/dantakli/chonk/somatic_training_set/platinumgenomes.dup.bed`

#### Format of SV callset.

This is a BED file and it is formatted as such,

| CHROM | START (0-base) | END | SVTYPE | SAMPLE | GENOTYPE |
| ----- | -------------- | --- | ------ | ------ | -------- |
| 1     |  11060111      |  11062310 | DEL | NA19238 | 1\|0 |
| 1     |  11230357      |  11232287 | DEL | NA19017 | 0\|1 |
| 1     |  11682872      |  11683189 | DEL | HG00096 | 0\|1 |
| 1     |  11682872      |  11683189 | DEL | HG00268 | 1\|1 |

You are likely familiar with genotypes formatted like this, `0/1`. The `|` indicates the variants are phased, meaning we know which chromosome (maternal or paternal) the variant is on. 

So `1|0` and `0|1` are both heterozygous mutations. While `1|1` is a homozygous mutation.

For the DUPs with `2` in the genotype (such as `0|2` or `2|2`), treat the 2 allele as a 1; so `0|2` is heterozygous and `2|2` is homozygous

--------------------------------

## Pipeline

Before we start mixing the alignments, we should make smaller subsets of alignments that contain the SV in question. This approach will make analysis much faster at the cost of having more files to manage. 

### Steps:

1. Get alignments within SV region (-/+ 10,000bp from start/end positions)
2. Calculate the coverage of the SV region and flanking regions

**more steps will be added but for now this should be enough**

### Detailed Steps

1. Use samtools to retrieve reads within the SV region.

```
$ samtools view -bh NA19238.mapped.ILLUMINA.bwa.YRI.high_coverage_pcr_free.20130924.bam \
	1:11050111-11072310 >NA19238.del.het.1.11060111.11062310.bam 

$ samtools index NA19238.del.het.1.11060111.11062310.bam
```
`-bh` means to write the output in BAM format and the keep the header

The format for the output file must be in this pattern:

`SAMPLE.SVTYPE.GENOTYPE.CHROM.START.END.bam`

  * The `SVTYPE` can be either 'del' or 'dup'.
  * The `GENOTYPE` can be either `het` or `hom` for heterozygous (`1|0` or `0|1`) and homozygous (`1|1`) genotypes
  * `CHROM.START.END` are the SV positions.

2. Use this script I wrote to calculate the coverage

```
$ python /home/dantakli/chonk/somatic_training_set/src/get_coverage.py NA19238.del.het.1.11060111.11062310.bam
```

Here's the output from the script, formatted as, `SV SAMPLE GENOTYPE SV-COVERAGE LEFT-FLANK-COVERAGE RIGHT-FLANK-COVERAGE`
```
1:11060111-11062310	NA19238	het	36.394270122783084	58.6712	58.0546
```

Notice the SV coverage is around half of the flanking regions. This is good! 

Please use this script to make a file that contains the coverage of all the SVs with a `het` or `hom` genotype in the BED files above. 
