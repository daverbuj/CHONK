# k-mer genotyping
----------------------

Goal: Make a simple machine learning classifier that can genotype SV using kmers

* [K-mer](http://compbiolwiki.plos.org/wiki/K-mer)
  * unique subsequences of a sequence of length *k*

For example:

![k-mer](http://compbiolwiki.plos.org/w/images/thumb/a/a1/K_mers.svg/406px-K_mers.svg.png)

**The two 3-mers contained within the sequence ATGG**

------------

## Methodology:

1. Make your training set
2. [Pick a classifier]https://scikit-learn.org/stable/tutorial/machine_learning_map/index.html)
3. Train your model and report it's performance
  * Train deletions and duplications separately 

#### Detailed Methodology.

##### 1. Training Set

* SV Callset:
  * Deletions: `/projects/ps-gleesonlab5/reference_panels/1kgp/platinumgenomes/platinumgenomes.del.bed`
  * Duplications: `/projects/ps-gleesonlab5/reference_panels/1kgp/platinumgenomes/platinumgenomes.dup.bed`

* BAM Files:
  * `/projects/ps-gleesonlab5/reference_panels/1kgp/platinumgenomes/*bam`

If a sample does not have a deletion/duplication variant in the bed file, it has a reference (copy number 2) genotype

For example:

```
1       766599  769112  DEL     HG03006 0|1
1       766599  769112  DEL     NA19239 0|1
```

Two samples (HG03006 and NA19239) have a deletion at `1:766599-769112` the remaining 25 samples do not have a deletion

##### Variables for Training Features

Please please experiment with these! Have fun :joy:

* k-mer length
  * use odd numbers
  * try `[5,7,9,11 ... 27,29,31,33]`
  * the bigger the *k* the more k-mers you will have (and thus more memory will be consumed)

* Hamming distance
  * this is a measure of how 

* `SA` tags
  * test reads with/without `SA` tags
  * use `Aln.has_tag('SA')` to find them

* Soft clips
  * test reads with/without soft clips
  * use `if 'S' in Aln.cigarstring` to find them


