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
2. [Pick a classifier](https://scikit-learn.org/stable/tutorial/machine_learning_map/index.html)
3. Train your model and report it's performance
  * **Train deletions and duplications separately**
  * You will make 2 classifiers one for deletions and one for duplications

#### Detailed Methodology.

##### 1. Training Set

Use this code: 

`/home/dantakli/chonk/kmers/src/kmer_training_set.py`

You will need an input file `-i`

* Deletions: `/home/dantakli/chonk/kmers/sv_list_subset_dels.txt`
* Duplications : `/home/dantakli/chonk/kmers/sv_list_dup.txt`

Example: 

```
$ python /home/dantakli/chonk/kmers/src/kmer_training_set.py -i /home/dantakli/chonk/kmers/sv_list_subset_dels.txt -k 9 -d 3
```


Please please experiment with these! Have fun :joy:

* k-mer length `-k INT` 
  * use odd numbers
  * try `[5,7,9,11 ... 27,29,31,33]`
  * the bigger the *k* the more k-mers you will have (and thus more memory will be consumed)

* Hamming distance `-d INT`
  * this is a measure of how two strings are different
  * For example: ATGC vs ATGA have a Hamming distance of 1.
  * We want to exclude k-mers that match very closely to the Reference
  * test different thresholds of Hamming distance `[4,5,6, ... 10,11,12]` 
  * Do not make this value greater than *k*

* `SA` tags `-sa`
  * use the `-sa` option to only use reads with SA tags
  * test reads with/without `SA` tags
  * use `Aln.has_tag('SA')` to find them

* Soft clips `-sc`
  * use the `-sc` option to only use reads with soft-clips
  * test reads with/without soft clips=
  * use `if 'S' in Aln.cigarstring` to find them

* Flanking sequence `-f INT`
  * this is the number of base pairs to search for reads from the start and end positions
  * default is 100bp
  * try `[50, 100, 150, 200, 250, 500, 1000]`


##### Output

Output is in the current working directory you are running you code in with this format:

`kmers_training_set_K{}_HDIST{}_FLANK{}_SA{}_SC{}.txt`
  * K{} = your value of K
  * HDIST = your value of the Hamming distance
  * FLANK{} = your value of the flanking base pairs
  * SA{} = 0 if you did NOT use the `-sa` option, 1 if you did
  * SC{} = 0 if you did NOT use the `-sc` option, 1 if you did

Output Columns:
  * CHROM START END SVTYPE ID
    * these are just informing you about the variant and the ID. Do NOT train on these values
  * GT 
    * this is the genotype. **Use this column as your training labels** 
  * START_RATIO
    * this is the ratio of unique k-mers to total k-mers for the start position
  * END_RATIO
    * this is the ratio of unique k-mers to total k-mers for the end position 

* **Use the START_RATIO and END_RATIO as your training features**

### Tips:

Pick a wide range of K and Hamming distance to find an area that has good performance. Then we will hone in to find the optimal parameter. 

For example: make tables with K= 5,9,13,21,27,33. If we notice that K<13 performs better than K>21, then we should do fine-scale parameter searching, only for K<13

If you are ready to train your model. Follow the examples on scikit learn and use a linear SVM for now. 

-----------------------


----------------------

