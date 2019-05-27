# CHONK Flowchart
-----------------

## Overview

1. Metadata & Obtain Informative Reads
2. Call Breakpoints and Merge them
3. Extract Features
4. Genotype: Germline & Somatic

-----------------

## 1. Metadata & Informative Reads

See `Metadata` under the `Output Contents` section for details.

This step will perform the following steps:

1. Find "inverse" exclude regions
2. Collect metadata on depth of coverage
3. Collect metadata on template lengths
4. Parse out SV informative reads

### Detailed Steps

1. Inverse Excluded Regions

> user will provide an exclude file, if not use the entire lengths of the chromosome found in the BAM header
> using the lengths of the chromosomes, make a temporary BED file of the starts and ends of the chromosomes and subtract the excluded regions from them
> with the new "un-excluded" regions, retain regions >10kb and window them 100bp and 1kbp sliding windows
> calculate the GC content of the bins and save them to a temporary file
> the "un-excluded" file (original and windowed) will be used in all steps but genotyping 

2. Coverage Metadata

> calculate the mean coverage and standard deviation for each chromosome
> calculate the mean coverage and standard deviation for each windowed region (100bp & 1kbp) with respect to GC content

3. Template Length Metadata

> calculate the median template length and median absolute difference for each chromosome

4. SV Informative Reads

> use flags to parse out the following types of reads: split-reads, discordant paired-ends, reads with an unmapped mate, and unmapped reads with a mapped mate.

### Output: metadata JSON file

The metadata JSON output file will contain all the relevant information used for breakpoint calling
and feature extraction and adjustment. It will also contain information on the reference used for alignment,
the sample name and paths to output files. 

The benefits of a JSON output file are the user can provide this file for subsequent CHONK commands and the user
can provide multiple JSON output files to genotype multiple samples.

#### Output Contents:
* File paths
  
  * FASTA File

  * Alignment Files

    * Original BAM file
      * Chromosome Names
      * Chromosome Lengths
      * Sample Name

    * Split-Read / Discordant Paired-end BAM file
    * Anchored / Unmapped Paired-end BAM file
  
  * SV Files
    * CHONK raw breakpoints
    * CHONK merged breakpoints

  * Metadata
    * Depth of Coverage (per chromosome)
      * mean coverage 
      * standard deviation
      * mean coverage of 100bp windows binned by GC content
      * standard deviation of 100bp windows binned by GC content
      * mean coverage of 1kbp windows binned by GC content
      * standard deviation of 1kbp windows binned by GC content

    * Template Length (per chromosome)
      * median template length
      * median absolute deviation 
      * template length thresholds

  * Exclude File
    * inverse exclude file

  * Output
    * path to output
    * path to temporary files
    * path to log files

-----------------

## 2. Breakpoints

This step will perform the following steps:

1. Parse 

-----------------

## 3. Feature Extraction

-----------------

## 4. Genotyping

-----------------
