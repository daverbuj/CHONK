# CHONK

**CHONK** is a machine learning–based tool for detecting and genotyping germline structural variants (GSVs) and calling somatic structural variants (SSVs) from short-read paired-end whole-genome sequencing data. It classifies deletions (DEL), tandem duplications (DUP), and inversions (INV), and is uniquely designed to identify SSVs in bulk (non-tumor) tissue at allele frequencies as low as 1%.

> Masters thesis — San Diego State University, 2020
> Author: Dan Averbuj

---

## Table of Contents

- [Background](#background)
- [Overview](#overview)
- [Pipeline](#pipeline)
- [Installation](#installation)
- [Usage](#usage)
  - [1. metadata](#1-metadata)
  - [2. breakpoints](#2-breakpoints)
  - [3. features](#3-features)
- [Feature Description](#feature-description)
- [Performance](#performance)
- [Training Data](#training-data)
- [Dependencies](#dependencies)
- [Project Structure](#project-structure)
- [References](#references)

---

## Background

Structural variants (SVs) — genomic mutations larger than 50 bp — alter copy number, break or fuse genes, and change gene regulation. They contribute to 0.5–1% of heritable differences between individuals and are implicated in autism, schizophrenia, rare developmental disorders, and cancer.

Somatic SVs (SSVs) arise post-zygotically and are present only in a subset of cells. Because their allele frequency is far below that of germline variants (typically below 25%), they produce a weaker signal and are far harder to detect reliably. Existing callers are limited to tumor cells, microarray resolution, or small variant classes. CHONK is the first tool designed for unbiased SSV detection in bulk WGS using short-read paired-end sequencing, operating at allele frequencies down to 1%.

---

## Overview

CHONK implements a five-stage pipeline:

```
BAM + FASTA
     │
     ▼
┌──────────┐     ┌────────────────┐     ┌──────────────────┐
│ metadata │────▶│  breakpoints   │────▶│    features      │
│          │     │                │     │                  │
│ DOC, GC  │     │ split-read     │     │ coverage, SR/DPE,│
│ insert   │     │ discordant PE  │     │ k-mers, context, │
│ lengths  │     │                │     │ repeat overlap   │
└──────────┘     └────────────────┘     └──────────────────┘
   .json            .bed (SVs)          .DEL.txt / .DUP.txt
                                        .INV.txt
                                             │
                                             ▼
                                    ┌─────────────────┐     ┌──────────────────┐
                                    │    genotype     │────▶│     somatic      │
                                    │                 │     │                  │
                                    │ RF classifier   │     │ RF classifier    │
                                    │ GSV genotyping  │     │ SSV calling      │
                                    │ (6 models)      │     │ (6 models)       │
                                    └─────────────────┘     └──────────────────┘
                                     .vcf (GSVs ≥80%)        .vcf (SSVs ≥50%)
```

Stages 1–3 (this repository) form the preprocessing pipeline: metadata extraction, breakpoint detection, and feature extraction. Stages 4–5 apply trained random forest classifiers to genotype GSVs and call SSVs; each stage uses six independent models stratified by SV type (DEL / DUP / INV) and size (small < 1 kb / large ≥ 1 kb).

The metadata step is expensive but only needs to be run once per sample.

---

## Pipeline

### Stage 1 — metadata

Scans the BAM file through copy-number-neutral regions to compute:

- Per-chromosome **depth of coverage** (DOC), mean read length, and insert-length statistics (mean ± std) via Welford's online algorithm
- A **GC-bin null model**: mean and std of read counts in 25 bp and 1 kb windows, binned by GC content at 4 percentage-point intervals, used to normalise coverage features in the absence of structural variation

Results are saved to a JSON file reused by all subsequent stages.

### Stage 2 — breakpoints

Traverses the BAM file and detects SV breakpoints using two complementary signals:

| Signal | Description |
|--------|-------------|
| **Split reads (SR)** | Reads with an SA tag spanning two positions; strand orientation determines SV type: same strand → DEL or DUP; opposite strands → INV. Provides single-base-pair breakpoint resolution. |
| **Discordant paired ends (DPE)** | Read pairs with insert length > μ + 3.5σ, or with unexpected strand orientation. Opposite-strand → DEL/DUP; same-strand → INV. |
| **Clipped reads (CR)** | Soft-clipped reads lacking a secondary alignment. Used as supporting evidence in feature extraction but not for breakpoint calling. |

Any single supporting alignment is sufficient to record a candidate SV, maximising sensitivity. Candidates are merged by 80% reciprocal overlap. Output is a BED file of candidate loci.

### Stage 3 — features

For each candidate SV, extracts four feature groups (63 features total):

| Group | Features |
|-------|----------|
| **Coverage** (9) | Fold-change DOC for SV body, left flank, right flank; GC-normalised mean/std for each |
| **Supporting fragments** (12) | Ratios of SR, DPE, and clipped fragments; mapping-quality and base-quality statistics for supporting and non-supporting reads |
| **K-mer junction** (10 for DEL/DUP, 12 for INV) | Overlap of read k-mers with reference/alt junction k-mer sets; soft-clip pseudo-alignment ratios |
| **Sequence context** (13) | GC content and sequence complexity for SV body, flanks, and CI regions; log SV length; CI width |
| **Overlap** *(optional, 15)* | Fraction of each region overlapping repeat-masker, segmental-duplication, and STR annotations |

Feature importance varies by SV type: k-mer features are most informative for DELs, coverage features for DUPs, and supporting-alignment features for INVs.

### Stage 4 — genotype *(classifier, not in this repository)*

Applies six random forest models (small/large × DEL/DUP/INV) trained on high-confidence GSV calls to predict genotype (0/0, 0/1, 1/1). SVs with prediction probability ≥ 0.8 are reported as GSVs in VCF format; those with probability < 0.8 are passed to stage 5 as potential SSVs.

### Stage 5 — somatic *(classifier, not in this repository)*

Applies six additional random forest models to classify poorly-genotyped SVs as somatic or not-somatic. Predicted SSVs with probability > 0.5 are reported in a VCF file.

---

## Installation

```bash
git clone https://github.com/daverbuj/chonk.git
cd chonk
pip install .
```

### System requirements

- Python ≥ 3.6
- [`samtools`](http://www.htslib.org/) (accessible in `$PATH`)
- [`bedtools`](https://bedtools.readthedocs.io/) (accessible in `$PATH`)

---

## Usage

### 1. metadata

```
chonk metadata -i BAM -f FASTA [-r CONTIGS] [-x EXCLUDE] [-s SEED] [-o OUTDIR]
```

| Argument | Description |
|----------|-------------|
| `-i` | BAM alignment file (must be indexed) |
| `-f` | Reference FASTA (must be indexed with `samtools faidx`) |
| `-r` | Comma-separated list of contigs to process (default: all) |
| `-x` | BED file of regions to exclude (e.g. repetitive / self-chaining regions) |
| `-s` | Random seed for region sampling (default: 42) |
| `-o` | Output directory (default: same directory as the BAM file) |

**Output:** `<sample>_metadata.chonk.json`

**Example:**
```bash
chonk metadata \
  -i HG00096.bam \
  -f human_g1k_v37.fasta \
  -x cn2.regions.grch37.masked.bed \
  -o metadata/
```

> **Tip:** Providing a copy-number-neutral mask (`-x`) substantially reduces runtime by restricting sampling to stable diploid regions.

---

### 2. breakpoints

```
chonk breakpoints -m METADATA [-r CONTIGS] [-o OUTDIR]
```

| Argument | Description |
|----------|-------------|
| `-m` | Metadata JSON from `chonk metadata` |
| `-r` | Comma-separated list of contigs (default: all in metadata) |
| `-o` | Output directory (default: same directory as metadata file) |

**Output:** `<sample>_breakpoints.chonk.bed`

**Example:**
```bash
chonk breakpoints \
  -m metadata/HG00096_metadata.chonk.json \
  -o breakpoints/
```

---

### 3. features

```
chonk features -bp BED -metadir DIR -fasta FASTA \
               -k K -rlen RLEN -pk PK -o OUT \
               [-rm RM_BED] [-sd SD_BED] [-str STR_BED]
```

| Argument | Description |
|----------|-------------|
| `-bp` | Breakpoint BED file (`chrom start end svtype cipos ciend iid gt`) |
| `-metadir` | Directory of per-sample metadata JSON files |
| `-fasta` | Reference FASTA file |
| `-k` | K-mer size for junction features (e.g. `20`) |
| `-rlen` | Read-length modifier, float (e.g. `1.0`) |
| `-pk` | K-mer size for pseudo-alignment features (e.g. `15`) |
| `-o` | Output file prefix (must end in `.txt`) |
| `-rm` | *(optional)* Repeat-masker merged BED |
| `-sd` | *(optional)* Segmental-duplication merged BED |
| `-str` | *(optional)* Short tandem repeat merged BED |

**Output:** `<prefix>.DEL.txt`, `<prefix>.DUP.txt`, `<prefix>.INV.txt`
(If all three overlap tracks are provided, final output is `<prefix>.DEL.final.txt`, etc.)

**Example:**
```bash
chonk features \
  -bp breakpoints/HG00096_breakpoints.chonk.bed \
  -metadir metadata/ \
  -fasta human_g1k_v37.fasta \
  -k 20 -rlen 1.0 -pk 15 \
  -o features/output.txt \
  -rm hg19_repeatmasker_merged.bed \
  -sd hg19_segdup_merged.bed \
  -str hg19_str_merged.bed
```

---

## Feature Description

### Coverage features (9)

| Feature | Description |
|---------|-------------|
| `sv_doc_fc` | SV body DOC / chromosome DOC |
| `sv_gc_mean` | Mean GC-normalised read-count fold-change across SV windows |
| `sv_gc_std` | Std of above |
| `lf_doc_fc` | Left-flank DOC fold-change (1 kb) |
| `lf_gc_mean` | Left-flank GC-normalised fold-change mean |
| `lf_gc_std` | Left-flank GC-normalised fold-change std |
| `rf_doc_fc` | Right-flank DOC fold-change (1 kb) |
| `rf_gc_mean` | Right-flank GC-normalised fold-change mean |
| `rf_gc_std` | Right-flank GC-normalised fold-change std |

### Supporting fragment features (12)

| Feature | Description |
|---------|-------------|
| `sf_ratio` | Fraction of all fragments that are supporting (SR + DPE + clip) |
| `split_ratio` | Fraction that are split-read supporting |
| `disc_ratio` | Fraction that are discordant PE supporting |
| `clip_ratio` | Fraction that are soft-clipped supporting |
| `sf_mapq_mean/median` | Mapping quality of supporting fragments |
| `nonsf_mapq_mean/median` | Mapping quality of non-supporting fragments |
| `sf_baseq_mean/median` | Base quality of supporting fragments |
| `nonsf_baseq_mean/median` | Base quality of non-supporting fragments |

### K-mer features (10 for DEL/DUP, 12 for INV)

| Feature | Description |
|---------|-------------|
| `ll`, `lr`, `la` | Left-window reads vs. left-ref, right-ref, alt k-mers |
| `rl`, `rr`, `ra` | Right-window reads vs. left-ref, right-ref, alt k-mers |
| `start_ratio_left/right` | Soft-clip k-mers at SV start vs. left/right ref |
| `end_ratio_left/right` | Soft-clip k-mers at SV end vs. left/right ref |

### Sequence context features (13)

| Feature | Description |
|---------|-------------|
| `sv_gc`, `lf_gc`, `rf_gc`, `lo_gc`, `ro_gc` | GC fraction of body, flanks, CI regions |
| `sv_comp`, `lf_comp`, `rf_comp`, `lo_comp`, `ro_comp` | Sequence complexity (zlib-compressed / raw length ratio) |
| `log_sv_len` | log₁₀(SV length) |
| `bp_start_ci`, `bp_end_ci` | Confidence-interval width / read length |

### Overlap features (15, optional)

Fraction of each region (body, left/right flank, left/right CI) overlapping:
`*_rm` (repeat masker), `*_sd` (segmental duplication), `*_str` (short tandem repeats)

---

## Performance

Preliminary results were evaluated on the Polaris cohort using the germline (GSV) and somatic (SSV) classifiers. Note that context-based features were not yet included in these runs, and results are expected to improve with full feature sets and larger training data.

### Germline genotyping (GSV)

Random forest outperformed all other classifiers tested on all three SV types:

| SV type | F1 score |
|---------|----------|
| DEL     | 82%      |
| DUP     | 75%      |
| INV     | 93%      |

### Somatic calling (SSV)

| Metric | Value |
|--------|-------|
| Sensitivity | 70% |
| False discovery rate | 5% |
| SVs correctly identified | 9 / 13 |
| Detectable allele frequency range | 1% – 25% |

The `sf_ratio` (supporting fragment ratio) feature is a strong predictor of allele frequency.

---

## Training Data

### Germline classifier

Trained on three cohorts processed through LUMPY and Manta for SV calling, then filtered by a **consensus + Mendelian expectation** approach:
- Genotypes were determined by consensus of four independent genotypers: DELLY, SVTyper, SV2, and Paragraph.
- Only SVs consistent with Mendelian inheritance within each trio were retained.

| Cohort | Description |
|--------|-------------|
| Polaris | Public WGS trios |
| HGSV | Human Genome Structural Variation Consortium trios (pre-existing callset) |
| PTCD | Pontine Tegmental Cap Dysplasia internal trio dataset |

### Somatic classifier

Trained using **digital and physical spike-in mixtures** of known GSVs to simulate SSVs at controlled allele frequencies. Labels (true positive / false positive) were assigned by checking overlap with the mixture's known variant set and re-genotyping with CHONK.

| Cohort / Mixture | Type | Coverage | AF range |
|-----------------|------|----------|----------|
| Polaris + PTCD (12 children downsampled and mixed) | Digital | 196× | 1–50% |
| GIAB Ashkenazi Jewish trio (two parents mixed) | Digital | 250× | 1.25–50% |
| BSMN physical spike-in (Platinum + HGSV + controls) | Physical | — | ~0.1–2% and 32–68% |

---

## Dependencies

| Package | Use |
|---------|-----|
| `pysam` | BAM file access |
| `numpy` | Numerical aggregation |
| `scikit-learn` | Random forest classifiers (genotype / somatic stages) |
| `samtools` | FASTA sequence extraction |
| `bedtools` | Window generation, GC annotation, BED intersection |

---

## Project Structure

```
chonk/                     # Installable Python package
├── __init__.py
├── cli.py                 # Entry point and argument parsing
├── bam.py                 # BAM file access (BamFile class)
├── alignment.py           # Split-read alignment data structures
├── metadata.py            # Metadata extraction (Metadata class + pipeline)
├── breakpoints.py         # SV breakpoint detection
├── utils.py               # Welford, GC utilities, shared helpers
└── features/
    ├── coverage.py        # Depth-of-coverage features
    ├── fragments.py       # Supporting fragment features (Fragment class)
    ├── kmers.py           # K-mer junction and pseudo-alignment features
    ├── context.py         # Sequence context features
    ├── overlap.py         # Repeat-element overlap features
    └── extract.py         # Feature extraction pipeline orchestrator
setup.py
```

---

## References

- Sudmant et al. (2015). An integrated map of structural variation in 2,504 human genomes. *Nature*. [doi:10.1038/nature15394](https://doi.org/10.1038/nature15394)
- Antaki, Brandler & Sebat (2018). SV2: Accurate structural variation genotyping and de novo mutation detection from whole genomes. *Bioinformatics*. [doi:10.1093/bioinformatics/bty054](https://doi.org/10.1093/bioinformatics/bty054)
- Layer et al. (2014). LUMPY: A probabilistic framework for structural variant discovery. *Genome Biology*. [doi:10.1186/gb-2014-15-6-r84](https://doi.org/10.1186/gb-2014-15-6-r84)
- Chen et al. (2016). Manta: Rapid detection of structural variants and indels for germline and cancer sequencing applications. *Bioinformatics*. [doi:10.1093/bioinformatics/btv710](https://doi.org/10.1093/bioinformatics/btv710)
- Rausch et al. (2012). DELLY: Structural variant discovery by integrated paired-end and split-read analysis. *Bioinformatics*. [doi:10.1093/bioinformatics/bts378](https://doi.org/10.1093/bioinformatics/bts378)
- Chiang et al. (2015). SpeedSeq: Ultra-fast personal genome analysis and interpretation. *Nature Methods*. [doi:10.1038/nmeth.3505](https://doi.org/10.1038/nmeth.3505)
- Pedersen & Quinlan (2019). Duphold: scalable, depth-based annotation and curation of high-confidence structural variant calls. *GigaScience*. [doi:10.1093/gigascience/giz040](https://doi.org/10.1093/gigascience/giz040)
- Chen et al. (2019). Paragraph: A graph-based structural variant genotyper for short-read sequence data. *bioRxiv*. [doi:10.1101/635011](https://doi.org/10.1101/635011)
- Chaisson et al. (2019). Multi-platform discovery of haplotype-resolved structural variation in human genomes. *Nature Communications*. [doi:10.1038/s41467-018-08148-z](https://doi.org/10.1038/s41467-018-08148-z)
- Varoquaux et al. (2015). Scikit-learn. *GetMobile*. [doi:10.1145/2786984.2786995](https://doi.org/10.1145/2786984.2786995)
- D'Gama & Walsh (2018). Somatic mosaicism and neurodevelopmental disease. *Nature Neuroscience*. [doi:10.1038/s41593-018-0257-3](https://doi.org/10.1038/s41593-018-0257-3)
- Medvedev et al. (2009). Computational methods for discovering structural variation. *Nature Methods*. [doi:10.1038/nmeth.1374](https://doi.org/10.1038/nmeth.1374)
