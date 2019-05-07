#!/usr/bin/env python
import pysam
import sys
import numpy as np
#####################################

class Metadata():
        def __init__(self,filename):
                meta_list = filename.split('/').pop().split('.')
                self.sample = meta_list[0]
                self.svtype = meta_list[1]
                self.genotype = meta_list[2]
                self.sv = '{}:{}-{}'.format(meta_list[3],meta_list[4],meta_list[5])
                self.chrom = meta_list[3]
                self.left_flank = '{}:{}-{}'.format(self.chrom,
                        int(meta_list[4])-10000, int(meta_list[4]))
                self.right_flank = '{}:{}-{}'.format(self.chrom,
                        int(meta_list[5]), int(meta_list[5])+10000)

def parse_region(region):
        chrom, pos = region.split(':')
        start, end = pos.split('-')
        return chrom,int(start),int(end)

def get_coverage(Bam,region):
        rlen=[]
        for Aln in Bam.fetch(region=region):
                if Aln.is_duplicate or Aln.is_qcfail or Aln.is_unmapped: continue
                rlen.append(Aln.reference_length)

        # coverage = average read length * number of reads / number of base pairs
        chrom,start,end = parse_region(region)
        G = end-start
        N = len(rlen)
        L = np.mean(rlen)
        return (N*L)/float(G)

if __name__ == '__main__':
        filename = sys.argv[1]
        meta = Metadata(filename)
        Bam = pysam.AlignmentFile(filename)
        left_flank_cov = get_coverage(Bam,meta.left_flank)
        right_flank_cov = get_coverage(Bam,meta.right_flank)
        sv_cov = get_coverage(Bam,meta.sv)

        # output format:
        # sv, sample, genotype, sv_coverage, left_flank_cov, right_flank_cov
        out = '\t'.join(map(str,(
                meta.sv,
                meta.sample,
                meta.genotype,
                sv_cov,
                left_flank_cov,
                right_flank_cov)))

        print(out)
