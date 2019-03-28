import sys
import pysam

bam_fh = sys.argv[1]
bam = pysam.AlignmentFile(bam_fh,'r')

for aln in bam.fetch(reference='chr21'):
    prim_chr = aln.reference_name
    mate_chr = aln.next_reference_name
    if aln.mate_is_reverse:
        mate_strand = '-'
    else:
        mate_strand = '+'
    if aln.is_reverse:
        prim_strand = '-'
    else:
        prim_strand = '+'
    if prim_chr == mate_chr:
        if prim_strand != mate_strand:
            temp_len = aln.template_length
            if temp_len > 0:
                print(str(temp_len) + ',')
