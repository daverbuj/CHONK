import sys
import pysam

bam_fh = sys.argv[1]
bam = pysam.AlignmentFile(bam_fh,'r')
correct_start_end=[]
wrong_start_end=[]

for aln in bam.fetch(reference='chr21'):
    if aln.has_tag("SA") == True:                       # Alignment is split (SA tag)
        sa_tag = aln.get_tag("SA").split(',')
        aln_two_strand = sa_tag[2]
        aln_two_ref = sa_tag[0]
        aln_two_pos = sa_tag[1]
        if aln.is_reverse:
            aln_one_strand = '-'
        else:
            aln_one_strand ='+'
        if aln_one_strand == aln_two_strand:            # Both alignments are on same strand
            if aln.reference_name == aln_two_ref:       # Both alignments are on the same chr
                cig=aln.cigartuples
                match_length=0
                for val in cig:
                    if val[0] == 0:
                        match_length+=val[1]
                if len(cig) % 2 == 0:                   # If cigarstring has even values
                    if cig[0][0] == 0:
                        aln_side = 'L'
                    else:
                        aln_side = 'X'
                elif len(cig) % 2:                      # If cigarstring has odd values     
                    if cig[0][1] < cig[-1][1]:
                        aln_side = 'L'
                    else:
                        aln_side = 'X'
                if aln_side == 'L':
                    sv_start = aln.pos + match_length
                    sv_end = aln_two_pos
                    if int(sv_start) < int(sv_end):
                        correct_start_end.append(1)
                    else:
                        wrong_start_end.append(1)
                    print(aln.query_name, aln.mapping_quality, sv_start, sv_end)

cor=len(correct_start_end)
wrn=len(wrong_start_end)
print(round(cor/(cor+wrn),3))                            # Prints percent of alignments that have start position befoe end position
