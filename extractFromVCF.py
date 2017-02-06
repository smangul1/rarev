

import vcf
import sys
import argparse

ap = argparse.ArgumentParser()
ap.add_argument('inVCF', help='vcf file of the individual')
ap.add_argument('out', help='list of HOMO SNPs')

args = ap.parse_args()


fOut=open(args.out,"w")

vcf_reader = vcf.Reader(open(args.inVCF, 'r'))
for record in vcf_reader:
	#print record.num_hom_alt, record.ALT,len(record.ALT)
        if len(record.ALT)==1 and len(record.ALT[0])==1:
            fOut.write(str(record.POS)+","+str(record.ALT[0]))
            fOut.write("\n")
	
	
fOut.close()