import pysam
import argparse
import csv



ap = argparse.ArgumentParser()
ap.add_argument('inbam', help='Mapped reads in bam format')
ap.add_argument('inVCF', help='Mapped reads in bam format')

args = ap.parse_args()




#------
print "Read homo variants"

homoSet=set()
homoDict={}


#13156,C


with open(args.inVCF) as csvfile:
    readCSV = csv.reader(csvfile, delimiter=',')
    for row in readCSV:
        pos=int(row[0])-1
        allele=row[1]
        homoSet.add(pos)
        homoDict[pos]=allele


samfile = pysam.AlignmentFile(args.inbam, "rb" )
for pileupcolumn in samfile.pileup("1", 1, 500000):
    
    if pileupcolumn.pos in homoSet:
    	if pileupcolumn.n>0:
		list=[]
		for pileupread in pileupcolumn.pileups:
			
			list[:]=[]
            		if not pileupread.is_del and not pileupread.is_refskip:
                	# query position is None if is_del or is_refskip is set.
               			list.append(pileupread.alignment.query_sequence[pileupread.query_position])

		print list
# print ('\tbase in read %s = %s' %
#                      (pileupread.alignment.query_name,
#                       pileupread.alignment.query_sequence[pileupread.query_position]))

samfile.close()


#462236 [T] C

#coverage at base 462235 = 3
#base in read ERR018420.29136049 = T
#base in read ERR018418.10629699 = C
#base in read ERR018418.17833472 = C


#462236