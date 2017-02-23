import pysam
import argparse
import csv
import sys


ap = argparse.ArgumentParser()
ap.add_argument('inBam', help='Mapped reads in bam format')
args = ap.parse_args()



samfile = pysam.AlignmentFile(args.inBam, "rb" )
for pileupcolumn in samfile.pileup("3"):
    
    if pileupcolumn.pos in snpSet:
    	if pileupcolumn.n>0:

            allelesSet=set()
            allelesSet.clear()
            
            alleleList=[]
            alleleList[:]=[]
            
            for pileupread in pileupcolumn.pileups:
			
			
            		if not pileupread.is_del and not pileupread.is_refskip:
                	# query position is None if is_del or is_refskip is set.
               			alleleList.append(pileupread.alignment.query_sequence[pileupread.query_position])
        
        
            allelesSet=set(alleleList)
            altCount=alleleList.count(snpDict[pileupcolumn.pos][1])
            
            
            if snpDict[pileupcolumn.pos][1] in allelesSet and snpDict[pileupcolumn.pos][0] not in allelesSet and len(allelesSet)==1:
                if altCount>10:
                            fileWES.write("WES,"+str(pileupcolumn.pos)+","+str(snpDict[pileupcolumn.pos][0])+","+str(snpDict[pileupcolumn.pos][1])+","+str(altCount))
                            fileWES.write("\n")
                            #print allelesSet, snpDict[pileupcolumn.pos], pileupcolumn.pos, altCount
                            homoSet.add(pileupcolumn.pos)
                            dictWES[pileupcolumn.pos]=altCount

#print ('\tbase in read %s = %s' %(pileupread.alignment.query_name,pileupread.alignment.query_sequence[pileupread.query_position]))

samfile.close()



