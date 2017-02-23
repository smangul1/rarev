import pysam
import argparse
import csv
import sys


ap = argparse.ArgumentParser()
ap.add_argument('WESbam', help='Mapped reads in bam format')
ap.add_argument('WGSbam', help='Mapped reads in bam format')
ap.add_argument('outWGS', help='Mapped reads in bam format')
ap.add_argument('outWES', help='Mapped reads in bam format')
ap.add_argument('inVCF', help='Mapped reads in bam format')
ap.add_argument('fpReads', help='Mapped reads in bam format')


args = ap.parse_args()




#------
print "Read variants from vcf", args.inVCF

snpSet=set()
snpDict={}


homoSet=set()

dictWES={}


#'3', '7215930', '2', '21280', 'G:21277', 'A:3'

with open(args.inVCF) as csvfile:
    readCSV = csv.reader(csvfile, delimiter='\t')
    for row in readCSV:
        if row!=[]:
            if "CHROM" not in row[0]:
                pos=int(row[1])-1
                allele_ref=row[4].split(":")[0]
                allele_alt=row[5].split(":")[0]
                snpSet.add(pos)
                snpDict[pos]=(allele_ref,allele_alt)





print "Number of variants from vcf", len(snpSet)




fileWES=open(args.outWES,"w")


print "Open",args.WESbam

samfile = pysam.AlignmentFile(args.WESbam, "rb" )
for pileupcolumn in samfile.pileup("3"):
    
    if pileupcolumn.pos in snpSet:
        
    	if pileupcolumn.n>30:
            
            

            allelesSet=set()
            allelesSet.clear()
            
            alleleList=[]
            alleleList[:]=[]
            
            for pileupread in pileupcolumn.pileups:
			
			
            		if not pileupread.is_del and not pileupread.is_refskip:
                	# query position is None if is_del or is_refskip is set.
               			alleleList.append(pileupread.alignment.query_sequence[pileupread.query_position])
        
        
            allelesSet=set(alleleList)
            mainCount=alleleList.count(snpDict[pileupcolumn.pos][0])

            altCount=alleleList.count(snpDict[pileupcolumn.pos][1])
            
            
            
            if snpDict[pileupcolumn.pos][1] in allelesSet and snpDict[pileupcolumn.pos][0] in allelesSet and len(allelesSet)==2:
                
                
                #print alleleList.count(snpDict[pileupcolumn.pos][0]), alleleList.count(snpDict[pileupcolumn.pos][1]), alleleList


                if altCount>=10 and mainCount>10:
                    #print alleleList.count(snpDict[pileupcolumn.pos][0]), alleleList.count(snpDict[pileupcolumn.pos][1]), alleleList

                            fileWES.write("WES,"+str(pileupcolumn.pos)+","+str(snpDict[pileupcolumn.pos][0])+","+str(snpDict[pileupcolumn.pos][1])+","+str(altCount))
                            fileWES.write("\n")
                            #print allelesSet, snpDict[pileupcolumn.pos], pileupcolumn.pos, altCount
                            homoSet.add(pileupcolumn.pos)
                            dictWES[pileupcolumn.pos]=altCount



#print ('\tbase in read %s = %s' %(pileupread.alignment.query_name,pileupread.alignment.query_sequence[pileupread.query_position]))

samfile.close()



fileWES.close()



print "Numbed of homo variantas supported by more than 10 reads according to WES is", len(homoSet)




fileWGS=open(args.outWGS,"w")

readSet=set()
readSetFP=set()
out=open(args.fpReads,"w")


#--------------------------------------
samfile = pysam.AlignmentFile(args.WGSbam, "rb" )
for pileupcolumn in samfile.pileup("3"):
    if pileupcolumn.pos in homoSet:
        if pileupcolumn.n>0:
            alleleList=[]
            alleleList[:]=[]
            
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    alleleList.append(pileupread.alignment.query_sequence[pileupread.query_position])
        
        
                    readSet.add(pileupread.alignment.query_name)
                    
                    
                    if pileupread.alignment.query_sequence[pileupread.query_position]!=snpDict[pileupcolumn.pos][1] and (pileupread.alignment.query_sequence[pileupread.query_position]!=snpDict[pileupcolumn.pos][0]) :
                        #print ('\tbase in read %s = %s' % (pileupread.alignment.query_name, pileupread.alignment.query_sequence[pileupread.query_position]))
                                
                        out.write(str(pileupcolumn.pos)+","+snpDict[pileupcolumn.pos][0]+","+snpDict[pileupcolumn.pos][1] + ","+str(pileupread.alignment)+","+str(pileupread.query_position)+","+str(pileupread.alignment.query_sequence[pileupread.query_position]) +","+ str(snpDict[pileupcolumn.pos][1])+","+ str(pileupread.alignment.query_qualities[pileupread.query_position]))
                        out.write("\n")
                        readSetFP.add(pileupread.alignment.query_name)

                
            allelesSet=set(alleleList)
            
            
            ALT_count=alleleList.count(snpDict[pileupcolumn.pos][1])
            MAIN_count=alleleList.count(snpDict[pileupcolumn.pos][0])


                      
            if len(alleleList)>0:
                
                
                r=ALT_count/float(len(alleleList))

                fileWGS.write(str(pileupcolumn.pos)+","+str(dictWES[pileupcolumn.pos])+","+str(snpDict[pileupcolumn.pos][1])+","+str(len(alleleList))+","+str(ALT_count)+","+str(r))
                fileWGS.write("\n")
            
            #print alleleList.count(snpDict[pileupcolumn.pos][1]),len(alleleList),allelesSet, pileupcolumn.pos, r

fileWGS.close()

print len(readSetFP), len(readSet), float(len(readSetFP))/len(readSet)




