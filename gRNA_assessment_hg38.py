#---------------------------------import---------------------------------------

import urllib2
import math
import time
import random
import socket
import re
import sys, csv, operator
import ssl
import pandas as pd


from threading import Timer

#------------------------------------------------------------------------------

######################################################################################################
### get sequence
def flanking_seq(chr_,snp_):#find the sequence for gRNA design
    try:
        start_= int(snp_)-30
        end_=int(snp_)+30 #gRNA is 20 bp long so get in range
    except:
        print "Wrong pos"
    else:
        url_="http://genome.ucsc.edu/cgi-bin/das/hg38/dna?segment="+chr_+"%3A"+str(start_)+","+str(end_)
        req_ = urllib2.Request(url_)
        resp_ = urllib2.urlopen(req_)
        respHtml_ = resp_.read()
        DNA_seq=respHtml_[respHtml_.find('</DNA>')-63:respHtml_.find('</DNA>')]
        DNA_=DNA_seq.replace('\n','')
    return DNA_
######################################################################################################


######################################################################################################
### find gRNAs
def search_gRNA(id_,seq,strand):#find gRNAs
    pos=[]
    guides=[]
    output=[]
    if strand == '-':
        seq=''.join(["atcg"["tagc".index(n)] for n in seq[::-1]])#reverse complementary  
    pos = [m.start() for m in re.finditer('(?=gg)', seq[20:len(seq)])]
    if len(pos) == 0:
        print "No matches in strand "+strand+" of "+seq
        return []
    else:
        for i in range (0,len(pos)):
#cut site range: xxx(here) xxx N[the position given]GG
            pos[i]+=20
            sub_seq=seq[pos[i]-20:pos[i]-1]
#filter GC content threshold:
            gc_rate=100*(sub_seq.count("g")+ sub_seq.count("c"))/len(sub_seq)
            print gc_rate
            if gc_rate >=20 and gc_rate <=80 and abs(31-(pos[i]-4))<=10: 
#gRNA within resonable range
              if pos[i] >= 35:
                output = [sub_seq,-31+(pos[i]-4),strand]
              else:
                output = [sub_seq,31-(pos[i]-4),strand]
              guides.append(output)

    return guides


#https://www.abmgood.com/marketing/knowledge_base/CRISPR_Cas9_gRNA_Design.php
#referred from http://stackoverflow.com/questions/4664850/find-all-occurrences-of-a-substring-in-python

######################################################################################################

######################################################################################################
#0 get file info
SNP_file="SNP.csv"

#1 Get the list of SNPs
SNP_list=[]
Cant_find_SNP=[]
output_SNP=[]
output_SNP_mm=[]#The ones do not fit the mm requirements or distance requirements (but maybe acceptable)
ID_list=[]

#read in file
SNP_list= pd.read_csv(SNP_file)
print SNP_list

####select columns
for i in range(0,len(SNP_list)):
    print i 
    SNP_list['ID']=('chr'+SNP_list['RSID Chr'].astype(str)+':'+SNP_list['RSID BP'].astype(str))
    SNP_list['chr_']=SNP_list['RSID Chr'].astype(str)
    SNP_list['snp_']=SNP_list['RSID BP'].astype(str)
    break

SNP_list1 = SNP_list[['ID','chr_', 'snp_']]
print SNP_list1

####file that contains no_gRNA SNP
file_object = open('no_gRNA_SNP.txt','w')
file_object.write("no_gRNA_SNP\n")
file_object.close( )
####

######outputfile######Focal output######
outfile = open('SNP_test.txt', 'w');
outfile.write("RSID    gRNA    SNP-gRNA Distance    Strand  On target    1mm    2mm    total_mm"+"\n");
outfile.close();
######outputfile######Focal output######

for i in range(0,len(SNP_list1)):
    print i
    print SNP_list1.iloc[i]
    result=[]
    result_filter=[]
               
    
    #Get sequence based on chr no. and base position 
    target_seq = flanking_seq(SNP_list1['chr_'].iloc[i], SNP_list1['snp_'].iloc[i])
    gRNA_seq_sense = search_gRNA(SNP_list1['ID'].iloc[i], target_seq,'+')
    gRNA_seq_anti = search_gRNA(SNP_list1['ID'].iloc[i], target_seq,'-')
    gRNA_seq=gRNA_seq_sense+gRNA_seq_anti
    print "All the gRNAs for "+SNP_list1['ID'].iloc[i]+":"
    print gRNA_seq
               
    try:
        len(gRNA_seq[0])
    except:
        Cant_find_SNP.append(SNP_list1['ID'].iloc[i])
        print "No matched gRNA found for " + SNP_list1['ID'].iloc[i]
        file_object = open('no_gRNA_SNP.txt','a')
        file_object.write(SNP_list1['ID'].iloc[i]+"\n")
        file_object.close( )
        pass
    else:
        for j in range(0,len(gRNA_seq)):
            print gRNA_seq[j][0]
            result.append([SNP_list1['ID'].iloc[i]]+gRNA_seq[j])
	    print result
        
        for j in range (0, len(result)):
            if result[j][2]<=10: #with in 10 basepairs
               result_filter.append(result[j])

        if len(result_filter)==0:#no suitable gRNA for this SNP.
            Cant_find_SNP.append(SNP_list1['ID'].iloc[i])
            print "No matched gRNA found for " + SNP_list1['ID'].iloc[i]
            file_object = open('no_gRNA_SNP.txt','a')
            file_object.write(SNP_list1['ID'].iloc[i]+"\n")
            file_object.close( )
        
        else:
            ID_list.append(SNP_list1['ID'].iloc[i])
            output_SNP=output_SNP+result_filter

print "SNP_list"
print SNP_list1

print "Cant_find_SNP"
print Cant_find_SNP

print "output_SNP is "
print output_SNP

#####################################################################################################
######outputfile
outfile = open("hg38_gRNA_to_process.txt", "w");
outfile.write("Focal_RSID	gRNA	SNP-gRNA Distance	Strand"+"\n");
for i in range(0,len(output_SNP)):
    for j in range(0,len(output_SNP[i])):
        outfile.write(str(output_SNP[i][j])+"	")
    outfile.write("\n")
outfile.close();
#####################################################################################################





