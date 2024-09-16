import numpy as np
import os 
import time 
import pandas as pd
import gzip
##import matplotlib.pyplot as plt
from Bio import pairwise2
##from Bio.pairwise2 import format_alignment
import sys

tick=time.time()
index = int(sys.argv[1])-1
library_name = sys.argv[2] # example: 3018__0__HTS_RAM_062

############ UPLOAD REFERENCE TABLE ################

ref_df = pd.read_excel('lRM001_TWIST_library_V10.xlsx')

###### FUNCTIONS ########

match_dic = {}
match_dic[('A', 'A')] = 1
match_dic[('A', 'C')] = -1
match_dic[('A', 'T')] = -1
match_dic[('A', 'G')] = -1
match_dic[('C', 'A')] = -1
match_dic[('C', 'C')] = 1
match_dic[('C', 'T')] = -1
match_dic[('C', 'G')] = -1
match_dic[('G', 'A')] = -1
match_dic[('G', 'C')] = -1
match_dic[('G', 'T')] = -1
match_dic[('G', 'G')] = 1
match_dic[('T', 'A')] = -1
match_dic[('T', 'C')] = -1
match_dic[('T', 'T')] = 1
match_dic[('T', 'G')] = -1
match_dic[('N', 'N')] = 0
match_dic[('N', 'A')] = 0
match_dic[('N', 'T')] = 0
match_dic[('N', 'G')] = 0
match_dic[('N', 'C')] = 0
match_dic[('A', 'N')] = 0
match_dic[('T', 'N')] = 0
match_dic[('G', 'N')] = 0
match_dic[('C', 'N')] = 0

def hamming_distance(seq1,seq2):

    hdist=0
    #print(seq1,seq2)

    for c,char in enumerate(seq1):
        if char.upper()!=seq2[c].upper() and char.upper()!="N" and seq2[c].upper()!="N":
            hdist+=1

    return(hdist)

def reverse_complement(seq):
    rc=[]
    for char in seq.upper():
        if char.upper()=='A':
            rc.append('T')
        elif char.upper()=='T':
            rc.append('A')
        elif char.upper()=='C':
            rc.append('G')
        elif char.upper()=='G':
            rc.append('C')
        elif char.upper()=="N":
            rc.append("N")
        else:
            print('ERROR',char)

    return("".join(rc)[::-1].upper())

def filter(seq,quals):
    
    convert = {'F':37,':':25,',':11,'#':2}
    final_seq=''
    
    for b,base in enumerate(quals):
        if convert[base] > 20:
            final_seq+=seq[b]
        else:
            final_seq+='N'
    
    return(final_seq.upper())

def get_extension(R1_read, match_dict):

    scaffold = 'GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC'
    riboswitch_polyT = 'CGCGGTTCTATCTAGTTACGCGTTAAACCAACTAGAATTTTTTT'
    polyA = 'AAAAAAATT'
    barcode=[]

    try:

        scaffold_alignment = pairwise2.align.localds(scaffold, R1_read[1].upper(), match_dict, -5, -1, one_alignment_only=True)
        scaffold_alignment_percentage = round(scaffold_alignment[0][2]/len(scaffold)*100,0)

        extension_start = scaffold_alignment[0][4]

        riboswitch_alignment = pairwise2.align.localds(riboswitch_polyT, R1_read[1].upper(), match_dict, -5, -1, one_alignment_only=True)
        riboswitch_alignment_percentage = round(riboswitch_alignment[0][2]/len(riboswitch_polyT)*100,0)

        extension_end = riboswitch_alignment[0][3]

        extension_seq = R1_read[1][extension_start:extension_end]
        extension_quals = R1_read[3][extension_start:extension_end]

        return([extension_seq,extension_quals,scaffold_alignment_percentage,riboswitch_alignment_percentage])

    except:

        return(None)


############################## MAIN #########################################

for r,row in enumerate(ref_df.iterrows()):
    if r==index:
   
        exceptions=0
        valid_reads=0
        WT_reads,PE_reads,other_reads,bad_reads,index_tracker=([],[],[],[],[])
        on=time.time()
        ID = row[1]['ID']
        spacer_ref = row[1]['Spacer'].upper()
        barcode_ref = row[1]['Barcode'].upper()
        pegRNA_end_ref = row[1]['pegRNA 3\' end'].upper()
        target_site_ref = row[1]['Target'].upper()
        group = row[1]['Group ID']
        
        R1_file = 'demuxed/'+library_name+'_'+str(index)+'_'+barcode_ref+'_R1.fq.gz'
        R3_file = 'demuxed/'+library_name+'_'+str(index)+'_'+barcode_ref+'_R3.fq.gz'
        
        log_file = 'output2/'+library_name+'_'+str(index)+'_'+barcode_ref+'_results.txt'

        with open(log_file,'w') as w:
            for l in ['initiate v2','\nlibrary',library_name,'index',index]:
                w.write(str(l)+'\n')

        control_alignment = pairwise2.align.localds(pegRNA_end_ref.upper(), reverse_complement(target_site_ref).upper(), match_dic, -3, -1, one_alignment_only=True)   
        #print(control_alignment[0][0])
        #print(control_alignment[0][1])
        prime_edited_sequence = ''
        switch=0
        for c,char in enumerate(control_alignment[0][0]):
            if char=='-' and switch==0:
                prime_edited_sequence+=control_alignment[0][1][c].upper()
            elif char!='-' and switch==0 and len(set(control_alignment[0][0][c+1:]))>1: #other bases ahead in seq
                switch=1
                prime_edited_sequence+=char.upper()
            elif char!='-' and switch==1 and len(set(control_alignment[0][0][c+1:]))==1: #other bases ahead in seq
                switch=0
                prime_edited_sequence+=char.upper()
            elif char!='-' and switch==1 and len(set(control_alignment[0][0][c+1:]))>1: #other bases ahead in seq
                prime_edited_sequence+=char.upper()
            elif char=='-' and switch==1:
                pass

        #print(prime_edited_sequence)

        print(ID,spacer_ref,pegRNA_end_ref,barcode_ref)
        print(target_site_ref,reverse_complement(prime_edited_sequence))

        with open(log_file,'a') as w:
            for l in ['pegRNA ID',ID,'barcode',barcode_ref,'spacer',spacer_ref,'pegRNA 3\' end',pegRNA_end_ref]:
                w.write(str(l)+'\n')
            for l in ['target site reference',target_site_ref,'precise edit reference',reverse_complement(prime_edited_sequence)]:
                w.write(str(l)+'\n')

        R1_read=[]
        R3_read=[]
        read_index=0

        with gzip.open(R1_file,'rt') as f1, gzip.open(R3_file,'rt') as f2:

            for x,y in zip(f1,f2):
                R1_read.append(x.rstrip())
                R3_read.append(y.rstrip())

                if len(R1_read)==4 and len(R3_read)==4 and R1_read[0].split()[0]==R3_read[0].split()[0]:
                   
                    extension = get_extension(R1_read,match_dic)

                    if extension==None:
                        bad_reads.append([R1_read,R3_read])
                    else:
                        #filtered_extension = filter(extension[0],extension[1]) # get_extension() reads [seq,quals]
                    
                        if len(extension[0])!=len(pegRNA_end_ref) or extension[2]<50 or extension[3]<50: #extension[0] is seq, extension[2] is scaffold_alignment_percentage, extension[3] is riboswitch_alignment_percentage
                            bad_reads.append([R1_read,R3_read])
                        
                        else:
                            
                            if hamming_distance(extension[0].upper(),pegRNA_end_ref.upper())==0:

                                valid_reads+=1

                                target_seq = filter(R3_read[1],R3_read[3]) #R3 starts w/ constant GGATCC

                                if hamming_distance(target_seq[6:6+50].upper(),reverse_complement(target_site_ref).upper())==0:
                                    WT_reads.append([R1_read,R3_read])
                                elif hamming_distance(target_seq[6:6+len(prime_edited_sequence)].upper(),prime_edited_sequence.upper())==0:
                                    PE_reads.append([R1_read,R3_read])
                                else:
                                    other_reads.append([R1_read,R3_read])

                            else:
                                bad_reads.append([R1_read,R3_read])

                    R1_read=[]
                    R3_read=[]
                    read_index+=1

                    #if read_index%1000==0 and read_index!=0:
                    #    print(read_index,time.time()-on,'s')
                    #    with open(log_file,'a') as w:
                    #        for l in [read_index,time.time()-on]:
                    #            w.write(str(l)+'\n')
                    #    on=time.time()
                        
        print(read_index)
        print(valid_reads,len(WT_reads),len(PE_reads),len(other_reads),len(bad_reads))

        with open(log_file,'a') as w:
            for l in ['total reads',read_index,'valid reads',valid_reads,'WT reads',len(WT_reads),'PE reads',len(PE_reads),'other reads',len(other_reads),'bad reads',len(bad_reads)]:
                w.write(str(l)+'\n')

        print('writing to file...')

        with gzip.open('output2/'+library_name+'_'+str(index)+'_'+barcode_ref+'_WT_R1.fq.gz','wt') as x1:
            for read in [z[0] for z in WT_reads]:
                for line in read:
                    x1.write(line+'\n')

        with gzip.open('output2/'+library_name+'_'+str(index)+'_'+barcode_ref+'_PE_R1.fq.gz','wt') as x1:
            for read in [z[0] for z in PE_reads]:
                for line in read:
                    x1.write(line+'\n')

        with gzip.open('output2/'+library_name+'_'+str(index)+'_'+barcode_ref+'_other_R1.fq.gz','wt') as x1:
            for read in [z[0] for z in other_reads]:
                for line in read:
                    x1.write(line+'\n')

        with gzip.open('output2/'+library_name+'_'+str(index)+'_'+barcode_ref+'_bad_R1.fq.gz','wt') as x1:
            for read in [z[0] for z in bad_reads]:
                for line in read:
                    x1.write(line+'\n')

        ######################################################################

        with gzip.open('output2/'+library_name+'_'+str(index)+'_'+barcode_ref+'_WT_R3.fq.gz','wt') as y1:
            for read in [z[1] for z in WT_reads]:
                for line in read:
                    y1.write(line+'\n')

        with gzip.open('output2/'+library_name+'_'+str(index)+'_'+barcode_ref+'_PE_R3.fq.gz','wt') as y1:
            for read in [z[1] for z in PE_reads]:
                for line in read:
                    y1.write(line+'\n')

        with gzip.open('output2/'+library_name+'_'+str(index)+'_'+barcode_ref+'_other_R3.fq.gz','wt') as y1:
            for read in [z[1] for z in other_reads]:
                for line in read:
                    y1.write(line+'\n')

        with gzip.open('output2/'+library_name+'_'+str(index)+'_'+barcode_ref+'_bad_R3.fq.gz','wt') as y1:
            for read in [z[1] for z in bad_reads]:
                for line in read:
                    y1.write(line+'\n')

        ######################################################################

        break

with open(log_file,'a') as w:
    for l in ['total time (min)',(time.time()-tick)/60,'\ncomplete']:
        w.write(str(l)+'\n')

