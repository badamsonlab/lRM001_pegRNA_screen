import numpy as np
import os 
import time 
import pandas as pd
import gzip
import sys

R1_file = sys.argv[1]
R3_file = sys.argv[2]

log_file = 'demuxed/'+R1_file.split('-')[0]+'_demux_log.txt'

with open(log_file,'w') as w:
    for l in ['initiate']:
        w.write(l+'\n')

############ UPLOAD REFERENCE TABLE ################

ref_df = pd.read_excel('PEmbryo_TWIST_library_V10.xlsx')

###### FUNCTIONS ########

def hamming_distance(seq1,seq2):

    hdist=0
    #print(seq1,seq2)

    for c,char in enumerate(seq1):
        if char.upper()!=seq2[c].upper() and char.upper()!="N" and seq2[c].upper()!="N":
            hdist+=1

    return(hdist)

def reverse_complement(seq):
    rc=[]
    for char in seq.lower():
        if char.lower()=='a':
            rc.append('t')
        elif char.lower()=='t':
            rc.append('a')
        elif char.lower()=='g':
            rc.append('c')
        elif char.lower()=='c':
            rc.append('g')
        elif char.lower()=="n":
            rc.append("n")
        else:
            print('ERROR',char)

    return("".join(rc)[::-1].upper())

def extract_barcode_R3(R3):
    
    potential_barcodes = []
    
    reverse_seq = R3[1][::-1].upper()
    
    for b in np.arange(0,len(reverse_seq)-10,1):
        
        query = reverse_seq[b:b+7]
        
        if hamming_distance(query,'AAAAAAA')<=1:
            
            potential_barcode_1 = reverse_complement(reverse_seq[b+8:b+8+12][::-1])
            potential_barcode_2 = reverse_complement(reverse_seq[b+7:b+7+12][::-1])
            potential_barcode_3 = reverse_complement(reverse_seq[b+6:b+6+12][::-1])
            
            potential_barcodes= [potential_barcode_1,potential_barcode_2,potential_barcode_3]
            
            break
            
    return(potential_barcodes)

def extract_barcode_R1(R1):

    potential_barcodes = []

    reverse_seq = R1[1][::-1].upper()

    for b in np.arange(0,len(reverse_seq)-10,1):

        query = reverse_seq[b:b+8]

        if hamming_distance(query,'TTTTTTTT')<=1:

            potential_barcode_1 = reverse_seq[b-12:b][::-1]
            potential_barcode_2 = reverse_seq[b-11:b+1][::-1]
            potential_barcode_3 = reverse_seq[b-10:b+2][::-1]
            potential_barcode_4 = reverse_seq[b-9:b+3][::-1]

            potential_barcodes= [potential_barcode_1,potential_barcode_2,potential_barcode_3,potential_barcode_4]

            break

    return(potential_barcodes)


###################### MAIN #################################

barcodes = ref_df['Barcode'].tolist()
spacers = ref_df['Spacer'].tolist()
read_index=0
on=time.time()
total_time=0
R1_read=[]
R3_read=[]
assigned_reads,unassigned_reads,recombined_reads=(0,0,0)

with open(log_file,'a') as w:
    for l in [R1_file,R3_file,len(barcodes),len(spacers)]:
        w.write(str(l)+'\n')

with gzip.open(R1_file,'rt') as f1, gzip.open(R3_file,'rt') as f2:

    for x,y in zip(f1,f2):
        R1_read.append(x.rstrip())
        R3_read.append(y.rstrip())

        if len(R1_read)==4 and len(R3_read)==4 and R1_read[0].split()[0]==R3_read[0].split()[0]:

            matches=0
            valid_barcode=''
            for barcode in extract_barcode_R3(R3_read):
                if barcode in barcodes:
                    matches+=1
                    valid_barcode=barcode

            if matches==0:
                for barcode in extract_barcode_R1(R1_read):
                    if barcode in barcodes:
                        matches+=1
                        valid_barcode=barcode

            if matches==1:

                i = barcodes.index(valid_barcode)
                h_dist = hamming_distance(R1_read[1][1:1+20].upper(),spacers[i].upper())

                if h_dist<=3:
                    assigned_reads+=1

                    with gzip.open('demuxed/'+R1_file.split('-')[0]+'_'+str(i)+'_'+valid_barcode+'_R1.fq.gz','at') as x1:
                        for line in R1_read:
                            x1.write(line+'\n')

                    with gzip.open('demuxed/'+R1_file.split('-')[0]+'_'+str(i)+'_'+valid_barcode+'_R3.fq.gz','at') as x1:
                        for line in R3_read:
                            x1.write(line+'\n')
        
                else:
                    recombined_reads+=1

                    with gzip.open('demuxed/'+R1_file.split('-')[0]+'_recombined_R1.fq.gz','at') as x1:
                        for line in R1_read:
                            x1.write(line+'\n')

                    with gzip.open('demuxed/'+R1_file.split('-')[0]+'_recombined_R3.fq.gz','at') as x1:
                        for line in R3_read:
                            x1.write(line+'\n')

            else:
                unassigned_reads+=1

                with gzip.open('demuxed/'+R1_file.split('-')[0]+'_unassigned_R1.fq.gz','at') as x1:
                    for line in R1_read:
                        x1.write(line+'\n')

                with gzip.open('demuxed/'+R1_file.split('-')[0]+'_unassigned_R3.fq.gz','at') as x1:
                    for line in R3_read:
                        x1.write(line+'\n')


            R1_read=[]
            R3_read=[]
            read_index+=1

            if read_index%100000==0 and read_index!=0:
                interval = time.time()-on
                print(read_index,interval,'s')
                total_time+=interval
                on=time.time()

                with open(log_file,'a') as w:
                    for l in [read_index,'time elapsed: '+str(total_time/3600)+ ' h']:
                        w.write(str(l)+'\n')
                

print(read_index,assigned_reads,recombined_reads,unassigned_reads)

with open(log_file,'a') as w:
    for l in [read_index,assigned_reads,recombined_reads,unassigned_reads,'complete']:
        w.write(str(l)+'\n')

