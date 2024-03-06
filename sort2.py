import csv
import sys
import re
import os
import argparse
from Bio import SeqIO

# This script sorts input reads based on their gene of origin; ambiguous reads are put in to "others" category

#### Functions ####
# "Main" function - main code for read sorting
def main1(print1=True):
    # Return float(x) if x is in the defined range
    def float1(x):
        x=float(x)
        if (x < 0.0) or (x > 1.0): raise argparse.ArgumentTypeError("%f not in range [0.0, 1.0]"%(x))
        return x
    
    # Gather input arguments    
    parser = argparse.ArgumentParser(description='ort',formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--fastq', action="store", dest="fastq",required=True)
    parser.add_argument('--pos', action="store", dest="pos",required=True)
    parser.add_argument('--nucl', action="store", dest="nucl",required=True)
    parser.add_argument('--refs', action="store", dest="refs",required=True)
    parser.add_argument('--ref', action="store", dest="ref",required=True)
    parser.add_argument('--dirOut', action="store", dest="dirOut",required=True)
    parser.add_argument('--title', action="store", dest="title",required=True)
    parser.add_argument('--match_percentage',action="store", dest="match_percentage",required=True,type=float1)
    parser.add_argument('--offset',action="store", dest="offset_range",required=True,type=int)
    p=parser.parse_args()
    print p
    
    # Read the input
    fq0=open(p.fastq,'r')
    fq1=filter(lambda row: row!="\n",fq0)
    temp1=[]
    fastq_handles={}
    
    # Iterate over input ".fastq" files to sort the reads
    for i,r in enumerate(fq1):
        r=r.rstrip()
        j=(i+1)%4
        if j==1: # Fastq header
            if len(temp1)!=0: sys.exit('Error: read array should be empty before reading new read')
            if r[0:1]!='@': sys.exit('Error in header format: doesn\'t start with @')
            temp1.append(r)
        elif j==2: # Fastq DNA sequence
            if len(temp1)!=1: sys.exit('Read array should contain only header before reading read sequence')
            temp1.append(r)
        elif j==3: # Fastq "+" tag
            if len(temp1)!=2: sys.exit('Read array should contain only header and sequence before reading fastq tag')
            if r[0:1]!='+': sys.exit('Err in fastq format: expected to find \'+\' tag before quality string')
            temp1.append(r)
        elif j==0: # Fastq quality string
            if len(temp1)!=3: sys.exit('Read array should contain only header, sequence and "+" tag before reading quality string')
            temp1.append(r)
            
            sort1(p,temp1,fastq_handles) # Sort the read and its associated quality by gene of origin
            temp1=[]
        
        if print1:
            if (i+1)%4000000==0: print str((i+1)/4)
        #if i>1000: break # Uncomment this line to perform test run on the first 1000 entries of the ".fastq" file
    
    # Close I/O after analysis completion
    fq0.close()
    for k in fastq_handles.keys():
        print "closing \'%s\'"%(k)
        fastq_handles[k].close()

# Sort input reads based on their gene of origin
def sort1(p,fastq_record,fastq_handles):
    seq1 = fastq_record[1] # Read sequence
    qual1 = fastq_record[3] # Read sequencing quality string
    header1 = fastq_record[0] # Read header
    
    # Gather sorting parameters - position and identity of sorting nucleotides
    # ";" - gene separator, here and in all parameters below
    nucl0 = p.nucl.split(";")
    pos0  = p.pos.split(";")
    refs = p.refs.split(";")
    ref = p.ref.split(";")
    match_percentage = p.match_percentage
    
    # Test correctness of the sorting parameters
    if ((len(nucl0)!=len(pos0))or(len(nucl0)!=len(refs))): sys.exit("Error: mismatch between number of sorting nucleotides and number of sorting positions")
    if (len(ref)!=1): sys.exit("\"ref\" field in factors_table can contain only single entry")
    if ((len(nucl0)==0) or (len(pos0)==0)): sys.exit("Error: either sorting nucleotides or sorting positions field is empty")
    
     # Initialize analysis variables
    sort_to_ref=ref[0]+".others" # Reads would be sorted to "others" file by default if gene of origin is not identified
    put_in_logfile=False
    test_sort_to_ref=[] #
    test_put_in_logfile=[] #
    
    # Iterate over different sorting genes to check if given read matches to any of them
    for i in range(len(pos0)):
        '''
        # Haplotypes separator - "&".
        # The rows below can be uncommented (instead of rows 104, 105) to run pipeline with haplotypes: 
        #nucl1 = nucl0[i].split("&") 
        #pos1  = pos0[i].split("&")
        '''
        nucl1 = [nucl0[i]]
        pos1  = [pos0[i]]  
        
        if (len(nucl1)!=len(pos1)) or (len(pos1)==0): sys.exit("Mismatch between number of haplotype nucleotides and haplotype positions!")
        
        z=0
        haplotypes_with_offsets=0
        
        # Visit a haplotype (here haplotype means one or more alleles with know bp distance between each other)
        for j in range(len(pos1)):
            nucl2 = nucl1[j].split(",")
            pos2  = map(int,pos1[j].split(","))
            if (len(nucl2)!=len(pos2)) or (len(pos2)==0): sys.exit("Mismatch between amount of sorting position and sorting nucleotides")
            
            # To account for potential indels in the reads, offset is added to sorting position when searching for sorting nucleotides
            offset_range = p.offset_range
            stop = 0;
            
            for offset0 in range(offset_range+1):
                if stop > 0: break
                
                #print "Absolute offset is %i"%(offset0)
                for modifier in [-1,1]:
                    offset = int(offset0*modifier)
                    #print "Checking offset %i"%(offset)
                    zz=0
                    # Check how many read nucleotides at sorting positions (while accounting for the offset) match sorting nucleotides
                    for k in range(len(pos2)):
                        pos3=pos2[k]
                        nucl3=nucl2[k]
                        s = seq1[(pos3-1+offset):(pos3+offset)]
                        if s==nucl3: zz+=1
                    
                    # If percent of matching nucleotides higher than user-defined cutoff, mark the read as sorted
                    if (float(zz)/float(len(pos2))) >= float(match_percentage):
                        #print "Match found!"
                        z+=1
                        if offset!=0: haplotypes_with_offsets+=1
                        stop = 1
                        break
        
        # Indicate potential origin genes of the read
        if z==len(pos1):
            test_sort_to_ref.append(refs[i])
            if haplotypes_with_offsets>0: test_put_in_logfile.append(True)
            else: test_put_in_logfile.append(False)
    
    # Only if read sequence matches a single reference gene we know where to sort it; Otherwise it will go to the default "others" file
    if len(test_sort_to_ref)==1:
        sort_to_ref = test_sort_to_ref[0]
        put_in_logfile=test_put_in_logfile[0]
    
    # Write sorted read to the output file and create new output file if file doesn't exist
    # Uncomment "withIndels.log" lines below to create a separate file with reads that were sorted with an offset   
    title2="%s.%s"%(p.title,sort_to_ref)
    
    if fastq_handles.get(title2,"")=="":
        fastq_handles[title2] = open("%s/%s.fastq"%(p.dirOut,title2),'w')
        #fastq_handles[title2+".withIndels.log"] = open("%s/%s.withIndels.log"%(p.dirOut,title2),'w')
        newline1=""
    else: newline1="\n"
    
    fastq_handles[title2].write("%s%s"%(newline1,"\n".join([header1,seq1,'+',qual1]))) # Write sorted reads
    
    #if put_in_logfile:
        #fastq_handles[title2+".withIndels.log"].write("%s%s"%(newline1,"\n".join([header1,seq1,'+',qual1])))

#### Main_run ####
main1()
