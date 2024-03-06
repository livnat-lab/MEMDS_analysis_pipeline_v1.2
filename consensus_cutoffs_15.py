
import sys
import re
import math
from numpy import double

# This script compares summary table of mutations per read family against a set of cutoff criteria to determine which mutations should be considered as "true" 

''' 
Input table format:

* = data corresponding to mutation positions in different BC3 groups within given read family

Barcode                   = BC5, a unique identifier of each family
* Consensus               = semicolon separated mutations
* Consensus_TSS           = semicolon separated mutations, with position relative to the translation start site        
* Mutation_freqs       = semicolon separated HQ mutation frequencies
* Mutation_counts         = semicolon separated HQ mutation read counts
total_count               = total read counts
* Positions               = number of mutations
* WT_freqs_in_mutations   = semicolon separated frequency of WT in positions of mutations
* barcode3_existing       = semicolon separated number of different types of BC3 linked to each mutation position (or among all reads, if Barcode is WT)
* barcode3mut_above1reads = semicolon separated number of different types of BC3 linked to each mutation position and represented by at least 2 reads
* barcode3WT_above1reads  = semicolon separated number of different types of BC3 linked to WT reads in each mutation position and represented by least 2 reads
* barcode3_above1reads    = semicolon separated number of different types of BC3 linked to any kind of read, represented by least 2 reads each
* barcode3_all_WTcount       = semicolon separated number of different types of BC3 linked to WT, in mutation position
unmutatedReads_counts       = count of reads with no mutations at all
unmutatedReads_freq       = frequency of reads with no mutations at all from all reads
'''



###### Functions #####
# Substitute 'N/A' values by '0'
def float1(x):
    if x=='NA': x=0
    return(float(x))

# Select mutations passing the cutoff criteria. Substitute mutations at ambiguous positions by 'N's.
'''
# mut_name_ - an array of mutation names; x - mutation information to update
# toMut,toN - mutated positions that pass cutoff criteria and WT positions that do not pass the former, respectively
# writeN - check if add 'N' at ambiguous positions, default: false
# add - check if multiple 'N's should be displayed at same position,default: true
# maxVal=False - Deprecated
'''

def select1(mut_name_,x,toMut,toN,writeN=False,add=True,maxVal=False): 
    mut_val=[]
    mut_name=[]
    pos1=positions(mut_name_)
    pos_accepted={} # Tells which positions have accepted mutations, since in such positions ambiguous mutations should not be reported
    
    # Create a dictionary of mutated positions passing the cutoff
    for i in range(len(mut_name_)):
        pos_accepted[pos1[i]]=False
    for i in range(len(mut_name_)):
        if toMut[i]: pos_accepted[pos1[i]]=True
        
    for i in range(len(x)):
        # Mutation is accepted, add its information
        if toMut[i]: 
            mut_val.append(x[i]) 
            mut_name.append(mut_name_[i])
        
        # Ambiguous position: both mutation and WT are not accepted at given position    
        elif toN[i] and not(pos_accepted[pos1[i]]): 
            mut_name0 = (re.sub('[A-Za-z\-]+','N',mut_name_[i])) # Change mutation name to 'N'
            if writeN: x0=re.sub('[A-Za-z\-]+','N',x[i]) # Change mutation info to 'N', if writeN=True
            else:      x0=x[i]
            
            # Add updated information to mutation information array, while collapsing multiple 'N's at same position
            if (len(mut_name)==0) or (mut_name[len(mut_name)-1] != mut_name0):
                mut_val.append(x0)
                mut_name.append(mut_name0)
            elif add:
                if maxVal:
                    pass
                else:
                    mut_val[len(mut_val)-1] += x0
    return mut_val

# Extract mutation position from variant name (e.g. 39CT -> 39)
def positions(mut_name_):
    positions1=[]
    for i in range(len(mut_name_)):
        m1=re.match('.*?(\d+).*?',mut_name_[i])
        if m1: positions1.append(int(m1.group(1)))
        else: sys.exit('Err: unclear mutation position in %s!'%(mut_name_[i]))
    return(positions1)

# Check if x is 'NA' and return its value as integer if not 'NA'
def intIfNotNA(x):
    if x=='NA': return('NA')
    else: return(int(x))

####### Main ########
# Collect input parameters
f1=      sys.argv[1]
outIdx=  sys.argv[2]
second_barcode_size = int(sys.argv[3])

# Organize I/O table headers
col = {'Barcode':0, 'Consensus':1, 'Consensus_TSS':2, 'Mutation_freqs':3, 'Mutation_counts':4,
       'total_count':5,'Positions':6,'WT_freqs_in_mutations':7,'barcode3_existing':8, 'barcode3mut_above1reads':9,
       'barcode3WT_above1reads':10,'barcode3_above1reads':11,'barcode3_all_WTcount':12,'unmutatedReads_counts':13,'unmutatedReads_freq':14}

col_out = col

col_names     = [z[0] for z in sorted(col.items(),     key=lambda x: x[1])]
col_names_out = [z[0] for z in sorted(col_out.items(), key=lambda x: x[1])]

# Define cutoff criteria
min_freq=map(float, sys.argv[4].split(","))
min_count=map(int, sys.argv[5].split(","))
if second_barcode_size > 0:
    minCountBarcode = [int(sys.argv[6])]
    bc3groupCountOK=int(sys.argv[7])
    bc3groupsWithReadsAbove1_OK=int(sys.argv[8])
    
    print("min mut freq = {}".format(min_freq))
    print("min read count = {}".format(min_count))    
    print ("min count 3'BC = {}".format(minCountBarcode))
    print "bc3groupCountOK=%i bc3groupsWithReadsAbove1_OK=%i"%(bc3groupCountOK,bc3groupsWithReadsAbove1_OK)
    
else:
    minCountBarcode = [0]
    bc3groupCountOK=0
    bc3groupsWithReadsAbove1_OK=0

    print("min mut freq = {}".format(min_freq))
    print("min read count = {}".format(min_count))    
    print ("min count 3'BC = {}".format(minCountBarcode))
    print "bc3groupCountOK=%i bc3groupsWithReadsAbove1_OK=%i"%(bc3groupCountOK,bc3groupsWithReadsAbove1_OK)

# Create output files
out_files={}
out_files_names={}
out_files2={}
rejected_cases={}
for f in range(len(min_freq)):
    for c in range(len(min_count)):
        for bc in range(len(minCountBarcode)):
            name1="mutFreq%s_readCount%i_BC3WithMut%i_BC3above%i"%(str(round(min_freq[f],2)),min_count[c],minCountBarcode[bc],bc3groupsWithReadsAbove1_OK)
            out_files_names[name1] = "%s.%s.txt"%(outIdx,name1)
            out_files[name1]=open("%s.%s.txt"%(outIdx,name1),'w')
            out_files[name1].write(("\t".join(col_names_out))+"\n")
            out_files2[name1]=open("%s.%s.cons-count.txt"%(outIdx,name1),'w')
            out_files2[name1].write("consensus\tbarcodes_count\n")
            rejected_cases["%s.low_count"%(name1)]=0
            rejected_cases["%s.other"%(name1)]=0
            rejected_cases["%s.from_WThaplotype"%(name1)]=0

# Iterate over input file to compare mutation profile of the read families against cutoff criteria
h1=open(f1,'r')
for i,r in enumerate(h1):
    r1=r.rstrip().split('\t')
    
    # Check that input header contains expected column names at expected positions
    if i==0:
        if col_names != r1: sys.exit('Unexpected column names in %s\n'%(f1))
        for y in range(len(r1)):
            if col.get(r1[y],"")=="": sys.exit('Wrong column header %s at column %s'%(r1[y],y+1))
            elif str(col[r1[y]]) != str(y): sys.exit('Column header at wrong position: %s instead of %s'%(y, col[r1[y]]))
    
    # Parse input data
    else:
        # Extract mutation statistics per read family
        freqs0=   r1[col['Mutation_freqs']] # # High quality mutation frequencies
        freqs1=   map(float1,freqs0.split(";"))
        
        counts0=  r1[col['Mutation_counts']] # High quality mutation counts
        counts1 = map(int,counts0.split(";"))
        
        count =   int(r1[col['total_count']]) # Family size
        
        wt_freq_in_mut0 = r1[col['WT_freqs_in_mutations']] # WT read frequency at mutated positions
        wt_freq_in_mut1 = map(float1,wt_freq_in_mut0.split(";"))
        
        barcode3_existing0 = r1[col['barcode3_existing']] # Counts of different BC3 types associated with each mutation
        barcode3_existing1 = map(int,barcode3_existing0.split(";"))
        
        barcode3mut_above1reads0 = r1[col['barcode3mut_above1reads']] # Counts of different BC3 types associated with each mutation and represented by amount of reads equal to or greater than the BC3_min cutoff
        barcode3mut_above1reads1 = map(intIfNotNA,barcode3mut_above1reads0.split(";"))
        
        barcode3WT_above1reads0 = r1[col['barcode3WT_above1reads']] # Count of different BC3 types associated with WT reads at each mutated position and represented by amount of reads equal to or greater than the BC3_min cutoff
        barcode3WT_above1reads1 = map(int,barcode3WT_above1reads0.split(";"))
        
        barcode3_above1reads0   = r1[col['barcode3_above1reads']] # Count of different BC3 types at each mutated position represented by amount of reads equal to or greater than the BC3_min cutoff
        barcode3_above1reads1   = map(int,barcode3_above1reads0.split(";"))
        
        barcode3_all_WTcount0 = r1[col['barcode3_all_WTcount']] ## Counts of different BC3 types associated with WT reads at each mutated position 
        barcode3_all_WTcount1 = map(intIfNotNA,barcode3_all_WTcount0.split(";"))
        
        ####
        aboveMin=bc3groupsWithReadsAbove1_OK
        groupsMin=bc3groupCountOK
        
        # Check that no information is missing from mutation profile of read family 
        if (len(freqs1)!=len(counts1)) or (len(freqs1)!=len(wt_freq_in_mut1)) or (len(freqs1)!=len(barcode3_existing1)) \
            or (len(freqs1)!=len(barcode3mut_above1reads1) or (len(freqs1)!=len(barcode3WT_above1reads1)) \
            or (len(freqs1)!=len(barcode3_above1reads1)) or (len(freqs1)!=len(barcode3_all_WTcount1))):
            sys.exit('Missing parameters in line %s %s'%(i+1, r1[0]))
        
        
        if(min(counts1)<0):
            print("excluded:")
            print(r1)
        else:
            allels1 = r1[col['Consensus']].split(';')
            allels2 = r1[col['Consensus_TSS']].split(';')
            if (len(allels1)==len(freqs1)):
                for f in range(len(min_freq)):
                    for c in range(len(min_count)):
                        for bc in range(len(minCountBarcode)):
                            name1 = "mutFreq%s_readCount%i_BC3WithMut%i_BC3above%i"%(str(round(min_freq[f],2)),min_count[c],minCountBarcode[bc],bc3groupsWithReadsAbove1_OK)
                            
                            # Check if read family size passes the cutoff criteria, reject low size families
                            if (count >= min_count[c]):
                                
                                # Case: read family contains mutations
                                if r1[col['Consensus']] != 'WT':
                                    is_mut_accepted0    = [ frq >= min_freq[f] for frq in freqs1] # Check if mutation frequency passes cut-off criteria at each mutated position
                                    is_WT_accepted0     = [ frq >=  min_freq[f] for frq in wt_freq_in_mut1] # Check if WT read frequency at each mutated position passes cut-off criteria 
                                    
                                    # Check if the number of BC3 types associated with mutation or WT reads at each mutated position passes the cutoff criteria
                                    is_secondBC_accepted_mut =    [ (barcode3_existing1[u]    >= minCountBarcode[bc]) and
                                                                    ((barcode3mut_above1reads1[u] >=  aboveMin) or (barcode3_existing1[u] >= groupsMin))
                                                                    for u in range(len(barcode3_existing1)) ]
                                    
                                    is_secondBC_accepted_wt  =    [ (barcode3_all_WTcount1[u]    >= minCountBarcode[bc]) and
                                                                    ((barcode3WT_above1reads1[u] >=  aboveMin) or (barcode3_all_WTcount1[u] >= groupsMin))
                                                                    for u in range(len(barcode3_all_WTcount1)) ]
                                    
                                    # Check for each mutated position if mutation- or WT-linked reads are passing the combined cutoff criteria
                                    is_mut_accepted = [is_mut_accepted0[zzz] and  is_secondBC_accepted_mut[zzz] for zzz in range(len(is_mut_accepted0))]
                                    is_WT_accepted  = [is_WT_accepted0[zzz]  and  is_secondBC_accepted_wt[zzz]  for zzz in range(len(is_WT_accepted0))]
                                    is_WT_not_accepted = [not(zzz) for zzz in is_WT_accepted]
                                    
                                    '''
                                    Output of the cutoff check:
                                    1) Mutation name in positions with mutation above the combined cutoff criteria (e.g. 39CT).
                                    2) 'N's in positions where neither a mutation, nor 'WT' nucleotide are above combined cutoff criteria.
                                    3) 'WT', if no mutations or 'N's are reported, and 'WT' nucleotides pass the combine cutoff criteria in all positions with suspected mutations.
                                    
                                    4) If multiple mutations are accepted at the same position, then all are shown. 
                                    5) If multiple 'N's are declared at the same position, they are collapsed together (e.g.: instead of '49NN' -> '49N' would be shown).
                                    '''
                                    
                                    # Update mutation summary data to return accepted mutations or 'N's for ambiguous positions
                                    if (sum(is_mut_accepted) > 0) or (sum(is_WT_not_accepted) > 0): #
                                        r2 = r1[:]
                                        r2[col['Consensus']] =             ";".join(             select1(allels1,allels1,                 is_mut_accepted,is_WT_not_accepted,writeN=True,add=False))
                                        r2[col['Consensus_TSS']] =         ";".join(             select1(allels1,allels2,                 is_mut_accepted,is_WT_not_accepted,writeN=True,add=False))
                                        r2[col['Mutation_freqs']] =        ";".join(map(str,     select1(allels1,freqs1,                  is_mut_accepted,is_WT_not_accepted)))
                                        r2[col['Mutation_counts']] =       ";".join(map(str,     select1(allels1,counts1,                 is_mut_accepted,is_WT_not_accepted)))
                                        r2[col['WT_freqs_in_mutations']] = ";".join(map(str,     select1(allels1,wt_freq_in_mut1,         is_mut_accepted,is_WT_not_accepted)))
                                        r2[col['Positions']] =                          str( len(select1(allels1,allels1,                 is_mut_accepted,is_WT_not_accepted)))
                                        r2[col['barcode3_existing']] =      ";".join(map(str,    select1(allels1,barcode3_existing1,      is_mut_accepted,is_WT_not_accepted)))
                                        r2[col['barcode3mut_above1reads']] =";".join(map(str,    select1(allels1,barcode3mut_above1reads1,is_mut_accepted,is_WT_not_accepted)))
                                        r2[col['barcode3_all_WTcount']] =   ";".join(map(str,    select1(allels1,barcode3_all_WTcount1,   is_mut_accepted,is_WT_not_accepted)))
                                        
                                        r2[col['barcode3mut_above1reads']] =   ";".join(map(str,    select1(allels1,barcode3mut_above1reads1,is_mut_accepted,is_WT_not_accepted)))
                                        r2[col['barcode3WT_above1reads']] =    ";".join(map(str,    select1(allels1,barcode3WT_above1reads1, is_mut_accepted,is_WT_not_accepted)))
                                        r2[col['barcode3_all_WTcount']] =      ";".join(map(str,    select1(allels1,barcode3_all_WTcount1,   is_mut_accepted,is_WT_not_accepted)))
                                        r2[col['barcode3_above1reads']] =      ";".join(map(str,    select1(allels1,barcode3_above1reads1,   is_mut_accepted,is_WT_not_accepted)))
                                        
                                        out_files[name1].write(("\t".join(r2))+"\n") # Write updated mutation data to the output file
                                    
                                    # Update mutation summary data to present read families as 'WT' if all the mutations in the family are rejected, but their associated WT reads are all accepted
                                    elif sum(is_WT_not_accepted) == 0:
                                        r2 = r1[:]
                                        r2[col['Consensus']] =            "WT"
                                        r2[col['Consensus_TSS']] =        "WT"
                                        r2[col['Mutation_freqs']] =       r2[col['unmutatedReads_freq']]
                                        r2[col['Mutation_counts']] =      r2[col['unmutatedReads_counts']]
                                        r2[col['WT_freqs_in_mutations']] = 'NA'
                                        r2[col['Positions']] =            str(0)
                                        r2[col['barcode3_existing']] =  r2[col['barcode3_all_WTcount']]
                                        r2[col['barcode3_above1reads']] = 'NA'
                                        out_files[name1].write(("\t".join(r2))+"\n") # Write updated mutation data to the output file
                                    else:
                                        rejected_cases["%s.other"%(name1)] +=1 # Reject and count cases not fitting to any of the criteria above 
                                
                                # Case: read family has no high quality mutations
                                else:
                                    # Check that input data for WT families is valid
                                    if(len(barcode3_existing1)!=1): sys.exit('Err: for WT family single entry in \'barcode3_existing1\' column is expected!')
                                    if(len(barcode3WT_above1reads1)!=1): sys.exit('Err: for WT family single entry in \'barcode3WT_above1reads1\' column is expected!')
                                    
                                    # Check and output WT read families passing the combined criteria cutoff 
                                    if (barcode3_existing1[0] >= minCountBarcode[bc]) and ((barcode3WT_above1reads1[0] >=  aboveMin) or (barcode3_existing1[0] >= groupsMin)):
                                        r2 = r1[:]
                                        r2[col['Consensus']] =            "WT"
                                        r2[col['Consensus_TSS']] =        "WT"
                                        r2[col['Mutation_freqs']] =       r2[col['unmutatedReads_freq']]
                                        r2[col['Mutation_counts']] =      r2[col['unmutatedReads_counts']]
                                        r2[col['WT_freqs_in_mutations']] = 'NA'
                                        r2[col['Positions']] =            str(1)
                                        r2[col['barcode3_above1reads']] = 'NA'
                                        out_files[name1].write(("\t".join(r2))+"\n") # Write updated mutation data to the output file
                                    else:
                                        rejected_cases["%s.from_WThaplotype"%(name1)] +=1 # Reject and count cases of WT families not passing any of the cutoff criteria 
                            else:
                                rejected_cases["%s.low_count"%(name1)] +=1 # Reject and count cases of low size read families
            else: sys.exit('Length of variant array does not match length of mut_freq array in %i: %s'%(i+1, r1[0]))
        # if i>3000: break # Uncomment to perform test run of the script on the first 3000 lines of the input

h1.close()

for k in out_files.keys(): out_files[k].close()

# Write output to consensus count files summing counts of rejected families and number of families containing each mutation set
for k in out_files_names.keys():
    hap_bcCount={}
    h2=open(out_files_names[k],'r')
    for i,r in enumerate(h2):
        if i > 0:
            r1=r.rstrip().split('\t')
            hap=r1[col['Consensus']]
            if hap_bcCount.get(hap,"") == "": hap_bcCount[hap]=1
            else: hap_bcCount[hap] += 1
    h2.close()
    hap_bcCount1=sorted(hap_bcCount.items(), reverse=True,key=lambda k: k[1])
    out_files2[k].write("%s\t%i\n"%('rejected.low_count',rejected_cases["%s.low_count"%(k)]))
    out_files2[k].write("%s\t%i\n"%('rejected.from_WThaplotype',rejected_cases["%s.from_WThaplotype"%(k)]))
    out_files2[k].write("%s\t%i\n"%('rejected.other',rejected_cases["%s.other"%(k)]))
    for hap in hap_bcCount1:
        out_files2[k].write("%s\t%i\n"%(hap[0],hap[1]))

print "Done\n" # Check for cluster runs, to validate that the run was completed and the script did not silently die mid-run
