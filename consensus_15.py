import re
import re
import sys
import string
from tb2 import tb

''' 
This script summarizes mutation counts per read family (group of reads sharing same primary barcode).
####
Input table format:

# Each BC5 family can be represented by one or more lines, where each line stores statistics of a single mutation.   
# WT families (not containing any mutation in any of the reads) are represented by a single row containing whole family statistics.
# If specific BC5 group includes both reads with and without mutations, then the WT statistics are stored in one of the lines, and the other line/s store mutation/s statistics.
    
# The input table contains the following columns
   
   * ID                             = BC5 string
   * Mutation                       = mutation string
   * Mutation_count                 = reads counts with mutations
   * HQ_mutation_count              = high-quality reads counts with mutations
   * Total_count                    = total reads counts
   * HQ_total_count                 = total high-quality reads counts
   * Mutation_freq                  = Mutation frequency calculated as HQ_mutation_count/Total_HQ_reads,
                                      where Total_HQ_reads denotes the total count of non-mutated reads in specific position.
   
   **  barcode3_IDs                 = BC3 strings types, linked to mutations, semicolon separated
   **  barcode3_HQ_mut_counts       = high quality reads counts (semicolon separated, corresponding to barcode3_IDs)
   **  barcode3_existing            = number of different types of BC3 strings, linked to mutations
   **  barcode3_missing             = Deprecated, showed amount of reads in the family with missing BC3
   **  barcode3_HQ_mut_freq         = frequency of mutations in each BC3 string, (semicolon separated, corresponding to barcode3_IDs)
   
   *** barcode3_all                 = all BC3 strings types
   *** barcode3_all_counts          = reads counts (semicolon separated, corresponding to barcode3_all)
   *** barcode3_all_WTcounts        = reads counts linked to WT positions, or reads without mutations, (semicolon separated, corresponding to barcode3_all)
   *** barcode3_all_WTcount         = number of different types of BC3 strings linked to WT positions, or reads without mutations.
   *** barcode3_missing_all_WTcount = Deprecated
   
   For lines that store WT statistics (read families with no mutations):
   barcode3_all_WTcounts, barcode3_all_WTcounts, barcode3_missing_all_WTcount are **NA**. Their equivalents are: 
   
   * barcode3_all_WTcounts         ->  barcode3_HQ_mut_counts
   * barcode3_all_WTcount          ->  barcode3_existing
   * barcode3_missing_all_WTcount  =  Deprecated
   
'''

#### Functions ####
# Convert 'NA' values to '0'
def float1(x):
    if x=='NA': x=0
    return(float(x))

# Calculate frequency of WT reads in read families w/ mutations
def wt_freq(positions,freqs):
    pos_freq={}
    if len(positions)==len(freqs):
        for i in range(len(positions)):
            if pos_freq.get(positions[i],"")=="": pos_freq[positions[i]]= float1(freqs[i])
            else: pos_freq[positions[i]] = pos_freq[positions[i]] + float1(freqs[i])
        wt_freq1=[ 1.0 - pos_freq[positions[i]] for i in range(len(positions))]
    else: sys.exit('Missing frequencies for some mutation positions')
    return wt_freq1

# Return counts of 3' barcodes associated with number of reads above the defined cutoff (default = 1)
def casesAboveMinimunCount(countsStr,minCount):
    def int1(x): return(int(x) if str(x)!="NA" else 0) # expecting "NA" when no HQ mutation is found
    counts1 = map(int1,countsStr.split(";"))
    return sum([c > minCount for c in counts1])    

# Summarize consensus statistics (mutation and WT counts) per BC5 read family
def consensus(m_,translationStartSite,BC3_min=1):
    t=tb(cols)
    m_1 = t.ne(m_,'Mutation',"WT")
    m_2 = t.eq(m_,'Mutation',"WT")
    HQmutations_found=False
    
    # Check if high-quality mutations are found in the read family
    if(len(m_1) > 0):
        if sum(map(int,t.sel(m_1,'HQ_mutation_count'))) > 0:
            HQmutations_found=True
    
    # Case: HQ mutations were found
    if HQmutations_found:
        
        # Extract positions and names of the HQ mutations
        positions_to_sort=[int(re.match(r'.*?(\d+).*?',x).group(1)) for x in t.sel(m_1,'Mutation')]
        m_1 = t.ins(m_1,'pos',positions_to_sort)
        m_1 = t.sort(m_1,'pos')
        mut = [(re.match(r'.*?([A-Za-z\-]+).*?',x).group(1)) for x in t.sel(m_1,'Mutation')]
        pos=t.sel(m_1,'pos')
        
        # Check and output cases of two different mutations at same position
        checkConflicts=[True if (int(pos[c])==int(pos[c-1])) else False for c in range(1,len(pos))] 
        if sum(checkConflicts) != 0:
             print "Two types of mutations at the same position in \"%s\""%(t.sel(m_1,'ID')[0])
        
        # Generate consensus statistics
        # In case we want to process only mutations meeting specific conditions - we can place the statement 'procede1=True' under required condition and return empty data frame if condition is not met 
        procede1=True
        if procede1:
            # Aggregate mutation names, their counts and frequencies; aggregate total read count and HQ read counts per mutated position
            h = ";".join(t.sel(m_1,'Mutation'))
            total_count = t.sel(m_1,'Total_count')[0]
            counts_hq = ";".join(map(str,t.sel(m_1,'HQ_mutation_count')))
            total_hq =  ";".join(map(str,t.sel(m_1,'HQ_total_count')))
            
            freqs_hq0 = t.sel(m_1,'Mutation_freq')
            freqs_hq =  ";".join(map(str,freqs_hq0))
            wt_freqs1 = ";".join(map(str,wt_freq(t.sel(m_1,'pos'),freqs_hq0)))
            
            if sum([True for zzz in freqs_hq0 if ((float(zzz) < 0.0) or (float(zzz) > 1.0))])>0:
                print "Frequency value outside [0,1] range for ID=\"%s\""%(t.sel(m_1,'ID')[0]) 
            
            # Aggregate BC3 statistics
            barcode3_existing =   ";".join(map(str,t.sel(m_1,'barcode3_existing'))) # Counts of different BC3 types linked to each mutation (or for all reads if the reads are all WT)
            barcode3_missing =    ";".join(map(str,t.sel(m_1,'barcode3_missing'))) # Deprecated
            barcode3_all_WTcount = ";".join(map(str,t.sel(m_1,'barcode3_all_WTcount'))) # Counts of different BC3 types linked to WT reads in each position where mutation was found
            barcode3mut_above1reads = ";".join([str(casesAboveMinimunCount(c,BC3_min)) for c in t.sel(m_1,'barcode3_HQ_mut_counts')]) # Counts of different BC3 types linked to each mutation and represented by amount of reads equal or greater than the BC3_min cutoff
            barcode3WT_above1reads  = ";".join([str(casesAboveMinimunCount(c,BC3_min)) for c in t.sel(m_1,'barcode3_all_WTcounts')]) # Counts of different BC3 types linked to WT reads in each position where mutation was found and represented by amount of reads equal or greater than the BC3_min cutoff
            barcode3_above1reads    = ";".join([str(casesAboveMinimunCount(c,BC3_min)) for c in t.sel(m_1,'barcode3_all_counts')]) # Counts of different BC3 types (either mutation or WT linked) in each position where mutation was found and represented by amount of reads equal or greater than the BC3_min cutoff
            
            # Add mutation position relative to TSS (translation start site):
            TSS = translationStartSite.split(",")
            pos_tss = [int(TSS[0]) + (int(TSS[1])* (int(p)-1)) for p in t.sel(m_1,'pos')]

            t.ins(m_1,'pos_TSS',pos_tss)
            haplotype_TSS=[]
            for x in range(len(mut)):
                pos0 = pos_tss[x]
                mut0 = mut[x]
                
                if int(TSS[0]) < 0 and pos0 >= 0: pos0 = pos0 + 1
                if int(TSS[0]) > 0 and pos0 <= 0: pos0 = pos0 - 1
                
                if int(TSS[1]) < 0: 
                    mut0 = mut0.translate(string.maketrans("ATGC","TACG"))
                
                # Uncomment the lines below to add i,d labels (insertions/deletions)    
                #if re.match(r'[A-Za-z]+-',mut0):   mut0=re.sub(r'-','d',mut0)
                #elif re.match(r'-[A-Za-z]+',mut0): mut0=re.sub(r'-','i',mut0)
                haplotype_TSS.append("%s%s"%(pos0,mut0))
            h_tss = ";".join(haplotype_TSS)
            
            # Add to data frame: 
            # 'Barcode', 'Consensus', 'Consensus_TSS', 'Mutation_freqs', 'Mutation_counts', 'total_count','Positions',
            # 'WT_freqs_in_mutations', 'barcode3_existing', 'barcode3mut_above1reads', 'barcode3WT_above1reads', 'barcode3_above1reads', 'barcode3_all_WTcount'
            df=[t.sel(m_1,'ID')[0],h,h_tss,freqs_hq,counts_hq,total_count,len(m_1),
                wt_freqs1,barcode3_existing,barcode3mut_above1reads,barcode3WT_above1reads,barcode3_above1reads,barcode3_all_WTcount]
        
        else:
            df=[]
            print "skipped ID=\"%s\" - ambiguous mutations"%(t.sel(m_1,'ID')[0])
    
    # Case: WT reads, but not HQ mutations were found in the read family
    elif(len(m_2) > 0):
        # if the data refers to a WT group, then in the input tables mutations statistics (frequency,count) should include zero values. The total count column should include the total count of reads.
        barcode3_existing =        ";".join(map(str,t.sel(m_2,'barcode3_existing'))) # Total count of different BC3 types in read family
        barcode3_missing =         ";".join(map(str,t.sel(m_2,'barcode3_missing'))) # Deprecated
        barcode3mut_above1reads =  ";".join([str(casesAboveMinimunCount(c,BC3_min)) for c in t.sel(m_2,'barcode3_HQ_mut_counts')]) # Count of different BC3 types linked to WT reads represented by amount of reads equal or greater than the BC3_min cutoff
        barcode3_above1reads =     ";".join([str(casesAboveMinimunCount(c,BC3_min)) for c in t.sel(m_2,'barcode3_all_counts')]) # Counts of different BC3 types represented by amount of reads equal or greater than the BC3_min cutoff
        
        total_count = t.sel(m_2,'Total_count')[0] 
        freqs_hq =  ";".join(map(str,t.sel(m_2,'Mutation_freq')))
        
        # Add to data frame: 
        # 'Barcode', 'Consensus', 'Consensus_TSS', 'Mutation_freqs' = 0, 'Mutation_counts' = 0, 'total_count','Positions' = 0,
        # 'WT_freqs_in_mutations', 'barcode3_existing', 'barcode3mut_above1reads' = NA, 'barcode3WT_above1reads' = 'barcode3_above1reads', 'barcode3_above1reads', 'barcode3_all_WTcount' = 'barcode3_existing'
        df=[t.sel(m_2,'ID')[0],'WT','WT',0,0,total_count,0,
            freqs_hq,barcode3_existing,'NA',barcode3_above1reads,barcode3_above1reads,barcode3_existing]
    # Case: no HQ mutations and no WT reads were found - skip read family
    else:
        df=[]
        print "skipped ID=\"%s\""%(t.sel(m_1,'ID')[0] if(len(m_1) > 0) else t.sel(m_2,'ID')[0])
    
    # Add columns: 'unmutatedReads_counts','unmutatedReads_freq'    
    if len(df) > 0:
        if(len(m_2) > 0):
            if len(m_2)!=1: sys.exit('Err: read family should contain a single WT entry for \'%s\''%(df[0]))
            
            freq_hq2 =  t.sel(m_2,'Mutation_freq')[0]
            count_hq2 = t.sel(m_2,'HQ_mutation_count')[0]
            total_count2 = t.sel(m_2,'Total_count')[0]
            
            if HQmutations_found:
                if (int(total_count2) != int(df[5])):
                    sys.exit('Err in \'%s\': Total number of reads in the family should not differ between WT and variant read data'%(df[0]))
            
            # Add WT read statistics
            df = df + [count_hq2,freq_hq2]
        else:
            df = df + [0,0] # Add '0' if no WT reads are present in the read family
    
    return(df)

#### Main ####
if __name__ == '__main__':
    # Input parameters
    in1 = sys.argv[1]
    out1 = sys.argv[2]
    translationStartSite = sys.argv[3]
    filter_pos = sys.argv[4]
    BC3_min = int(sys.argv[5])
    
    # Input column names
    cols={'ID':0,'Mutation':1,'Mutation_count':2,'HQ_mutation_count':3,
     'Total_count':4,'HQ_total_count':5,'Mutation_freq':6,
     'barcode3_IDs':7, 'barcode3_HQ_mut_counts':8,'barcode3_existing':9,'barcode3_missing':10,'barcode3_HQ_mut_freq':11,
     'barcode3_all':12,'barcode3_all_counts':13,'barcode3_all_WTcounts':14,'barcode3_all_WTcount':15,'barcode3_missing_all_WTcount':16}
    
    # Output column names
    out_cols=['Barcode', 'Consensus', 'Consensus_TSS', 'Mutation_freqs', 'Mutation_counts', 'total_count','Positions',
              'WT_freqs_in_mutations','barcode3_existing','barcode3mut_above%ireads'%(BC3_min),'barcode3WT_above%ireads'%(BC3_min),'barcode3_above%ireads'%(BC3_min),'barcode3_all_WTcount', 'unmutatedReads_counts','unmutatedReads_freq']    
    
    # Open I/O files
    t1=open(in1,'r')
    o1=open(out1,'w')
    
    prev_id=''
    barcode_group=[]
    visited={}
    o1.write("\t".join(out_cols)+"\n")
    
    # Iterate over mutation table to collect variants per read family
    for i,r in enumerate(t1):
        # Test that input table contains expected column names in expected order
        if i==0:
            x=r.rstrip().split("\t")
            for y in range(len(x)):
                if cols.get(x[y],"")=="": sys.exit('Unexpected column header %s'%(x[y]))
                elif str(cols[x[y]]) != str(y): sys.exit('Wrong column header position: %s instead of %s'%(y, cols[x[y]]))
        
        else:
            x=r.rstrip().split("\t")
            
            # Mutation at filtering position is used only to identify PCR artifacts, and has no research value, hence it is not added to mutation summary table
            if (re.match(filter_pos, x[cols['Mutation']])):
                print "%s: %s"%(i,r)
                continue
            
            # Find mutation profile for read families
            if len(x)!=len(cols): sys.exit('Missing data in line %s'%(i))
            if (prev_id != x[cols['ID']]) and (len(barcode_group)>0):
                c = consensus(barcode_group,translationStartSite,BC3_min)
                if len(c) > 0:
                    o1.write("\t".join(map(str,c))+"\n")
                else:
                    pass # print "ID=\"%s\" was skipped"%(prev_id)
                barcode_group=[]
                if visited.get(x[cols['ID']],"")!="": sys.exit('Unexpected barcode %s in line %s, input should be sorted'%(x[cols['ID']], i))
                visited[x[cols['ID']]]=True

            barcode_group.append(x)
            prev_id = x[cols['ID']]
        #if i > 10000: break # Uncomment this line to perform test run on first 10,000 read families in the input
        
    # Find mutation profile for the last read family in the input    
    c = consensus(barcode_group,translationStartSite,BC3_min)
    if len(c) > 0:
        o1.write("\t".join(map(str,c))+"\n")
        barcode_group=[]
        t1.close()
        o1.close()
    else:
        pass #print "ID=\"%s\" was skipped"%(prev_id)

print "Done" # Check for cluster runs, to validate that the run was completed and the script did not silently die mid-run
