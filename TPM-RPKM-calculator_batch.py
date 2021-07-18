### Script to calculate TPM (Transcripts Per Kilobase Million) and RPKM (Reads Per Kilobase Million)
### values of genes to normalize for sequencing depth and gene length
### written by Feng JU (richieju@eawag.ch)

# convert coverage to read abundance
# column i is average coverage while column j the reference length
def summarize_table(fn,i,j):
    from collections import OrderedDict
    assert type(i) is int, "name is not an int: %r" % i
    assert type(j) is int, "name is not an int: %r" % i

    k = 0
    d_cov, d_RL, d_rabu, d_RPK = OrderedDict(),{},{},{}
    
    for l in open(fn,'r'):
        k+=1
        try:
            lis = l.strip().split(',')
            ID  = lis[0]
            cov = float(lis[i-1])  # average coverage
            RL  = int(lis[j-1])/1000.0  # reference sequence length (kbp)
            rl  = 150 # read length, bps
            rabu = 1000*cov*RL/rl
            RPK = rabu/RL
            
            d_cov[ID] = cov
            d_RL[ID] = RL
            d_rabu[ID] = rabu
            d_RPK[ID] = RPK 
        except IndexError:
            print "Ignore line %s: not enough elements" % str(k)
        except ValueError:
            print "Ignore line %s: not float or int type" % str(k)

    return d_cov, d_RL, d_rabu, d_RPK
        

#1. Divide the read counts by the length of each gene in kilobases. This gives you reads per kilobase (RPK).
#2. "per million" scaling factor is calculated as the sum of all the RPK values in a sample divided by 1,000,000 
#3. Divide the RPK values by the "per million" scaling factor. This gives you TPM.
def TPM(d_RPK):
    #per million scaling factor for calculating TPM
    SF_TPM = sum(d_RPK.values())/1000000
    d_TPM = {k: v/SF_TPM for k, v in d_RPK.items()}
    return d_TPM


#1."per million" scaling factor is calculated as total number of mapped reads divided by 1,000,000
#2. Divide the read counts by the "per million" scaling factor. This normalizes for sequencing depth, giving you reads per million (RPM)
#3. Divide the RPM values by the length of the gene, in kilobases. This gives you RPKM.
def RPKM(d_rabu, d_RL):
    #per million scaling factor for calculating RPKM
    SF_RPKM = sum(d_rabu.values())/1000000
    d_RPM = {k: v/SF_RPKM for k, v in d_rabu.items()}
    d_RPKM = {k: (v/d_RL[k]) for k, v in d_rabu.items()}
    return d_RPKM

def making_folder(foldername):
    if os.path.exists(foldername):
        for root, dirs, files in os.walk(foldername):
            for name in files:
                os.remove(os.path.join(root,name))
    else:
        os.mkdir(foldername)

if __name__ == '__main__':
    import os
    foldername = 'ORFs_RNA_coverage'
    making_folder(foldername + '-TPM-RPKM')
    for root,dirs,files in os.walk(foldername):
        for file in files:
            print '------','Processing',file,'in prograss','------'
            d_cov, d_RL, d_rabu, d_RPK = summarize_table(os.path.join(root, file),2,3)
            d_TPM = TPM(d_RPK)
            d_RPKM = RPKM(d_rabu, d_RL)

            f=open(foldername + '-TPM-RPKM'+'/'+file,'w')
            f.write(','.join(['ID', 'TPM','RPKM','Coverage','Ref-length','Read-abu'])+'\n')
            for key in d_cov.keys():
                out_lis = map(str,[key,d_TPM[key],d_RPKM[key],d_cov[key],d_RL[key],d_rabu[key]])
                f.write(','.join(out_lis)+'\n')
        
print 'DONE!'
    


        
        
        
    
        
        
