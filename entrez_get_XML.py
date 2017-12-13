import Bio
import time
from Bio import Entrez
import os
    
if __name__ == '__main__':
    start_time = time.time()
    
    ##################################
    #INPUTS:
    #this file contains all Ref Seq IDs to be downloaded
    target_file = "seq_input/M20160921ACFE4208EAFA842A78A1B3BA7138A93DB4D25CK.list"
    Entrez.email = "plesa@ucla.edu"     # Always tell NCBI who you are
    entrez_dir = "folA_entrez/"
    download_on = True
    
    ##################################
    f = open(target_file, 'r')
    
    down_count = 1
    for line in f:
        #get AccID
        accID = line.rstrip('\n')
        
        if download_on:
            filename = entrez_dir+accID+".xml"
            if not os.path.isfile(filename):
                # Downloading...
                net_handle = Entrez.efetch(db="protein",id=accID,retmode="xml")
                out_handle = open(filename, "w")
                out_handle.write(net_handle.read())
                out_handle.close()
                net_handle.close()
                print("Saved seq # "+str(down_count))
                time.sleep(0.25)
            down_count += 1
    print(str(down_count)+" seqs total")
    print("--- %s seconds ---" % (time.time() - start_time))
