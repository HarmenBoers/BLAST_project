# -*- coding: utf-8 -*-

# This script downloads all the sequences in the list provided by the user, and puts them into a single file to be used as a database.
# This script also stores the fetched fasta sequences into individual files which will be used as queries for (PSI-)BLAST searches.

import sys
import urllib2
# This function fetch the fasta formatted sequence of the uniProtID supplied in the argument.
def fetchFasta(uniProtID):
    ##########################
    ### START CODING HERE ####
    ##########################
    # Define the variable 'url' as a string with the URL to the Fasta formatted sequence of the uniProtID in the uniprot website (www.uniprot.org).
    url="http://www.uniprot.org/uniprot/" + str(uniProtID.rstrip()) + ".fasta"
    
    ##########################
    ###  END CODING HERE  ####
    ##########################
    fh = urllib2.urlopen(url)
    result = fh.read()
    fh.close()
    return result

def main():
    uniProtID_list = open(sys.argv[1])
    db_file = open(sys.argv[2], "w")
    query_folder = sys.argv[3]
    if not query_folder.endswith("/"):
        query_folder = query_folder + "/" 
    # Caution! Please make sure that the directory exists before you run this script.


    for line in uniProtID_list:
        ##########################
        ### START CODING HERE ####
        ##########################
        # Fetch the fasta formatted sequence for each uniProtID.
        # Store the fasta sequences as individual fasta file in your query directory.
        # Store all the fasta sequences in one single fasta file as well. These individual files will be used as (PSI-)BLAST queries later on.

            line = line.rstrip()
            seq = fetchFasta(line)
            fh = open(str(query_folder)+str(line) + ".fasta", "w")
            fh.write(seq)
            fh.close()
            db_file.write(seq)
        ##########################
        ###  END CODING HERE  ####
        ##########################
    uniProtID_list.close()
    db_file.close()
if __name__ == "__main__":
    main()
