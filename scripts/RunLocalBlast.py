#!/usr/bin/python

# This script runs (PSI-)BLAST locally.

import subprocess
import math
import pylab
import numpy

# For Python versions prior to 3.2, use optparse:
import optparse

# For Python version 3.2 and later, use argparse:
# import argparse


# This function executes blast or psi-blast  for the given query and db.
# The variable psiblast is a binary variable that defines which program to be executed: 0 for BLAST and 1 for PSI-BLAST.
def Blast(query, db, psiblast):
    if (psiblast == 0):
        ##########################
        ### START CODING HERE ####
        ##########################
        # Define the variable 'cmd' as a string with the command for BLASTing 'query' against the specified database 'db'.
        # Note that it is is easier to parse the output if the output is in tabular format (add option -m 9 or -m 8).
        cmd = "blastall -p blastp -i " + query + " -d" + db + " -m 9 -e 10000"
        
        ##########################
        ###  END CODING HERE  ####
        ##########################
    else:
        ##########################
        ### START CODING HERE ####
        ##########################
        # Define the variable 'cmd' as a string with the command for PSI-BLASTing 'query' against the specified database 'db'.
        # Note that it is is easier to parse the output if the output is in tabular format (add option -m 9 or -m 8).
        cmd = "blastpgp -j 3 -i " + query + " -d" + db + " -m 9 -e 10000"
        
        ##########################
        ###  END CODING HERE  ####
        ##########################
    
    p = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
    blast_result = p.stdout.read()
    return blast_result

# This function parses the output of (PSI-)BLAST and stores the result in blast_dict (defined in main()).
# Note that it is is easier to parse the output if the output is in tabular format (add option -m 9 or -m 8).
def parseBlastResult(blast_result, blast_dict):
    for line in blast_result.split("\n"):
        if len(line) >1 and line[0] != "#" and line[0] != "[":
            query = line.split()[0].split("|")[1]
            subject = line.split()[1].split("|")[1]
            ##########################
            ### START CODING HERE ####
            ##########################
            # Parse the e-score corresponding to this line's (query, subject) pair and store it in blast_dict.
            escore = (float(line.split()[10]))
            blast_dict[(query, subject)] = escore

            ##########################
            ###  END CODING HERE  ####
            ##########################
            
# This function writes the scores of all-against-all protein pairs to the output file.
def writeOutput(uniprot_ids, output_filename, blast_dict):
    with open(output_filename, "w") as f:
        for a in uniprot_ids:
            for b in uniprot_ids:
                if a != b:
                    if (a, b) in blast_dict:
                        f.write(a+"\t"+b+"\t"+ str(blast_dict[(a, b)])+"\n")
                    else:
                        f.write(a+"\t"+b+"\t"+ "NA\n")

# This function plots the distribution of the log(e-value). pseudocount is added to avoid log(0).
# The pseudocount in this case is the smallest nonzero e-value divided by 1000.
def plotEValueDistribution(blast_dict):
    print(sorted(blast_dict.values()))
    sorted_e_val = sorted(blast_dict.values())
    nonzero_indices = numpy.nonzero(sorted_e_val)[0]
    pseudocount = sorted_e_val[nonzero_indices[0]]/1000.0
    pylab.hist(map(lambda x: math.log10(x+pseudocount), blast_dict.values()))
    pylab.xlabel("log(e-value)")
    pylab.ylabel("Frequency")
    pylab.show()
    ##########################
    ### START CODING HERE ####
    ##########################
    # Calculate the number of e-values lower than threshold.
    
    
    ##########################
    ###  END CODING HERE  ####
    ##########################
  
def main():
# Since our script needs many arguments from the user, it is advisable to create a nice command-line interface for the user. 
# Depending on your python version, this should be done with the optparse or argparse module as below.

    # For Python versions prior to 3.2:
    parser = optparse.OptionParser()
    parser.add_option("-l", help="the list of UniProt IDs")
    parser.add_option("-q", help="the query folder")
    parser.add_option("-d", help="the fasta file of the database")
    parser.add_option("-o", help="output file")
    parser.add_option("-t", default=0, type=int, help="Type of BLAST to be run: 0 = BLASTP (default); 1 = PSI-BLAST")
    (options, args) = parser.parse_args()

    # Assign the parsed arguments to the corresponding variables.
    if options.l == None:
        print "error: option -l is required"
        exit(1)
    if options.q == None:
        print "error: option -q is required"
        exit(1)
    if options.d == None:
        print "error: option -d is required"
        exit(1)
    if options.o == None:
        print "error: option -o is required"
        exit(1)
    uniProtID_list   = options.l
    query_folder     = options.q
    if not query_folder.endswith("/"):
        query_folder = query_folder + "/" 
    db               = options.d
    psiblast         = options.t
    output_filename  = options.o
    
    # For Python version 3.2 and later:
    """
    parser = argparse.ArgumentParser(description="Automatically running BLAST and PSI-BLAST")
    parser.add_argument("-l", help="the list of UniProt IDs", required=True)
    parser.add_argument("-q", help="the query folder", required=True)
    parser.add_argument("-d", help="the fasta file of the database", required=True)
    parser.add_argument("-o", help="output file", required=True)
    parser.add_argument("-t", choices=[0,1], default=0, type=int, help="Type of BLAST to be run: 0 = BLASTP (default); 1 = PSI-BLAST")
    args = parser.parse_args()
    
    # Assign the parsed arguments to the corresponding variables.
    uniProtID_list = vars(args)["l"]
    query_folder = vars(args)["q"]
    db = vars(args)["d"]
    psiblast = vars(args)["t"]
    output_filename = vars(args)["o"]
    """ 
    
    # The blast_dict dictionary will be used to store protein pair and the corresponding e-value.
    # Keys for blast_dict are the combination of query and subject/hit, e.g.:
    # key             = (query, subject)
    # blast_dict[key] = e_value
    blast_dict = dict()
    # uniprot_ids is a list to store all UniProt IDs contained in uniProtID_list.
    uniprot_ids = list()


    for line in open(uniProtID_list):   
        qry = line.strip()
        ##########################
        ### START CODING HERE ####
        ##########################
        # Run (PSI-)BLAST for all query proteins.
        # Store all the uniprot IDs in the uniprot_ids.
        # Parse and store the blast result in the blast_dict.
        #print(qry)
        uniprot_ids.append(qry)
        qry = query_folder+qry+".fasta"

        parseBlastResult(Blast(qry, db, psiblast), blast_dict)
        ##########################
        ###  END CODING HERE  ####
        ##########################

    writeOutput(uniprot_ids, output_filename, blast_dict)
    plotEValueDistribution(blast_dict)
    
if __name__ == "__main__":
    main()
