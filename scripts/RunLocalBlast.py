#!/usr/bin/python

# This script runs (PSI-)BLAST locally.

import subprocess
import math
import pylab
import numpy
import optparse

def Blast(query, db, psiblast):
    
    '''
    This function executes blast or psi-blast for the given query and db. 
    The variable psiblast is a binary variable that defines which program to be 
    executed: 0 for BLAST and 1 for PSI-BLAST.
    '''
    
    if (psiblast == 0):
        cmd = "blastall -p blastp -i " + query + " -d" + db + " -m 9"

    else:
        cmd = "blastpgp -j 3 -i " + query + " -d" + db + " -m 9"
    
    p = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, 
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
    blast_result = p.stdout.read()
    return blast_result

def parseBlastResult(blast_result, blast_dict):
    
    '''
    This function parses the output of (PSI-)BLAST and stores the result in 
    blast_dict (defined in main()).
    '''
    
    for line in blast_result.split("\n"):
        if len(line) >1 and line[0] != "#" and line[0] != "[":
            query = line.split()[0].split("|")[1]
            subject = line.split()[1].split("|")[1]
            escore = (float(line.split()[10]))
            blast_dict[(query, subject)] = escore
            
def writeOutput(uniprot_ids, output_filename, blast_dict):
    
    '''
    This function writes the scores of all-against-all protein pairs to the 
    output file.
    '''
    
    with open(output_filename, "w") as f:
        for a in uniprot_ids:
            for b in uniprot_ids:
                if a != b:
                    if (a, b) in blast_dict:
                        f.write(a+"\t"+b+"\t"+ str(blast_dict[(a, b)])+"\n")
                    else:
                        f.write(a+"\t"+b+"\t"+ "NA\n")

def plotEValueDistribution(blast_dict):
    
    '''
    This function plots the distribution of the log(e-value). 
    To avoid log(0), pseudocount is added. 
    In this case pseudocount is the smallest nonzero e-value divided by 1000.
    '''
    
    sorted_e_val = sorted(blast_dict.values())
    nonzero_indices = numpy.nonzero(sorted_e_val)[0]
    pseudocount = sorted_e_val[nonzero_indices[0]]/1000.0
    pylab.hist(map(lambda x: math.log10(x+pseudocount), blast_dict.values()))
    pylab.xlabel("log(e-value)")
    pylab.ylabel("Frequency")
    pylab.show()
  
def main():
    
    # Parse arguments
    parser = optparse.OptionParser()
    parser.add_option("-l", help="the list of UniProt IDs")
    parser.add_option("-q", help="the query folder")
    parser.add_option("-d", help="the fasta file of the database")
    parser.add_option("-o", help="output file")
    parser.add_option("-t", default=0, type=int, help="type of BLAST: 0 = BLASTP (default); 1 = PSI-BLAST")
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
        
    # The blast_dict dictionary will be used to store protein pair and the corresponding e-value.
    # Keys for blast_dict are the combination of query and subject/hit, e.g.:
    # key             = (query, subject)
    # blast_dict[key] = e_value
    blast_dict = dict()
    
    # uniprot_ids is a list to store all UniProt IDs contained in uniProtID_list.
    uniprot_ids = list()


    for line in open(uniProtID_list):   
        qry = line.strip()

        # Run (PSI-)BLAST for all query proteins.
        # Store all the uniprot IDs in the uniprot_ids.
        # Parse and store the blast result in the blast_dict.
        
        uniprot_ids.append(qry)
        qry = query_folder + qry + ".fasta"
        parseBlastResult(Blast(qry, db, psiblast), blast_dict)

    writeOutput(uniprot_ids, output_filename, blast_dict)
    plotEValueDistribution(blast_dict)
    
if __name__ == "__main__":
    main()
