#!/usr/bin/python

import os
import sys
import urllib2

def retrieve_go_terms(protein_list, cachefile):
    
    '''
    The following function returns a dictionary, in which each item has a 
    UniProt ID as its key, and the corresponding set of GO terms as its value. 
    The code: s = go_dict["Q9BS40"] would store the set of GO terms associated 
    to protein Q9BS40 in 's'.
    
    If a cachefile is specified as the script's [optional] 5th input argument, 
    the function checks to see if this file already exists. If it does, the GO 
    terms are read from the specified file. If it does not, the GO terms are 
    retrieved through EBI's QuickGO webservice.
    
    Finally, if no cachefile is specified, the GO terms will always be retrieved
    from EBI, and saved to the /tmp/default_go_cache.tmp file.
    '''
    
    # Get the GO data, from EBI's webservice or from a cache (depending on the [optionally] specified cache file).
    if not cachefile or not os.path.isfile(cachefile):
        urlformat = "http://www.ebi.ac.uk/QuickGO/GAnnotation?format=tsv&limit=-1&col=proteinID,goID&protein=%s"
        response  = urllib2.urlopen(urlformat % ','.join(protein_list))
        go_data   = response.read()
        if cachefile:
            f = open(cachefile, "w")
        else:
            f = open("/tmp/default_go_cache.tmp", "w")
        f.write(go_data)
        f.close()
    else:
        with open(cachefile) as f:
            go_data = f.read()
    
    # Define the output dictionary and its keys.
    go_dict = dict()
    for protein_id in protein_list:
        go_dict[protein_id] = set()

    # Loop over the lines in 'go_data' and add every GO annotation to the appropriate set.
    for line in go_data.split('\n'):
        if not line.startswith('ID\t'):
            terms = line.split("\t")
            if len(terms) == 2:
                protein_id    = terms[0]
                go_annotation = terms[1]
                go_dict[protein_id].add(go_annotation)

    return go_dict

def parse(filename):
    
    '''
    The following function returns the list of UniProt IDs contained in the file
    '''
    
    proteins_list = []
    with(open(filename, "rU")) as f:
        for row in f:
            proteins_list.append(row.rstrip())
    return proteins_list

def compute_score(go_set1, go_set2):
    
    '''
    Given two sets of GO terms this function computes the score function.
    '''
    
    intersec = go_set1.intersection(go_set2)
    union = go_set1.union(go_set2)
    if len(union) != 0:
        return float(len(intersec))/len(union)
    return 0.0

def compute_similarity(score, threshold1, threshold2):
    
    '''
    Given a score and two thresholds (threshold1 < threshold2), the following
    function returns the string "different" if the score is less than 
    threshold1; it returns the string "ambiguous" if the score is greater than 
    threshold1, but less than threshold2; it returns the string "similar" 
    if the score is greater than threshold2.
    '''

    if score < threshold1:
        return "different"
    elif threshold1 <= score < threshold2:
        return "ambiguous"
    elif score >= threshold2:
        return "similar"

def compute_pairs(proteins_list):
    
    '''
    The following function returns a list containing all unique protein pairs.
    You can add a pair of proteins to the list using the following code:
    pairs_list.append((protein1, protein2))
    '''
    
    pairs_list = list()
    
    for i in range(len(proteins_list)):
        for j in range(len(proteins_list)):
            
            # Prevent double selection of protein pairs (A-B == B-A) and
            # identical proteins (A-A)
            if i > j:
                pairs_list.append((proteins_list[i], proteins_list[j]))
                
    return pairs_list

# The following function writes the pairs, similarity and, score to the output file.
def write_results(filename, go_dictionary, pairs_list, threshold1, threshold2):
    with open(filename, "w") as f:
        for (p1,p2) in pairs_list:
            score      = compute_score(go_dictionary[p1], go_dictionary[p2])
            similarity = compute_similarity(score, threshold1, threshold2)
            f.write("%s\t%s\t%s\t%s\n" % (p1, p2, similarity, score))

def main():
    
    # Check the number of arguments
    narg = len(sys.argv)
    if narg < 5:
        print "usage: %s inputfile outputfile threshold1 threshold2 [cachefile]" % sys.argv[0]
        return 1
        
    # Parse arguments
    inputfile   = sys.argv[1]
    outputfile  = sys.argv[2]
    threshold1  = float(sys.argv[3])
    threshold2  = float(sys.argv[4])
    
    if narg == 6:
        cachefile = sys.argv[5]
    else:
        cachefile = False

    if threshold1 >= threshold2:
        print "threshold1 must be less than threshold2"
        return 2

    # Parse the input file and retrieve the GO terms associated with each protein.
    proteins_list = parse(inputfile)
    go_dictionary = retrieve_go_terms(proteins_list, cachefile)
    
    # Compute and prints the scores and the similarity string of each protein pair.
    pairs_list = compute_pairs(proteins_list)
    write_results(outputfile, go_dictionary, pairs_list, threshold1, threshold2)
    return 0

if __name__ == "__main__":
    sys.exit(main())