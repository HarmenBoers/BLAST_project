#!/usr/bin/python

# This script reads and parses your previously obtained results out of a blast output file, and a benchmark output file.
# It then creates the corresponding ROC plot, and calculates the AUC.

from operator import itemgetter
#import pylab
import sys
import os


def parse_blast_results(filename):
    blast_evalues = dict()
    f = open(filename)
    for line in f:
        line = line.rstrip()
        arr = line.split("\t")
        if len(arr) == 3:
            if (arr[0] != arr[1]):
                key = arr[0] + "_" + arr[1]
                if(arr[2] == "NA"):
                    value = 1e6
                else:
                    value = float(arr[2])
                if (key not in blast_evalues):
                    blast_evalues[key] = value
            else:
                print "Warning: Comparing protein to itself:", arr[0]
        else:
            print "Warning: the following line does not have three elements separated by a tab:", line
    f.close()
    return(blast_evalues)


def parse_benchmark_results(filename):
    benchmark_dict = dict()
    f = open(filename)
    for line in f:
        line = line.rstrip()
        arr = line.split("\t")
        if len(arr) == 4:
            if (arr[0] != arr[1]):
                key = arr[0] + "_" + arr[1]
                if (key not in benchmark_dict):
                    benchmark_dict[key] = arr[2]
                    #predictions benchmark are symmetric, so add both possible keys:
                    key = arr[1] + "_" + arr[0]
                    if (key not in benchmark_dict):
                      benchmark_dict[key] = arr[2]
            else:
                print "Warning: Comparing a protein to itself:", arr[0]
        else:
            print "Warning: the following line does not have three elements separated by a tab:", line
    f.close()
    return(benchmark_dict)


def roc_plots(blast_evalues, benchmark_dict):
    
    # Initialize returned variables
    coverage  = [0]
    error = [0]
    #lowest_pair = dict()
        
    # Reverse blast_evalues dictionnary to allow iteration on unique evalues
    inv_evalues = {}
    for k, v in blast_evalues.iteritems():
        v = float(v)
        inv_evalues[v] = inv_evalues.get(v, [])
        inv_evalues[v].append(k)
    sorted_evalues = sorted(inv_evalues.keys())
        
    # First E-value threshold is lowest observed E-value
    # Every pair is considered as being negative according to BLAST
    # True negatives are the number of similar pairs according to benchmark
    # False negatives are the number of different pairs according to benchmark
    tn = benchmark_dict.values().count('different') # True negative
    fn = benchmark_dict.values().count('similar') # False negative
    tp = 0 # True positive
    fp = 0 # False positive

    # Iterate through every unique E-value observed
    for evalue in sorted_evalues:
        
        proteins = inv_evalues[evalue]            
        benchmark = [benchmark_dict[p] for p in proteins]
        similar = benchmark.count('similar')
        different = benchmark.count('different')

        # Update true positive/negatives and false positive/negatives
        tp += similar
        fn -= similar
        tn -= different
        fp += different
        
        # Compute coverage and error
        coverage.append(float(tp) / (tp + fn))
        error.append(float(fp) / (fp + tn))

    return error, coverage
    
def integrate(x,y):
    #########################
    ### START CODING HERE ###
    #########################
    '''
    Description
        This functions integrates a function using the trapezoid method
    Arguments
        x: list. x-coordinates
        y: list. y-coordinates
    Value
        auc: float. Area under the curve
    '''
    # Check arguments
    if len(x) != len(y):
        print 'Error: Integrate arguments are invalid'
        return None

    auc = 0.0 # Area under the curve
    for i in xrange(0, len(x)-1): # Sum each trapezoid
        auc += (x[i+1] - x[i]) * ((y[i] + y[i+1]) / 2)
    return auc
    #########################
    ###  END CODING HERE  ###
    #########################

def main():
    blast_results_file     = sys.argv[1]
    benchmark_results_file = sys.argv[2]

    blast_evalues     = parse_blast_results(blast_results_file)
    benchmark_results = parse_benchmark_results(benchmark_results_file)
    x, y              = roc_plots(blast_evalues, benchmark_results)
    auc = integrate(x,y)
    print auc
    return auc

    #Always write a rocplot table in the current working directory to create a plot using R for example
    # fh = open(os.getcwd()+"/blast_rocplot.txt", "w")
    # fh.write("x\ty\n")
    # for i in range(len(x)):
    #     fh.write(str(x[i]) + "\t" + str(y[i]) + "\n")
    # fh.close()
    # pylab.plot(x,y)
    # pylab.show()

if __name__ == "__main__":
    main()