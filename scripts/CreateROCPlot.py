#!/usr/bin/python

# This script reads and parses your previously obtained results out of a blast output file, and a benchmark output file.
# It then creates the corresponding ROC plot, and calculates the AUC.

from operator import itemgetter
import pylab
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
    sorted_blast = sorted(blast_evalues.items(), key=itemgetter(1))

    # For evalue threshold=0, true and false positive rate are both 0
    true_positive_rate  = [0] # TP/(TP+FN), y-axis
    false_positive_rate = [0] # FP/(FP/TN), x-axis
    lowest_pair = dict()

    unique_evalue = sorted(set(blast_evalues.values()))
    print(len(unique_evalue))
    for threshold in unique_evalue:
        false_negatives = 0
        true_negatives = 0
        true_positives = 0
        false_positives = 0
        threshold = float(threshold)

        for prot in sorted_blast:

            #IF GO is similar and BLAST is above threshold = FN
            if benchmark_dict[prot[0]] == "similar" and float(prot[1]) > threshold:
                false_negatives += 1

            #IF GO different and BLAST is below threshold = TN

            elif benchmark_dict[prot[0]] == "different" and float(prot[1]) > threshold:
                true_negatives += 1

            #IF GO is similar and BLAST is above threshold = TP
            elif benchmark_dict[prot[0]] == "similar" and float(prot[1]) <= threshold:
                true_positives += 1

            #If GO is different and BLAST is above threshold = FP
            elif benchmark_dict[prot[0]] == "different" and float(prot[1]) <= threshold:
                false_positives += 1
                lowest_pair[prot[0]] = prot[1]

        true_positive_rate.append(true_positives / float(true_positives + false_negatives))
        false_positive_rate.append(false_positives / float(false_positives + true_negatives))

    #Always write a false_positives.txt sorted list to check whether there are any false false positives.
    # As asked in Q4.3
    sorted_pair = sorted(lowest_pair.items(), key=itemgetter(1))
    fh = open(os.getcwd() + "/false_positives.txt", "w")
    for i in range(len(sorted_pair)):
        fh.write(str(sorted_pair[i][0]) + "\t" + str(sorted_pair[i][1]) + "\n")
    fh.close()

    return false_positive_rate, true_positive_rate


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

    f = {k:v for k, v in zip(x,y)} # Function reprensentation
    auc = 0.0 # Area under the curve
    xs = sorted(x)

    # Sum each trapezoid
    for i in xrange(len(x) - 1):
        a = xs[i]
        b = xs[i+1]
        auc += (b - a) * ((f[a] + f[b]) / 2)

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
    print integrate(x,y)

    #Always write a rocplot table in the current working directory to create a plot using R for example
    fh = open(os.getcwd()+"/blast_rocplot.txt", "w")
    fh.write("x\ty\n")
    for i in range(len(x)):
        fh.write(str(x[i]) + "\t" + str(y[i]) + "\n")
    fh.close()
    # pylab.plot(x,y)
    # pylab.show()

if __name__ == "__main__":
    main()