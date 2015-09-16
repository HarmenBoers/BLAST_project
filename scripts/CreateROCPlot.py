#!/usr/bin/python

# This script reads and parses your previously obtained results out of a blast output file, and a benchmark output file.
# It then creates the corresponding ROC plot, and calculates the AUC.

from operator import itemgetter
import pylab
import sys
from sklearn import metrics

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
    false_negatives = 0
    true_negatives = 0
    true_positives = 0
    false_positives = 0
    blast_threshold = 1e-5
    for key, value in benchmark_dict.iteritems():        
        #########################
        ### START CODING HERE ###
        #########################
        # Ignore entries in the benchmark_dict classified as "ambiguous"
        #print(key,value)
        #THIS IS NOT VERY USEFUL... right?!
        if value == "ambiguous":
            continue
            #del benchmark_dict[key]
        #########################
        ###  END CODING HERE  ###
        #########################
    sorted_blast = sorted(blast_evalues.items(), key=itemgetter(1))
    # For evalue threshold=0, true and false positive rate are both 0
    true_positive_rate  = [0] # TP/(TP+FN), y-axis
    false_positive_rate = [0] # FP/(FP/TN), x-axis
    prev_evalue = -1.
    first = True
    for item in sorted_blast:
        # We have not counted TP, FP, TN and FN yet, so at first iteration, don't update TPR and FPR
        if prev_evalue != blast_evalues[item[0]] and not first:
            true_positive_rate.append(true_positives / float(true_positives + false_negatives))
            false_positive_rate.append(false_positives / float(false_positives + true_negatives))            
        first = False
        #########################
        ### START CODING HERE ###
        #########################
        #print(benchmark_dict[item[0]])
        #print(item[1])
        #TP FN FP TN

        #If GO is similar but BLAST is below threshold
        if benchmark_dict[item[0]] == "similar" and float(item[1]) > blast_threshold:
            false_negatives += 1

        #IF GO different/ambiguous and BLAST is below threshold
        if benchmark_dict[item[0]] == "different" and float(item[1]) > blast_threshold:
            true_negatives += 1

        #IF GO is similar and BLAST is above threshold
        if benchmark_dict[item[0]] == "similar" and float(item[1]) < blast_threshold:
            true_positives += 1

        #If GO is diffenrt and BLAST is above threshold
        if benchmark_dict[item[0]] == "different" and float(item[1]) < blast_threshold:
            false_positives += 1

        #########################
        ###  END CODING HERE  ###
        #########################
        # If the evalue of a protein pair is the same as the previous one, the ordering doesn't mean anything.
        # Therefore, all these values are taken together.

        #This is the exact same as the value of item[1] gg
        prev_evalue = blast_evalues[item[0]]
    # Append the very last item to the list, since the last item will be the same as the previous item.

    true_positive_rate.append(true_positives/float(true_positives + false_negatives))
    false_positive_rate.append(false_positives/float(false_positives + true_negatives))
    return false_positive_rate, true_positive_rate

def integrate(x,y):
    #########################
    ### START CODING HERE ###
    #########################
    # auc = 0
    # for i in xrange(len(x)-1):
    #     auc += (x[i+1] - x[i]) * ((y[i] + y[i + 1]) / 2)
    # #print auc
    auc = metrics.auc(x, y, reorder=True)
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
    #print(x,y)
    pylab.plot(x,y)
    pylab.show()

if __name__ == "__main__":
    main()