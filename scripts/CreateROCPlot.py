#!/usr/bin/python

# This script reads and parses your previously obtained results out of a blast output file, and a benchmark output file.
# It then creates the corresponding ROC plot, and calculates the AUC.

from operator import itemgetter
import pylab
import sys

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
                    
                    # Predictions benchmark are symmetric
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

    sorted_blast = sorted(blast_evalues.items(), key=itemgetter(1))
    
    # For evalue threshold=0, true and false positive rate are both 0
    true_positive_rate  = [0] # TP/(TP+FN), y-axis
    false_positive_rate = [0] # FP/(FP/TN), x-axis
    prev_evalue = -1.
    first = True
        
    for item in sorted_blast:
                
        # We have not counted TP, FP, TN and FN yet, so at first iteration, don't update TPR and FPR
        if prev_evalue != blast_evalues[item[0]] and not first:
            
            try:
                tpr = true_positives / float(true_positives + false_negatives)
            except ZeroDivisionError:
                     tpr = 0
                
            try:
                fpr = false_positives / float(false_positives + true_negatives)
            except ZeroDivisionError:
                     fpr = 0
            
            true_positive_rate.append(tpr)
            false_positive_rate.append(fpr)
              
        first = False

        # Count true/false positives/negatives
        if benchmark_dict[item[0]] == "similar" and float(item[1]) > blast_threshold:
            false_negatives += 1

        elif benchmark_dict[item[0]] == "different" and float(item[1]) > blast_threshold:
            true_negatives += 1

        elif benchmark_dict[item[0]] == "similar" and float(item[1]) < blast_threshold:
            true_positives += 1

        elif benchmark_dict[item[0]] == "different" and float(item[1]) < blast_threshold:
            false_positives += 1

        # If the evalue of a protein pair is the same as the previous one, the ordering doesn't mean anything. Therefore, all these values are taken together.
        prev_evalue = blast_evalues[item[0]]
        
    # Append the very last item to the list, since the last item will be the same as the previous item.
    true_positive_rate.append(true_positives/float(true_positives + false_negatives))
    false_positive_rate.append(false_positives/float(false_positives + true_negatives))
    
    return false_positive_rate, true_positive_rate
    
def integrate(x,y):
    
    '''
    This functions integrates a function using the trapezoid method
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

def main():
    blast_results_file     = sys.argv[1]
    benchmark_results_file = sys.argv[2]
    
    blast_evalues     = parse_blast_results(blast_results_file)
    benchmark_results = parse_benchmark_results(benchmark_results_file)
    x, y              = roc_plots(blast_evalues, benchmark_results)
    auc = integrate(x, y)
    
    print(auc)
    pylab.plot(x,y)
    pylab.show()

if __name__ == "__main__":
    main()