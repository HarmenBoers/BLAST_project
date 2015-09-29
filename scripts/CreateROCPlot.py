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
    sorted_blast = sorted(blast_evalues.items(), key=itemgetter(1))

    # For evalue threshold=0, true and false positive rate are both 0
    true_positive_rate  = [0] # TP/(TP+FN), y-axis
    false_positive_rate = [0] # FP/(FP/TN), x-axis
    lowest_pair = dict()
    prev_evalue = -1
    first = True
    unique_evalue = sorted(set(blast_evalues.values()))
    print(len(unique_evalue))
    for threshold in unique_evalue:
        false_negatives = 0
        true_negatives = 0
        true_positives = 0
        false_positives = 0
        threshold = float(threshold)
        # if prev_evalue != first:
        #     true_positive_rate.append(true_positives / float(true_positives + false_negatives))
        #     false_positive_rate.append(false_positives / float(false_positives + true_negatives))
        # first = False
        for prot in sorted_blast:
            #print(prot)


            if benchmark_dict[prot[0]] == "similar" and float(prot[1]) > threshold:
                false_negatives += 1

            #IF GO different/ambiguous and BLAST is below threshold

            elif benchmark_dict[prot[0]] == "different" and float(prot[1]) > threshold:
                true_negatives += 1

            #IF GO is similar and BLAST is above threshold
            elif benchmark_dict[prot[0]] == "similar" and float(prot[1]) <= threshold:
                true_positives += 1

            #If GO is diffenrt and BLAST is above threshold
            elif benchmark_dict[prot[0]] == "different" and float(prot[1]) <= threshold:
                false_positives += 1
        true_positive_rate.append(true_positives / float(true_positives + false_negatives))
        false_positive_rate.append(false_positives / float(false_positives + true_negatives))
        # prev_evalue = blast_evalues[threshold[0]]






    # for thresh in sorted_blast:
    #     # We have not counted TP, FP, TN and FN yet, so at first iteration, don't update TPR and FPR
    #     if prev_evalue != blast_evalues[thresh[0]] and not first:
    #         true_positive_rate.append(true_positives / float(true_positives + false_negatives))
    #         false_positive_rate.append(false_positives / float(false_positives + true_negatives))
    #     first = False
    #
    #     #########################
    #     ### START CODING HERE ###
    #     #########################
    #     #print(benchmark_dict[item[0]])
    #     #print(item[1])
    #     #If GO is similar but BLAST is below threshold
    #     for prot in sorted_blast:
    #         if benchmark_dict[prot[0]] == "similar" and float(prot[1]) > thresh:
    #             false_negatives += 1
    #
    #         #IF GO different/ambiguous and BLAST is below threshold
    #         if benchmark_dict[prot[0]] == "different" and float(prot[1]) > thresh:
    #             true_negatives += 1
    #
    #         #IF GO is similar and BLAST is above threshold
    #         if benchmark_dict[prot[0]] == "similar" and float(prot[1]) < thresh:
    #             true_positives += 1
    #
    #         #If GO is diffenrt and BLAST is above threshold
    #         if benchmark_dict[prot[0]] == "different" and float(prot[1]) < thresh:
    #             false_positives += 1
    #             #lowest_pair[item[0]] = item[1]
    #     #########################
    #     ###  END CODING HERE  ###
    #     #########################
    #     # If the evalue of a protein pair is the same as the previous one, the ordering doesn't mean anything.
    #     # Therefore, all these values are taken together.
    #
    #     #This is the exact same as the value of item[1] gg
    #     prev_evalue = blast_evalues[thresh[0]]
    # Append the very last item to the list, since the last item will be the same as the previous item.

    # true_positive_rate.append(true_positives/float(true_positives + false_negatives))
    # false_positive_rate.append(false_positives/float(false_positives + true_negatives))


    sorted_pair = sorted(lowest_pair.items(), key=itemgetter(1))
    fh = open("/Users/harmen/PycharmProjects/blast_project/false_positives.txt", "w")
    for i in range(len(sorted_pair)):
        fh.write(str(sorted_pair[i][0]) + "\t" + str(sorted_pair[i][1]) + "\n")
    #print(sorted_pair[0:])
    fh.close()
    return false_positive_rate, true_positive_rate

def integrate(x,y):
    #########################
    ### START CODING HERE ###
    #########################
    # auc = 0
    # for i in xrange(len(x)-1):
    #     auc += (x[i+1] - x[i]) * ((y[i] + y[i + 1]) / 2)
    auc = metrics.auc(x, y, reorder=True)
    print auc
    return auc
    #########################
    ###  END CODING HERE  ###
    #########################


def drange(start, stop, step):
    r = start
    while r < stop:
        yield r
        r += step

def main():
    blast_results_file     = sys.argv[1]
    benchmark_results_file = sys.argv[2]

    blast_evalues     = parse_blast_results(blast_results_file)
    benchmark_results = parse_benchmark_results(benchmark_results_file)
    x, y              = roc_plots(blast_evalues, benchmark_results)

    print integrate(x,y)
    #print(x,y)
    fh = open("/Users/harmen/PycharmProjects/blast_project/blast_rocplot.txt", "w")
    fh.write("x\ty\n")
    for i in range(len(x)):
        fh.write(str(x[i]) + "\t" + str(y[i]) + "\n")
    fh.close()
    #
    pylab.plot(x,y)
    pylab.show()

if __name__ == "__main__":
    main()