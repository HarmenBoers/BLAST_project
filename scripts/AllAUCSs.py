#/usr/bin/python

import argparse 
import subprocess

# Define arguments
parser = argparse.ArgumentParser(description='Computes ROC-plot AUC for all threshold combinations')
parser.add_argument('blast', help='BLAST  file')
parser.add_argument('uniprot', help='UniprotID list file')
parser.add_argument('-g', '--classify_go', help='ClassifyGO.py script')
parser.add_argument('-r', '--create_roc_plot', help='CreateROCPlot.py script')
parser.add_argument('-o', '--output', help='Output file')
parser.add_argument('-c', '--cache', help='Cache file')

args = parser.parse_args()

# Parse arguments
blast_file = args.blast
uniprot_file = args.uniprot

if not args.output:
    output_file = './results/all_aucs.out'
else:
    output_file = args.output
    
if not args.create_roc_plot:
    classify_go = './scripts/ClassifyGO.py'
else:
    classify_go = args.classify_go
    
if not args.create_roc_plot:
    create_roc_plot = './scripts/CreateROCPlot.py'
else:
    create_roc_plot = args.create_roc_plot
    
if not args.cache:
    cache = "/tmp/default_go_cache.tmp"
else:
    cache = args.cache

# Define commands
go_file = './results/go_results.out'
open(go_file, 'w').close()
gen_go_cmd = ['python', classify_go, uniprot_file, go_file]
roc_cmd = ['python', create_roc_plot, blast_file, go_file]

# Main
out = open(output_file, 'w')
auc_list = []
for i in xrange(1,101,5):
    for j in xrange(1,101,5):
        if i >= j: continue
        t1 = str(float(i)/100)
        t2 = str(float(j)/100)
        go_cmd = gen_go_cmd + [t1, t2, cache]
        go = subprocess.Popen(go_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        go.wait()
        roc = subprocess.Popen(roc_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        auc = roc.communicate()[0].strip()
        print(t1, t2, auc)
        out.write(','.join([t1, t2, str(auc)]) + '\n')
out.close()