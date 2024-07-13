from primer_opt.primer_problem import primer_problem
#from desdeo_emo.EAs.RVEA import RVEA
from desdeo_emo.EAs.OfflineRVEA import ProbRVEA_MC
from desdeo_emo.utilities.plotlyanimate import animate_init_, animate_next_
import pandas as pd
import pickle
import fastaparser


path_to_data = 'primer_data'
path_to_results = 'primer_opt_results' 

max_seq = 3
fasta = open(path_to_data+'/rep_set.fasta', 'r')
reader = fastaparser.Reader(fasta, parse_method='quick')
count = 0
for seq in reader:
    dna_sequence = seq.sequence
    run = seq.header[1:]
    path_to_file = path_to_results + '/Run_' + str(run)
    infile = open(path_to_file, 'rb')
    results_data=pickle.load(infile)
    infile.close()
    