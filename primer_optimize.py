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
    #dna_sequence = "AGAGTTTGATCCTGGCTCAGATTGAACGCTGGCGGCACGCCTAACACATGCAAGTCGAACGGCAGCGGGGGAAAGCTTGCTTTCCTGCCGGCGAGTGGCGGACGGGTGAGTAATGCGTAGGAATTTGCCATTAAGAGGGGGACAACTCGGGGAAACTCGAGCTAATACCACATAATCTCTTCGGAGCAAAGAAGGGGATTCTTCGGAACCTTTCGCTTAATGAGAAGCCTACGTTGGATTAGCTTGTTGGTGGGGTAAAGGCTCACCAAGGCGATGATCTATAGCTGGTCTGAGAGGATGATCAGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTGAGGAATTTTGGACAATGGGGGAAACCCTGATCCAGCGATGCCGCGTGTGTGAAGAAGGCCTAAGGGTTGTAAAGCACTTTTAGTGAGGAAGAGAGTAAGTCGGTTAATACCCGGCTTGCAAGACGTTACTCACAGAAAAAGCGCCGGCTAACTCTGTGCCAGCAGCCGCGGTAATACAGAGGGTGCAAGCGTTAATCGGATTGACTGGGCGTAAAGGGCGCGTAGGCGGTAAGATAAGTCAGATGTTAAAAACCCGAGCTCAACTTGGGGACTGCATTTGAAACTATCTCACTAGAGTACAGTAGAGGAGAGCGGAATTTCCGGTGTAGCGGTGAAATGCGTAGATATCGGAAGGAACACCAGTGGCGAAGGCGGCTCTCTGGACTGACACTGACGCTGAGGCGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCTGTAAACGATGAGAACTAGCTGTTGGTACGTTTAGTATCAGTAGCGCAGCTAACGCGTTAAGTTCTCCGCCTGGGGATTACGGTCGCAAGACTAAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGCGGTTTAATTCGATGCAACCCGAAAAACCTTACCTACCCTTGACATCCCGCGAAGCCTGTAGAGATACGGGCGTGCTCGAAAGAGAACGCGGTGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGTAACGAGCGCAACCCTTGTCCTTAGTTGCCATCTACATTAGTAGGGAACTCTAAGGAGACTGCCGGCGATAAGTCGGAGGAAGGTGGGGACGATGTCAAGTCATCATGGCCTTTATGGGTAGGGCTACACGCGTGCTACAATGGGCAGTACAAAGGGAAGCGAAGCTGTGAAGTGGAGCAAACCTCAGAAAGCTGCTCGTAATCCGGATTGAAGTCTGCAACTCGACTTCATGAGGTTGGAATCGCTAGTAATCGCAGATCAGCATGCTGCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGAAGTGGGTTGTACCAGAAGTAGGAGAGCTAACCTTCGGGAGGCATCTTACCACGGTATGATTCATGACTGGGGTGAAGTCGTAACAAGGTA"
    pp = primer_problem(s=dna_sequence, 
                    US_GC_percent=50,
                    US_Tm_f=59,
                    US_Tm_r=59,
                    Ct=50,
                    US_Tm_diff=0,
                    US_len_diff=0,
                    US_hairpin=0)

    evolver = ProbRVEA_MC(problem=pp, n_gen_per_iter=10, n_iterations=5, keep_archive=True)
    individual, solutions, archive = evolver.end()

    figure = animate_init_(solutions, filename=path_to_results + "/primer_plot_"+str(run)+".html")
    while evolver.continue_evolution():
        print(f"Running iteration {evolver._iteration_counter+1}")
        #print(evolver.population.objectives)
        evolver.iterate()
        #non_dominated = evolver.population.non_dominated_fitness()
        figure = animate_next_(
            evolver.population.objectives, #[non_dominated],
            figure,
            filename=path_to_results + "/primer_plot_"+str(run)+".html",
            generation=evolver._iteration_counter,
        )

    results_dict = {
            'archive': evolver.archive,
            'individuals_solutions': evolver.population.individuals,
            'obj_solutions': evolver.population.objectives,
            'uncertainty_solutions': evolver.population.uncertainity
        }
    path_to_file = path_to_results + '/Run_' + str(run)
    outfile = open(path_to_file, 'wb')
    pickle.dump(results_dict, outfile)
    outfile.close()
    count += 1
    if count >= max_seq:
        break

