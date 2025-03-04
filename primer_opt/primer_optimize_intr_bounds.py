from primer_opt.primer_problem_interactive import primer_problem
#from desdeo_emo.EAs.RVEA import RVEA
from desdeo_emo.EAs.OfflineRVEA import ProbRVEA_MC
from desdeo_emo.utilities.plotlyanimate import animate_init_, animate_next_
import pandas as pd
import numpy as np
import fastaparser
import pickle

#fasta = open('rep_set.fasta', 'r')
#reader = fastaparser.Reader(fasta, parse_method='quick')


path_to_results = 'primer_opt_results/interactive_results' 
#dna_sequence = "AGAGTTTGATCCTGGCTCAGATTGAACGCTGGCGGCACGCCTAACACATGCAAGTCGAACGGCAGCGGGGGAAAGCTTGCTTTCCTGCCGGCGAGTGGCGGACGGGTGAGTAATGCGTAGGAATTTGCCATTAAGAGGGGGACAACTCGGGGAAACTCGAGCTAATACCACATAATCTCTTCGGAGCAAAGAAGGGGATTCTTCGGAACCTTTCGCTTAATGAGAAGCCTACGTTGGATTAGCTTGTTGGTGGGGTAAAGGCTCACCAAGGCGATGATCTATAGCTGGTCTGAGAGGATGATCAGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTGAGGAATTTTGGACAATGGGGGAAACCCTGATCCAGCGATGCCGCGTGTGTGAAGAAGGCCTAAGGGTTGTAAAGCACTTTTAGTGAGGAAGAGAGTAAGTCGGTTAATACCCGGCTTGCAAGACGTTACTCACAGAAAAAGCGCCGGCTAACTCTGTGCCAGCAGCCGCGGTAATACAGAGGGTGCAAGCGTTAATCGGATTGACTGGGCGTAAAGGGCGCGTAGGCGGTAAGATAAGTCAGATGTTAAAAACCCGAGCTCAACTTGGGGACTGCATTTGAAACTATCTCACTAGAGTACAGTAGAGGAGAGCGGAATTTCCGGTGTAGCGGTGAAATGCGTAGATATCGGAAGGAACACCAGTGGCGAAGGCGGCTCTCTGGACTGACACTGACGCTGAGGCGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCTGTAAACGATGAGAACTAGCTGTTGGTACGTTTAGTATCAGTAGCGCAGCTAACGCGTTAAGTTCTCCGCCTGGGGATTACGGTCGCAAGACTAAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGCGGTTTAATTCGATGCAACCCGAAAAACCTTACCTACCCTTGACATCCCGCGAAGCCTGTAGAGATACGGGCGTGCTCGAAAGAGAACGCGGTGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGTAACGAGCGCAACCCTTGTCCTTAGTTGCCATCTACATTAGTAGGGAACTCTAAGGAGACTGCCGGCGATAAGTCGGAGGAAGGTGGGGACGATGTCAAGTCATCATGGCCTTTATGGGTAGGGCTACACGCGTGCTACAATGGGCAGTACAAAGGGAAGCGAAGCTGTGAAGTGGAGCAAACCTCAGAAAGCTGCTCGTAATCCGGATTGAAGTCTGCAACTCGACTTCATGAGGTTGGAATCGCTAGTAATCGCAGATCAGCATGCTGCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGAAGTGGGTTGTACCAGAAGTAGGAGAGCTAACCTTCGGGAGGCATCTTACCACGGTATGATTCATGACTGGGGTGAAGTCGTAACAAGGTA"
dna_sequence = "GAGTAACGCGTAGGAACCAACCTTAGAGAGTGGAATAACCTTGGGAAACTAAGGCTAATACCGCATATACCTCGAGAGGGAAAGGAGAGTAATCTCTGCTCTAGGACGGGCCTGCGCCCGATTAGCTTGTTGGTAAGGTAATGGCTTACCAAGGCATCGATCGGTAGCTGGTCTGAGAGGACGATCAGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGGACAATGGGCGCAAGCTTGATCCAGCCATGCCGCGTGAGTGAAGAAGGCCTTCGGGTTGTAAAGCTCTTTCACACGCGACGATGATGACGGTAGCGTGAGAAGAAGCCCCGGCTAACTTCGTGCCAGCAGCCGCGGTAATACGAAGGGGGCTAGCGTTGTTCGGAATTACTGGGCGTAAAGGGCGCGTAGGCGGTCTAATTTGTCGGGGGTGAAATCCCAGGGCTTAACCTTGGAAGTGCCTTCGGGACAATTAGGCTTGAGACCGGGAGAGGATGGCGGAATTCCCAGTGTAGAGGTGAAATTCGTAGATATTGGGAAGAACACCGGTGGCGAAAGCGGCCATCTGGTCCGGTTCTGACGCTAAAGCGCGAAAGCGTGGGGGAGCGAACAGGATTAGATACCCTGGTAGCCACGCCGTAAACGATGTGTGCTGGATGTCGGGGGGCATGCTCTTCGGTGTCGTAGCTAACGCGTGAAGCACACCGTCTGGGGAGTACGGCCGCAAGGTTAAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGCAGAACCTTACCAGCCTTTGACATGCCCTTTATATCCTAAAGAGACTTGGGAGTCGGTTCGGCCGGAAGGGACACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCCTATTCTCAGTTGCCATCGGGTCATGCCGGGCACTCTGAGGGGACTGCCGGTGACAAGCCGGAGGAAGGTGGGGATGACGTCAAGTCCTCATGGCCCTTACAGGCTGGGCTACACACGTGCTACAATGGCGGTGACAATGGGTTATCAGGCGACTCTGCGAAGAGGAGCGAATCCTAAAAGACCGTCTTAGTTCGGATTGCACTCTGCAACCCGGGTGCATGAAGTTGGAATCGCTAGTAATCGCGGATCAGCACGCCGCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGTTGTCTTTACTCGAAGACAGTGTGCCAACCTTAA"

n_objs = 8
pp = primer_problem(s=dna_sequence,Ct=50)

preference_input = np.zeros((8,2))
evolver = ProbRVEA_MC(problem=pp, n_gen_per_iter=100, n_iterations=10, interact=True, lattice_resolution=4)
evolver.set_interaction_type('Preferred ranges')
individual, solutions, archive = evolver.end()

for i in range(1):
    evolver.iterate()

figure = animate_init_(solutions, filename=path_to_results + "/primer_plot_bounds.html")
while evolver.continue_evolution():
    results_dict = {
    'archive': evolver.archive,
    'individuals_solutions': evolver.population.individuals,
    'obj_solutions': evolver.population.objectives,
    'uncertainty_solutions': evolver.population.uncertainity,
    'ref_point': preference_input
    }
    path_to_file = path_to_results + '/Run_bounds_' + str(evolver._iteration_counter)
    outfile = open(path_to_file, 'wb')
    pickle.dump(results_dict, outfile)
    outfile.close()
    print(f"Running iteration {evolver._iteration_counter+1}")
    pref, plot = evolver.start()
    print(pref.content['message'])
    #for i in range(7):
    #    preference_input[i] = input("Reference point obj "+str(i+1)+" : ")
    for i in range(n_objs):
        preference_input[i,0] = input("Lower bound obj "+str(i+1)+" : ")
        preference_input[i,1] = input("Upper bound obj "+str(i+1)+" : ")
    
    response = preference_input
    #pref.response = pd.DataFrame(response, columns=pref.content['dimensions_data'].columns)
    pref.response = response
    evolver.iterate(pref)

    #non_dominated = evolver.population.non_dominated_fitness()
    figure = animate_next_(
        evolver.population.objectives, #[non_dominated],
        figure,
        filename=path_to_results + "/primer_plot_bounds.html",
        generation=evolver._iteration_counter,
    )
individual, solutions, archive = evolver.end()