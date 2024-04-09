from primer_opt.primer_problem_interactive import primer_problem
import pickle
import matplotlib.pyplot as plt
from plot_parcoord import plot_vals
import numpy as np

path_to_results = 'primer_opt_results/interactive_results_final'
preference_type = 'bounds'
#preference_type = 'refpnt'
n_objs = 8

for gen in range(1,5):
    path_to_file = path_to_results + '/Run_'+preference_type+'_' + str(gen)
    infile = open(path_to_file, 'rb')
    results_data=pickle.load(infile)
    infile.close()


    objs = results_data['obj_solutions']
    uncs = results_data['uncertainty_solutions']
    if gen == 1:
        pref = None
        sorted = objs[:,1].argsort()[::-1]
        objs = objs[sorted][1:]
        uncs = uncs[sorted][1:]
    else:
        pref = results_data['ref_point'].transpose()
    unc_avg_all = np.mean(uncs, axis=1)
    unc_avg_all_max = np.max(unc_avg_all)
    unc_avg_all_min = np.min(unc_avg_all)
    plot_vals(objs=objs,
            unc=uncs, 
            preference=pref,
            iteration = preference_type, 
            interaction_count=gen, 
            ideal=np.min(objs-0.1,axis=0),
            nadir=np.max(objs,axis=0),
            min=unc_avg_all_min,
            max=unc_avg_all_max,
            path=path_to_results)


x=[1,2,3,4,5,6,7,8]
fig, (ax1,ax2,ax3,ax4,ax5,ax6,ax7) = plt.subplots(1, 7, sharey=False)
ax = (ax1,ax2,ax3,ax4,ax5,ax6,ax7)
y= results_data['obj_solutions'][0] 
err = results_data['uncertainty_solutions'][0]
err_y = [err,err]
print(y)
print(err)
# plot subplots and set xlimit
for i in range(7):
    ax[i].plot(x,y,'g')
    ax[i].set_xlim([ x[i],x[i+1]])
for i in range(7):
    ax[i].errorbar(x, y, xerr=0, yerr=err, fmt ='o')

plt.subplots_adjust(wspace=0)
plt.show()