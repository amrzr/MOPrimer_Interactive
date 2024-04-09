from plotly.offline.offline import plot
import plotly_express as ex
import pandas as pd
import numpy as np
import matplotlib

def plot_vals(objs, unc, preference, iteration, interaction_count, ideal, nadir, min, max, path):
    objs_orig = objs
    columns = ["f_"+str(i+1) for i in range(np.shape(objs)[1])]
    range_plot = np.vstack((ideal,nadir))
    range_plot = np.hstack((range_plot,[[3],[3]]))
    if np.shape(objs)[0] > 0:
        unc_avg = np.mean(unc, axis=1)

        #unc_avg = (unc_avg-np.min(unc_avg))/(np.max(unc_avg)-np.min(unc_avg))
        unc_avg = (unc_avg - min) / (max-min)
        objs_col = unc_avg.reshape(-1, 1)
        objs = np.hstack((objs, objs_col))
    objs = np.vstack((objs, range_plot))
    objs = pd.DataFrame(objs, columns=columns + ["color"])

    if preference is not None:
        if np.shape(preference)[0] > 2:
            col_prefs = [[2]]
        else:
            col_prefs = np.ones((np.shape(preference)[0],1))*2
        if np.shape(preference)[0] > 2:
            preference = preference.reshape(1,-1)
        pref = pd.DataFrame(np.hstack((preference, col_prefs)), columns=columns + ["color"])
        data_final = pd.concat([objs, pref])
    else:
        data_final = objs

    color_scale_custom = [(0.0, 'rgb(69,2,86)'), (0.083, 'rgb(59,28,140)'), (0.167, 'rgb(33,144,141)'),
                          (0.25, 'rgb(90,200,101)'), (0.334, 'rgb(249,231,33)'),
                          (0.334, 'red'), (0.7, 'red'), (0.7, 'white'),
                          (1.0, 'white')]


    fig = ex.parallel_coordinates(
            data_final,
            dimensions=columns,
            color="color", color_continuous_scale=color_scale_custom, range_color=(0,3),width=600, height=400)
    #plot(fig, filename= path + "/solutions_" + str(iteration) + "_" + str(interaction_count) + ".html")
    fig.update_layout()
    fig.write_image(path + "/solutions_" + str(iteration) + "_" + str(interaction_count) + "2.pdf")
    """
    color_scale_custom2 = [(0.0, (69,2,86,1)), (0.083, (59,28,140,1)), 
                        (0.167, (33,144,141,1)),(0.25, (90,200,101,1)), (0.334, (249,231,33,1)),
                        (0.334, (237, 9, 9,1)), (0.7, (237, 9, 9,1)), (0.7, (247, 247, 242,1)),(1.0, (247, 247, 242,1))]
    cm = matplotlib.colors.LinearSegmentedColormap.from_list("",color_scale_custom2)
    data_final['color']=data_final['color']/3
    data_final.plot(legend=False)
    pd.plotting.parallel_coordinates(data_final, "color", colormap=cm).get_figure().savefig(path + "/solutions_" + str(iteration) + "_" + str(interaction_count) + ".pdf")
    """