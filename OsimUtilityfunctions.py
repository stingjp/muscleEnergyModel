import os
os.add_dll_directory('C:/Users/jonstingel/opensim-core-4.5.1-2024-08-23-cf3ef35/bin')
import opensim as osim
import matplotlib.pyplot as plt
import simulationSetups as simset
import pandas as pd
import numpy as np
import scipy.io as sio
import time


# activate the metabolics probe for use in simulations
def probeActivate(model):
    # get the probeset
    probeset = model.getProbeSet()
    numProbes = probeset.getSize()
    # need to loop through and set them all to be enabled hopefully
    for p in range(numProbes):
        probe = probeset.get(p)
        probe.setEnabled(True)
    # update the model to be returned. 
    model.updProbeSet()
    return model

# create an ID helper function that runs the analysis and plots the moments. 
def IDplotter(solution, tag, showornot):
    # get the column labels
    labels = solution.getColumnLabels()
    # get the time
    time = solution.getIndependentColumn()
    # get the number of columns
    numCols = solution.getNumColumns()
    # create a figure plotting all of the moments. 
    fig, axs = plt.subplots(1, 3, figsize=(10, 6))
    for i in range(numCols):
        temp = labels[i]
        if 'knee_angle' in temp and 'beta' not in temp:
            if '_l_' in temp:
                axs[0].plot(time, solution.getDependentColumn(temp).to_numpy(), label=temp, linestyle=':', color='red')
            else:
                axs[0].plot(time, solution.getDependentColumn(temp).to_numpy(), label=temp, linestyle='--', color='red')
        elif 'ankle_angle' in temp:
            if '_l_' in temp:
                axs[0].plot(time, solution.getDependentColumn(temp).to_numpy(), label=temp, linestyle=':', color='blue')
            else:
                axs[0].plot(time, solution.getDependentColumn(temp).to_numpy(), label=temp, linestyle='--', color='blue')
        elif 'pelvis' in temp: 
            axs[1].plot(time, solution.getDependentColumn(temp).to_numpy(), label=temp, linestyle=':')   
        else:
            axs[1].plot(time, solution.getDependentColumn(temp).to_numpy(), label=temp)
    
    axs[0].set_title('Knee/Ankle - ' + tag)
    axs[0].legend()
    axs[1].set_title('Other coordinates - ' + tag)
    # axs[1].legend()
    axs[2].axis('off')
    handles, heads = axs[1].get_legend_handles_labels()
    axs[2].legend(handles, heads, loc='center')
    plt.savefig(tag + '_ID.png')
    if showornot:
        plt.show()
    return

# create a dataframe for all the gait timings
def subjectGaitTimings():
    cycles = {}
    cycles['welk001_welknatural_i'] = [0.341, 1.017, 50.255, 50.907]
    cycles['welk001_welknatural_f'] = [1.017, 1.675, 50.907, 51.578]
    cycles['welk001_welkexo_i'] = [0.511, 1.147, 48.41, 49.04]
    cycles['welk001_welkexo_f'] = [1.147, 1.775, 49.04, 49.673]

    cycles['welk002_welknatural_i'] = [0.933, 1.629, 57.547, 58.249]
    cycles['welk002_welknatural_f'] = [1.629, 2.326, 58.249, 58.947]
    cycles['welk002_welkexo_i'] = [0.707, 1.363, 57.079, 57.745]
    cycles['welk002_welkexo_f'] = [1.363, 2.021, 57.745, 58.398]

    cycles['welk003_welknatural_i'] = [0.782, 1.502, 58.123, 58.841]
    cycles['welk003_welknatural_f'] = [1.502, 2.217, 58.841, 59.551]
    cycles['welk003_welkexo_i'] = [3.006, 3.696, 58.043, 58.718]
    cycles['welk003_welkexo_f'] = [3.696, 4.378, 58.718, 59.388]

    cycles['welk004_welknatural_i'] = [3.572, 4.277, 49.2, 49.893]
    cycles['welk004_welknatural_f'] = [4.277, 4.965, 49.893, 50.597]
    cycles['welk004_welkexo_i'] = [9.281, 9.945, 50.286, 50.931]
    cycles['welk004_welkexo_f'] = [9.945, 10.586, 50.931, 51.582]

    cycles['welk005_welknatural_i'] = [5.595, 6.279, 37.859, 38.536]
    cycles['welk005_welknatural_f'] = [6.279, 6.964, 38.536, 39.227]
    cycles['welk005_welkexo_i'] = [7.205, 7.854, 51.819, 52.475]
    cycles['welk005_welkexo_f'] = [7.854, 8.495, 52.475, 53.119]

    cycles['welk007_welknatural_i'] = [0.713, 1.426, 55.597, 56.31]
    cycles['welk007_welknatural_f'] = [1.426, 2.142, 56.31, 57.016]
    cycles['welk007_welkexo_i'] = [19.294, 19.965, 44.973, 45.651]
    cycles['welk007_welkexo_f'] = [19.965, 20.636, 45.651, 46.324]

    cycles['welk008_welknatural_i'] = [67.197, 67.965, 113.761, 114.53]
    cycles['welk008_welknatural_f'] = [67.965, 68.726, 114.53, 115.299]
    cycles['welk008_welkexo_i'] = [41.961, 42.661, 20.208, 20.915]
    cycles['welk008_welkexo_f'] = [42.661, 43.359, 20.915, 21.619]

    cycles['welk009_welknatural_i'] = [44.199, 44.908, 100.343, 101.058]
    cycles['welk009_welknatural_f'] = [44.908, 45.617, 101.058, 101.77]
    cycles['welk009_welkexo_i'] = [45.963, 46.577, 110.378, 111.007]
    cycles['welk009_welkexo_f'] = [46.577, 47.201, 111.007, 111.622]

    cycles['welk010_welknatural_i'] = [0.316, 1.082, 6.436, 7.198]
    cycles['welk010_welknatural_f'] = [1.082, 1.853, 7.198, 7.963]
    cycles['welk010_welkexo_i'] = [1.222, 1.938, 11.217, 11.915]
    cycles['welk010_welkexo_f'] = [1.938, 2.645, 11.915, 12.624]

    cycles['welk013_welknatural_i'] = [2.071, 2.754, 7.579, 8.261]
    cycles['welk013_welknatural_f'] = [2.754, 3.452, 8.261, 8.945]
    cycles['welk013_welkexo_i'] = [2.599, 3.216, 6.876, 7.487]
    cycles['welk013_welkexo_f'] = [3.216, 3.826, 7.487, 8.098]

    return cycles


    
    
    # edit the experimental data for simulations
    # renameExperimentalData()
    # run the simulations of the subject, and get metabolic cost of motion
    # close all
    # metabolicsModelSetup()
    # close all
    # torqueStateTrackGRFPrescribe()
    # close all
    # torqueStatePrescribeGRFPrescribe()
    # close all
    # torqueMarkerTrackGRFPrescribe()
    # close all
    # muscleStatePrescribeGRFPrescribe_2()
    # close all
    # torqueStateTrackGRFTrack()
    # close all
    # torqueMarkerTrackGRFTrack()
    # close all
    # Issues = muscleStatePrescribeGRFPrescribe(Issues)
    # close all
    return
