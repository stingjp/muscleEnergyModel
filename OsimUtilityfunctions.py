import os
os.add_dll_directory('C:/Users/jonstingel/opensim-core-4.5.1-2024-08-23-cf3ef35/bin')
import opensim as osim
import matplotlib.pyplot as plt
import simulationSetups as simset
import pandas as pd
import numpy as np
import scipy.io as sio
import time
import pdb


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
def IDplotter(solution, tag, showornot, trialinfo):
    subjectname = trialinfo[0]
    conditionname = trialinfo[1]
    trialname = trialinfo[2]
    IDmoments = {}
    trialdir = os.getcwd()
    model = osim.Model(os.path.join(trialdir, 'simple_model_all_the_probes.osim'))
    modelmass = model.getTotalMass(model.initSystem())
    # first gather and set the ID moments
    IDmoments = getIDMoments(trialdir, IDmoments, modelmass)

    # get the simulation column labels
    labels = solution.getColumnLabels()
    # get the time
    time = solution.getIndependentColumn()
    # get the number of columns
    numCols = solution.getNumColumns()
    # create a figure plotting all of the moments. 
    fig, axs = plt.subplots(1, 4, figsize=(12, 5))
    for i in range(numCols):
        temp = labels[i]
        if temp not in IDmoments and 'force' not in temp:
            print('moment not found in IDmoments')
            print(temp)
            raise Exception('Expected moment was missing')
            break
        if 'knee_angle' in temp and 'beta' not in temp:
            if '_l_' in temp:
                axs[0].plot(time, solution.getDependentColumn(temp).to_numpy(), label=temp, linestyle=':', color='red')
                axs[0].plot(np.linspace(time[0], time[-1], 100), IDmoments[temp].flatten()*modelmass, label='ID_'+temp, linestyle='-', color='red')
            else:
                axs[0].plot(time, solution.getDependentColumn(temp).to_numpy(), label=temp, linestyle='--', color='red')
                axs[0].plot(np.linspace(time[0], time[-1], 100), IDmoments[temp].flatten()*modelmass, label='ID_'+temp, linestyle='-', color='red')
        elif 'ankle_angle' in temp:
            if '_l_' in temp:
                axs[0].plot(time, solution.getDependentColumn(temp).to_numpy(), label=temp, linestyle=':', color='blue')
                axs[0].plot(np.linspace(time[0], time[-1], 100), IDmoments[temp].flatten()*modelmass, label='ID_'+temp, linestyle='-', color='blue')
            else:
                axs[0].plot(time, solution.getDependentColumn(temp).to_numpy(), label=temp, linestyle='--', color='blue')
                axs[0].plot(np.linspace(time[0], time[-1], 100), IDmoments[temp].flatten()*modelmass, label='ID_'+temp, linestyle='-', color='blue')
        elif 'pelvis' in temp: 
            axs[2].plot(time, solution.getDependentColumn(temp).to_numpy(), label=temp, linestyle=':')   
        else:
            axs[2].plot(time, solution.getDependentColumn(temp).to_numpy(), label=temp)
    
    axs[0].set_title('Knee/Ankle - ' + tag)
    # axs[0].legend()
    handles, heads = axs[0].get_legend_handles_labels()
    axs[1].axis('off')
    axs[1].legend(handles, heads, loc='center')
    axs[1].set_title('Other coordinates - ' + tag)
    # axs[1].legend()
    axs[3].axis('off')
    handles2, heads2 = axs[2].get_legend_handles_labels()
    axs[3].legend(handles2, heads2, loc='center')
    
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

# function for getting the muscle activations. 
def getMuscleActivations(trialdir, muscleacts):
    # first need to establish the file that we are going to load
    trialactfile = osim.TimeSeriesTable(os.path.join(trialdir, 'muscletrack_redo_states_py.sto'))
    trialtime = trialactfile.getIndependentColumn()
    time100 = np.linspace(trialtime[0], trialtime[-1], 100).reshape(100,1)
    # get the column labels
    labels = trialactfile.getColumnLabels()
    # get the number of columns
    numCols = trialactfile.getNumColumns()
    # loop through the columns and store the muscle activations
    for i in range(numCols):
        temp = labels[i]
        if '_r/activation' in temp:
            tempact = trialactfile.getDependentColumn(temp).to_numpy()
            tempact = tempact.reshape(len(tempact),1)
            tempact = np.interp(time100.flatten(), trialtime, tempact.flatten())
            # check if the key exists in the dictionary
            if temp not in muscleacts:
                muscleacts[temp] = tempact.reshape(len(tempact),1)
            else: # add the data as a new column in the same key
                muscleacts[temp] = np.column_stack((muscleacts[temp], tempact))
    return muscleacts

# function for getting the muscle moments. 
def getJointMoments(trialdir, moments, modelmass):
    # load the right file for the given trial 
    trialmomfile = osim.TimeSeriesTable(os.path.join(trialdir, 'muscletrack_redo_moments_py.sto'))
    trialtime = trialmomfile.getIndependentColumn()
    time100 = np.linspace(trialtime[0], trialtime[-1], 100).reshape(100,1)
    # get the column labels
    labels = trialmomfile.getColumnLabels()
    # get the number of columns
    numCols = trialmomfile.getNumColumns()
    # loop through the columns and store the muscle moments
    for i in range(numCols):
        temp = labels[i]
        if '_l_moment' not in temp:
            tempmom = trialmomfile.getDependentColumn(temp).to_numpy()
            tempmom = tempmom / modelmass
            tempmom = tempmom.reshape(len(tempmom),1)
            tempmom = np.interp(time100.flatten(), trialtime, tempmom.flatten())
            # check if the key exists in the dictionary
            if temp not in moments:
                moments[temp] = tempmom.reshape(len(tempmom),1)
            else: # add the data as a new column in the same key
                moments[temp] = np.column_stack((moments[temp], tempmom))
    return moments

# function to gather and plot exotendon tension from simulations. 
def getExotendonTension(subject, trial, trialdir, exoTensions):
    subjecttrial = subject + '_' + trial
    
    # load the right file for the given trial 
    trialtensionfile = osim.TimeSeriesTable(os.path.join(trialdir, 'muscletrack_redo_exotendon_py.sto'))
    trialtime = trialtensionfile.getIndependentColumn()
    time100 = np.linspace(trialtime[0], trialtime[-1], 100).reshape(100,1)
    # get the column labels
    labels = trialtensionfile.getColumnLabels()
    # get the number of columns
    numCols = trialtensionfile.getNumColumns()
    # loop through the columns and store the exotendon tensions
    for i in range(numCols):
        temp = labels[i]
        if 'tension' in temp:
            temptension = trialtensionfile.getDependentColumn(temp).to_numpy()
            temptension = temptension
            temptension = temptension.reshape(len(temptension),1)
            temptension = np.interp(time100.flatten(), trialtime, temptension.flatten())
            # check if the key exists in the dictionary
            if subjecttrial not in exoTensions:
                exoTensions[subjecttrial] = temptension.reshape(len(temptension),1)
            else: # add the data as a new column in the same key
                exoTensions[subjecttrial] = np.column_stack((exoTensions[subjecttrial], temptension))
    return exoTensions

# get the moments from ID in the same way that we get the simulation ones
def getIDMoments(trialdir, moments, modelmass):
    # load the right file for the given trial
    trialmomfile = osim.TimeSeriesTable(os.path.join(trialdir, 'IDactual', 'inverse_dynamics.sto'))
    # get the start and stop times for this trial - same as the simulation script
    cycles = subjectGaitTimings()
    # get the subject name
    pathnames = trialdir.split('\\')
    subjectname = pathnames[-3]
    conditionname = pathnames[-2]
    trialname = pathnames[-1]
    # get the trial key
    trialkey_i = subjectname + '_' + conditionname + '_i'
    trialkey_f = subjectname + '_' + conditionname + '_f'
    if '01' in trialname:
        gait_start = cycles[trialkey_i][0]
        gait_end = cycles[trialkey_f][0]
    elif '02' in trialname:
        gait_start = cycles[trialkey_i][1]
        gait_end = cycles[trialkey_f][1]
    elif '03' in trialname:
        gait_start = cycles[trialkey_i][2]
        gait_end = cycles[trialkey_f][2]
    elif '04' in trialname:
        gait_start = cycles[trialkey_i][3]
        gait_end = cycles[trialkey_f][3]
    else:
        print('trial name not recognized')
        return

    ## ISSUE: welk002 has ID for an early gait cycle, 
    # and then we solve for much later in the trial
    # need to adjust the gait start and end times for this trial, 
    # or handle the timings separately. 

    # okay so we have the file, now we need to get the data from it
    trialtime = trialmomfile.getIndependentColumn()
    trialtime = np.array(trialtime)
    # get the column labels
    labels = trialmomfile.getColumnLabels()
    # get the number of columns
    numCols = trialmomfile.getNumColumns()
    # loop through the columns. and store the data that we want
    for i in range(numCols):
        # first get the column name
        temp = labels[i]
        # check if the column is a moment
        if 'moment' in temp or 'force' in temp:
            tempmom = trialmomfile.getDependentColumn(temp).to_numpy()
            # okay now we have to shorten to the right time frame and get it all situated. 
            # Create a mask for the time values that fall between gait_start and gait_end
            try:
                mask = (trialtime >= gait_start) & (trialtime <= gait_end)
                # Apply the mask to trialtime and tempmom
                trialtime_shortened = trialtime[mask]
                tempmom_shortened = tempmom[mask]
                tempmom_interpolated = np.interp(np.linspace(gait_start, gait_end, 100), trialtime_shortened, tempmom_shortened).flatten()
            except:
                trialtime_shortened = trialtime
                tempmom_shortened = tempmom
                tempmom_interpolated = np.interp(np.linspace(trialtime[0], trialtime[-1], 100), trialtime_shortened, tempmom_shortened).flatten()
            
            tempmom_interpolated = tempmom_interpolated / modelmass
            # Check if the key exists in the dictionary
            if temp not in moments:
                moments[temp] = tempmom_interpolated.reshape(len(tempmom_interpolated), 1)
            else:  # Add the data as a new column in the same key
                moments[temp] = np.column_stack((moments[temp], tempmom_interpolated))
    return moments

# function for gathering the muscle forces
def getMuscleForces(trialdir, activeForces, passiveForces, totalForces, modelmass):
    # load in all the right files
    forcesdata = osim.TimeSeriesTable(os.path.join(trialdir, 'muscletrack_redo_muscleforces_py.sto'))
    forcetime = forcesdata.getIndependentColumn()
    time100 = np.linspace(forcetime[0], forcetime[-1], 100).reshape(100,1)
    # get the column labels
    labels = forcesdata.getColumnLabels()
    # get the number of columns
    numCols = forcesdata.getNumColumns()
    # loop through and store the forces data
    for i in range(numCols):
        temp = labels[i]
        if '_r|active_fiber_force_along_tendon' in temp:
            tempact = forcesdata.getDependentColumn(temp).to_numpy()
            tempact = tempact / modelmass
            tempact = tempact.reshape(len(tempact),1)
            tempact = np.interp(time100.flatten(), forcetime, tempact.flatten())
            # check if the key exists in the dictionary
            if temp not in activeForces:
                activeForces[temp] = tempact.reshape(len(tempact),1)
            else: # add the data as a new column in the same key
                activeForces[temp] = np.column_stack((activeForces[temp], tempact))
        elif '_r|passive_fiber_force_along_tendon' in temp:
            temppass = forcesdata.getDependentColumn(temp).to_numpy()
            temppass = temppass / modelmass
            temppass = temppass.reshape(len(temppass),1)
            temppass = np.interp(time100.flatten(), forcetime, temppass.flatten())
            # check if the key exists in the dictionary
            if temp not in passiveForces:
                passiveForces[temp] = temppass.reshape(len(temppass),1)
            else: # add the data as a new column in the same key
                passiveForces[temp] = np.column_stack((passiveForces[temp], temppass))
        elif '_r|fiber_force_along_tendon' in temp:
            temptotal = forcesdata.getDependentColumn(temp).to_numpy()
            temptotal = temptotal / modelmass
            temptotal = temptotal.reshape(len(temptotal),1)
            temptotal = np.interp(time100.flatten(), forcetime, temptotal.flatten())
            # check if the key exists in the dictionary
            if temp not in totalForces:
                totalForces[temp] = temptotal.reshape(len(temptotal),1)
            else: # add the data as a new column in the same key
                totalForces[temp] = np.column_stack((totalForces[temp], temptotal))
    return activeForces, passiveForces, totalForces
