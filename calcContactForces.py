import os
# os.add_dll_directory('C:/Users/jonstingel/opensim/opensim-core-4.5-2024-05-15-a1a2282/bin')
os.add_dll_directory('C:/Users/jonstingel/opensim-core-4.5.1-2024-08-23-cf3ef35/bin')

import opensim as osim
import numpy as np
import pdb
import time
import matplotlib.pyplot as plt
import scipy
import sys
import OsimUtilityfunctions as ouf
import pandas as pd
from scipy.interpolate import interp1d
from matplotlib.ticker import FormatStrFormatter

# naturalcolor = '#fdb863'
ncolor = '#e66101'
ncolorlight = '#fee0b6'

# ncolor1 = '#fff5eb'
# ncolor2 = '#fee6ce'
# ncolor3 = '#fdd0a2'
# ncolor4 = '#fdae6b'
# ncolor5 = '#fd8d3c'
# ncolor6 = '#f16913'
# ncolor7 = '#d94801'
# ncolor8 = '#a63603'
# ncolor9 = '#7f2704'

ncolor1 = '#7f2704'
ncolor2 = '#a63603'
ncolor3 = '#d94801'
ncolor4 = '#f16913'
ncolor5 = '#fd8d3c'
ncolor6 = '#fdae6b'
ncolor7 = '#fdd0a2'
ncolor8 = '#fee6ce'
ncolor9 = '#fff5eb'


# plt.rcParams['font.family'] = 'Times New Roman'
# exotendoncolor = '#f1a340'
ecolor = '#5e3c99'
ecolormid = '#b2abd2'
ecolormiddark = '#807dba'
ecolorlight = '#d8daeb'

# function from Nick Bianco - not used in script, but used as reference for moco
def calc_knee_reaction_force(self, root_dir, solution):
    modelProcessor = osim.ModelProcessor()
    model = modelProcessor.process()
    jr = osim.analyzeSpatialVec(model, solution,
                                ['.*walker_knee.*reaction_on_parent.*'])
    jr = jr.flatten(['_mx', '_my', '_mz', '_fx', '_fy', '_fz'])
    traj = np.empty(jr.getNumRows())
    max = -np.inf
    for itime in range(jr.getNumRows()):
        for irxn in range(int(jr.getNumColumns() / 6)):
            fx = jr.getDependentColumnAtIndex(6 * irxn + 3)[itime]
            fy = jr.getDependentColumnAtIndex(6 * irxn + 4)[itime]
            fz = jr.getDependentColumnAtIndex(6 * irxn + 5)[itime]
            norm = np.sqrt(fx**2 + fy**2 + fz**2)
            traj[itime] = norm
            max = np.max([norm, max])
    time = jr.getIndependentColumn()
    avg = np.trapz(traj, x=time) / (time[-1] - time[0])
    # fig = plt.figure()
    # ax = fig.add_subplot(1, 1, 1)
    # ax.plot(time, traj)
    # plt.show()
    g = np.abs(model.get_gravity()[1])
    state = model.initSystem()
    mass = model.getTotalMass(state)
    weight = mass * g
    return max / weight, avg / weight
# simple helper to get the model mass
def get_model_total_mass(wkdir, filename):
    model = osim.Model(os.path.join(wkdir, filename))
    modelmass = model.getTotalMass(model.initSystem())
    grav = np.abs(model.get_gravity()[1])
    return modelmass #*grav
# function that actually runs the analysis to compute knee contact force, and
# transforms it to the tibia frame
def computeKneeContact(trimmodel, initTime, finalTime, trialdir, tag, whichleg, runtool):
    # '''
    # intersegmental forces - method 2
    # try the analysis
    jr_tool = osim.AnalyzeTool()
    jr_tool.setName('jr_analysis_100con')
    # jr_tool.setModelFilename(os.path.join(trialdir, 'post_simple_model_all_the_probes_muscletrack.osim'))
    
    # I don't think this is going to work. 
    # jr_tool.setStatesStorage(statesStorage)
    # jr_tool.setStatesFileName('testfibsolution.sto')
    trimmingstates = osim.Storage('trimmingStates_' + tag + '.sto')
    jr_tool.setStatesStorage(trimmingstates)    
    
    # jr_tool.setExternalLoadsFileName('grf_walk.xml')
    # jr_tool.updControllerSet().cloneAndAppend(osim.PrescribedController(os.path.join(trialdir, 'muscletrack_controls_100con.sto')))
    jr_tool.updControllerSet().cloneAndAppend(osim.PrescribedController(os.path.join(trialdir, 'trimmingControls_' + tag + '.sto')))
    
    jra = osim.JointReaction()
    jra.setName('jra_' + tag)
    wherestr = osim.ArrayStr(); wherestr.append('child')
    jra.setInFrame(wherestr)
    
    jr_tool.updAnalysisSet().cloneAndAppend(jra)
    jr_tool.setInitialTime(initTime)
    jr_tool.setFinalTime(finalTime)
    jr_tool.setResultsDir(trialdir)
    
    # jr_tool.setModelFilename('jratestingmodel.osim')
    trimmodel.addAnalysis(jra)
    jr_tool.setModel(trimmodel)

    ## uncomment to rerun the analysis
    jr_tool.printToXML(os.path.join(trialdir, 'jr_setup.xml'))
    # time.sleep(0.5)
    # jr_tool = osim.AnalyzeTool(os.path.join(trialdir, 'jr_setup.xml'))
    
    if runtool:
        jr_tool.run()
        time.sleep(0.5)
    else:
        print('\n TOOL NOT RUNNING - JUST COLLECTING OLD RESULTS')

    # figure out how to do an intersegmental with proper value of forces showing up.
    trimjra = osim.TimeSeriesTable('jr_analysis_100con_jra_' + tag + '_ReactionLoads.sto')
    trimjralabels = trimjra.getColumnLabels()
    tiby = trimjra.getDependentColumn('walker_knee_r_on_tibia_r_in_tibia_r_fy').to_numpy()

    # import matplotlib.pyplot as plt
    # plt.figure()
    # plt.plot(np.array(trimjra.getIndependentColumn()), tiby)
    return tiby
# this is the same script ^ but for the prescribed mocoinverse solution
def computeKneeContactRedo(trimmodel, initTime, finalTime, trialdir, tag, whichleg, runtool):
    print('Joint reaction analysis for ' + tag)

    # intersegmental forces - method 2
    # try the analysis
    jr_tool = osim.AnalyzeTool()
    jr_tool.setName('jr_analysis_redo')
    # jr_tool.setModelFilename(os.path.join(trialdir, 'post_simple_model_all_the_probes_muscletrack.osim'))
    
    # I don't think this is going to work. 
    # jr_tool.setStatesStorage(statesStorage)
    # jr_tool.setStatesFileName('testfibsolution.sto')
    trimmingstates = osim.Storage('trimmingStates_redo_' + tag + '.sto')
    jr_tool.setStatesStorage(trimmingstates)    
    
    # jr_tool.setExternalLoadsFileName('grf_walk.xml')
    # jr_tool.updControllerSet().cloneAndAppend(osim.PrescribedController(os.path.join(trialdir, 'muscletrack_controls_100con.sto')))
    jr_tool.updControllerSet().cloneAndAppend(osim.PrescribedController(os.path.join(trialdir, 'trimmingControls_redo_' + tag + '.sto')))
    
    jra = osim.JointReaction()
    jra.setName('jra_redo_' + tag)
    wherestr = osim.ArrayStr(); wherestr.append('child')
    jra.setInFrame(wherestr)
    
    jr_tool.updAnalysisSet().cloneAndAppend(jra)
    jr_tool.setInitialTime(initTime)
    jr_tool.setFinalTime(finalTime)
    jr_tool.setResultsDir(trialdir)
    
    # jr_tool.setModelFilename('jratestingmodel.osim')
    trimmodel.addAnalysis(jra)
    jr_tool.setModel(trimmodel)

    ## uncomment to rerun the analysis
    jr_tool.printToXML(os.path.join(trialdir, 'jr_setup_redo.xml'))
    # time.sleep(0.5)
    # jr_tool = osim.AnalyzeTool(os.path.join(trialdir, 'jr_setup.xml'))
    
    # commented out since I have all the results currently processed
    # unmcomment in order to actually reproduce the files with new results. 
    
    if runtool:
        jr_tool.run()
        time.sleep(0.5)
    else:
        print('\n TOOL NOT RUNNING - JUST COLLECTING OLD RESULTS')
    
    # '''
    # figure out how to do an intersegmental with proper value of forces showing up.
    trimjra = osim.TimeSeriesTable('jr_analysis_redo_jra_redo_' + tag + '_ReactionLoads.sto')
    trimjralabels = trimjra.getColumnLabels()
    if whichleg == 'right':
        tibx = trimjra.getDependentColumn('walker_knee_r_on_tibia_r_in_tibia_r_fx').to_numpy()
        tiby = trimjra.getDependentColumn('walker_knee_r_on_tibia_r_in_tibia_r_fy').to_numpy()
        tibz = trimjra.getDependentColumn('walker_knee_r_on_tibia_r_in_tibia_r_fz').to_numpy()
    elif whichleg == 'left':
        tibx = trimjra.getDependentColumn('walker_knee_l_on_tibia_l_in_tibia_l_fx').to_numpy()
        tiby = trimjra.getDependentColumn('walker_knee_l_on_tibia_l_in_tibia_l_fy').to_numpy()
        tibz = trimjra.getDependentColumn('walker_knee_l_on_tibia_l_in_tibia_l_fz').to_numpy()
    elif whichleg == 'both':
        tibxr = trimjra.getDependentColumn('walker_knee_r_on_tibia_r_in_tibia_r_fx').to_numpy()
        tibxl = trimjra.getDependentColumn('walker_knee_l_on_tibia_l_in_tibia_l_fx').to_numpy()
        tibyr = trimjra.getDependentColumn('walker_knee_r_on_tibia_r_in_tibia_r_fy').to_numpy()
        tibyl = trimjra.getDependentColumn('walker_knee_l_on_tibia_l_in_tibia_l_fy').to_numpy()
        tibzr = trimjra.getDependentColumn('walker_knee_r_on_tibia_r_in_tibia_r_fz').to_numpy()
        tibzl = trimjra.getDependentColumn('walker_knee_l_on_tibia_l_in_tibia_l_fz').to_numpy()
        # tiby = np.array([tibyr, tibyl])

        # do a little peak finding, then get the index, and orient the index to match the peak in the right
        # then add the two together.
        # Find the peak value and index for tibyl
        peak_tibyl = np.min(tibyl)
        peak_index_tibyl = np.argmin(tibyl)

        # Find the peak value and index for tibyr
        peak_tibyr = np.min(tibyr)
        peak_index_tibyr = np.argmin(tibyr)

        # Rotate tibyl to match the peak index of tibyr
        tibyl_rotated = np.roll(tibyl, peak_index_tibyr - peak_index_tibyl)
        tibxl_rotated = np.roll(tibxl, peak_index_tibyr - peak_index_tibyl)
        tibzl_rotated = np.roll(tibzl, peak_index_tibyr - peak_index_tibyl)
        
        # print('Peak index tibyl: ' + str(peak_index_tibyl))
        # print('Peak index tibyr: ' + str(peak_index_tibyr))
        # print('Peak value tibyl: ' + str(peak_tibyl))
        # print('Peak value tibyr: ' + str(peak_tibyr))
        # pdb.set_trace()
        # plt.plot(tibyr, label='tibyr')
        # plt.plot(tibyl, label='tibyl')
        # plt.plot(tibyl_rotated, label='tibyl_rotated')
        # plt.legend()
        # plt.show()
        # pdb.set_trace()
        # plt.figure()
        # plt.plot(tibxr, label='tibxr')
        # plt.plot(tibxl, label='tibxl')
        # plt.plot(tibxl_rotated, label='tibxl_rotated')
        # plt.legend()
        # plt.show()
        # pdb.set_trace()
        # plt.figure()
        # plt.plot(tibzr, label='tibzr')
        # plt.plot(tibzl, label='tibzl')
        # plt.plot(-tibzl_rotated, label='tibzl_rotated')
        # plt.legend()
        # plt.show()
        # pdb.set_trace()

        # Combine the right and left data
        tibx = (tibxr+tibxl_rotated)/2
        tiby = (tibyr+tibyl_rotated)/2
        tibz = (tibzr-tibzl_rotated)/2 # z is negated since it is in the opposite side of the model. 

    # pdb.set_trace()
    # import matplotlib.pyplot as plt
    # plt.figure()
    # plt.plot(np.array(trimjra.getIndependentColumn()), tiby)
    return tibx, tiby, tibz
# this is the same script ^ but for the prescribed mocoinverse solution
def computeKneeContactRedoPoly(trimmodel, initTime, finalTime, trialdir, tag, whichleg, runtool):
    print('Joint reaction analysis for ' + tag)

    # intersegmental forces - method 2
    # try the analysis
    jr_tool = osim.AnalyzeTool()
    jr_tool.setName('jr_analysis_redo_poly')
    # jr_tool.setModelFilename(os.path.join(trialdir, 'post_simple_model_all_the_probes_muscletrack.osim'))
    
    # I don't think this is going to work. 
    # jr_tool.setStatesStorage(statesStorage)
    # jr_tool.setStatesFileName('testfibsolution.sto')
    trimmingstates = osim.Storage('trimmingStates_redo_poly_' + tag + '.sto')
    jr_tool.setStatesStorage(trimmingstates)    
    
    # jr_tool.setExternalLoadsFileName('grf_walk.xml')
    # jr_tool.updControllerSet().cloneAndAppend(osim.PrescribedController(os.path.join(trialdir, 'muscletrack_controls_100con.sto')))
    jr_tool.updControllerSet().cloneAndAppend(osim.PrescribedController(os.path.join(trialdir, 'trimmingControls_redo_poly_' + tag + '.sto')))
    
    jra = osim.JointReaction()
    jra.setName('jra_redo_poly' + tag)
    wherestr = osim.ArrayStr(); wherestr.append('child')
    jra.setInFrame(wherestr)
    
    jr_tool.updAnalysisSet().cloneAndAppend(jra)
    jr_tool.setInitialTime(initTime)
    jr_tool.setFinalTime(finalTime)
    jr_tool.setResultsDir(trialdir)
    
    # jr_tool.setModelFilename('jratestingmodel.osim')
    trimmodel.addAnalysis(jra)
    jr_tool.setModel(trimmodel)

    ## uncomment to rerun the analysis
    jr_tool.printToXML(os.path.join(trialdir, 'jr_setup_redo_poly.xml'))
    # time.sleep(0.5)
    # jr_tool = osim.AnalyzeTool(os.path.join(trialdir, 'jr_setup.xml'))
    
    # commented out since I have all the results currently processed
    # unmcomment in order to actually reproduce the files with new results. 
    
    if runtool:
        jr_tool.run()
        time.sleep(0.5)
    else:
        print('\n TOOL NOT RUNNING - JUST COLLECTING OLD RESULTS')
    
    # '''
    # figure out how to do an intersegmental with proper value of forces showing up.
    trimjra = osim.TimeSeriesTable('jr_analysis_redo_poly_jra_redo_' + tag + '_ReactionLoads.sto')
    trimjralabels = trimjra.getColumnLabels()
    if whichleg == 'right':
        tiby = trimjra.getDependentColumn('walker_knee_r_on_tibia_r_in_tibia_r_fy').to_numpy()
    elif whichleg == 'left':
        tiby = trimjra.getDependentColumn('walker_knee_l_on_tibia_l_in_tibia_l_fy').to_numpy()
    elif whichleg == 'both':
        tibyr = trimjra.getDependentColumn('walker_knee_r_on_tibia_r_in_tibia_r_fy').to_numpy()
        tibyl = trimjra.getDependentColumn('walker_knee_l_on_tibia_l_in_tibia_l_fy').to_numpy()
        # tiby = np.array([tibyr, tibyl])

        # do a little peak finding, then get the index, and orient the index to match the peak in the right
        # then add the two together.
        # Find the peak value and index for tibyl
        peak_tibyl = np.min(tibyl)
        peak_index_tibyl = np.argmin(tibyl)

        # Find the peak value and index for tibyr
        peak_tibyr = np.min(tibyr)
        peak_index_tibyr = np.argmin(tibyr)

        # Rotate tibyl to match the peak index of tibyr
        tibyl_rotated = np.roll(tibyl, peak_index_tibyr - peak_index_tibyl)
        
        # print('Peak index tibyl: ' + str(peak_index_tibyl))
        # print('Peak index tibyr: ' + str(peak_index_tibyr))
        # print('Peak value tibyl: ' + str(peak_tibyl))
        # print('Peak value tibyr: ' + str(peak_tibyr))
        # pdb.set_trace()
        # plt.plot(tibyr, label='tibyr')
        # plt.plot(tibyl, label='tibyl')
        # plt.plot(tibyl_rotated, label='tibyl_rotated')
        # plt.legend()
        # plt.show()
        # pdb.set_trace()

        # Combine the right and left data
        tiby = (tibyr+tibyl_rotated)/2
         

    # pdb.set_trace()
    # import matplotlib.pyplot as plt
    # plt.figure()
    # plt.plot(np.array(trimjra.getIndependentColumn()), tiby)
    return tiby
# this is the same script ^ but for the prescribed mocoinverse solution
def computeKneeContactPrescribe(trimmodel, initTime, finalTime, trialdir, tag):
    # '''
    # intersegmental forces - method 2
    # try the analysis
    jr_tool = osim.AnalyzeTool()
    jr_tool.setName('jr_analysis_prescribe')
    # jr_tool.setModelFilename(os.path.join(trialdir, 'post_simple_model_all_the_probes_muscletrack.osim'))
    
    # I don't think this is going to work. 
    # jr_tool.setStatesStorage(statesStorage)
    # jr_tool.setStatesFileName('testfibsolution.sto')
    trimmingstates = osim.Storage('trimmingStates_prescribe_' + tag + '.sto')
    jr_tool.setStatesStorage(trimmingstates)    
    
    # jr_tool.setExternalLoadsFileName('grf_walk.xml')
    # jr_tool.updControllerSet().cloneAndAppend(osim.PrescribedController(os.path.join(trialdir, 'muscletrack_controls_100con.sto')))
    jr_tool.updControllerSet().cloneAndAppend(osim.PrescribedController(os.path.join(trialdir, 'trimmingControls_prescribe_' + tag + '.sto')))
    
    jra = osim.JointReaction()
    jra.setName('jra_prescribe_' + tag)
    wherestr = osim.ArrayStr(); wherestr.append('child')
    jra.setInFrame(wherestr)
    
    jr_tool.updAnalysisSet().cloneAndAppend(jra)
    jr_tool.setInitialTime(initTime)
    jr_tool.setFinalTime(finalTime)
    jr_tool.setResultsDir(trialdir)
    
    # jr_tool.setModelFilename('jratestingmodel.osim')
    trimmodel.addAnalysis(jra)
    jr_tool.setModel(trimmodel)

    ## uncomment to rerun the analysis
    jr_tool.printToXML(os.path.join(trialdir, 'jr_setup_prescribe.xml'))
    # time.sleep(0.5)
    # jr_tool = osim.AnalyzeTool(os.path.join(trialdir, 'jr_setup.xml'))
    jr_tool.run()
    time.sleep(0.5)
    # '''
    # figure out how to do an intersegmental with proper value of forces showing up.
    trimjra = osim.TimeSeriesTable('jr_analysis_prescribe_jra_prescribe_' + tag + '_ReactionLoads.sto')
    trimjralabels = trimjra.getColumnLabels()
    tiby = trimjra.getDependentColumn('walker_knee_r_on_tibia_r_in_tibia_r_fy').to_numpy()
    # pdb.set_trace()
    # import matplotlib.pyplot as plt
    # plt.figure()
    # plt.plot(np.array(trimjra.getIndependentColumn()), tiby)
    return tiby
# method for computing the individual muscle contributions to knee contact force
def getKneeContactributions(trialdir, musclesWanted_split, tag, whichleg, runtool):
    # os.chdir(trialdir)
    # solution = osim.MocoTrajectory('muscle_statetrack_grfprescribe_solution_100con.sto')
    # model = osim.Model('post_simple_model_all_the_probes_muscletrack.osim')
    # # attempting with just setting muscles to not apply force
    # muscles = model.getMuscles()
    # for m in range(muscles.getSize()):
    #     musc = muscles.get(m)
    #     muscname = musc.getName()
    #     if muscname not in musclesWanted_split:
    #         musc.setMaxIsometricForce(0.0)
    #     # else:
    #     #     print(muscname)
    
    # model.initSystem()
    if musclesWanted_split == []: # this is the inter or intersegmental condition
        statesStorage = osim.Storage('muscletrack_states_100con.sto')
        statesTable = osim.TimeSeriesTable('muscletrack_states_100con.sto')
        stateslabels = statesTable.getColumnLabels()
        statesTableTrim = osim.TimeSeriesTable('muscletrack_states_100con.sto')
        statesTableTrim2 = osim.TimeSeriesTable('muscletrack_states_100con.sto')
        
        # get a trimmed set of states that is just the joint angles speeds, and whatever muscles you want. 
        for stat in stateslabels: 
            if 'forceset' in stat:
                statesTableTrim.removeColumn(stat)
        osim.STOFileAdapter.write(statesTableTrim, 'trimmingStates_' + tag + '.sto')

        for stat in stateslabels:
            if 'forceset' in stat or 'speed' in stat: 
                statesTableTrim2.removeColumn(stat)
        osim.STOFileAdapter.write(statesTableTrim2, 'trimmingStates2_' + tag + '.sto')

        # get a version of the controls that matches. 
        controlsTable = osim.TimeSeriesTable('muscletrack_controls_100con.sto')
        controlslabels = controlsTable.getColumnLabels()
        controlsTableTrim = osim.TimeSeriesTable('muscletrack_controls_100con.sto')
        
        # get a version of the controls that is trimmed down. 
        for con in controlslabels:
            if 'reserve' not in con and 'lumbar' not in con:
                controlsTableTrim.removeColumn(con)
        osim.STOFileAdapter.write(controlsTableTrim, 'trimmingControls_' + tag + '.sto')
        
        # get a version of the model with no muscles in it
        # muscmodel = osim.Model('post_simple_model_all_the_probes_muscletrack.osim')
        trimmodelprocessor = osim.ModelProcessor('post_simple_model_all_the_probes_muscletrack.osim')
        trimmodelprocessor.append(osim.ModOpRemoveMuscles())
        trimmodel = trimmodelprocessor.process()
        trimmodel.initSystem()
        trimmodel.printToXML('trimmingmodel_' + tag + '.osim')
        
        initTime = np.array(statesTable.getIndependentColumn())[0]
        finalTime = np.array(statesTable.getIndependentColumn())[-1]
        
        # now try the positionMotion
        # pomostorage = osim.Storage('trimmingStates2_inter.sto')
        pomostorage = osim.Storage('trimmingStates_' + tag + '.sto')
        pomotraj = osim.StatesTrajectory.createFromStatesStorage(trimmodel, pomostorage)
        pomo = osim.PositionMotion.createFromStatesTrajectory(trimmodel, pomotraj)
        
        trimmodel.addComponent(pomo)
        trimmodel.initSystem()
        trimmodel.printToXML('trimmingmodel2' + tag + '.osim')
        
        # call the analyze tool to actually do the analysis and get the values. 
        jray = computeKneeContact(trimmodel, initTime, finalTime, trialdir, tag, whichleg, runtool)
    
    
    elif 'all' in musclesWanted_split:
        statesStorage = osim.Storage('muscletrack_states_100con.sto')
        statesTable = osim.TimeSeriesTable('muscletrack_states_100con.sto')
        stateslabels = statesTable.getColumnLabels()
        statesTableTrim = osim.TimeSeriesTable('muscletrack_states_100con.sto')
        statesTableTrim2 = osim.TimeSeriesTable('muscletrack_states_100con.sto')
        # in this case, we want all the muscles in the model
        osim.STOFileAdapter.write(statesTableTrim, 'trimmingStates_' + tag + '.sto')
        
        for stat in stateslabels:
            if 'speed' in stat:
                statesTableTrim2.removeColumn(stat)
        osim.STOFileAdapter.write(statesTableTrim2, 'trimmingStates2_' + tag + '.sto')


        # get a version of the controls that matches. 
        controlsTable = osim.TimeSeriesTable('muscletrack_controls_100con.sto')
        controlslabels = controlsTable.getColumnLabels()
        controlsTableTrim = osim.TimeSeriesTable('muscletrack_controls_100con.sto')
        
        # in this case we want all the controls, not getting rid of any muscles
        osim.STOFileAdapter.write(controlsTableTrim, 'trimmingControls_' + tag + '.sto')
        
        
        # get a version of the model with no muscles in it
        # muscmodel = osim.Model('post_simple_model_all_the_probes_muscletrack.osim')
        trimmodel = osim.Model('post_simple_model_all_the_probes_muscletrack.osim')
        
        # model keeping all the muscles again for this one
        trimmodel.initSystem()
        trimmodel.printToXML('trimmingmodel_' + tag + '.osim')
        
        initTime = np.array(statesTable.getIndependentColumn())[0]
        finalTime = np.array(statesTable.getIndependentColumn())[-1]
        
        # now try the positionMotion
        pomostorage = osim.Storage('trimmingStates_' + tag + '.sto')
        pomotraj = osim.StatesTrajectory.createFromStatesStorage(trimmodel, pomostorage)
        pomo = osim.PositionMotion.createFromStatesTrajectory(trimmodel, pomotraj)
        
        trimmodel.addComponent(pomo)
        trimmodel.initSystem()
        trimmodel.printToXML('trimmingmodel2' + tag + '.osim')
    
        ###
        jray = computeKneeContact(trimmodel, initTime, finalTime, trialdir, tag, whichleg, runtool)
    
    elif 'reserve' in musclesWanted_split:
        statesStorage = osim.Storage('muscletrack_states_100con.sto')
        statesTable = osim.TimeSeriesTable('muscletrack_states_100con.sto')
        stateslabels = statesTable.getColumnLabels()
        statesTableTrim = osim.TimeSeriesTable('muscletrack_states_100con.sto')
        statesTableTrim2 = osim.TimeSeriesTable('muscletrack_states_100con.sto')
        
        # get a trimmed set of states that is just the joint angles speeds, and whatever muscles you want. 
        for stat in stateslabels: 
            if 'forceset' in stat:
                statesTableTrim.removeColumn(stat)
        osim.STOFileAdapter.write(statesTableTrim, 'trimmingStates_' + tag + '.sto')

        for stat in stateslabels:
            if 'forceset' in stat or 'speed' in stat: 
                statesTableTrim2.removeColumn(stat)
        osim.STOFileAdapter.write(statesTableTrim2, 'trimmingStates2_' + tag + '.sto')

        # get a version of the controls that matches. 
        controlsTable = osim.TimeSeriesTable('muscletrack_controls_100con.sto')
        controlslabels = controlsTable.getColumnLabels()
        controlsTableTrim = osim.TimeSeriesTable('muscletrack_controls_100con.sto')
        
        # get a version of the controls that is trimmed down. 
        for con in controlslabels:
            if 'lumbar' not in con:
                controlsTableTrim.removeColumn(con)
        osim.STOFileAdapter.write(controlsTableTrim, 'trimmingControls_' + tag + '.sto')
        
        # get a version of the model with no muscles in it (or reserves?)
        # muscmodel = osim.Model('post_simple_model_all_the_probes_muscletrack.osim')
        trimmodelprocessor = osim.ModelProcessor('post_simple_model_all_the_probes_muscletrack.osim')
        trimmodelprocessor.append(osim.ModOpRemoveMuscles())
        trimmodel = trimmodelprocessor.process()
        trimmodel.initSystem()
        
        # now have to do it for the forceset too?
        trimforces = trimmodel.getForceSet()
        numForces = trimforces.getSize()  
        count = 0
        for f in range(numForces):
            fo = trimforces.get(f-count)
            # print(fo.getName())
            if 'lumbar' not in fo.getName() and 'HOBL' not in fo.getName():
                getrid = True
                trimforces.remove(f-count)
                count += 1
        
        trimmodel.initSystem()
        trimmodel.printToXML('trimmingmodel_' + tag + '.osim')
        
        initTime = np.array(statesTable.getIndependentColumn())[0]
        finalTime = np.array(statesTable.getIndependentColumn())[-1]
        
        # now try the positionMotion
        # pomostorage = osim.Storage('trimmingStates2_inter.sto')
        pomostorage = osim.Storage('trimmingStates_' + tag + '.sto')
        pomotraj = osim.StatesTrajectory.createFromStatesStorage(trimmodel, pomostorage)
        pomo = osim.PositionMotion.createFromStatesTrajectory(trimmodel, pomotraj)
        
        trimmodel.addComponent(pomo)
        trimmodel.initSystem()
        trimmodel.printToXML('trimmingmodel2' + tag + '.osim')
        
        # call the analyze tool to actually do the analysis and get the values. 
        jray = computeKneeContact(trimmodel, initTime, finalTime, trialdir, tag, whichleg, runtool)
    
    elif 'none' in musclesWanted_split:
        statesStorage = osim.Storage('muscletrack_states_100con.sto')
        statesTable = osim.TimeSeriesTable('muscletrack_states_100con.sto')
        stateslabels = statesTable.getColumnLabels()
        statesTableTrim = osim.TimeSeriesTable('muscletrack_states_100con.sto')
        statesTableTrim2 = osim.TimeSeriesTable('muscletrack_states_100con.sto')
        
        # get a trimmed set of states that is just the joint angles speeds, and whatever muscles you want. 
        for stat in stateslabels: 
            if 'forceset' in stat:
                statesTableTrim.removeColumn(stat)
        osim.STOFileAdapter.write(statesTableTrim, 'trimmingStates_' + tag + '.sto')

        for stat in stateslabels:
            if 'forceset' in stat or 'speed' in stat: 
                statesTableTrim2.removeColumn(stat)
        osim.STOFileAdapter.write(statesTableTrim2, 'trimmingStates2_' + tag + '.sto')

        # get a version of the controls that matches. 
        controlsTable = osim.TimeSeriesTable('muscletrack_controls_100con.sto')
        controlslabels = controlsTable.getColumnLabels()
        controlsTableTrim = osim.TimeSeriesTable('muscletrack_controls_100con.sto')
        
        # get a version of the controls that is trimmed down. 
        for con in controlslabels:
            controlsTableTrim.removeColumn(con)
        osim.STOFileAdapter.write(controlsTableTrim, 'trimmingControls_' + tag + '.sto')
        
        # get a version of the model with no muscles in it (or reserves?)
        # muscmodel = osim.Model('post_simple_model_all_the_probes_muscletrack.osim')
        trimmodelprocessor = osim.ModelProcessor('post_simple_model_all_the_probes_muscletrack.osim')
        trimmodelprocessor.append(osim.ModOpRemoveMuscles())
        trimmodel = trimmodelprocessor.process()
        trimmodel.initSystem()
        
        # now have to do it for the forceset too?
        trimforces = trimmodel.getForceSet()
        numForces = trimforces.getSize()
        count = 0
        for f in range(numForces):
            fo = trimforces.get(f-count)
            if 'HOBL' not in fo.getName():
                getrid = True
                trimforces.remove(f-count)
                count += 1
        
        trimmodel.initSystem()
        trimmodel.printToXML('trimmingmodel_' + tag + '.osim')
        
        initTime = np.array(statesTable.getIndependentColumn())[0]
        finalTime = np.array(statesTable.getIndependentColumn())[-1]
        
        # now try the positionMotion
        # pomostorage = osim.Storage('trimmingStates2_inter.sto')
        pomostorage = osim.Storage('trimmingStates_' + tag + '.sto')
        pomotraj = osim.StatesTrajectory.createFromStatesStorage(trimmodel, pomostorage)
        pomo = osim.PositionMotion.createFromStatesTrajectory(trimmodel, pomotraj)
        
        trimmodel.addComponent(pomo)
        trimmodel.initSystem()
        trimmodel.printToXML('trimmingmodel2' + tag + '.osim')
        
        # call the analyze tool to actually do the analysis and get the values. 
        jray = computeKneeContact(trimmodel, initTime, finalTime, trialdir, tag, whichleg, runtool)
        
    
    else:
        statesStorage = osim.Storage('muscletrack_states_100con.sto')
        statesTable = osim.TimeSeriesTable('muscletrack_states_100con.sto')
        stateslabels = statesTable.getColumnLabels()
        statesTableTrim = osim.TimeSeriesTable('muscletrack_states_100con.sto')
        statesTableTrim2 = osim.TimeSeriesTable('muscletrack_states_100con.sto')
        
        # musclesWanted_split = ['1'2'3'4']
        # get a trimmed set of states that is just the joint angles speeds, and whatever muscles you want. 
        for stat in stateslabels:
            if 'forceset' in stat:
                # print(stat)
                getrid = True
                for want in musclesWanted_split:
                    if want in stat:
                        # want this one
                        # print('want this one')
                        # print(stat)
                        getrid = False
                if getrid:
                    statesTableTrim.removeColumn(stat)
        osim.STOFileAdapter.write(statesTableTrim, 'trimmingStates_' + tag + '.sto')

        for stat in stateslabels:
            if 'forceset' in stat or 'speed' in stat:
                getrid = True
                for want in musclesWanted_split:
                    if want in stat:
                        # want to keep this one
                        # print(stat)
                        getrid = False
                if getrid:
                    statesTableTrim2.removeColumn(stat)
        osim.STOFileAdapter.write(statesTableTrim2, 'trimmingStates2_' + tag + '.sto')


        # get a version of the controls that matches. 
        controlsTable = osim.TimeSeriesTable('muscletrack_controls_100con.sto')
        controlslabels = controlsTable.getColumnLabels()
        controlsTableTrim = osim.TimeSeriesTable('muscletrack_controls_100con.sto')
        
        # get a version of the controls that is trimmed down. 
        for con in controlslabels:
            # print(con)
            if'lumbar' not in con:
                getrid = True
                for want in musclesWanted_split:
                    if want in con:
                        # want this
                        # print('want this one')
                        # print(con)
                        getrid = False
                if getrid:
                    controlsTableTrim.removeColumn(con)
        osim.STOFileAdapter.write(controlsTableTrim, 'trimmingControls_' + tag + '.sto')
        
        
        # get a version of the model with no muscles in it
        # muscmodel = osim.Model('post_simple_model_all_the_probes_muscletrack.osim')
        trimmodel = osim.Model('post_simple_model_all_the_probes_muscletrack.osim')
        # get rid of muscles that we don't want
        trimmuscles = trimmodel.getMuscles()
        numMuscles = trimmuscles.getSize()
        count = 0
        for m in range(numMuscles):
            musc = trimmuscles.get(m-count)
            # print(musc.getName())
            getrid = True
            for mu in musclesWanted_split:
                # print(mu)
                if mu == musc.getName():
                    # print(musc.getName())
                    # print(mu)
                    getrid = False
            if getrid:
                trimmuscles.remove(musc)
                count +=1
                
        # now have to do it for the forceset too?
        trimforces = trimmodel.getForceSet()
        numForces = trimforces.getSize()  
        count = 0
        for f in range(numForces):
            fo = trimforces.get(f-count)
            # print(fo.getName())
            if 'lumbar' not in fo.getName() and 'HOBL' not in fo.getName():
                getrid = True
                for mu in musclesWanted_split:
                    if mu == fo.getName():
                        getrid = False
                if getrid:
                    trimforces.remove(f-count)
                    count += 1
        
        # model should only have muscles that we want now. 
        trimmodel.initSystem()
        trimmodel.printToXML('trimmingmodel_' + tag + '.osim')
        
        
        initTime = np.array(statesTable.getIndependentColumn())[0]
        finalTime = np.array(statesTable.getIndependentColumn())[-1]
        
        # now try the positionMotion
        pomostorage = osim.Storage('trimmingStates_' + tag + '.sto')
        pomotraj = osim.StatesTrajectory.createFromStatesStorage(trimmodel, pomostorage)
        pomo = osim.PositionMotion.createFromStatesTrajectory(trimmodel, pomotraj)
        
        trimmodel.addComponent(pomo)
        trimmodel.initSystem()
        trimmodel.printToXML('trimmingmodel2' + tag + '.osim')
    
        ###
        jray = computeKneeContact(trimmodel, initTime, finalTime, trialdir, tag, whichleg, runtool)
    
    return jray
# method for computing the individual muscle contributions to knee contact force
def getKneeContactributionsRedo(trialdir, musclesWanted_split, tag, whichleg, runtool):
    if not runtool:
        print('\n NOT RUNNING TOOL - JUST COLLECTING OLD RESULTS')
        if musclesWanted_split == []:
            # load the model
            trimmodel = osim.Model('trimmingmodel2_redo_' + tag + '.osim')
            # timings - these don't matter since we are not rerunning, just need filler for the function
            initTime = 0.0
            finalTime = 1.0
            jrax, jray, jraz = computeKneeContactRedo(trimmodel, initTime, finalTime, trialdir, tag, whichleg, runtool)
        elif 'all' in musclesWanted_split:
            # load the model
            trimmodel = osim.Model('trimmingmodel2_redo_' + tag + '.osim')
            # timings - these don't matter since we are not rerunning, just need filler for the function
            initTime = 0.0
            finalTime = 1.0
            jrax, jray, jraz = computeKneeContactRedo(trimmodel, initTime, finalTime, trialdir, tag, whichleg, runtool)
        elif 'reserve' in musclesWanted_split:
            # load the model
            trimmodel = osim.Model('trimmingmodel2_redo_' + tag + '.osim')
            # timings - these don't matter since we are not rerunning, just need filler for the function
            initTime = 0.0
            finalTime = 1.0
            jrax, jray, jraz = computeKneeContactRedo(trimmodel, initTime, finalTime, trialdir, tag, whichleg, runtool)
        elif 'none' in musclesWanted_split:
            # load the model
            trimmodel = osim.Model('trimmingmodel2_redo_' + tag + '.osim')
            # timings - these don't matter since we are not rerunning, just need filler for the function
            initTime = 0.0
            finalTime = 1.0
            jrax, jray, jraz = computeKneeContactRedo(trimmodel, initTime, finalTime, trialdir, tag, whichleg, runtool)
        else:
            # load the model
            trimmodel = osim.Model('trimmingmodel2_redo_' + tag + '.osim')
            # timings - these don't matter since we are not rerunning, just need filler for the function
            initTime = 0.0
            finalTime = 1.0
            jrax, jray, jraz = computeKneeContactRedo(trimmodel, initTime, finalTime, trialdir, tag, whichleg, runtool)

    else:
        print('\n RUNNING TOOL')

        # os.chdir(trialdir)
        # solution = osim.MocoTrajectory('muscle_statetrack_grfprescribe_solution_100con.sto')
        # model = osim.Model('post_simple_model_all_the_probes_muscletrack.osim')
        # # attempting with just setting muscles to not apply force
        # muscles = model.getMuscles()
        # for m in range(muscles.getSize()):
        #     musc = muscles.get(m)
        #     muscname = musc.getName()
        #     if muscname not in musclesWanted_split:
        #         musc.setMaxIsometricForce(0.0)
        #     # else:
        #     #     print(muscname)
        
        # model.initSystem()
        if musclesWanted_split == []: # this is the inter or intersegmental condition
            statesStorage = osim.Storage('muscletrack_redo_states_py.sto')
            statesTable = osim.TimeSeriesTable('muscletrack_redo_states_py.sto')
            stateslabels = statesTable.getColumnLabels()
            statesTableTrim = osim.TimeSeriesTable('muscletrack_redo_states_py.sto')
            statesTableTrim2 = osim.TimeSeriesTable('muscletrack_redo_states_py.sto')
            
            # get a trimmed set of states that is just the joint angles speeds, and whatever muscles you want. 
            for stat in stateslabels: 
                if 'forceset' in stat:
                    statesTableTrim.removeColumn(stat)
            osim.STOFileAdapter.write(statesTableTrim, 'trimmingStates_redo_' + tag + '.sto')

            for stat in stateslabels:
                if 'forceset' in stat or 'speed' in stat: 
                    statesTableTrim2.removeColumn(stat)
            osim.STOFileAdapter.write(statesTableTrim2, 'trimmingStates2_redo_' + tag + '.sto')

            # get a version of the controls that matches. 
            controlsTable = osim.TimeSeriesTable('muscletrack_redo_controls_py.sto')
            controlslabels = controlsTable.getColumnLabels()
            controlsTableTrim = osim.TimeSeriesTable('muscletrack_redo_controls_py.sto')
            
            # get a version of the controls that is trimmed down. 
            for con in controlslabels:
                if 'reserve' not in con and 'lumbar' not in con:
                    controlsTableTrim.removeColumn(con)
            osim.STOFileAdapter.write(controlsTableTrim, 'trimmingControls_redo_' + tag + '.sto')
            
            # get a version of the model with no muscles in it
            # muscmodel = osim.Model('post_simple_model_all_the_probes_muscletrack_redo.osim')
            trimmodelprocessor = osim.ModelProcessor('post_simple_model_all_the_probes_muscletrack_redo.osim')
            trimmodelprocessor.append(osim.ModOpRemoveMuscles())
            trimmodel = trimmodelprocessor.process()
            trimmodel.initSystem()
            trimmodel.printToXML('trimmingmodel_redo_' + tag + '.osim')
            
            initTime = np.array(statesTable.getIndependentColumn())[0]
            finalTime = np.array(statesTable.getIndependentColumn())[-1]
            
            # now try the positionMotion
            # pomostorage = osim.Storage('trimmingStates2_inter.sto')
            pomostorage = osim.Storage('trimmingStates_redo_' + tag + '.sto')
            pomotraj = osim.StatesTrajectory.createFromStatesStorage(trimmodel, pomostorage)
            pomo = osim.PositionMotion.createFromStatesTrajectory(trimmodel, pomotraj)
            
            trimmodel.addComponent(pomo)
            trimmodel.initSystem()
            trimmodel.printToXML('trimmingmodel2_redo_' + tag + '.osim')
            
            # call the analyze tool to actually do the analysis and get the values. 
            jrax, jray, jraz = computeKneeContactRedo(trimmodel, initTime, finalTime, trialdir, tag, whichleg, runtool)
        elif 'all' in musclesWanted_split:
            statesStorage = osim.Storage('muscletrack_redo_states_py.sto')
            statesTable = osim.TimeSeriesTable('muscletrack_redo_states_py.sto')
            stateslabels = statesTable.getColumnLabels()
            statesTableTrim = osim.TimeSeriesTable('muscletrack_redo_states_py.sto')
            statesTableTrim2 = osim.TimeSeriesTable('muscletrack_redo_states_py.sto')
            # in this case, we want all the muscles in the model
            osim.STOFileAdapter.write(statesTableTrim, 'trimmingStates_redo_' + tag + '.sto')
            
            for stat in stateslabels:
                if 'speed' in stat:
                    statesTableTrim2.removeColumn(stat)
            osim.STOFileAdapter.write(statesTableTrim2, 'trimmingStates2_redo_' + tag + '.sto')


            # get a version of the controls that matches. 
            controlsTable = osim.TimeSeriesTable('muscletrack_redo_controls_py.sto')
            controlslabels = controlsTable.getColumnLabels()
            controlsTableTrim = osim.TimeSeriesTable('muscletrack_redo_controls_py.sto')
            
            # in this case we want all the controls, not getting rid of any muscles
            osim.STOFileAdapter.write(controlsTableTrim, 'trimmingControls_redo_' + tag + '.sto')
            
            
            # get a version of the model with no muscles in it
            # muscmodel = osim.Model('post_simple_model_all_the_probes_muscletrack_redo.osim')
            trimmodel = osim.Model('post_simple_model_all_the_probes_muscletrack_redo.osim')
            
            # model keeping all the muscles again for this one
            trimmodel.initSystem()
            trimmodel.printToXML('trimmingmodel_redo_' + tag + '.osim')
            
            initTime = np.array(statesTable.getIndependentColumn())[0]
            finalTime = np.array(statesTable.getIndependentColumn())[-1]
            
            # now try the positionMotion
            pomostorage = osim.Storage('trimmingStates_redo_' + tag + '.sto')
            pomotraj = osim.StatesTrajectory.createFromStatesStorage(trimmodel, pomostorage)
            pomo = osim.PositionMotion.createFromStatesTrajectory(trimmodel, pomotraj)
            
            trimmodel.addComponent(pomo)
            trimmodel.initSystem()
            trimmodel.printToXML('trimmingmodel2_redo_' + tag + '.osim')
        
            ###
            jrax, jray, jraz = computeKneeContactRedo(trimmodel, initTime, finalTime, trialdir, tag, whichleg, runtool)
        elif 'reserve' in musclesWanted_split:
            statesStorage = osim.Storage('muscletrack_redo_states_py.sto')
            statesTable = osim.TimeSeriesTable('muscletrack_redo_states_py.sto')
            stateslabels = statesTable.getColumnLabels()
            statesTableTrim = osim.TimeSeriesTable('muscletrack_redo_states_py.sto')
            statesTableTrim2 = osim.TimeSeriesTable('muscletrack_redo_states_py.sto')
            
            # get a trimmed set of states that is just the joint angles speeds, and whatever muscles you want. 
            for stat in stateslabels: 
                if 'forceset' in stat:
                    statesTableTrim.removeColumn(stat)
            osim.STOFileAdapter.write(statesTableTrim, 'trimmingStates_redo_' + tag + '.sto')

            for stat in stateslabels:
                if 'forceset' in stat or 'speed' in stat: 
                    statesTableTrim2.removeColumn(stat)
            osim.STOFileAdapter.write(statesTableTrim2, 'trimmingStates2_redo_' + tag + '.sto')

            # get a version of the controls that matches. 
            controlsTable = osim.TimeSeriesTable('muscletrack_redo_controls_py.sto')
            controlslabels = controlsTable.getColumnLabels()
            controlsTableTrim = osim.TimeSeriesTable('muscletrack_redo_controls_py.sto')
            
            # get a version of the controls that is trimmed down. 
            for con in controlslabels:
                if 'lumbar' not in con:
                    controlsTableTrim.removeColumn(con)
            osim.STOFileAdapter.write(controlsTableTrim, 'trimmingControls_redo_' + tag + '.sto')
            
            # get a version of the model with no muscles in it (or reserves?)
            # muscmodel = osim.Model('post_simple_model_all_the_probes_muscletrack_redo.osim')
            trimmodelprocessor = osim.ModelProcessor('post_simple_model_all_the_probes_muscletrack_redo.osim')
            trimmodelprocessor.append(osim.ModOpRemoveMuscles())
            trimmodel = trimmodelprocessor.process()
            trimmodel.initSystem()
            
            # now have to do it for the forceset too?
            trimforces = trimmodel.getForceSet()
            numForces = trimforces.getSize()  
            count = 0
            for f in range(numForces):
                fo = trimforces.get(f-count)
                # print(fo.getName())
                if 'lumbar' not in fo.getName() and 'HOBL' not in fo.getName():
                    getrid = True
                    trimforces.remove(f-count)
                    count += 1
            
            trimmodel.initSystem()
            trimmodel.printToXML('trimmingmodel_redo_' + tag + '.osim')
            
            initTime = np.array(statesTable.getIndependentColumn())[0]
            finalTime = np.array(statesTable.getIndependentColumn())[-1]
            
            # now try the positionMotion
            # pomostorage = osim.Storage('trimmingStates2_inter.sto')
            pomostorage = osim.Storage('trimmingStates_redo_' + tag + '.sto')
            pomotraj = osim.StatesTrajectory.createFromStatesStorage(trimmodel, pomostorage)
            pomo = osim.PositionMotion.createFromStatesTrajectory(trimmodel, pomotraj)
            
            trimmodel.addComponent(pomo)
            trimmodel.initSystem()
            trimmodel.printToXML('trimmingmodel2_redo_' + tag + '.osim')
            
            # call the analyze tool to actually do the analysis and get the values. 
            jrax, jray, jraz = computeKneeContactRedo(trimmodel, initTime, finalTime, trialdir, tag, whichleg, runtool)  
        elif 'none' in musclesWanted_split:
            statesStorage = osim.Storage('muscletrack_redo_states_py.sto')
            statesTable = osim.TimeSeriesTable('muscletrack_redo_states_py.sto')
            stateslabels = statesTable.getColumnLabels()
            statesTableTrim = osim.TimeSeriesTable('muscletrack_redo_states_py.sto')
            statesTableTrim2 = osim.TimeSeriesTable('muscletrack_redo_states_py.sto')
            
            # get a trimmed set of states that is just the joint angles speeds, and whatever muscles you want. 
            for stat in stateslabels: 
                if 'forceset' in stat:
                    statesTableTrim.removeColumn(stat)
            osim.STOFileAdapter.write(statesTableTrim, 'trimmingStates_redo_' + tag + '.sto')

            for stat in stateslabels:
                if 'forceset' in stat or 'speed' in stat: 
                    statesTableTrim2.removeColumn(stat)
            osim.STOFileAdapter.write(statesTableTrim2, 'trimmingStates2_redo_' + tag + '.sto')

            # get a version of the controls that matches. 
            controlsTable = osim.TimeSeriesTable('muscletrack_redo_controls_py.sto')
            controlslabels = controlsTable.getColumnLabels()
            controlsTableTrim = osim.TimeSeriesTable('muscletrack_redo_controls_py.sto')
            
            # get a version of the controls that is trimmed down. 
            for con in controlslabels:
                controlsTableTrim.removeColumn(con)
            osim.STOFileAdapter.write(controlsTableTrim, 'trimmingControls_redo_' + tag + '.sto')
            
            # get a version of the model with no muscles in it (or reserves?)
            # muscmodel = osim.Model('post_simple_model_all_the_probes_muscletrack_redo.osim')
            trimmodelprocessor = osim.ModelProcessor('post_simple_model_all_the_probes_muscletrack_redo.osim')
            trimmodelprocessor.append(osim.ModOpRemoveMuscles())
            trimmodel = trimmodelprocessor.process()
            trimmodel.initSystem()
            
            # now have to do it for the forceset too?
            trimforces = trimmodel.getForceSet()
            numForces = trimforces.getSize()
            count = 0
            for f in range(numForces):
                fo = trimforces.get(f-count)
                if 'HOBL' not in fo.getName():
                    getrid = True
                    trimforces.remove(f-count)
                    count += 1
            
            trimmodel.initSystem()
            trimmodel.printToXML('trimmingmodel_redo_' + tag + '.osim')
            
            initTime = np.array(statesTable.getIndependentColumn())[0]
            finalTime = np.array(statesTable.getIndependentColumn())[-1]
            
            # now try the positionMotion
            # pomostorage = osim.Storage('trimmingStates2_inter.sto')
            pomostorage = osim.Storage('trimmingStates_redo_' + tag + '.sto')
            pomotraj = osim.StatesTrajectory.createFromStatesStorage(trimmodel, pomostorage)
            pomo = osim.PositionMotion.createFromStatesTrajectory(trimmodel, pomotraj)
            
            trimmodel.addComponent(pomo)
            trimmodel.initSystem()
            trimmodel.printToXML('trimmingmodel2_redo_' + tag + '.osim')
            
            # call the analyze tool to actually do the analysis and get the values. 
            jrax, jray, jraz = computeKneeContactRedo(trimmodel, initTime, finalTime, trialdir, tag, whichleg, runtool)
        else:
            statesStorage = osim.Storage('muscletrack_redo_states_py.sto')
            statesTable = osim.TimeSeriesTable('muscletrack_redo_states_py.sto')
            stateslabels = statesTable.getColumnLabels()
            statesTableTrim = osim.TimeSeriesTable('muscletrack_redo_states_py.sto')
            statesTableTrim2 = osim.TimeSeriesTable('muscletrack_redo_states_py.sto')
            
            # musclesWanted_split = ['1'2'3'4']
            # get a trimmed set of states that is just the joint angles speeds, and whatever muscles you want. 
            for stat in stateslabels:
                if 'forceset' in stat:
                    # print(stat)
                    getrid = True
                    for want in musclesWanted_split:
                        if want in stat:
                            # want this one
                            # print('want this one')
                            # print(stat)
                            getrid = False
                    if getrid:
                        statesTableTrim.removeColumn(stat)
            osim.STOFileAdapter.write(statesTableTrim, 'trimmingStates_redo_' + tag + '.sto')

            for stat in stateslabels:
                if 'forceset' in stat or 'speed' in stat:
                    getrid = True
                    for want in musclesWanted_split:
                        if want in stat:
                            # want to keep this one
                            # print(stat)
                            getrid = False
                    if getrid:
                        statesTableTrim2.removeColumn(stat)
            osim.STOFileAdapter.write(statesTableTrim2, 'trimmingStates2_redo_' + tag + '.sto')


            # get a version of the controls that matches. 
            controlsTable = osim.TimeSeriesTable('muscletrack_redo_controls_py.sto')
            controlslabels = controlsTable.getColumnLabels()
            controlsTableTrim = osim.TimeSeriesTable('muscletrack_redo_controls_py.sto')
            
            # get a version of the controls that is trimmed down. 
            for con in controlslabels:
                # print(con)
                if'lumbar' not in con:
                    getrid = True
                    for want in musclesWanted_split:
                        if want in con:
                            # want this
                            # print('want this one')
                            # print(con)
                            getrid = False
                    if getrid:
                        controlsTableTrim.removeColumn(con)
            osim.STOFileAdapter.write(controlsTableTrim, 'trimmingControls_redo_' + tag + '.sto')
            
            # get a version of the model with no muscles in it
            # muscmodel = osim.Model('post_simple_model_all_the_probes_muscletrack_redo.osim')
            trimmodel = osim.Model('post_simple_model_all_the_probes_muscletrack_redo.osim')
            # get rid of muscles that we don't want
            trimmuscles = trimmodel.getMuscles()
            numMuscles = trimmuscles.getSize()
            count = 0
            for m in range(numMuscles):
                musc = trimmuscles.get(m-count)
                # print(musc.getName())
                getrid = True
                for mu in musclesWanted_split:
                    # print(mu)
                    if mu == musc.getName():
                        # print(musc.getName())
                        # print(mu)
                        getrid = False
                if getrid:
                    trimmuscles.remove(musc)
                    count +=1
                    
            # now have to do it for the forceset too?
            trimforces = trimmodel.getForceSet()
            numForces = trimforces.getSize()  
            count = 0
            for f in range(numForces):
                fo = trimforces.get(f-count)
                # print(fo.getName())
                if 'lumbar' not in fo.getName() and 'HOBL' not in fo.getName():
                    getrid = True
                    for mu in musclesWanted_split:
                        if mu == fo.getName():
                            getrid = False
                    if getrid:
                        trimforces.remove(f-count)
                        count += 1
            
            # model should only have muscles that we want now. 
            trimmodel.initSystem()
            trimmodel.printToXML('trimmingmodel_redo_' + tag + '.osim')
            
            
            initTime = np.array(statesTable.getIndependentColumn())[0]
            finalTime = np.array(statesTable.getIndependentColumn())[-1]
            
            # now try the positionMotion
            pomostorage = osim.Storage('trimmingStates_redo_' + tag + '.sto')
            pomotraj = osim.StatesTrajectory.createFromStatesStorage(trimmodel, pomostorage)
            pomo = osim.PositionMotion.createFromStatesTrajectory(trimmodel, pomotraj)
            # pdb.set_trace()
            trimmodel.addComponent(pomo)
            trimmodel.initSystem()
            trimmodel.printToXML('trimmingmodel2_redo_' + tag + '.osim')
            jrax, jray, jraz = computeKneeContactRedo(trimmodel, initTime, finalTime, trialdir, tag, whichleg, runtool)

    return jrax, jray, jraz
# method for computing the individual muscle contributions to knee contact force
def getKneeContactributionsRedoPoly(trialdir, musclesWanted_split, tag, whichleg, runtool):
    if not runtool:
        print('\n NOT RUNNING TOOL - JUST COLLECTING OLD RESULTS')
        if musclesWanted_split == []:
            # load the model
            trimmodel = osim.Model('trimmingmodel2_redo_poly_' + tag + '.osim')
            # timings - these don't matter since we are not rerunning, just need filler for the function
            initTime = 0.0
            finalTime = 1.0
            jray = computeKneeContactRedo(trimmodel, initTime, finalTime, trialdir, tag, whichleg, runtool)
        elif 'all' in musclesWanted_split:
            # load the model
            trimmodel = osim.Model('trimmingmodel2_redo_poly_' + tag + '.osim')
            # timings - these don't matter since we are not rerunning, just need filler for the function
            initTime = 0.0
            finalTime = 1.0
            jray = computeKneeContactRedo(trimmodel, initTime, finalTime, trialdir, tag, whichleg, runtool)
        elif 'reserve' in musclesWanted_split:
            # load the model
            trimmodel = osim.Model('trimmingmodel2_redo_poly_' + tag + '.osim')
            # timings - these don't matter since we are not rerunning, just need filler for the function
            initTime = 0.0
            finalTime = 1.0
            jray = computeKneeContactRedo(trimmodel, initTime, finalTime, trialdir, tag, whichleg, runtool)
        elif 'none' in musclesWanted_split:
            # load the model
            trimmodel = osim.Model('trimmingmodel2_redo_poly_' + tag + '.osim')
            # timings - these don't matter since we are not rerunning, just need filler for the function
            initTime = 0.0
            finalTime = 1.0
            jray = computeKneeContactRedo(trimmodel, initTime, finalTime, trialdir, tag, whichleg, runtool)
        else:
            # load the model
            trimmodel = osim.Model('trimmingmodel2_redo_poly_' + tag + '.osim')
            # timings - these don't matter since we are not rerunning, just need filler for the function
            initTime = 0.0
            finalTime = 1.0
            jray = computeKneeContactRedo(trimmodel, initTime, finalTime, trialdir, tag, whichleg, runtool)

    else:
        print('\n RUNNING TOOL')

        # os.chdir(trialdir)
        # solution = osim.MocoTrajectory('muscle_statetrack_grfprescribe_solution_100con.sto')
        # model = osim.Model('post_simple_model_all_the_probes_muscletrack.osim')
        # # attempting with just setting muscles to not apply force
        # muscles = model.getMuscles()
        # for m in range(muscles.getSize()):
        #     musc = muscles.get(m)
        #     muscname = musc.getName()
        #     if muscname not in musclesWanted_split:
        #         musc.setMaxIsometricForce(0.0)
        #     # else:
        #     #     print(muscname)
        
        # model.initSystem()
        if musclesWanted_split == []: # this is the inter or intersegmental condition
            statesStorage = osim.Storage('muscletrack_redo_states_poly_py.sto')
            statesTable = osim.TimeSeriesTable('muscletrack_redo_states_poly_py.sto')
            stateslabels = statesTable.getColumnLabels()
            statesTableTrim = osim.TimeSeriesTable('muscletrack_redo_states_poly_py.sto')
            statesTableTrim2 = osim.TimeSeriesTable('muscletrack_redo_states_poly_py.sto')
            
            # get a trimmed set of states that is just the joint angles speeds, and whatever muscles you want. 
            for stat in stateslabels: 
                if 'forceset' in stat:
                    statesTableTrim.removeColumn(stat)
            osim.STOFileAdapter.write(statesTableTrim, 'trimmingStates_redo_poly_' + tag + '.sto')

            for stat in stateslabels:
                if 'forceset' in stat or 'speed' in stat: 
                    statesTableTrim2.removeColumn(stat)
            osim.STOFileAdapter.write(statesTableTrim2, 'trimmingStates2_redo_poly_' + tag + '.sto')

            # get a version of the controls that matches. 
            controlsTable = osim.TimeSeriesTable('muscletrack_redo_controls_poly_py.sto')
            controlslabels = controlsTable.getColumnLabels()
            controlsTableTrim = osim.TimeSeriesTable('muscletrack_redo_controls_poly_py.sto')
            
            # get a version of the controls that is trimmed down. 
            for con in controlslabels:
                if 'reserve' not in con and 'lumbar' not in con:
                    controlsTableTrim.removeColumn(con)
            osim.STOFileAdapter.write(controlsTableTrim, 'trimmingControls_redo_poly_' + tag + '.sto')
            
            # get a version of the model with no muscles in it
            # muscmodel = osim.Model('post_simple_model_all_the_probes_muscletrack_redo.osim')
            trimmodelprocessor = osim.ModelProcessor('post_simple_model_all_the_probes_muscletrack_redo.osim')
            trimmodelprocessor.append(osim.ModOpRemoveMuscles())
            trimmodel = trimmodelprocessor.process()
            trimmodel.initSystem()
            trimmodel.printToXML('trimmingmodel_redo_poly_' + tag + '.osim')
            
            initTime = np.array(statesTable.getIndependentColumn())[0]
            finalTime = np.array(statesTable.getIndependentColumn())[-1]
            
            # now try the positionMotion
            # pomostorage = osim.Storage('trimmingStates2_inter.sto')
            pomostorage = osim.Storage('trimmingStates_redo_poly_' + tag + '.sto')
            pomotraj = osim.StatesTrajectory.createFromStatesStorage(trimmodel, pomostorage)
            pomo = osim.PositionMotion.createFromStatesTrajectory(trimmodel, pomotraj)
            
            trimmodel.addComponent(pomo)
            trimmodel.initSystem()
            trimmodel.printToXML('trimmingmodel2_redo_poly_' + tag + '.osim')
            
            # call the analyze tool to actually do the analysis and get the values. 
            jray = computeKneeContactRedo(trimmodel, initTime, finalTime, trialdir, tag, whichleg, runtool)
        elif 'all' in musclesWanted_split:
            statesStorage = osim.Storage('muscletrack_redo_states_poly_py.sto')
            statesTable = osim.TimeSeriesTable('muscletrack_redo_states_poly_py.sto')
            stateslabels = statesTable.getColumnLabels()
            statesTableTrim = osim.TimeSeriesTable('muscletrack_redo_states_poly_py.sto')
            statesTableTrim2 = osim.TimeSeriesTable('muscletrack_redo_states_poly_py.sto')
            # in this case, we want all the muscles in the model
            osim.STOFileAdapter.write(statesTableTrim, 'trimmingStates_redo_poly_' + tag + '.sto')
            
            for stat in stateslabels:
                if 'speed' in stat:
                    statesTableTrim2.removeColumn(stat)
            osim.STOFileAdapter.write(statesTableTrim2, 'trimmingStates2_redo_poly_' + tag + '.sto')


            # get a version of the controls that matches. 
            controlsTable = osim.TimeSeriesTable('muscletrack_redo_controls_poly_py.sto')
            controlslabels = controlsTable.getColumnLabels()
            controlsTableTrim = osim.TimeSeriesTable('muscletrack_redo_controls_poly_py.sto')
            
            # in this case we want all the controls, not getting rid of any muscles
            osim.STOFileAdapter.write(controlsTableTrim, 'trimmingControls_redo_poly_' + tag + '.sto')
            
            
            # get a version of the model with no muscles in it
            # muscmodel = osim.Model('post_simple_model_all_the_probes_muscletrack_redo.osim')
            trimmodel = osim.Model('post_simple_model_all_the_probes_muscletrack_redo.osim')
            
            # model keeping all the muscles again for this one
            trimmodel.initSystem()
            trimmodel.printToXML('trimmingmodel_redo_poly_' + tag + '.osim')
            
            initTime = np.array(statesTable.getIndependentColumn())[0]
            finalTime = np.array(statesTable.getIndependentColumn())[-1]
            
            # now try the positionMotion
            pomostorage = osim.Storage('trimmingStates_redo_poly_' + tag + '.sto')
            pomotraj = osim.StatesTrajectory.createFromStatesStorage(trimmodel, pomostorage)
            pomo = osim.PositionMotion.createFromStatesTrajectory(trimmodel, pomotraj)
            
            trimmodel.addComponent(pomo)
            trimmodel.initSystem()
            trimmodel.printToXML('trimmingmodel2_redo_poly_' + tag + '.osim')
        
            ###
            jray = computeKneeContactRedo(trimmodel, initTime, finalTime, trialdir, tag, whichleg, runtool)
        elif 'reserve' in musclesWanted_split:
            statesStorage = osim.Storage('muscletrack_redo_states_poly_py.sto')
            statesTable = osim.TimeSeriesTable('muscletrack_redo_states_poly_py.sto')
            stateslabels = statesTable.getColumnLabels()
            statesTableTrim = osim.TimeSeriesTable('muscletrack_redo_states_poly_py.sto')
            statesTableTrim2 = osim.TimeSeriesTable('muscletrack_redo_states_poly_py.sto')
            
            # get a trimmed set of states that is just the joint angles speeds, and whatever muscles you want. 
            for stat in stateslabels: 
                if 'forceset' in stat:
                    statesTableTrim.removeColumn(stat)
            osim.STOFileAdapter.write(statesTableTrim, 'trimmingStates_redo_poly_' + tag + '.sto')

            for stat in stateslabels:
                if 'forceset' in stat or 'speed' in stat: 
                    statesTableTrim2.removeColumn(stat)
            osim.STOFileAdapter.write(statesTableTrim2, 'trimmingStates2_redo_poly_' + tag + '.sto')

            # get a version of the controls that matches. 
            controlsTable = osim.TimeSeriesTable('muscletrack_redo_controls_poly_py.sto')
            controlslabels = controlsTable.getColumnLabels()
            controlsTableTrim = osim.TimeSeriesTable('muscletrack_redo_controls_poly_py.sto')
            
            # get a version of the controls that is trimmed down. 
            for con in controlslabels:
                if 'lumbar' not in con:
                    controlsTableTrim.removeColumn(con)
            osim.STOFileAdapter.write(controlsTableTrim, 'trimmingControls_redo_poly_' + tag + '.sto')
            
            # get a version of the model with no muscles in it (or reserves?)
            # muscmodel = osim.Model('post_simple_model_all_the_probes_muscletrack_redo.osim')
            trimmodelprocessor = osim.ModelProcessor('post_simple_model_all_the_probes_muscletrack_redo.osim')
            trimmodelprocessor.append(osim.ModOpRemoveMuscles())
            trimmodel = trimmodelprocessor.process()
            trimmodel.initSystem()
            
            # now have to do it for the forceset too?
            trimforces = trimmodel.getForceSet()
            numForces = trimforces.getSize()  
            count = 0
            for f in range(numForces):
                fo = trimforces.get(f-count)
                # print(fo.getName())
                if 'lumbar' not in fo.getName() and 'HOBL' not in fo.getName():
                    getrid = True
                    trimforces.remove(f-count)
                    count += 1
            
            trimmodel.initSystem()
            trimmodel.printToXML('trimmingmodel_redo_poly_' + tag + '.osim')
            
            initTime = np.array(statesTable.getIndependentColumn())[0]
            finalTime = np.array(statesTable.getIndependentColumn())[-1]
            
            # now try the positionMotion
            # pomostorage = osim.Storage('trimmingStates2_inter.sto')
            pomostorage = osim.Storage('trimmingStates_redo_poly_' + tag + '.sto')
            pomotraj = osim.StatesTrajectory.createFromStatesStorage(trimmodel, pomostorage)
            pomo = osim.PositionMotion.createFromStatesTrajectory(trimmodel, pomotraj)
            
            trimmodel.addComponent(pomo)
            trimmodel.initSystem()
            trimmodel.printToXML('trimmingmodel2_redo_poly_' + tag + '.osim')
            
            # call the analyze tool to actually do the analysis and get the values. 
            jray = computeKneeContactRedo(trimmodel, initTime, finalTime, trialdir, tag, whichleg, runtool)  
        elif 'none' in musclesWanted_split:
            statesStorage = osim.Storage('muscletrack_redo_states_poly_py.sto')
            statesTable = osim.TimeSeriesTable('muscletrack_redo_states_poly_py.sto')
            stateslabels = statesTable.getColumnLabels()
            statesTableTrim = osim.TimeSeriesTable('muscletrack_redo_states_poly_py.sto')
            statesTableTrim2 = osim.TimeSeriesTable('muscletrack_redo_states_poly_py.sto')
            
            # get a trimmed set of states that is just the joint angles speeds, and whatever muscles you want. 
            for stat in stateslabels: 
                if 'forceset' in stat:
                    statesTableTrim.removeColumn(stat)
            osim.STOFileAdapter.write(statesTableTrim, 'trimmingStates_redo_poly_' + tag + '.sto')

            for stat in stateslabels:
                if 'forceset' in stat or 'speed' in stat: 
                    statesTableTrim2.removeColumn(stat)
            osim.STOFileAdapter.write(statesTableTrim2, 'trimmingStates2_redo_poly_' + tag + '.sto')

            # get a version of the controls that matches. 
            controlsTable = osim.TimeSeriesTable('muscletrack_redo_controls_poly_py.sto')
            controlslabels = controlsTable.getColumnLabels()
            controlsTableTrim = osim.TimeSeriesTable('muscletrack_redo_controls_poly_py.sto')
            
            # get a version of the controls that is trimmed down. 
            for con in controlslabels:
                controlsTableTrim.removeColumn(con)
            osim.STOFileAdapter.write(controlsTableTrim, 'trimmingControls_redo_poly_' + tag + '.sto')
            
            # get a version of the model with no muscles in it (or reserves?)
            # muscmodel = osim.Model('post_simple_model_all_the_probes_muscletrack_redo.osim')
            trimmodelprocessor = osim.ModelProcessor('post_simple_model_all_the_probes_muscletrack_redo.osim')
            trimmodelprocessor.append(osim.ModOpRemoveMuscles())
            trimmodel = trimmodelprocessor.process()
            trimmodel.initSystem()
            
            # now have to do it for the forceset too?
            trimforces = trimmodel.getForceSet()
            numForces = trimforces.getSize()
            count = 0
            for f in range(numForces):
                fo = trimforces.get(f-count)
                if 'HOBL' not in fo.getName():
                    getrid = True
                    trimforces.remove(f-count)
                    count += 1
            
            trimmodel.initSystem()
            trimmodel.printToXML('trimmingmodel_redo_poly_' + tag + '.osim')
            
            initTime = np.array(statesTable.getIndependentColumn())[0]
            finalTime = np.array(statesTable.getIndependentColumn())[-1]
            
            # now try the positionMotion
            # pomostorage = osim.Storage('trimmingStates2_inter.sto')
            pomostorage = osim.Storage('trimmingStates_redo_poly_' + tag + '.sto')
            pomotraj = osim.StatesTrajectory.createFromStatesStorage(trimmodel, pomostorage)
            pomo = osim.PositionMotion.createFromStatesTrajectory(trimmodel, pomotraj)
            
            trimmodel.addComponent(pomo)
            trimmodel.initSystem()
            trimmodel.printToXML('trimmingmodel2_redo_poly_' + tag + '.osim')
            
            # call the analyze tool to actually do the analysis and get the values. 
            jray = computeKneeContactRedo(trimmodel, initTime, finalTime, trialdir, tag, whichleg, runtool)
        else:
            statesStorage = osim.Storage('muscletrack_redo_states_poly_py.sto')
            statesTable = osim.TimeSeriesTable('muscletrack_redo_states_poly_py.sto')
            stateslabels = statesTable.getColumnLabels()
            statesTableTrim = osim.TimeSeriesTable('muscletrack_redo_states_poly_py.sto')
            statesTableTrim2 = osim.TimeSeriesTable('muscletrack_redo_states_poly_py.sto')
            
            # musclesWanted_split = ['1'2'3'4']
            # get a trimmed set of states that is just the joint angles speeds, and whatever muscles you want. 
            for stat in stateslabels:
                if 'forceset' in stat:
                    # print(stat)
                    getrid = True
                    for want in musclesWanted_split:
                        if want in stat:
                            # want this one
                            # print('want this one')
                            # print(stat)
                            getrid = False
                    if getrid:
                        statesTableTrim.removeColumn(stat)
            osim.STOFileAdapter.write(statesTableTrim, 'trimmingStates_redo_poly_' + tag + '.sto')

            for stat in stateslabels:
                if 'forceset' in stat or 'speed' in stat:
                    getrid = True
                    for want in musclesWanted_split:
                        if want in stat:
                            # want to keep this one
                            # print(stat)
                            getrid = False
                    if getrid:
                        statesTableTrim2.removeColumn(stat)
            osim.STOFileAdapter.write(statesTableTrim2, 'trimmingStates2_redo_poly_' + tag + '.sto')


            # get a version of the controls that matches. 
            controlsTable = osim.TimeSeriesTable('muscletrack_redo_controls_poly_py.sto')
            controlslabels = controlsTable.getColumnLabels()
            controlsTableTrim = osim.TimeSeriesTable('muscletrack_redo_controls_poly_py.sto')
            
            # get a version of the controls that is trimmed down. 
            for con in controlslabels:
                # print(con)
                if'lumbar' not in con:
                    getrid = True
                    for want in musclesWanted_split:
                        if want in con:
                            # want this
                            # print('want this one')
                            # print(con)
                            getrid = False
                    if getrid:
                        controlsTableTrim.removeColumn(con)
            osim.STOFileAdapter.write(controlsTableTrim, 'trimmingControls_redo_poly_' + tag + '.sto')
            
            # get a version of the model with no muscles in it
            # muscmodel = osim.Model('post_simple_model_all_the_probes_muscletrack_redo.osim')
            trimmodel = osim.Model('post_simple_model_all_the_probes_muscletrack_redo.osim')
            # get rid of muscles that we don't want
            trimmuscles = trimmodel.getMuscles()
            numMuscles = trimmuscles.getSize()
            count = 0
            for m in range(numMuscles):
                musc = trimmuscles.get(m-count)
                # print(musc.getName())
                getrid = True
                for mu in musclesWanted_split:
                    # print(mu)
                    if mu == musc.getName():
                        # print(musc.getName())
                        # print(mu)
                        getrid = False
                if getrid:
                    trimmuscles.remove(musc)
                    count +=1
                    
            # now have to do it for the forceset too?
            trimforces = trimmodel.getForceSet()
            numForces = trimforces.getSize()  
            count = 0
            for f in range(numForces):
                fo = trimforces.get(f-count)
                # print(fo.getName())
                if 'lumbar' not in fo.getName() and 'HOBL' not in fo.getName():
                    getrid = True
                    for mu in musclesWanted_split:
                        if mu == fo.getName():
                            getrid = False
                    if getrid:
                        trimforces.remove(f-count)
                        count += 1
            
            # model should only have muscles that we want now. 
            trimmodel.initSystem()
            trimmodel.printToXML('trimmingmodel_redo_poly_' + tag + '.osim')
            
            
            initTime = np.array(statesTable.getIndependentColumn())[0]
            finalTime = np.array(statesTable.getIndependentColumn())[-1]
            
            # now try the positionMotion
            pomostorage = osim.Storage('trimmingStates_redo_poly_' + tag + '.sto')
            pomotraj = osim.StatesTrajectory.createFromStatesStorage(trimmodel, pomostorage)
            pomo = osim.PositionMotion.createFromStatesTrajectory(trimmodel, pomotraj)
            # pdb.set_trace()
            trimmodel.addComponent(pomo)
            trimmodel.initSystem()
            trimmodel.printToXML('trimmingmodel2_redo_poly_' + tag + '.osim')
            jray = computeKneeContactRedo(trimmodel, initTime, finalTime, trialdir, tag, whichleg, runtool)

    return jray
# method for computing the individual muscle contributions to knee contact force
def getKneeContactributionsPrescribe(trialdir, musclesWanted_split, tag):
    # os.chdir(trialdir)
    # solution = osim.MocoTrajectory('muscle_statetrack_grfprescribe_solution_100con.sto')
    # model = osim.Model('post_simple_model_all_the_probes_muscletrack.osim')
    # # attempting with just setting muscles to not apply force
    # muscles = model.getMuscles()
    # for m in range(muscles.getSize()):
    #     musc = muscles.get(m)
    #     muscname = musc.getName()
    #     if muscname not in musclesWanted_split:
    #         musc.setMaxIsometricForce(0.0)
    #     # else:
    #     #     print(muscname)
    
    # model.initSystem()
    if musclesWanted_split == []: # this is the inter or intersegmental condition
        statesStorage = osim.Storage('muscleprescribe_states_redoarms.sto')
        statesTable = osim.TimeSeriesTable('muscleprescribe_states_redoarms.sto')
        stateslabels = statesTable.getColumnLabels()
        statesTableTrim = osim.TimeSeriesTable('muscleprescribe_states_redoarms.sto')
        statesTableTrim2 = osim.TimeSeriesTable('muscleprescribe_states_redoarms.sto')
        
        # get a trimmed set of states that is just the joint angles speeds, and whatever muscles you want. 
        for stat in stateslabels: 
            if 'forceset' in stat:
                statesTableTrim.removeColumn(stat)
        osim.STOFileAdapter.write(statesTableTrim, 'trimmingStates_prescribe_' + tag + '.sto')

        for stat in stateslabels:
            if 'forceset' in stat or 'speed' in stat: 
                statesTableTrim2.removeColumn(stat)
        osim.STOFileAdapter.write(statesTableTrim2, 'trimmingStates2_prescribe_' + tag + '.sto')

        # get a version of the controls that matches. 
        controlsTable = osim.TimeSeriesTable('muscleprescribe_controls_redoarms.sto')
        controlslabels = controlsTable.getColumnLabels()
        controlsTableTrim = osim.TimeSeriesTable('muscleprescribe_controls_redoarms.sto')
        
        # get a version of the controls that is trimmed down. 
        for con in controlslabels:
            if 'reserve' not in con and 'lumbar' not in con:
                controlsTableTrim.removeColumn(con)
        osim.STOFileAdapter.write(controlsTableTrim, 'trimmingControls_prescribe_' + tag + '.sto')
        
        # get a version of the model with no muscles in it
        # muscmodel = osim.Model('post_simple_model_all_the_probes_muscletrack.osim')
        trimmodelprocessor = osim.ModelProcessor('post_simple_model_all_the_probes_muscleprescribe.osim')
        trimmodelprocessor.append(osim.ModOpRemoveMuscles())
        trimmodel = trimmodelprocessor.process()
        trimmodel.initSystem()
        trimmodel.printToXML('trimmingmodel_prescribe_' + tag + '.osim')
        
        initTime = np.array(statesTable.getIndependentColumn())[0]
        finalTime = np.array(statesTable.getIndependentColumn())[-1]
        
        # now try the positionMotion
        # pomostorage = osim.Storage('trimmingStates2_inter.sto')
        pomostorage = osim.Storage('trimmingStates_prescribe_' + tag + '.sto')
        pomotraj = osim.StatesTrajectory.createFromStatesStorage(trimmodel, pomostorage)
        pomo = osim.PositionMotion.createFromStatesTrajectory(trimmodel, pomotraj)
        
        trimmodel.addComponent(pomo)
        trimmodel.initSystem()
        trimmodel.printToXML('trimmingmodel2_prescribe_' + tag + '.osim')
        
        # call the analyze tool to actually do the analysis and get the values. 
        jray = computeKneeContactPrescribe(trimmodel, initTime, finalTime, trialdir, tag)
    
    
    elif 'all' in musclesWanted_split:
        statesStorage = osim.Storage('muscleprescribe_states_redoarms.sto')
        statesTable = osim.TimeSeriesTable('muscleprescribe_states_redoarms.sto')
        stateslabels = statesTable.getColumnLabels()
        statesTableTrim = osim.TimeSeriesTable('muscleprescribe_states_redoarms.sto')
        statesTableTrim2 = osim.TimeSeriesTable('muscleprescribe_states_redoarms.sto')
        # in this case, we want all the muscles in the model
        osim.STOFileAdapter.write(statesTableTrim, 'trimmingStates_prescribe_' + tag + '.sto')
        
        for stat in stateslabels:
            if 'speed' in stat:
                statesTableTrim2.removeColumn(stat)
        osim.STOFileAdapter.write(statesTableTrim2, 'trimmingStates2_prescribe_' + tag + '.sto')


        # get a version of the controls that matches. 
        controlsTable = osim.TimeSeriesTable('muscleprescribe_controls_redoarms.sto')
        controlslabels = controlsTable.getColumnLabels()
        controlsTableTrim = osim.TimeSeriesTable('muscleprescribe_controls_redoarms.sto')
        
        # in this case we want all the controls, not getting rid of any muscles
        osim.STOFileAdapter.write(controlsTableTrim, 'trimmingControls_prescribe_' + tag + '.sto')
        
        
        # get a version of the model with no muscles in it
        # muscmodel = osim.Model('post_simple_model_all_the_probes_muscletrack.osim')
        trimmodel = osim.Model('post_simple_model_all_the_probes_muscleprescribe.osim')
        
        # model keeping all the muscles again for this one
        trimmodel.initSystem()
        trimmodel.printToXML('trimmingmodel_prescribe_' + tag + '.osim')
        
        initTime = np.array(statesTable.getIndependentColumn())[0]
        finalTime = np.array(statesTable.getIndependentColumn())[-1]
        
        # now try the positionMotion
        pomostorage = osim.Storage('trimmingStates_prescribe_' + tag + '.sto')
        pomotraj = osim.StatesTrajectory.createFromStatesStorage(trimmodel, pomostorage)
        pomo = osim.PositionMotion.createFromStatesTrajectory(trimmodel, pomotraj)
        
        trimmodel.addComponent(pomo)
        trimmodel.initSystem()
        trimmodel.printToXML('trimmingmodel2_prescribe_' + tag + '.osim')
    
        ###
        jray = computeKneeContactPrescribe(trimmodel, initTime, finalTime, trialdir, tag)
    
    elif 'reserve' in musclesWanted_split:
        statesStorage = osim.Storage('muscleprescribe_states_redoarms.sto')
        statesTable = osim.TimeSeriesTable('muscleprescribe_states_redoarms.sto')
        stateslabels = statesTable.getColumnLabels()
        statesTableTrim = osim.TimeSeriesTable('muscleprescribe_states_redoarms.sto')
        statesTableTrim2 = osim.TimeSeriesTable('muscleprescribe_states_redoarms.sto')
        
        # get a trimmed set of states that is just the joint angles speeds, and whatever muscles you want. 
        for stat in stateslabels: 
            if 'forceset' in stat:
                statesTableTrim.removeColumn(stat)
        osim.STOFileAdapter.write(statesTableTrim, 'trimmingStates_prescribe_' + tag + '.sto')

        for stat in stateslabels:
            if 'forceset' in stat or 'speed' in stat: 
                statesTableTrim2.removeColumn(stat)
        osim.STOFileAdapter.write(statesTableTrim2, 'trimmingStates2_prescribe_' + tag + '.sto')

        # get a version of the controls that matches. 
        controlsTable = osim.TimeSeriesTable('muscleprescribe_controls_redoarms.sto')
        controlslabels = controlsTable.getColumnLabels()
        controlsTableTrim = osim.TimeSeriesTable('muscleprescribe_controls_redoarms.sto')
        
        # get a version of the controls that is trimmed down. 
        for con in controlslabels:
            if 'lumbar' not in con:
                controlsTableTrim.removeColumn(con)
        osim.STOFileAdapter.write(controlsTableTrim, 'trimmingControls_prescribe_' + tag + '.sto')
        
        # get a version of the model with no muscles in it (or reserves?)
        # muscmodel = osim.Model('post_simple_model_all_the_probes_muscletrack.osim')
        trimmodelprocessor = osim.ModelProcessor('post_simple_model_all_the_probes_muscleprescribe.osim')
        trimmodelprocessor.append(osim.ModOpRemoveMuscles())
        trimmodel = trimmodelprocessor.process()
        trimmodel.initSystem()
        
        # now have to do it for the forceset too?
        trimforces = trimmodel.getForceSet()
        numForces = trimforces.getSize()  
        count = 0
        for f in range(numForces):
            fo = trimforces.get(f-count)
            # print(fo.getName())
            if 'lumbar' not in fo.getName() and 'HOBL' not in fo.getName():
                getrid = True
                trimforces.remove(f-count)
                count += 1
        
        trimmodel.initSystem()
        trimmodel.printToXML('trimmingmodel_prescribe_' + tag + '.osim')
        
        initTime = np.array(statesTable.getIndependentColumn())[0]
        finalTime = np.array(statesTable.getIndependentColumn())[-1]
        
        # now try the positionMotion
        # pomostorage = osim.Storage('trimmingStates2_inter.sto')
        pomostorage = osim.Storage('trimmingStates_prescribe_' + tag + '.sto')
        pomotraj = osim.StatesTrajectory.createFromStatesStorage(trimmodel, pomostorage)
        pomo = osim.PositionMotion.createFromStatesTrajectory(trimmodel, pomotraj)
        
        trimmodel.addComponent(pomo)
        trimmodel.initSystem()
        trimmodel.printToXML('trimmingmodel2_prescribe_' + tag + '.osim')
        
        # call the analyze tool to actually do the analysis and get the values. 
        jray = computeKneeContactPrescribe(trimmodel, initTime, finalTime, trialdir, tag)
    
    elif 'none' in musclesWanted_split:
        statesStorage = osim.Storage('muscleprescribe_states_redoarms.sto')
        statesTable = osim.TimeSeriesTable('muscleprescribe_states_redoarms.sto')
        stateslabels = statesTable.getColumnLabels()
        statesTableTrim = osim.TimeSeriesTable('muscleprescribe_states_redoarms.sto')
        statesTableTrim2 = osim.TimeSeriesTable('muscleprescribe_states_redoarms.sto')
        
        # get a trimmed set of states that is just the joint angles speeds, and whatever muscles you want. 
        for stat in stateslabels: 
            if 'forceset' in stat:
                statesTableTrim.removeColumn(stat)
        osim.STOFileAdapter.write(statesTableTrim, 'trimmingStates_prescribe_' + tag + '.sto')

        for stat in stateslabels:
            if 'forceset' in stat or 'speed' in stat: 
                statesTableTrim2.removeColumn(stat)
        osim.STOFileAdapter.write(statesTableTrim2, 'trimmingStates2_prescribe_' + tag + '.sto')

        # get a version of the controls that matches. 
        controlsTable = osim.TimeSeriesTable('muscleprescribe_controls_redoarms.sto')
        controlslabels = controlsTable.getColumnLabels()
        controlsTableTrim = osim.TimeSeriesTable('muscleprescribe_controls_redoarms.sto')
        
        # get a version of the controls that is trimmed down. 
        for con in controlslabels:
            controlsTableTrim.removeColumn(con)
        osim.STOFileAdapter.write(controlsTableTrim, 'trimmingControls_prescribe_' + tag + '.sto')
        
        # get a version of the model with no muscles in it (or reserves?)
        # muscmodel = osim.Model('post_simple_model_all_the_probes_muscletrack.osim')
        trimmodelprocessor = osim.ModelProcessor('post_simple_model_all_the_probes_muscleprescribe.osim')
        trimmodelprocessor.append(osim.ModOpRemoveMuscles())
        trimmodel = trimmodelprocessor.process()
        trimmodel.initSystem()
        
        # now have to do it for the forceset too?
        trimforces = trimmodel.getForceSet()
        numForces = trimforces.getSize()
        count = 0
        for f in range(numForces):
            fo = trimforces.get(f-count)
            if 'HOBL' not in fo.getName():
                getrid = True
                trimforces.remove(f-count)
                count += 1
        
        trimmodel.initSystem()
        trimmodel.printToXML('trimmingmodel_prescribe_' + tag + '.osim')
        
        initTime = np.array(statesTable.getIndependentColumn())[0]
        finalTime = np.array(statesTable.getIndependentColumn())[-1]
        
        # now try the positionMotion
        # pomostorage = osim.Storage('trimmingStates2_inter.sto')
        pomostorage = osim.Storage('trimmingStates_prescribe_' + tag + '.sto')
        pomotraj = osim.StatesTrajectory.createFromStatesStorage(trimmodel, pomostorage)
        pomo = osim.PositionMotion.createFromStatesTrajectory(trimmodel, pomotraj)
        
        trimmodel.addComponent(pomo)
        trimmodel.initSystem()
        trimmodel.printToXML('trimmingmodel2_prescribe_' + tag + '.osim')
        
        # call the analyze tool to actually do the analysis and get the values. 
        jray = computeKneeContactPrescribe(trimmodel, initTime, finalTime, trialdir, tag)
        
    
    else:
        statesStorage = osim.Storage('muscleprescribe_states_redoarms.sto')
        statesTable = osim.TimeSeriesTable('muscleprescribe_states_redoarms.sto')
        stateslabels = statesTable.getColumnLabels()
        statesTableTrim = osim.TimeSeriesTable('muscleprescribe_states_redoarms.sto')
        statesTableTrim2 = osim.TimeSeriesTable('muscleprescribe_states_redoarms.sto')
        
        # musclesWanted_split = ['1'2'3'4']
        # get a trimmed set of states that is just the joint angles speeds, and whatever muscles you want. 
        for stat in stateslabels:
            if 'forceset' in stat:
                # print(stat)
                getrid = True
                for want in musclesWanted_split:
                    if want in stat:
                        # want this one
                        # print('want this one')
                        # print(stat)
                        getrid = False
                if getrid:
                    statesTableTrim.removeColumn(stat)
        osim.STOFileAdapter.write(statesTableTrim, 'trimmingStates_prescribe_' + tag + '.sto')

        for stat in stateslabels:
            if 'forceset' in stat or 'speed' in stat:
                getrid = True
                for want in musclesWanted_split:
                    if want in stat:
                        # want to keep this one
                        # print(stat)
                        getrid = False
                if getrid:
                    statesTableTrim2.removeColumn(stat)
        osim.STOFileAdapter.write(statesTableTrim2, 'trimmingStates2_prescribe_' + tag + '.sto')


        # get a version of the controls that matches. 
        controlsTable = osim.TimeSeriesTable('muscleprescribe_controls_redoarms.sto')
        controlslabels = controlsTable.getColumnLabels()
        controlsTableTrim = osim.TimeSeriesTable('muscleprescribe_controls_redoarms.sto')
        
        # get a version of the controls that is trimmed down. 
        for con in controlslabels:
            # print(con)
            if 'lumbar' not in con:
                getrid = True
                for want in musclesWanted_split:
                    if want in con:
                        # want this
                        # print('want this one')
                        # print(con)
                        getrid = False
                if getrid:
                    controlsTableTrim.removeColumn(con)
        osim.STOFileAdapter.write(controlsTableTrim, 'trimmingControls_prescribe_' + tag + '.sto')
        
        # get a version of the model with no muscles in it
        # muscmodel = osim.Model('post_simple_model_all_the_probes_muscletrack.osim')
        trimmodel = osim.Model('post_simple_model_all_the_probes_muscleprescribe.osim')
        # get rid of muscles that we don't want
        trimmuscles = trimmodel.getMuscles()
        numMuscles = trimmuscles.getSize()
        count = 0
        for m in range(numMuscles):
            musc = trimmuscles.get(m-count)
            # print(musc.getName())
            getrid = True
            for mu in musclesWanted_split:
                # print(mu)
                if mu == musc.getName():
                    # print(musc.getName())
                    # print(mu)
                    getrid = False
            if getrid:
                trimmuscles.remove(musc)
                count +=1
        # now have to do it for the forceset too?
        trimforces = trimmodel.getForceSet()
        numForces = trimforces.getSize()  
        count = 0
        for f in range(numForces):
            fo = trimforces.get(f-count)
            if 'lumbar' not in fo.getName() and 'HOBL' not in fo.getName() and 'shoulder' not in fo.getName() and 'elbow' not in fo.getName() and 'pro_sup' not in fo.getName():
            # if 'HOBL' not in fo.getName():
                # print(fo.getName())
                getrid = True
                for mu in musclesWanted_split:
                    if mu == fo.getName():
                        getrid = False
                if getrid:
                    trimforces.remove(f-count)
                    count += 1
        
        # model should only have muscles that we want now. 
        trimmodel.initSystem()
        trimmodel.printToXML('trimmingmodel_prescribe_' + tag + '.osim')
        
        initTime = np.array(statesTable.getIndependentColumn())[0]
        finalTime = np.array(statesTable.getIndependentColumn())[-1]
        
        # now try the positionMotion
        pomostorage = osim.Storage('trimmingStates_prescribe_' + tag + '.sto')
        pomotraj = osim.StatesTrajectory.createFromStatesStorage(trimmodel, pomostorage)
        pomo = osim.PositionMotion.createFromStatesTrajectory(trimmodel, pomotraj)
        
        trimmodel.addComponent(pomo)
        trimmodel.initSystem()
        trimmodel.printToXML('trimmingmodel2_prescribe_' + tag + '.osim')
    
        ###
        jray = computeKneeContactPrescribe(trimmodel, initTime, finalTime, trialdir, tag)
    
    return jray
# plotting method that works for each component of the knee contact force
def plotKneeContactForce(tagcomponent, analyzedir, welksubjects, ncolor, ecolor, n_timespercent101, e_timespercent101, datacombined):
    # start by separating the data...
    ninterseg_combine = datacombined['ninterseg_combine']
    einterseg_combine = datacombined['einterseg_combine']
    ntfl_combine = datacombined['ntfl_combine']
    etfl_combine = datacombined['etfl_combine']
    ngas_combine = datacombined['ngas_combine']
    egas_combine = datacombined['egas_combine']
    nhams_combine = datacombined['nhams_combine']
    ehams_combine = datacombined['ehams_combine']
    nquads_combine = datacombined['nquads_combine']
    equads_combine = datacombined['equads_combine']
    nall_combine = datacombined['nall_combine']
    eall_combine = datacombined['eall_combine']


    ###########################################################################
    # plotting for the joint contacts - natural and exotendon
    if len(welksubjects) == 1:
        ncolors = ['#fee0b6', '#fdae6b', '#fd8d3c', '#e66101']
        ecolors = ['#c7eae5', '#80cdc1', '#35978f', '#01665e']
    

    ###########################################################################
    # figure: segmenting all the muscles between exo and nat 
    ## really nice figure for seeing what is going on, but likely not going to 
    ## be in the paper...
    fig9, ax9 = plt.subplots(2, 4, figsize=(14, 6))# , dpi=300)
    ax9 = ax9.flatten()
    # intersegmental forces average
    for i, curve in enumerate(ninterseg_combine):
         ax9[0].plot(n_timespercent101, curve, label='natural'+str(i), color=ncolors[i] if len(welksubjects) == 1 else ncolor)
    for i, curve in enumerate(einterseg_combine):
        ax9[0].plot(e_timespercent101, curve, label='exotendon'+str(i), color=ecolors[i] if len(welksubjects) == 1 else ecolor)
    ax9[0].set_xlabel('% Gait cycle')
    ax9[0].set_ylabel('Force (BW)')
    ax9[0].set_title('intersegmental')
    # ax9[0].legend()
    
    # tfl forces
    for i, curve in enumerate(ntfl_combine):
        ax9[1].plot(n_timespercent101, curve, label='natural'+str(i), color=ncolors[i] if len(welksubjects) == 1 else ncolor)
    for i, curve in enumerate(etfl_combine):
        ax9[1].plot(e_timespercent101, curve, label='exotendon'+str(i), color=ecolors[i] if len(welksubjects) == 1 else ecolor)
    # ax9[1].plot(n_timespercent101, ntfl_combine, label='natural', color=ncolor)
    # ax9[1].plot(e_timespercent101, etfl_combine, label='exotendon', color=ecolor)
    ax9[1].set_xlabel('% Gait cycle')
    # ax9[1].set_ylabel('Force (BW)')
    # ax9[1].legend()
    ax9[1].set_title('tfl')

    # gastroc forces
    for i, curve in enumerate(ngas_combine):
        ax9[2].plot(n_timespercent101, curve, label='natural'+str(i), color=ncolors[i] if len(welksubjects) == 1 else ncolor)
    for i, curve in enumerate(egas_combine):
        ax9[2].plot(e_timespercent101, curve, label='exotendon'+str(i), color=ecolors[i] if len(welksubjects) == 1 else ecolor)
    # ax9[2].plot(n_timespercent101, ngas_combine, label='natural', color=ncolor)
    # ax9[2].plot(e_timespercent101, egas_combine, label='exotendon', color=ecolor)
    ax9[2].set_xlabel('% Gait cycle')
    # ax9[2].set_ylabel('Force (BW)')
    # ax9[2].legend()
    ax9[2].set_title('gastroc')
    
    # hamstring forces
    for i, curve in enumerate(nhams_combine):
        ax9[3].plot(n_timespercent101, curve, label='natural'+str(i), color=ncolors[i] if len(welksubjects) == 1 else ncolor)
    for i, curve in enumerate(ehams_combine):
        ax9[3].plot(e_timespercent101, curve, label='exotendon'+str(i), color=ecolors[i] if len(welksubjects) == 1 else ecolor)
    # ax9[3].plot(n_timespercent101, nhams_combine, label='natural', color=ncolor)
    # ax9[3].plot(e_timespercent101, ehams_combine, label='exotendon', color=ecolor)
    ax9[3].set_xlabel('% Gait cycle')
    # ax9[3].set_ylabel('Force (BW)')
    # ax9[3].legend()
    ax9[3].set_title('hamstrings')

    # quads forces
    for i, curve in enumerate(nquads_combine):
        ax9[4].plot(n_timespercent101, curve, label='natural'+str(i), color=ncolors[i] if len(welksubjects) == 1 else ncolor)
    for i, curve in enumerate(equads_combine):
        ax9[4].plot(e_timespercent101, curve, label='exotendon'+str(i), color=ecolors[i] if len(welksubjects) == 1 else ecolor)
    # ax9[4].plot(n_timespercent101, nquads_combine, label='natural', color=ncolor)
    # ax9[4].plot(e_timespercent101, equads_combine, label='exotendon', color=ecolor)
    ax9[4].set_xlabel('% Gait cycle')
    # ax9[4].set_ylabel('Force (BW)')
    # ax9[4].legend()
    ax9[4].set_title('quadriceps')

    # # reserve forces
    # ax9[5].plot(n_timespercent101, nreserve_combine, label='natural', color=ncolor)
    # ax9[5].plot(e_timespercent101, ereserve_combine, label='exotendon', color=ecolor)
    # ax9[5].set_xlabel('% Gait cycle')
    # # ax9[5].set_ylabel('Force (BW)')
    # # ax9[5].legend()
    # ax9[5].set_title('reserves')
    
    # # added all forces
    # ax9[5].plot(n_timespercent101, nquads_combine+ nhams_combine+ ngas_combine+ ntfl_combine+ ninterseg_combine+ nreserve_combine, label='natural', color=ncolor)
    # ax9[5].plot(e_timespercent101, equads_combine+ ehams_combine+ egas_combine+ etfl_combine+ einterseg_combine+ ereserve_combine, label='exotendon', color=ecolor)
    # ax9[5].set_xlabel('% Gait cycle')
    # # ax9[6].set_ylabel('Force (BW)')
    # # ax9[6].legend()
    # ax9[5].set_title('Total Vertical Contact')

    # all forces from whole analysis
    for i, curve in enumerate(nall_combine):
        ax9[5].plot(n_timespercent101, curve, label='natural'+str(i), color=ncolors[i] if len(welksubjects) == 1 else ncolor)
    for i, curve in enumerate(eall_combine):
        ax9[5].plot(e_timespercent101, curve, label='exotendon'+str(i), color=ecolors[i] if len(welksubjects) == 1 else ecolor)
    # ax9[5].plot(n_timespercent101, nall_combine, label='natural', color=ncolor)
    # ax9[5].plot(e_timespercent101, eall_combine, label='exotendon', color=ecolor)
    ax9[5].set_xlabel('% Gait cycle')
    # ax9[5].set_ylabel('Force (BW)')
    ax9[5].set_title('Total ' + tagcomponent + ' contact')
    # ax9[5].legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10)
    # Hide the last subplot and use it to display the legend   
    ax9[6].axis('off')
    handles, labels = ax9[5].get_legend_handles_labels()
    if len(welksubjects) == 1:
        ax9[6].legend(handles, labels, loc='center', fontsize=14)

    fig9.tight_layout()
    plt.savefig(analyzedir + '\\contact1_' + whichleg + '_' + tagcomponent + '.png')

    ###########################################################################
    # figure: segmenting all the muscles between exo and nat 
    ## really nice figure for seeing what is going on, but likely not going to 
    ## be in the paper...
    fig10, ax10 = plt.subplots(1,7, figsize=(18,3), dpi=500)
    fontz = 16
    font_properties = {'fontsize': 16, 'fontfamily': 'serif', 'fontname': 'Times New Roman'}
    # tick_font_properties = {'fontfamily': 'serif', 'fontname': 'Times New Roman'}

    # all forces from whole analysis
    ax10[0].fill_between(n_timespercent101, np.mean(nall_combine, 0) - np.std(nall_combine, 0), np.mean(nall_combine, 0) + np.std(nall_combine, 0), color=ncolor, alpha=0.2)
    ax10[0].fill_between(e_timespercent101, np.mean(eall_combine, 0) - np.std(eall_combine, 0), np.mean(eall_combine, 0) + np.std(eall_combine, 0), color=ecolor, alpha=0.2)
    ax10[0].plot(n_timespercent101, np.mean(nall_combine,0), label='natural', color=ncolor)
    ax10[0].plot(e_timespercent101, np.mean(eall_combine,0), label='exotendon', color=ecolor)
    ax10[0].set_xlabel('% Gait cycle', **font_properties)
    ax10[0].set_ylabel('Force (BW)', **font_properties)
    # ax10[0].legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10)
    ax10[0].set_title('Total\n' + tagcomponent + ' Knee Contact', **font_properties)
    ax10[0].tick_params(axis='both', which='major', labelsize=fontz)

    # quads forces
    ax10[1].fill_between(n_timespercent101, np.mean(nquads_combine, 0) - np.std(nquads_combine, 0), np.mean(nquads_combine, 0) + np.std(nquads_combine, 0), color=ncolor, alpha=0.2)#, label='Natural St. Dev.')
    ax10[1].fill_between(e_timespercent101, np.mean(equads_combine, 0) - np.std(equads_combine, 0), np.mean(equads_combine, 0) + np.std(equads_combine, 0), color=ecolor, alpha=0.2)#, label='Exotendon St. Dev.')
    ax10[1].plot(n_timespercent101, np.mean(nquads_combine, 0), label='Natural (Mean \u00B1 Std.)', color=ncolor)
    ax10[1].plot(e_timespercent101, np.mean(equads_combine, 0), label='Exotendon (Mean \u00B1 Std.)', color=ecolor)
    ax10[1].set_xlabel('% Gait cycle', **font_properties)
    ax10[1].tick_params(axis='both', which='major', labelsize=fontz)
    # ax10[1].set_ylabel('Force (BW)')
    # ax10[1].legend()
    ax10[1].set_title('Contribution of\nQuadriceps', **font_properties)

    # gastroc forces
    ax10[2].fill_between(n_timespercent101, np.mean(ngas_combine, 0) - np.std(ngas_combine, 0), np.mean(ngas_combine, 0) + np.std(ngas_combine, 0), color=ncolor, alpha=0.2)
    ax10[2].fill_between(e_timespercent101, np.mean(egas_combine, 0) - np.std(egas_combine, 0), np.mean(egas_combine, 0) + np.std(egas_combine, 0), color=ecolor, alpha=0.2)
    ax10[2].plot(n_timespercent101, np.mean(ngas_combine, 0), label='natural', color=ncolor)
    ax10[2].plot(e_timespercent101, np.mean(egas_combine, 0), label='exotendon', color=ecolor)
    ax10[2].set_xlabel('% Gait cycle', **font_properties)
    ax10[2].tick_params(axis='both', which='major', labelsize=fontz)
    # ax10[2].set_ylabel('Force (BW)')
    # ax10[2].legend()
    ax10[2].set_title('Contribution of\nGastrocnemius', **font_properties)

    # hamstring forces
    ax10[3].fill_between(n_timespercent101, np.mean(nhams_combine, 0) - np.std(nhams_combine, 0), np.mean(nhams_combine, 0) + np.std(nhams_combine, 0), color=ncolor, alpha=0.2)
    ax10[3].fill_between(e_timespercent101, np.mean(ehams_combine, 0) - np.std(ehams_combine, 0), np.mean(ehams_combine, 0) + np.std(ehams_combine, 0), color=ecolor, alpha=0.2)
    ax10[3].plot(n_timespercent101, np.mean(nhams_combine, 0), label='natural', color=ncolor)
    ax10[3].plot(e_timespercent101, np.mean(ehams_combine, 0), label='exotendon', color=ecolor)
    ax10[3].set_xlabel('% Gait cycle', **font_properties)
    ax10[3].tick_params(axis='both', which='major', labelsize=fontz)
    # ax10[3].set_ylabel('Force (BW)')
    # ax10[3].legend()
    ax10[3].set_title('Contribution of\nHamstrings', **font_properties)

    # tfl forces
    ax10[4].fill_between(n_timespercent101, np.mean(ntfl_combine, 0) - np.std(ntfl_combine, 0), np.mean(ntfl_combine, 0) + np.std(ntfl_combine, 0), color=ncolor, alpha=0.2)
    ax10[4].fill_between(e_timespercent101, np.mean(etfl_combine, 0) - np.std(etfl_combine, 0), np.mean(etfl_combine, 0) + np.std(etfl_combine, 0), color=ecolor, alpha=0.2)
    ax10[4].plot(n_timespercent101, np.mean(ntfl_combine, 0), label='natural', color=ncolor)
    ax10[4].plot(e_timespercent101, np.mean(etfl_combine, 0), label='exotendon', color=ecolor)
    ax10[4].set_xlabel('% Gait cycle', **font_properties)
    ax10[4].tick_params(axis='both', which='major', labelsize=fontz)
    # ax10[4].set_ylabel('Force (BW)')
    # ax10[4].legend()
    ax10[4].set_title('Contribution of\nTensor Fascia Latae', **font_properties)

    # intersegmental forces average
    ax10[5].fill_between(n_timespercent101, np.mean(ninterseg_combine, 0) - np.std(ninterseg_combine, 0), np.mean(ninterseg_combine, 0) + np.std(ninterseg_combine, 0), color=ncolor, alpha=0.2)
    ax10[5].fill_between(e_timespercent101, np.mean(einterseg_combine, 0) - np.std(einterseg_combine, 0), np.mean(einterseg_combine, 0) + np.std(einterseg_combine, 0), color=ecolor, alpha=0.2)
    ax10[5].plot(n_timespercent101, np.mean(ninterseg_combine, 0), label='natural', color=ncolor)
    ax10[5].plot(e_timespercent101, np.mean(einterseg_combine, 0), label='exotendon', color=ecolor)
    ax10[5].set_xlabel('% Gait cycle', **font_properties)
    # ax10[5].set_ylabel('Force (BW)', **font_properties)
    ax10[5].set_title('Contribution of\nIntersegmental Forces', **font_properties)
    ax10[5].tick_params(axis='both', which='major', labelsize=fontz)    

    # # reserve forces
    # ax10[5].plot(n_timespercent101, np.mean(nreserve_combine, 0), label='natural', color=ncolor)
    # ax10[5].plot(e_timespercent101, np.mean(ereserve_combine, 0), label='exotendon', color=ecolor)
    # ax10[5].set_xlabel('% Gait cycle', **font_properties)
    # # ax10[5].set_ylabel('Force (BW)')
    # # ax10[5].legend()
    # ax10[5].set_title('reserves')

    # # added all forces
    # ax10[5].plot(n_timespercent101, np.mean(nquads_combine,0) + np.mean(nhams_combine,0) + np.mean(ngas_combine,0) + np.mean(ntfl_combine,0) + np.mean(ninterseg_combine,0) + np.mean(nreserve_combine,0), label='natural', color=ncolor)
    # ax10[5].plot(e_timespercent101, np.mean(equads_combine,0) + np.mean(ehams_combine,0) + np.mean(egas_combine,0) + np.mean(etfl_combine,0) + np.mean(einterseg_combine,0) + np.mean(ereserve_combine,0), label='exotendon', color=ecolor)
    # ax10[5].set_xlabel('% Gait cycle', **font_properties)
    # # ax10[6].set_ylabel('Force (BW)')
    # # ax10[6].legend()
    # ax10[5].set_title('Total Vertical Contact')

    # Hide the last subplot and use it to display the legend
    ax10[6].axis('off')
    handles, labels = ax10[1].get_legend_handles_labels()
    # ax10[6].legend(handles, labels, loc='center', fontsize=fontz)
    fig10.tight_layout(pad=2.0, w_pad=0.5, h_pad=1.0)
    plt.savefig(analyzedir + '\\contact2_' + whichleg + '_' + tagcomponent + '.png')

    ###########################################################################
    ### figure: differences between conditions for each and all
    #### this is a simplified look at the same as above. We can see nice stuff, 
    #### but likely not going to be in the paper. 
    fig11, ax11 = plt.subplots(1,7, figsize=(14,3))#, dpi=300)
    # intersegmental forces average
    ax11[0].plot(n_timespercent101, np.mean(einterseg_combine, 0) - np.mean(ninterseg_combine, 0))
    ax11[0].set_xlabel('% Gait cycle')
    ax11[0].set_ylabel('Force (BW)')
    ax11[0].set_title('intersegmental')
    # ax11[0].legend()
    
    # tfl forces
    ax11[1].plot(n_timespercent101, np.mean(etfl_combine, 0) - np.mean(ntfl_combine, 0))
    ax11[1].set_xlabel('% Gait cycle')
    # ax11[1].set_ylabel('Force (BW)')
    # ax11[1].legend()
    ax11[1].set_title('tfl')

    # gastroc forces
    ax11[2].plot(n_timespercent101, np.mean(egas_combine, 0) - np.mean(ngas_combine, 0))
    ax11[2].set_xlabel('% Gait cycle')
    # ax11[2].set_ylabel('Force (BW)')
    # ax11[2].legend()
    ax11[2].set_title('gastroc')

    # hamstring forces
    ax11[3].plot(n_timespercent101, np.mean(ehams_combine, 0) - np.mean(nhams_combine, 0))
    ax11[3].set_xlabel('% Gait cycle')
    # ax11[3].set_ylabel('Force (BW)')
    # ax11[3].legend()
    ax11[3].set_title('hamstrings')

    # quads forces
    ax11[4].plot(n_timespercent101, np.mean(equads_combine, 0) - np.mean(nquads_combine, 0))
    ax11[4].set_xlabel('% Gait cycle')
    # ax11[4].set_ylabel('Force (BW)')
    # ax11[4].legend()
    ax11[4].set_title('quadriceps')

    # # added all forces
    # ax11[5].plot(n_timespercent101, (np.mean(equads_combine,0) + np.mean(ehams_combine,0) + np.mean(egas_combine,0) + np.mean(etfl_combine,0) + np.mean(einterseg_combine,0)) - (np.mean(nquads_combine,0) + np.mean(nhams_combine,0) + np.mean(ngas_combine,0) + np.mean(ntfl_combine,0) + np.mean(ninterseg_combine,0)))
    # ax11[5].set_xlabel('% Gait cycle')
    # # ax11[5].set_ylabel('Force (BW)')
    # # ax11[5].legend()
    # ax11[5].set_title('Total Vertical Contact')

    # all forces from whole analysis
    ax11[5].plot(n_timespercent101, np.mean(eall_combine,0) - np.mean(nall_combine,0), label='Exotendon diff from Natural')
    ax11[5].set_xlabel('% Gait cycle')
    # ax11[5].set_ylabel('Force (BW)')
    # ax11[5].legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10)
    ax11[5].set_title('Total ' + tagcomponent + ' contact')

    # Hide the last subplot and use it to display the legend
    ax11[6].axis('off')
    handles, labels = ax11[5].get_legend_handles_labels()
    ax11[6].legend(handles, labels, loc='center', fontsize=10)
    fig11.tight_layout()
    plt.savefig(analyzedir + '\\contact3_' + whichleg + '_' + tagcomponent + '.png')
    

    ###########################################################################
    ### figure: changes in force segmented together on plot
    #### Possible paper figure for R3.    
    colors4 = ['#648FFF','#785EF0','#DC267F','#FE6100', '#FFB000']
    fig12, ax12 = plt.subplots(1, 2, figsize=(12,4.55), dpi=500)
    ax12[0].axhline(0, color='black', linewidth=1)
    # intersegmental forces average
    ax12[0].plot(n_timespercent101, np.mean(einterseg_combine, 0) - np.mean(ninterseg_combine, 0), label='Intersegmental', color=ecolormid, linewidth=3)
    # # tfl forces
    ax12[0].plot(n_timespercent101, np.mean(etfl_combine, 0) - np.mean(ntfl_combine, 0), label='Tensor Fasciae Latae', color=ncolor1, linewidth=3) #colors4[0]
    # # gastroc forces
    ax12[0].plot(n_timespercent101, np.mean(egas_combine, 0) - np.mean(ngas_combine, 0), label='Gastrocnemius', color=ncolor3, linewidth=3) #colors4[1]
    # # hamstring forces
    ax12[0].plot(n_timespercent101, np.mean(ehams_combine, 0) - np.mean(nhams_combine, 0), label='Hamstrings', color=ncolor5, linewidth=3) #colors4[2]
    # quads forces
    ax12[0].plot(n_timespercent101, np.mean(equads_combine, 0) - np.mean(nquads_combine, 0), label='Quadriceps', color=ncolor7, linewidth=3) #colors4[3]
    # # reserve forces
    # # ax12[0].plot(n_timespercent101, np.mean(ereserve_combine, 0) - np.mean(nreserve_combine, 0), label='reserves', color=colors4[4])
    # added all forces
    # ax12[0].plot(n_timespercent101, (np.mean(equads_combine,0) + np.mean(ehams_combine,0) + np.mean(egas_combine,0) + np.mean(etfl_combine,0) + np.mean(einterseg_combine,0) + np.mean(ereserve_combine,0)) - (np.mean(nquads_combine,0) + np.mean(nhams_combine,0) + np.mean(ngas_combine,0) + np.mean(ntfl_combine,0) + np.mean(ninterseg_combine,0) + np.mean(nreserve_combine,0)), label='Total')
    # all forces from whole analysis
    ax12[0].plot(n_timespercent101, np.mean(eall_combine,0) - np.mean(nall_combine,0), label='Total Compressive force', linestyle='dashed', color='black', linewidth=3)

    ax12[0].set_xlabel('% Gait cycle', fontsize=16, fontweight='bold')
    ax12[0].set_ylabel('Compressive knee\n contact force difference (BW)', fontsize=16, fontweight='bold')
    # ax12[0].set_title('Exotendon change in contact force', fontsize=16)
    ax12[0].tick_params(axis='both', which='major', labelsize=16)
    # hide the second subplots and use it for the legend
    ax12[1].axis('off')
    handles, labels = ax12[0].get_legend_handles_labels()
    ax12[1].legend(handles, labels, loc='center', fontsize=16)
    fig12.tight_layout()
    fig12.savefig(analyzedir + '\\contact4_' + whichleg + '_' + tagcomponent + '.png')

    
    ###########################################################################
    ### figure: changes in force segmented together on plot
    # combine muscles other than quads.
    eother_combine = etfl_combine + egas_combine + ehams_combine
    nother_combine = ntfl_combine + ngas_combine + nhams_combine
    eother_mean = np.mean(eother_combine, 0)
    nother_mean = np.mean(nother_combine, 0)
    other_diff = eother_mean - nother_mean

    fig22, ax22 = plt.subplots(1, 2, figsize=(12,4.55), dpi=500)
    ax22[0].axhline(0, color='black', linewidth=1)
    ax22[0].plot(n_timespercent101, np.mean(eall_combine,0) - np.mean(nall_combine,0), label='Total Compressive force', linestyle='dashed', color='black', linewidth=3)
    ax22[0].plot(n_timespercent101, np.mean(einterseg_combine, 0) - np.mean(ninterseg_combine, 0), label='Intersegmental', color=ecolormiddark, linewidth=3)    
    ax22[0].plot(n_timespercent101, np.mean(equads_combine, 0) - np.mean(nquads_combine, 0), label='Quadriceps', color=ncolor5, linewidth=3) #colors4[3]
    ax22[0].plot(n_timespercent101, other_diff, label="Hamstrings, Gastrocnemius\n& Tensor Fasciae Latae", color=ncolor3, linewidth=3) #colors4[0]

    ax22[0].set_xlabel('% Gait cycle', fontsize=16, fontweight='bold')
    ax22[0].set_ylabel('Compressive knee\n contact force difference (BW)', fontsize=16, fontweight='bold')
    # ax22[0].set_title('Exotendon change in contact force', fontsize=16)
    ax22[0].tick_params(axis='both', which='major', labelsize=16)
    # hide the second subplots and use it for the legend
    ax22[1].axis('off')
    handles, labels = ax22[0].get_legend_handles_labels()
    ax22[1].legend(handles, labels, loc='center', fontsize=16)
    fig22.tight_layout()
    fig22.savefig(analyzedir + '\\contact14_' + whichleg + '_' + tagcomponent + '.png')




    # TODO: figure out why the difference in total and all added together. 
    # okay so not in how I am adding/averaging. has to be something in how the analysis is done between them.... am I missing something??
    
    fig13, ax13 = plt.subplots(1,3, figsize=(14,5))#, dpi=300)
    # intersegmental forces average - natural
    ax13[0].plot(n_timespercent101, np.mean(ninterseg_combine, 0), label='natural_interseg')
    ax13[0].plot(n_timespercent101, np.mean(ninterseg_combine + ntfl_combine, 0), label='natural_interseg + tfl')
    ax13[0].plot(n_timespercent101, np.mean(ninterseg_combine+ntfl_combine+nhams_combine, 0), label='natural_interseg+tfl+hams')
    ax13[0].plot(n_timespercent101, np.mean(ninterseg_combine+ntfl_combine+nhams_combine+ngas_combine, 0), label='natural_interseg+tfl+hams+gas')
    ax13[0].plot(n_timespercent101, np.mean(ninterseg_combine+ntfl_combine+nhams_combine+ngas_combine+nquads_combine, 0), label='natural_interseg+tfl+hams+gas+quads')
    # ax13[0].plot(n_timespercent101, np.mean(ninterseg_combine+ntfl_combine+nhams_combine+ngas_combine+nquads_combine+nreserve_combine, 0), label='natural_interseg+tfl+hams+gas+quads+reserve')
    ax13[0].plot(n_timespercent101, np.mean(nall_combine, 0), label='nat_all', linestyle='dotted')
    ax13[0].set_xlabel('% Gait cycle')
    ax13[0].set_ylabel('Force (BW)')
    ax13[0].set_title('natural')
    # ax13[0].legend()
    # intersegmental forces average - exotendon
    ax13[1].plot(e_timespercent101, np.mean(einterseg_combine, 0), label='exo_interseg')
    ax13[1].plot(e_timespercent101, np.mean(einterseg_combine+etfl_combine, 0), label='exo_interseg + tfl')
    ax13[1].plot(e_timespercent101, np.mean(einterseg_combine+etfl_combine+ehams_combine, 0), label='exo_interseg+tfl+hams')
    ax13[1].plot(e_timespercent101, np.mean(einterseg_combine+etfl_combine+ehams_combine+egas_combine, 0), label='exo_interseg+tfl+hams+gas')
    ax13[1].plot(e_timespercent101, np.mean(einterseg_combine+etfl_combine+ehams_combine+egas_combine+equads_combine, 0), label='exo_interseg+tfl+hams+gas+quads')
    # ax13[1].plot(e_timespercent101, np.mean(einterseg_combine+etfl_combine+ehams_combine+egas_combine+equads_combine+ereserve_combine, 0), label='exo_interseg+tfl+hams+gas+quads+reserve')
    ax13[1].plot(e_timespercent101, np.mean(eall_combine, 0), label='exo_all', linestyle='dotted')
    ax13[1].set_xlabel('% Gait cycle')
    ax13[1].set_ylabel('Force (BW)')
    ax13[1].set_title('exotendon')
    # ax13[1].legend()

    # Hide the last subplot and use it to display the legend
    ax13[2].axis('off')
    handles, labels = ax13[1].get_legend_handles_labels()
    ax13[2].legend(handles, labels, loc='center', fontsize=14)
    fig13.tight_layout()
    plt.savefig(analyzedir + '\\contact5_' + whichleg + '_' + tagcomponent + '.png')
    


    ###########################################################################
    ### Figure: total pop average for right leg between nat and exo
    #### Figure in the paper R1... 
    figcon6, axcon6 = plt.subplots(1, 2, figsize=(12,4.55), dpi=500)
    axcon6[0].fill_between(n_timespercent101, np.mean(nall_combine,0)-np.std(nall_combine,0), np.mean(nall_combine,0)+np.std(nall_combine,0), color=ncolorlight, alpha=0.75)
    axcon6[0].fill_between(e_timespercent101, np.mean(eall_combine,0)-np.std(eall_combine,0), np.mean(eall_combine,0)+np.std(eall_combine,0), color=ecolorlight, alpha=0.75)
    axcon6[0].plot(n_timespercent101, np.mean(nall_combine,0), color=ncolor, label='Natural (Mean \u00B1 Std.)', linewidth=3)
    axcon6[0].plot(e_timespercent101, np.mean(eall_combine,0), color=ecolor, label='Exotendon (Mean \u00B1 Std.)', linewidth=3)
    axcon6[0].set_xlabel('% Gait cycle', fontsize=16)
    axcon6[0].set_ylabel(tagcomponent + ' knee contact force (BW)', fontsize=16)
    axcon6[0].set_title('Total ' + tagcomponent + ' contact force', fontsize=16)
    axcon6[0].tick_params(axis='both', which='major', labelsize=16)
    # axcon6[0].xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    axcon6[0].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    # hide the second subplots and use it for the legend
    axcon6[1].axis('off')
    handles, labels = axcon6[0].get_legend_handles_labels()
    axcon6[1].legend(handles, labels, loc='center', fontsize=16)
    figcon6.tight_layout()
    figcon6.savefig(analyzedir + '\\contact6_' + whichleg + '_' + tagcomponent + '.png')


    ###########################################################################
    # output the data in an excel sheet
    # Create a dictionary to hold all the data
    data_dict = {
        'ninterseg_combine': ninterseg_combine,
        'einterseg_combine': einterseg_combine,
        'nquads_combine': nquads_combine,
        'equads_combine': equads_combine,
        'nhams_combine': nhams_combine,
        'ehams_combine': ehams_combine,
        'ngas_combine': ngas_combine,
        'egas_combine': egas_combine,
        'ntfl_combine': ntfl_combine,
        'etfl_combine': etfl_combine,
        'nall_combine': nall_combine,
        'eall_combine': eall_combine,
        'nreserve_combine': nreserve_combine,
        'ereserve_combine': ereserve_combine,
        'nnone_combine': nnone_combine,
        'enone_combine': enone_combine
    }

        # 'muscleacts_nat': muscleacts_nat,
        # 'muscleacts_exo': muscleacts_exo,
        # 'moments_nat': moments_nat,
        # 'moments_exo': moments_exo,
        # 'IDmoments_nat': IDmoments_nat,
        # 'IDmoments_exo': IDmoments_exo,
        # 'activeforces_nat': activeforces_nat,
        # 'activeforces_exo': activeforces_exo,
        # 'passiveforces_nat': passiveforces_nat,
        # 'passiveforces_exo': passiveforces_exo,
        # 'totalforces_nat': totalforces_nat,
        # 'totalforces_exo': totalforces_exo

    os.chdir(analyzedir)
    # Create a Pandas Excel writer using XlsxWriter as the engine
    with pd.ExcelWriter('analysis_results'+'_'+whichleg+'_'+tagcomponent+'.xlsx', engine='xlsxwriter') as writer:
        for key, value in data_dict.items():
            print(key)
            if isinstance(value, dict):
                print('found a dict')
                pdb.set_trace()
                for sub_key, sub_value in value.items():
                    df = pd.DataFrame(sub_value)
                    df.to_excel(writer, sheet_name=f'{key}_{sub_key}')
            else:
                df = pd.DataFrame(value)
                df.to_excel(writer, sheet_name=key)


    pdb.set_trace()
    plt.show()
    print('\n\nbeyond this is the breakdown paper figures.')
    # sys.exit()

    # ###########################################################################
    # # stats for the data
    # # here is the stats for R1 - differences in the peak vertical JCF
    # mean_nall_combine = np.mean(nall_combine,0)
    # std_nall_combine = np.std(nall_combine,0)
    # mean_eall_combine = np.mean(eall_combine,0)
    # std_eall_combine = np.std(eall_combine,0)

    # # get the peaks first and then do the stats on them...
    # peaks_nall_combine = np.max(nall_combine,1)
    # idx_peaks_nall_combine = nall_combine.argmax(1)
    # peaks_eall_combine = np.max(eall_combine,1)
    # idx_peaks_eall_combine = eall_combine.argmax(1)

    # mean_peaks_nall_combine = np.mean(peaks_nall_combine)
    # mean_peaks_eall_combine = np.mean(peaks_eall_combine)
    # std_peaks_nall_combine = np.std(peaks_nall_combine)
    # std_peaks_eall_combine = np.std(peaks_eall_combine)

    # # differences in peaks
    # diff_peaks_all_combine = peaks_nall_combine - peaks_eall_combine
    # mean_diff_peaks_all_combine = np.mean(peaks_nall_combine - peaks_eall_combine)
    # std_diff_peaks_all_combine = np.std(peaks_nall_combine - peaks_eall_combine)

    # # shapiro test for normal on the peaks differences
    # res = scipy.stats.shapiro(diff_peaks_all_combine)
    
    # # t - test for peaks differences is in the excel sheet...
    
    # # percent change in peak
    # perc_diff_peaks_all_combine = (peaks_eall_combine - peaks_nall_combine) / peaks_nall_combine * 100
    # mean_perc_diff_peaks_all_combine = np.mean(perc_diff_peaks_all_combine)
    # std_perc_diff_peaks_all_combine = np.std(perc_diff_peaks_all_combine)
    

    # ###########################################################################
    # more polished figures
    # pdb.set_trace()
    
    # ### Figure: total pop average for right leg for nat - segmented shaded
    # plt.figure(dpi=300)
    # plt.plot(n_timespercent101, np.mean(ninterseg_combine,0), label='intersegmental', color=ecolor)
    # plt.fill_between(n_timespercent101, np.mean(ninterseg_combine,0), color=ecolorlight)
    # plt.plot(n_timespercent101, np.mean(ntfl_combine,0) + np.mean(ninterseg_combine,0), label='inter + tfl', color=ncolor3)
    # plt.fill_between(n_timespercent101, np.mean(ninterseg_combine,0), np.mean(ntfl_combine,0) + np.mean(ninterseg_combine,0), color=ncolor2)

    # plt.plot(n_timespercent101, np.mean(ngas_combine,0) + np.mean(ntfl_combine,0) + np.mean(ninterseg_combine,0), label='inter + tfl + gas', color=ncolor5)
    # plt.fill_between(n_timespercent101, np.mean(ntfl_combine,0) + np.mean(ninterseg_combine,0), np.mean(ngas_combine,0) + np.mean(ntfl_combine,0) + np.mean(ninterseg_combine,0), color=ncolor4)

    # plt.plot(n_timespercent101, np.mean(nhams_combine,0) + np.mean(ngas_combine,0) + np.mean(ntfl_combine,0) + np.mean(ninterseg_combine,0), label='inter + tfl + gas + hams', color=ncolor7)
    # plt.fill_between(n_timespercent101, np.mean(ngas_combine,0) + np.mean(ntfl_combine,0) + np.mean(ninterseg_combine,0), np.mean(nhams_combine,0) + np.mean(ngas_combine,0) + np.mean(ntfl_combine,0) + np.mean(ninterseg_combine,0), color=ncolor6)

    # plt.plot(n_timespercent101, np.mean(nquads_combine,0) + np.mean(nhams_combine,0) + np.mean(ngas_combine,0) + np.mean(ntfl_combine,0) + np.mean(ninterseg_combine,0), label='inter + tfl + gas + hams + quads', color=ncolor9)
    # plt.fill_between(n_timespercent101, np.mean(nhams_combine,0) + np.mean(ngas_combine,0) + np.mean(ntfl_combine,0) + np.mean(ninterseg_combine,0), np.mean(nquads_combine,0) + np.mean(nhams_combine,0) + np.mean(ngas_combine,0) + np.mean(ntfl_combine,0) + np.mean(ninterseg_combine,0), color=ncolor8)

    # plt.legend()
    # plt.ylabel(tagcomponent + ' knee contact force (BW)', fontsize=16)
    # plt.xlabel('% Gait cycle', fontsize=16)
    # plt.title(tagcomponent + ' knee contact force - Natural', fontsize=16)    
    # plt.tick_params(axis='both', which='major', labelsize=16)
    # plt.tight_layout()
    # plt.savefig(analyzedir + '\\contact7_Natural_' + whichleg + '_' + tagcomponent + '.png')
    # plt.show()

    # pdb.set_trace()
    

    # original x (should be 0..100)
    x_old = np.array(n_timespercent101, dtype=float)
    # high-resolution x for smooth plotting
    x_fine = np.linspace(x_old.min(), x_old.max(), 2000)

    # mean component curves (vertical/contact values already in BW)
    mn_inter = np.mean(ninterseg_combine, 0)
    mn_tfl   = np.mean(ntfl_combine, 0)
    mn_gas   = np.mean(ngas_combine, 0)
    mn_hams  = np.mean(nhams_combine, 0)
    mn_quads = np.mean(nquads_combine, 0)

    # cumulative contributions used in the stacked plotting
    c1 = mn_inter
    c2 = c1 + mn_tfl
    c3 = c2 + mn_gas
    c4 = c3 + mn_hams
    c5 = c4 + mn_quads

    # use cubic spline interpolation for a smoother curve (requires scipy.interpolate.interp1d)
    # ensure x_old is strictly increasing (it should be 0..100)
    f_c1 = interp1d(x_old, c1, kind='cubic', axis=0)
    f_c2 = interp1d(x_old, c2, kind='cubic', axis=0)
    f_c3 = interp1d(x_old, c3, kind='cubic', axis=0)
    f_c4 = interp1d(x_old, c4, kind='cubic', axis=0)
    f_c5 = interp1d(x_old, c5, kind='cubic', axis=0)

    # evaluate on fine grid
    y_c1 = f_c1(x_fine)
    y_c2 = f_c2(x_fine)
    y_c3 = f_c3(x_fine)
    y_c4 = f_c4(x_fine)
    y_c5 = f_c5(x_fine)

    # plot using the high-resolution (smoothed) vectors
    plt.figure(dpi=500)
    # plt.plot(x_fine, y_c1, label='intersegmental', color=ecolor)
    plt.fill_between(x_fine, y_c1, color=ecolorlight)

    # plt.plot(x_fine, y_c2, label='inter + tfl', color=ncolor3)
    plt.fill_between(x_fine, y_c1, y_c2, color=ncolor1)

    # plt.plot(x_fine, y_c3, label='inter + tfl + gas', color=ncolor5)
    plt.fill_between(x_fine, y_c2, y_c3, color=ncolor3)

    # plt.plot(x_fine, y_c4, label='inter + tfl + gas + hams', color=ncolor7)
    plt.fill_between(x_fine, y_c3, y_c4, color=ncolor5)

    # plt.plot(x_fine, y_c5, label='inter + tfl + gas + hams + quads', color=ncolor9)
    plt.fill_between(x_fine, y_c4, y_c5, color=ncolor7)

    plt.legend()
    plt.ylabel('Compressive\n knee contact force (BW)', fontsize=16, fontweight='bold')
    plt.xlabel('% Gait cycle', fontsize=16, fontweight='bold')
    # plt.title(tagcomponent + ' knee contact force - Natural', fontsize=16)
    plt.tick_params(axis='both', which='major', labelsize=16)
    plt.tight_layout()

    # capture current axis limits and persist for later reuse
    ax = plt.gca()
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

    plt.savefig(os.path.join(analyzedir, f'contact7_Natural_{whichleg}_{tagcomponent}.png'))
    plt.show()


    
    # ###########################################################################
    # ### Figure: total pop average for right leg for exo - segmented shaded
    # plt.figure(dpi=300)
    # plt.plot(e_timespercent101, np.mean(einterseg_combine,0), label='intersegmental', color=ecolor)
    # plt.fill_between(e_timespercent101, np.mean(einterseg_combine,0), color=ecolorlight)

    # plt.plot(e_timespercent101, np.mean(etfl_combine,0) + np.mean(einterseg_combine,0), label='inter + tfl', color=ncolor3)
    # plt.fill_between(e_timespercent101, np.mean(einterseg_combine,0), np.mean(etfl_combine,0) + np.mean(einterseg_combine,0), color=ncolor2)

    # plt.plot(e_timespercent101, np.mean(egas_combine,0) + np.mean(etfl_combine,0) + np.mean(einterseg_combine,0), label='inter + tfl + gas', color=ncolor5)
    # plt.fill_between(e_timespercent101, np.mean(etfl_combine,0) + np.mean(einterseg_combine,0), np.mean(egas_combine,0) + np.mean(etfl_combine,0) + np.mean(einterseg_combine,0), color=ncolor4)

    # plt.plot(e_timespercent101, np.mean(ehams_combine,0) + np.mean(egas_combine,0) + np.mean(etfl_combine,0) + np.mean(einterseg_combine,0), label='inter + tfl + gas + hams', color=ncolor7)
    # plt.fill_between(e_timespercent101, np.mean(egas_combine,0) + np.mean(etfl_combine,0) + np.mean(einterseg_combine,0), np.mean(ehams_combine,0) + np.mean(egas_combine,0) + np.mean(etfl_combine,0) + np.mean(einterseg_combine,0),color=ncolor6)

    # plt.plot(e_timespercent101, np.mean(equads_combine,0) + np.mean(ehams_combine,0) + np.mean(egas_combine,0) + np.mean(etfl_combine,0) + np.mean(einterseg_combine,0), label='inter + tfl + gas + hams + quads', color=ncolor9)
    # plt.fill_between(e_timespercent101, np.mean(ehams_combine,0) + np.mean(egas_combine,0) + np.mean(etfl_combine,0) + np.mean(einterseg_combine,0), np.mean(equads_combine,0) + np.mean(ehams_combine,0) + np.mean(egas_combine,0) + np.mean(etfl_combine,0) + np.mean(einterseg_combine,0), color=ncolor8)

    # plt.legend()
    # plt.ylabel(tagcomponent + ' knee contact force (BW)', fontsize=16)
    # plt.xlabel('% Gait cycle', fontsize=16)
    # plt.title(tagcomponent + ' knee contact force - Exotendon', fontsize=16)
    # plt.tick_params(axis='both', which='major', labelsize=16)
    # plt.tight_layout()
    # plt.savefig(analyzedir + '\\contact8_exo_' + whichleg + '_' + tagcomponent + '.png')
    # plt.show()

    # Smoothed (cubic-spline) stacked plot for the Exotendon case (analogous to the Natural plot above)
    # prepares high-res x and cumulative component means, interpolates with cubic splines, and saves png.
    x_old_e = np.array(e_timespercent101, dtype=float)

    # high-resolution x for smooth plotting
    x_fine_e = np.linspace(x_old_e.min(), x_old_e.max(), 2000)

    # mean component curves (already in BW in the combined arrays)
    mn_inter_e = np.mean(einterseg_combine, 0)
    mn_tfl_e   = np.mean(etfl_combine, 0)
    mn_gas_e   = np.mean(egas_combine, 0)
    mn_hams_e  = np.mean(ehams_combine, 0)
    mn_quads_e = np.mean(equads_combine, 0)

    # cumulative contributions for stacked plotting
    c1_e = mn_inter_e
    c2_e = c1_e + mn_tfl_e
    c3_e = c2_e + mn_gas_e
    c4_e = c3_e + mn_hams_e
    c5_e = c4_e + mn_quads_e

    # cubic spline interpolation
    f_c1_e = interp1d(x_old_e, c1_e, kind='cubic', axis=0)
    f_c2_e = interp1d(x_old_e, c2_e, kind='cubic', axis=0)
    f_c3_e = interp1d(x_old_e, c3_e, kind='cubic', axis=0)
    f_c4_e = interp1d(x_old_e, c4_e, kind='cubic', axis=0)
    f_c5_e = interp1d(x_old_e, c5_e, kind='cubic', axis=0)

    # evaluate on fine grid
    y_c1_e = f_c1_e(x_fine_e)
    y_c2_e = f_c2_e(x_fine_e)
    y_c3_e = f_c3_e(x_fine_e)
    y_c4_e = f_c4_e(x_fine_e)
    y_c5_e = f_c5_e(x_fine_e)

    # plot using the smoothed vectors
    plt.figure(dpi=500)
    # plt.plot(x_fine_e, y_c1_e, label='intersegmental', color=ecolor)
    plt.fill_between(x_fine_e, y_c1_e, color=ecolorlight)

    # plt.plot(x_fine_e, y_c2_e, label='inter + tfl', color=ncolor3)
    plt.fill_between(x_fine_e, y_c1_e, y_c2_e, color=ncolor1)

    # plt.plot(x_fine_e, y_c3_e, label='inter + tfl + gas', color=ncolor5)
    plt.fill_between(x_fine_e, y_c2_e, y_c3_e, color=ncolor3)

    # plt.plot(x_fine_e, y_c4_e, label='inter + tfl + gas + hams', color=ncolor7)
    plt.fill_between(x_fine_e, y_c3_e, y_c4_e, color=ncolor5)

    # plt.plot(x_fine_e, y_c5_e, label='inter + tfl + gas + hams + quads', color=ncolor9)
    plt.fill_between(x_fine_e, y_c4_e, y_c5_e, color=ncolor7)

    # plt.legend()
    ax2 = plt.gca()
    ax2.set_xlim(xlim)  # reuse x-axis limits from Natural plot
    ax2.set_ylim(ylim)  # reuse y-axis limits from Natural plot

    plt.ylabel('Compressive\n knee contact force (BW)', fontsize=16, fontweight='bold')
    plt.xlabel('% Gait cycle', fontsize=16, fontweight='bold')
    # plt.title(tagcomponent + ' knee contact force - Exotendon', fontsize=16)
    plt.tick_params(axis='both', which='major', labelsize=16)
    plt.tight_layout()
    plt.savefig(os.path.join(analyzedir, f'contact7_Exotendon_{whichleg}_{tagcomponent}.png'))
    plt.show()



    
    ###########################################################################
    # stats for the individual muscle segmented contributions.... check out the excel sheet
    peaks_quads_natural = np.max(nquads_combine,1)
    peaks_gas_natural = np.max(ngas_combine,1)
    peaks_hams_natural = np.max(nhams_combine,1)
    peaks_tfl_natural = np.max(ntfl_combine,1)
    peaks_inter_natural = np.max(ninterseg_combine,1)
    peaks_reserve_natural = np.max(nreserve_combine,1)

    peaks_quads_exo = np.max(equads_combine,1)
    peaks_gas_exo = np.max(egas_combine,1)
    peaks_hams_exo = np.max(ehams_combine,1)
    peaks_tfl_exo = np.max(etfl_combine,1)
    peaks_inter_exo = np.max(einterseg_combine,1)
    peaks_reserve_exo = np.max(ereserve_combine,1)
    
    print('\n\nHere is a good spot to explore the peaks if you want...\n')
    pdb.set_trace()
    return

    

if __name__ == '__main__':
    # now to define all the setup that he has and is required
    basedir = os.getcwd()
    # repodir = 'G:\\Shared drives\\Exotendon\\muscleModel\\muscleEnergyModel';
    repodir = 'C:\\Users\\jonstingel\\code\\musclemodel\\muscleEnergyModel';
    
    # current results directory
    # resultsdir = os.path.join(repodir, '..\\results');
    # analyzedir = os.path.join(repodir, '..\\analysis');



    # combo results intermediate
    # resultsdir = 'C:\\Users\\jonstingel\\code\\musclemodel\\testresults\\results\\';
    # analyzedir = 'C:\\Users\\jonstingel\\code\\musclemodel\\testresults\\analysis\\';


    # combo results intermediate - copy
    # resultsdir = 'C:\\Users\\jonstingel\\code\\musclemodel\\testresults - Copy\\results\\';
    # analyzedir = 'C:\\Users\\jonstingel\\code\\musclemodel\\testresults - Copy\\analysis\\';

    # updated and consistent results of best simulation setup
    # resultsdir = 'C:\\Users\\jonstingel\\code\\musclemodel\\testresults - Copy - Copy\\results\\';
    # analyzedir = 'C:\\Users\\jonstingel\\code\\musclemodel\\testresults - Copy - Copy\\analysis\\';
    
    resultsdir = 'C:\\Users\\jonstingel\\code\\musclemodel\\testresults - Copy - Copy - Copy\\results\\';
    analyzedir = 'C:\\Users\\jonstingel\\code\\musclemodel\\testresults - Copy - Copy - Copy\\analysis\\';

    welkexoconditions = ['welkexo']
    welknaturalconditions = ['welknatural']
    welksubjects = ['welk002','welk010','welk003','welk005','welk008','welk009','welk013'];
    thingstoplot = ['contactForces']
    trials = ['trial01','trial02','trial03','trial04']
    whichleg = 'left'
    oldnotredo = False
    runtool = False
    indresults = False
    polycalc = False

    # get some results structures going
    welknaturalstruct_combine = {}
    welkexostruct_combine = {}
    naturalstruct_combine = {}
    exostruct_combine = {}
    naturalstruct_avg = {}
    exostruct_avg = {}

    muscleacts_nat = {}
    muscleacts_exo = {}
    moments_nat = {}
    moments_exo = {}
    IDmoments_nat = {}
    IDmoments_exo = {}
    activeforces_nat = {}
    activeforces_exo = {}
    passiveforces_nat = {}
    passiveforces_exo = {}
    totalforces_nat = {}
    totalforces_exo = {}

    # all of this was commented out... need to remember what all I was doing...
    # I think a lot of this got moved into functions above...
   
    ################################################################
    # do more analysis! yay 
    ################################################################
    # think about segmenting between the muscles and their individual contributions
    # figure out how to get each of the muscles and do the analysis for knee contact of just that muscle
    # need to figure out how to do joint reaction and isolate just the muscles
    
    # 1) load the model, load a verson of the model, set passives to zero
        # then set all controls but one we want to zero
        # then all activation states to zero that we don't want
    # data for y: vertical forces
    einterseg_combine = np.zeros((len(welksubjects)*len(trials), 101))
    ninterseg_combine = np.zeros((len(welksubjects)*len(trials), 101))
    equads_combine    = np.zeros((len(welksubjects)*len(trials), 101))
    nquads_combine    = np.zeros((len(welksubjects)*len(trials), 101))
    ehams_combine     = np.zeros((len(welksubjects)*len(trials), 101))
    nhams_combine     = np.zeros((len(welksubjects)*len(trials), 101))
    egas_combine      = np.zeros((len(welksubjects)*len(trials), 101))
    ngas_combine      = np.zeros((len(welksubjects)*len(trials), 101))
    etfl_combine      = np.zeros((len(welksubjects)*len(trials), 101))
    ntfl_combine      = np.zeros((len(welksubjects)*len(trials), 101))
    eall_combine      = np.zeros((len(welksubjects)*len(trials), 101))
    nall_combine      = np.zeros((len(welksubjects)*len(trials), 101))
    ereserve_combine  = np.zeros((len(welksubjects)*len(trials), 101))
    nreserve_combine  = np.zeros((len(welksubjects)*len(trials), 101))
    enone_combine     = np.zeros((len(welksubjects)*len(trials), 101))
    nnone_combine     = np.zeros((len(welksubjects)*len(trials), 101))
    
    # data for x: anterior-posterior forces
    einterseg_combine_x = np.zeros((len(welksubjects)*len(trials), 101))
    ninterseg_combine_x = np.zeros((len(welksubjects)*len(trials), 101))
    equads_combine_x    = np.zeros((len(welksubjects)*len(trials), 101))
    nquads_combine_x    = np.zeros((len(welksubjects)*len(trials), 101))
    ehams_combine_x     = np.zeros((len(welksubjects)*len(trials), 101))
    nhams_combine_x     = np.zeros((len(welksubjects)*len(trials), 101))
    egas_combine_x      = np.zeros((len(welksubjects)*len(trials), 101))
    ngas_combine_x      = np.zeros((len(welksubjects)*len(trials), 101))
    etfl_combine_x      = np.zeros((len(welksubjects)*len(trials), 101))
    ntfl_combine_x      = np.zeros((len(welksubjects)*len(trials), 101))
    eall_combine_x      = np.zeros((len(welksubjects)*len(trials), 101))
    nall_combine_x      = np.zeros((len(welksubjects)*len(trials), 101))
    ereserve_combine_x  = np.zeros((len(welksubjects)*len(trials), 101))
    nreserve_combine_x  = np.zeros((len(welksubjects)*len(trials), 101))
    enone_combine_x     = np.zeros((len(welksubjects)*len(trials), 101))
    nnone_combine_x     = np.zeros((len(welksubjects)*len(trials), 101))

    # data for z: medial-lateral forces
    einterseg_combine_z = np.zeros((len(welksubjects)*len(trials), 101))
    ninterseg_combine_z = np.zeros((len(welksubjects)*len(trials), 101))
    equads_combine_z    = np.zeros((len(welksubjects)*len(trials), 101))
    nquads_combine_z    = np.zeros((len(welksubjects)*len(trials), 101))
    ehams_combine_z     = np.zeros((len(welksubjects)*len(trials), 101))
    nhams_combine_z     = np.zeros((len(welksubjects)*len(trials), 101))
    egas_combine_z      = np.zeros((len(welksubjects)*len(trials), 101))
    ngas_combine_z      = np.zeros((len(welksubjects)*len(trials), 101))
    etfl_combine_z      = np.zeros((len(welksubjects)*len(trials), 101))
    ntfl_combine_z      = np.zeros((len(welksubjects)*len(trials), 101))
    eall_combine_z      = np.zeros((len(welksubjects)*len(trials), 101))
    nall_combine_z      = np.zeros((len(welksubjects)*len(trials), 101))
    ereserve_combine_z  = np.zeros((len(welksubjects)*len(trials), 101))
    nreserve_combine_z  = np.zeros((len(welksubjects)*len(trials), 101))
    enone_combine_z     = np.zeros((len(welksubjects)*len(trials), 101))
    nnone_combine_z     = np.zeros((len(welksubjects)*len(trials), 101))
    
    # define the muscles that we want in each of the splits
    musclesWanted = {}
    if whichleg == 'right':
        musclesWanted['inter'] = []
        musclesWanted['quads'] = ['vaslat_r', 'vasmed_r', 'vasint_r', 'recfem_r']
        musclesWanted['hams']  = ['bflh_r', 'bfsh_r', 'grac_r', 'sart_r', 'semimem_r', 'semiten_r']
        musclesWanted['gas']   = ['gaslat_r', 'gasmed_r']
        musclesWanted['tfl']   = ['tfl_r']
        musclesWanted['all'] = ['all']
        musclesWanted['reserve'] = ['reserve']
        musclesWanted['none'] = ['none']
    elif whichleg == 'left':
        musclesWanted['inter'] = []
        musclesWanted['quads'] = ['vaslat_l', 'vasmed_l', 'vasint_l', 'recfem_l']
        musclesWanted['hams']  = ['bflh_l', 'bfsh_l', 'grac_l', 'sart_l', 'semimem_l', 'semiten_l']
        musclesWanted['gas']   = ['gaslat_l', 'gasmed_l']
        musclesWanted['tfl']   = ['tfl_l']
        musclesWanted['all'] = ['all']
        musclesWanted['reserve'] = ['reserve']
        musclesWanted['none'] = ['none']
    elif whichleg == 'both':
        musclesWanted['inter'] = []
        musclesWanted['quads'] = ['vaslat_r', 'vasmed_r', 'vasint_r', 'recfem_r','vaslat_l', 'vasmed_l', 'vasint_l', 'recfem_l']
        musclesWanted['hams']  = ['bflh_r', 'bfsh_r', 'grac_r', 'sart_r', 'semimem_r', 'semiten_r','bflh_l', 'bfsh_l', 'grac_l', 'sart_l', 'semimem_l', 'semiten_l']
        musclesWanted['gas']   = ['gaslat_r', 'gasmed_r','gaslat_l', 'gasmed_l']
        musclesWanted['tfl']   = ['tfl_r','tfl_l']
        musclesWanted['all'] = ['all']
        musclesWanted['reserve'] = ['reserve']
        musclesWanted['none'] = ['none']
    # TODO
    
    # plot the stuff
    # figure out subtracting the intersegmental from each of the others?

    # # have the structures, now loop through and figure out how to fill them in
    # # loop the subjects
    # '''
    spot = 0
    for subj in range(len(welksubjects)):
        subject = welksubjects[subj]
        subjdir = os.path.join(resultsdir, subject)
        # create a structure for individual subject stuff

        # loop through conditions
        for cond in range(len(welkexoconditions)):
            condition = welkexoconditions[cond]
            condir = os.path.join(subjdir, condition)
            # loop the trials
            for tr in range(len(trials)):
                trial = trials[tr]
                trialdir = os.path.join(condir, trial)
                print(trialdir)
                # print(spot)
                
                # grab the model weight
                # modelmass = get_model_total_mass(trialdir, 'simple_model_all_the_probes_adjusted.osim')
                # modelmass = get_model_total_mass(trialdir, 'simple_model_all_the_probes.osim')
                modelmass = get_model_total_mass(trialdir, 'simple_model_all_the_probes.osim')
                naturalstruct_combine[subject] = modelmass
                
                ### now what do we want to do at each of the trials
                # write a method that gives the knee contact forces, 
                # pass in model, and the muscles that you want
                os.chdir(trialdir)                
                
                # # # quads and intersegmental forces - method 1
                # jrasrquads01 = getKneeContactributions(trialdir, musclesWanted['quads'], 'quads01')
                # jrasrquads001 = getKneeContactributions001(trialdir, musclesWanted['quads'], 'quads001')
                # jrasrquads0001 = getKneeContactributions0001(trialdir, musclesWanted['quads'], 'quads0001')
                # jrasrquads01_wt1 = getKneeContactributions(trialdir, musclesWanted['quads'], 'quads01_wt01')
                
                # # # intersegmental only
                # jrasrinter01 = getKneeContactributions(trialdir, musclesWanted['inter'], 'inter01')
                # jrasrinter001 = getKneeContactributions001(trialdir, musclesWanted['inter'], 'inter001')
                # jrasrinter0001 = getKneeContactributions0001(trialdir, musclesWanted['inter'], 'inter0001')
                # jrasrinter01_wt1 = getKneeContactributions(trialdir, musclesWanted['inter'], 'inter01_wt1')
                
                # plt.figure()
                # plt.plot(jrasrinter01, color='red', label='convergeTol = .01')
                # plt.plot(jrasrinter001, color='orange', label='convergeTol = .001')
                # plt.plot(jrasrinter0001, color='yellow', label='convergeTol = .0001')
                # plt.plot(jrasrinter01_wt1, color='pink', label='upped effort weight')
                # plt.legend(loc='lower right')
                # plt.ylabel('Tibia intersegmental vertical contact (N)')
                # plt.title(subject + condition + trial)
                
                # # # subtract intersegmental from quads
                # jrasrquadsonly01 = jrasrquads01 - jrasrinter01
                # jrasrquadsonly001 = jrasrquads001 - jrasrinter001
                # jrasrquadsonly0001 = jrasrquads0001 - jrasrinter0001
                # jrasrquadsonly01_wt1 = jrasrquads01_wt1 - jrasrinter01_wt1
                
                # plt.figure()
                # plt.plot(jrasrquadsonly01, color='red', label='convergeTol = .01')
                # plt.plot(jrasrquadsonly001, color='orange', label='convergeTol = .001')
                # plt.plot(jrasrquadsonly0001, color='yellow', label='convergeTol = .0001')
                # plt.plot(jrasrquadsonly01_wt1, color='pink', label='upped effort a bit')
                # plt.legend(loc='lower right')
                # plt.ylabel('Tibia quads vertical contact (N)')
                # plt.title(subject + condition + trial)
                
                
                # test = jrasr0001 - jrasr01
                # plt.figure(); plt.plot(test)




                try:
                    if oldnotredo:
                        ## okay now going to focus on the figures that I actually wanted 
                        jrasrquads = getKneeContactributions(trialdir, musclesWanted['quads'], 'quads', whichleg, runtool)
                        jrasrhams = getKneeContactributions(trialdir, musclesWanted['hams'], 'hams', whichleg, runtool)
                        jrasrgas = getKneeContactributions(trialdir, musclesWanted['gas'], 'gas', whichleg, runtool)
                        jrasrtfl = getKneeContactributions(trialdir, musclesWanted['tfl'], 'tfl', whichleg, runtool)
                        jrasrinter = getKneeContactributions(trialdir, musclesWanted['inter'], 'inter', whichleg, runtool)
                        jrasrall = getKneeContactributions(trialdir, musclesWanted['all'], 'all', whichleg, runtool)
                        jrasrinterreserve = getKneeContactributions(trialdir, musclesWanted['reserve'], 'reserve', whichleg, runtool)
                        jrasrnone = getKneeContactributions(trialdir, musclesWanted['none'], 'none', whichleg, runtool)
                    elif polycalc:
                        ## okay now going to focus on the figures that I actually wanted 
                        jrasrquadsx, jrasrquads, jrasrquadsz = getKneeContactributionsRedoPoly(trialdir, musclesWanted['quads'], 'quads', whichleg, runtool)
                        jrasrhamsx, jrasrhams, jrasrhamsz = getKneeContactributionsRedoPoly(trialdir, musclesWanted['hams'], 'hams', whichleg, runtool)
                        jrasrgasx, jrasrgas, jrasrgasz = getKneeContactributionsRedoPoly(trialdir, musclesWanted['gas'], 'gas', whichleg, runtool)
                        jrasrtflx, jrasrtfl, jrasrtflz = getKneeContactributionsRedoPoly(trialdir, musclesWanted['tfl'], 'tfl', whichleg, runtool)
                        jrasrinterx, jrasrinter, jrasrinterz = getKneeContactributionsRedoPoly(trialdir, musclesWanted['inter'], 'inter', whichleg, runtool)
                        jrasrallx, jrasrall, jrasrallz = getKneeContactributionsRedoPoly(trialdir, musclesWanted['all'], 'all', whichleg, runtool)
                        jrasrinterreservex, jrasrinterreserve, jrasrinterreservez = getKneeContactributionsRedoPoly(trialdir, musclesWanted['reserve'], 'reserve', whichleg, runtool)
                        jrasrnonex, jrasrnone, jrasrnonez = getKneeContactributionsRedoPoly(trialdir, musclesWanted['none'], 'none', whichleg, runtool)
                    else:
                        ## okay now going to focus on the figures that I actually wanted 
                        jrasrquadsx, jrasrquads, jrasrquadsz = getKneeContactributionsRedo(trialdir, musclesWanted['quads'], 'quads', whichleg, runtool)
                        jrasrhamsx, jrasrhams, jrasrhamsz = getKneeContactributionsRedo(trialdir, musclesWanted['hams'], 'hams', whichleg, runtool)
                        jrasrgasx, jrasrgas, jrasrgasz = getKneeContactributionsRedo(trialdir, musclesWanted['gas'], 'gas', whichleg, runtool)
                        jrasrtflx, jrasrtfl, jrasrtflz = getKneeContactributionsRedo(trialdir, musclesWanted['tfl'], 'tfl', whichleg, runtool)
                        jrasrinterx, jrasrinter, jrasrinterz = getKneeContactributionsRedo(trialdir, musclesWanted['inter'], 'inter', whichleg, runtool)
                        jrasrallx, jrasrall, jrasrallz = getKneeContactributionsRedo(trialdir, musclesWanted['all'], 'all', whichleg, runtool)
                        jrasrinterreservex, jrasrinterreserve, jrasrinterreservez = getKneeContactributionsRedo(trialdir, musclesWanted['reserve'], 'reserve', whichleg, runtool)
                        jrasrnonex, jrasrnone, jrasrnonez = getKneeContactributionsRedo(trialdir, musclesWanted['none'], 'none', whichleg, runtool)
                except:
                    print('Error with: ' + trialdir)
                    pdb.set_trace()
                    continue

                #### do some other data grabs here for the other data that we care about in each trial. 

                # start with muscle activations
                muscleacts_exo = ouf.getMuscleActivations(trialdir, muscleacts_exo)
                # and the joint moments
                moments_exo = ouf.getJointMoments(trialdir, moments_exo, modelmass)
                # compare to ID moments
                IDmoments_exo = ouf.getIDMoments(trialdir, IDmoments_exo, modelmass)
                # and the muscle forces
                activeforces_exo, passiveforces_exo, totalforces_exo = ouf.getMuscleForces(trialdir, activeforces_exo, passiveforces_exo, totalforces_exo, modelmass)

                # plt.figure(figsize=(11,8), dpi=300); 
                # plt.plot(jrasrall, label='all'); 
                # plt.plot(jrasrnone, label='none');
                # plt.plot(jrasrinterreserve, label='interreserve',linestyle='dashed'); 
                # plt.plot(jrasrtfl, label='tfl'); 
                # plt.plot(jrasrgas, label='gas'); 
                # plt.plot(jrasrhams, label='hams'); 
                # plt.plot(jrasrquads, label='quads');
                # plt.plot(jrasrinter, label='inter',linestyle='dotted')
                # plt.plot(jrasrnone+(jrasrtfl-jrasrnone)+(jrasrgas-jrasrnone)+(jrasrhams-jrasrnone)+(jrasrquads-jrasrnone), label='added', linestyle='dashed'); 
                
                # plt.legend()
                # important: interreserve has the reserves removed, where inter includes them still. interreserves is the only one that removes the reserves...
                

                ### here is where we have to do more work to get the other components.... 
                ### Y: vertical forces
                # subtract out the inter segmental
                jrasrinteronly = jrasrinterreserve
                jrasrreserveonly = jrasrinter - jrasrinterreserve
                jrasrnoneonly = jrasrnone

                jrasrquadsonly = jrasrquads - jrasrnone
                jrasrhamsonly = jrasrhams - jrasrnone
                jrasrgasonly = jrasrgas - jrasrnone
                jrasrtflonly = jrasrtfl - jrasrnone                
                
                # get percentages
                e_timespercent101 = np.arange(0,101,1)
                e_times = np.arange(0,len(jrasrquadsonly),1)
                e_timesinterp = np.linspace(0,len(e_times), 103)

                # get something in BW and interp to 100% gait cycle points. 
                jrasrquadsonly101 = -1*(np.interp(e_timesinterp, e_times, jrasrquadsonly)) / (modelmass * 9.81)
                jrasrhamsonly101 = -1*(np.interp(e_timesinterp, e_times, jrasrhamsonly)) / (modelmass * 9.81)
                jrasrgasonly101 = -1*(np.interp(e_timesinterp, e_times, jrasrgasonly)) / (modelmass * 9.81)
                jrasrtflonly101 = -1*(np.interp(e_timesinterp, e_times, jrasrtflonly)) / (modelmass * 9.81)
                jrasrinteronly101 = -1*(np.interp(e_timesinterp, e_times, jrasrinteronly)) / (modelmass * 9.81)
                jrasrallonly101 = -1*(np.interp(e_timesinterp, e_times, jrasrall)) / (modelmass * 9.81)
                jrasrreserveonly101 = -1*(np.interp(e_timesinterp, e_times, jrasrreserveonly)) / (modelmass * 9.81)
                jrasrnoneonly101 = -1*(np.interp(e_timesinterp, e_times, jrasrnone)) / (modelmass * 9.81)
                
                # data for y: vertical forces
                einterseg_combine[spot, :] = jrasrinteronly101[:-2]
                equads_combine[spot, :] = jrasrquadsonly101[:-2]
                ehams_combine[spot,:] = jrasrhamsonly101[:-2]
                egas_combine[spot,:] = jrasrgasonly101[:-2]
                etfl_combine[spot,:] = jrasrtflonly101[:-2]
                eall_combine[spot,:] = jrasrallonly101[:-2]
                ereserve_combine[spot,:] = jrasrreserveonly101[:-2]
                enone_combine[spot,:] = jrasrnoneonly101[:-2]
                

                ### X: anterior-posterior forces
                # subtract out the inter segmental
                jrasrinteronlyx = jrasrinterreservex
                jrasrreserveonlyx = jrasrinterx - jrasrinterreservex
                jrasrnoneonlyx = jrasrnonex

                jrasrquadsonlyx = jrasrquadsx - jrasrnonex
                jrasrhamsonlyx = jrasrhamsx - jrasrnonex
                jrasrgasonlyx = jrasrgasx - jrasrnonex
                jrasrtflonlyx = jrasrtflx - jrasrnonex        
                
                # get percentages
                e_timespercent101 = np.arange(0,101,1)
                e_timesx = np.arange(0,len(jrasrquadsonlyx),1)
                e_timesinterpx = np.linspace(0,len(e_timesx), 103)

                # get something in BW and interp to 100% gait cycle points. 
                jrasrquadsonly101x = -1*(np.interp(e_timesinterpx, e_timesx, jrasrquadsonlyx)) / (modelmass * 9.81)
                jrasrhamsonly101x = -1*(np.interp(e_timesinterpx, e_timesx, jrasrhamsonlyx)) / (modelmass * 9.81)
                jrasrgasonly101x = -1*(np.interp(e_timesinterpx, e_timesx, jrasrgasonlyx)) / (modelmass * 9.81)
                jrasrtflonly101x = -1*(np.interp(e_timesinterpx, e_timesx, jrasrtflonlyx)) / (modelmass * 9.81)
                jrasrinteronly101x = -1*(np.interp(e_timesinterpx, e_timesx, jrasrinteronlyx)) / (modelmass * 9.81)
                jrasrallonly101x = -1*(np.interp(e_timesinterpx, e_timesx, jrasrallx)) / (modelmass * 9.81)
                jrasrreserveonly101x = -1*(np.interp(e_timesinterpx, e_timesx, jrasrreserveonlyx)) / (modelmass * 9.81)
                jrasrnoneonly101x = -1*(np.interp(e_timesinterpx, e_timesx, jrasrnonex)) / (modelmass * 9.81)

                # data for x: anterior-posterior forces
                einterseg_combine_x[spot, :] = jrasrinteronly101x[:-2]
                equads_combine_x[spot, :] = jrasrquadsonly101x[:-2]
                ehams_combine_x[spot,:] = jrasrhamsonly101x[:-2]
                egas_combine_x[spot,:] = jrasrgasonly101x[:-2]
                etfl_combine_x[spot,:] = jrasrtflonly101x[:-2]
                eall_combine_x[spot,:] = jrasrallonly101x[:-2]
                ereserve_combine_x[spot,:] = jrasrreserveonly101x[:-2]
                enone_combine_x[spot,:] = jrasrnoneonly101x[:-2]
                

                ### Z: medial-lateral forces
                # subtract out the inter segmental
                jrasrinteronlyz = jrasrinterreservez
                jrasrreserveonlyz = jrasrinterz - jrasrinterreservez
                jrasrnoneonlyz = jrasrnonez

                jrasrquadsonlyz = jrasrquadsz - jrasrnonez
                jrasrhamsonlyz = jrasrhamsz - jrasrnonez
                jrasrgasonlyz = jrasrgasz - jrasrnonez
                jrasrtflonlyz = jrasrtflz - jrasrnonez        
                
                # get percentages
                e_timespercent101 = np.arange(0,101,1)
                e_timesz = np.arange(0,len(jrasrquadsonlyz),1)
                e_timesinterpz = np.linspace(0,len(e_timesz), 103)

                # get something in BW and interp to 100% gait cycle points. 
                jrasrquadsonly101z = -1*(np.interp(e_timesinterpz, e_timesz, jrasrquadsonlyz)) / (modelmass * 9.81)
                jrasrhamsonly101z = -1*(np.interp(e_timesinterpz, e_timesz, jrasrhamsonlyz)) / (modelmass * 9.81)
                jrasrgasonly101z = -1*(np.interp(e_timesinterpz, e_timesz, jrasrgasonlyz)) / (modelmass * 9.81)
                jrasrtflonly101z = -1*(np.interp(e_timesinterpz, e_timesz, jrasrtflonlyz)) / (modelmass * 9.81)
                jrasrinteronly101z = -1*(np.interp(e_timesinterpz, e_timesz, jrasrinteronlyz)) / (modelmass * 9.81)
                jrasrallonly101z = -1*(np.interp(e_timesinterpz, e_timesz, jrasrallz)) / (modelmass * 9.81)
                jrasrreserveonly101z = -1*(np.interp(e_timesinterpz, e_timesz, jrasrreserveonlyz)) / (modelmass * 9.81)
                jrasrnoneonly101z = -1*(np.interp(e_timesinterpz, e_timesz, jrasrnonez)) / (modelmass * 9.81)

                # data for z: medial-lateral forces
                einterseg_combine_z[spot, :] = jrasrinteronly101z[:-2]
                equads_combine_z[spot, :] = jrasrquadsonly101z[:-2]
                ehams_combine_z[spot,:] = jrasrhamsonly101z[:-2]
                egas_combine_z[spot,:] = jrasrgasonly101z[:-2]
                etfl_combine_z[spot,:] = jrasrtflonly101z[:-2]
                eall_combine_z[spot,:] = jrasrallonly101z[:-2]
                ereserve_combine_z[spot,:] = jrasrreserveonly101z[:-2]
                enone_combine_z[spot,:] = jrasrnoneonly101z[:-2]
                
                ## increase the spot - count of trials                
                spot += 1

    # '''        
    # now the natural condition
    spot = 0
    for subj in range(len(welksubjects)):
        subject = welksubjects[subj]
        subjdir = os.path.join(resultsdir, subject)
        # create a structure for individual subject stuff

        # loop through conditions
        for cond in range(len(welknaturalconditions)):
            condition = welknaturalconditions[cond]
            condir = os.path.join(subjdir, condition)
            # loop the trials
            for tr in range(len(trials)):
                trial = trials[tr]
                trialdir = os.path.join(condir, trial)
                print(trialdir)
                # print(spot)
                
                # grab the model weight
                # modelmass = get_model_total_mass(trialdir, 'simple_model_all_the_probes_adjusted.osim')
                modelmass = get_model_total_mass(trialdir, 'simple_model_all_the_probes.osim')
                naturalstruct_combine[subject] = modelmass
                
                ### now what do we want to do at each of the trials
                # write a method that gives the knee contact forces, 
                # pass in model, and the muscles that you want
                os.chdir(trialdir)                
                
                # # # quads and intersegmental forces - method 1
                # jrasrquads01 = getKneeContactributions(trialdir, musclesWanted['quads'], 'quads01')
                # jrasrquads001 = getKneeContactributions001(trialdir, musclesWanted['quads'], 'quads001')
                # jrasrquads0001 = getKneeContactributions0001(trialdir, musclesWanted['quads'], 'quads0001')
                # jrasrquads01wt1 = getKneeContactributions01_wt1(trialdir, musclesWanted['quads'], 'quads01wt1')
                # jrasrquads01wt2 = getKneeContactributions01_wt2(trialdir, musclesWanted['quads'], 'quads01wt2')
                
                # # intersegmental only
                # jrasrinter01 = getKneeContactributions(trialdir, musclesWanted['inter'], 'inter01')
                # jrasrinter001 = getKneeContactributions001(trialdir, musclesWanted['inter'], 'inter001')
                # jrasrinter0001 = getKneeContactributions0001(trialdir, musclesWanted['inter'], 'inter0001')
                # jrasrinter01wt1 = getKneeContactributions01_wt1(trialdir, musclesWanted['inter'], 'inter01wt1')
                # jrasrinter01wt2 = getKneeContactributions01_wt2(trialdir, musclesWanted['inter'], 'inter01wt2')


                # plt.figure()
                # plt.plot(jrasrinter01, color='red', label='convergeTol = .01')
                # plt.plot(jrasrinter001, color='orange', label='convergeTol = .001')
                # plt.plot(jrasrinter0001, color='yellow', label='convergeTol = .0001')
                # plt.plot(jrasrinter01wt1, color='pink', label='increase effort 1')
                # plt.plot(jrasrinter01wt2, color='purple', label='increase effort 2')
                
                # plt.legend(loc='lower right')
                # plt.ylabel('Tibia intersegmental vertical contact (N)')
                # plt.title(subject + condition + trial)
                
                # # subtract intersegmental from quads
                # jrasrquadsonly01 = jrasrquads01 - jrasrinter01
                # jrasrquadsonly001 = jrasrquads001 - jrasrinter001
                # jrasrquadsonly0001 = jrasrquads0001 - jrasrinter0001
                # jrasrquadsonly01wt1 = jrasrquads01wt1 - jrasrinter01wt1
                # jrasrquadsonly01wt2 = jrasrquads01wt2 - jrasrinter01wt2
                
                # plt.figure()
                # plt.plot(jrasrquadsonly01, color='red', label='convergeTol = .01')
                # plt.plot(jrasrquadsonly001, color='orange', label='convergeTol = .001')
                # plt.plot(jrasrquadsonly0001, color='yellow', label='convergeTol = .0001')
                # plt.plot(jrasrquadsonly01wt1, color='pink', label='increase effort 1')
                # plt.plot(jrasrquadsonly01wt2, color='purple', label='increase effort 2')
                
                
                # plt.legend(loc='lower right')
                # plt.ylabel('Tibia quads vertical contact (N)')
                # plt.title(subject + condition + trial)
                
                
                # test = jrasr0001 - jrasr01
                # plt.figure(); plt.plot(test)


                
                try:
                    if oldnotredo:
                        jrasrquads = getKneeContactributions(trialdir, musclesWanted['quads'], 'quads', whichleg, runtool)
                        jrasrhams = getKneeContactributions(trialdir, musclesWanted['hams'], 'hams', whichleg, runtool)
                        jrasrgas = getKneeContactributions(trialdir, musclesWanted['gas'], 'gas', whichleg, runtool)
                        jrasrtfl = getKneeContactributions(trialdir, musclesWanted['tfl'], 'tfl', whichleg, runtool)
                        jrasrinter = getKneeContactributions(trialdir, musclesWanted['inter'], 'inter', whichleg, runtool)
                        jrasrall = getKneeContactributions(trialdir, musclesWanted['all'], 'all', whichleg, runtool)
                        jrasrinterreserve = getKneeContactributions(trialdir, musclesWanted['reserve'], 'reserve', whichleg, runtool)
                        jrasrnone = getKneeContactributions(trialdir, musclesWanted['none'], 'none', whichleg, runtool)
                    elif polycalc:
                        jrasrquadsx, jrasrquads, jrasrquadsz = getKneeContactributionsRedoPoly(trialdir, musclesWanted['quads'], 'quads', whichleg, runtool)
                        jrasrhamsx, jrasrhams, jrasrhamsz = getKneeContactributionsRedoPoly(trialdir, musclesWanted['hams'], 'hams', whichleg, runtool)
                        jrasrgasx, jrasrgas, jrasrgasz = getKneeContactributionsRedoPoly(trialdir, musclesWanted['gas'], 'gas', whichleg, runtool)
                        jrasrtflx, jrasrtfl, jrasrtflz = getKneeContactributionsRedoPoly(trialdir, musclesWanted['tfl'], 'tfl', whichleg, runtool)
                        jrasrinterx, jrasrinter, jrasrinterz = getKneeContactributionsRedoPoly(trialdir, musclesWanted['inter'], 'inter', whichleg, runtool)
                        jrasrallx, jrasrall, jrasrallz = getKneeContactributionsRedoPoly(trialdir, musclesWanted['all'], 'all', whichleg, runtool)
                        jrasrinterreservex, jrasrinterreserve, jrasrinterreservez = getKneeContactributionsRedoPoly(trialdir, musclesWanted['reserve'], 'reserve', whichleg, runtool)
                        jrasrnonex, jrasrnone, jrasrnonez = getKneeContactributionsRedoPoly(trialdir, musclesWanted['none'], 'none', whichleg, runtool)
                    else:
                        ## okay now going to focus on the figures that I actually wanted 
                        jrasrquadsx, jrasrquads, jrasrquadsz = getKneeContactributionsRedo(trialdir, musclesWanted['quads'], 'quads', whichleg, runtool)
                        jrasrhamsx, jrasrhams, jrasrhamsz = getKneeContactributionsRedo(trialdir, musclesWanted['hams'], 'hams', whichleg, runtool)
                        jrasrgasx, jrasrgas, jrasrgasz = getKneeContactributionsRedo(trialdir, musclesWanted['gas'], 'gas', whichleg, runtool)
                        jrasrtflx, jrasrtfl, jrasrtflz = getKneeContactributionsRedo(trialdir, musclesWanted['tfl'], 'tfl', whichleg, runtool)
                        jrasrinterx, jrasrinter, jrasrinterz = getKneeContactributionsRedo(trialdir, musclesWanted['inter'], 'inter', whichleg, runtool)
                        jrasrallx, jrasrall, jrasrallz = getKneeContactributionsRedo(trialdir, musclesWanted['all'], 'all', whichleg, runtool)
                        jrasrinterreservex, jrasrinterreserve, jrasrinterreservez = getKneeContactributionsRedo(trialdir, musclesWanted['reserve'], 'reserve', whichleg, runtool)
                        jrasrnonex, jrasrnone, jrasrnonez = getKneeContactributionsRedo(trialdir, musclesWanted['none'], 'none', whichleg, runtool)
                except:
                    print('Error with: ' + trialdir)
                    pdb.set_trace()
                    continue
                # important: interreserve has the reserves removed, where inter includes them still. interreserves is the only one that removes the reserves...                

                # start with muscle activations
                muscleacts_nat = ouf.getMuscleActivations(trialdir, muscleacts_nat)
                # and the joint moments
                moments_nat = ouf.getJointMoments(trialdir, moments_nat, modelmass)
                # compare to ID moments
                IDmoments_nat = ouf.getIDMoments(trialdir, IDmoments_nat, modelmass)
                # and the muscle forces
                activeforces_nat, passiveforces_nat, totalforces_nat = ouf.getMuscleForces(trialdir, activeforces_nat, passiveforces_nat, totalforces_nat, modelmass)
                # now metabolics would be good as well

                # and do the kinematics as well

                # and make sure to look at the residuals too

                
                ### here is where we have to do more work to get the other components.... 
                ### Y: vertical forces
                # subtract out the inter segmental
                jrasrinteronly = jrasrinterreserve
                jrasrreserveonly = jrasrinter - jrasrinterreserve
                jrasrnoneonly = jrasrnone

                jrasrquadsonly = jrasrquads - jrasrnone
                jrasrhamsonly = jrasrhams - jrasrnone
                jrasrgasonly = jrasrgas - jrasrnone
                jrasrtflonly = jrasrtfl - jrasrnone                
                
                # get percentages
                n_timespercent101 = np.arange(0,101,1)
                n_times = np.arange(0,len(jrasrquadsonly),1)
                n_timesinterp = np.linspace(0,len(n_times), 103)

                # get something in BW and interp to 100% gait cycle points. 
                jrasrquadsonly101 = -1*(np.interp(n_timesinterp, n_times, jrasrquadsonly)) / (modelmass * 9.81)
                jrasrhamsonly101 = -1*(np.interp(n_timesinterp, n_times, jrasrhamsonly)) / (modelmass * 9.81)
                jrasrgasonly101 = -1*(np.interp(n_timesinterp, n_times, jrasrgasonly)) / (modelmass * 9.81)
                jrasrtflonly101 = -1*(np.interp(n_timesinterp, n_times, jrasrtflonly)) / (modelmass * 9.81)
                jrasrinteronly101 = -1*(np.interp(n_timesinterp, n_times, jrasrinteronly)) / (modelmass * 9.81)
                jrasrallonly101 = -1*(np.interp(n_timesinterp, n_times, jrasrall)) / (modelmass * 9.81)
                jrasrreserveonly101 = -1*(np.interp(n_timesinterp, n_times, jrasrreserveonly)) / (modelmass * 9.81)
                jrasrnoneonly101 = -1*(np.interp(n_timesinterp, n_times, jrasrnone)) / (modelmass * 9.81)
                
                # data for y: vertical forces
                ninterseg_combine[spot, :] = jrasrinteronly101[:-2]
                nquads_combine[spot, :] = jrasrquadsonly101[:-2]
                nhams_combine[spot,:] = jrasrhamsonly101[:-2]
                ngas_combine[spot,:] = jrasrgasonly101[:-2]
                ntfl_combine[spot,:] = jrasrtflonly101[:-2]
                nall_combine[spot,:] = jrasrallonly101[:-2]
                nreserve_combine[spot,:] = jrasrreserveonly101[:-2]
                nnone_combine[spot,:] = jrasrnoneonly101[:-2]
                

                ### X: anterior-posterior forces
                # subtract out the inter segmental
                jrasrinteronlyx = jrasrinterreservex
                jrasrreserveonlyx = jrasrinterx - jrasrinterreservex
                jrasrnoneonlyx = jrasrnonex

                jrasrquadsonlyx = jrasrquadsx - jrasrnonex
                jrasrhamsonlyx = jrasrhamsx - jrasrnonex
                jrasrgasonlyx = jrasrgasx - jrasrnonex
                jrasrtflonlyx = jrasrtflx - jrasrnonex        
                
                # get percentages
                n_timespercent101 = np.arange(0,101,1)
                n_timesx = np.arange(0,len(jrasrquadsonlyx),1)
                n_timesinterpx = np.linspace(0,len(n_timesx), 103)

                # get something in BW and interp to 100% gait cycle points. 
                jrasrquadsonly101x = -1*(np.interp(n_timesinterpx, n_timesx, jrasrquadsonlyx)) / (modelmass * 9.81)
                jrasrhamsonly101x = -1*(np.interp(n_timesinterpx, n_timesx, jrasrhamsonlyx)) / (modelmass * 9.81)
                jrasrgasonly101x = -1*(np.interp(n_timesinterpx, n_timesx, jrasrgasonlyx)) / (modelmass * 9.81)
                jrasrtflonly101x = -1*(np.interp(n_timesinterpx, n_timesx, jrasrtflonlyx)) / (modelmass * 9.81)
                jrasrinteronly101x = -1*(np.interp(n_timesinterpx, n_timesx, jrasrinteronlyx)) / (modelmass * 9.81)
                jrasrallonly101x = -1*(np.interp(n_timesinterpx, n_timesx, jrasrallx)) / (modelmass * 9.81)
                jrasrreserveonly101x = -1*(np.interp(n_timesinterpx, n_timesx, jrasrreserveonlyx)) / (modelmass * 9.81)
                jrasrnoneonly101x = -1*(np.interp(n_timesinterpx, n_timesx, jrasrnonex)) / (modelmass * 9.81)

                # data for x: anterior-posterior forces
                ninterseg_combine_x[spot, :] = jrasrinteronly101x[:-2]
                nquads_combine_x[spot, :] = jrasrquadsonly101x[:-2]
                nhams_combine_x[spot,:] = jrasrhamsonly101x[:-2]
                ngas_combine_x[spot,:] = jrasrgasonly101x[:-2]
                ntfl_combine_x[spot,:] = jrasrtflonly101x[:-2]
                nall_combine_x[spot,:] = jrasrallonly101x[:-2]
                nreserve_combine_x[spot,:] = jrasrreserveonly101x[:-2]
                nnone_combine_x[spot,:] = jrasrnoneonly101x[:-2]
                

                ### Z: medial-lateral forces
                # subtract out the inter segmental
                jrasrinteronlyz = jrasrinterreservez
                jrasrreserveonlyz = jrasrinterz - jrasrinterreservez
                jrasrnoneonlyz = jrasrnonez

                jrasrquadsonlyz = jrasrquadsz - jrasrnonez
                jrasrhamsonlyz = jrasrhamsz - jrasrnonez
                jrasrgasonlyz = jrasrgasz - jrasrnonez
                jrasrtflonlyz = jrasrtflz - jrasrnonez        
                
                # get percentages
                n_timespercent101 = np.arange(0,101,1)
                n_timesz = np.arange(0,len(jrasrquadsonlyz),1)
                n_timesinterpz = np.linspace(0,len(n_timesz), 103)

                # get something in BW and interp to 100% gait cycle points. 
                jrasrquadsonly101z = -1*(np.interp(n_timesinterpz, n_timesz, jrasrquadsonlyz)) / (modelmass * 9.81)
                jrasrhamsonly101z = -1*(np.interp(n_timesinterpz, n_timesz, jrasrhamsonlyz)) / (modelmass * 9.81)
                jrasrgasonly101z = -1*(np.interp(n_timesinterpz, n_timesz, jrasrgasonlyz)) / (modelmass * 9.81)
                jrasrtflonly101z = -1*(np.interp(n_timesinterpz, n_timesz, jrasrtflonlyz)) / (modelmass * 9.81)
                jrasrinteronly101z = -1*(np.interp(n_timesinterpz, n_timesz, jrasrinteronlyz)) / (modelmass * 9.81)
                jrasrallonly101z = -1*(np.interp(n_timesinterpz, n_timesz, jrasrallz)) / (modelmass * 9.81)
                jrasrreserveonly101z = -1*(np.interp(n_timesinterpz, n_timesz, jrasrreserveonlyz)) / (modelmass * 9.81)
                jrasrnoneonly101z = -1*(np.interp(n_timesinterpz, n_timesz, jrasrnonez)) / (modelmass * 9.81)

                # data for z: medial-lateral forces
                ninterseg_combine_z[spot, :] = jrasrinteronly101z[:-2]
                nquads_combine_z[spot, :] = jrasrquadsonly101z[:-2]
                nhams_combine_z[spot,:] = jrasrhamsonly101z[:-2]
                ngas_combine_z[spot,:] = jrasrgasonly101z[:-2]
                ntfl_combine_z[spot,:] = jrasrtflonly101z[:-2]
                nall_combine_z[spot,:] = jrasrallonly101z[:-2]
                nreserve_combine_z[spot,:] = jrasrreserveonly101z[:-2]
                nnone_combine_z[spot,:] = jrasrnoneonly101z[:-2]
                ## increase the spot - count of trials
                spot += 1



    ########################################################################################
    # plotting for activations, moments, etc. muscle insights for natural and exotendon
    ########################################################################################

    # if indresults: 
    #     # tweak the results directory to print out in an individual folder for subject. 
    #     analyzedir = os.path.join(analyzedir, welksubjects[0])
    #     print(analyzedir)
    
    # # create a figure for the muscle activations for natural and exotendon
    # fig1, ax1 = plt.subplots(5, 8, figsize=(20, 12), dpi=500)
    # muscles = list(muscleacts_nat.keys())
    # for i, muscle in enumerate(muscles):
    #     row = i // 8
    #     col = i % 8
    #     xnat = np.linspace(0, 100, len(muscleacts_nat[muscle]))
    #     xexo = np.linspace(0, 100, len(muscleacts_exo[muscle]))
        
    #     ax1[row, col].plot(xnat, muscleacts_nat[muscle], label='natural', color=ncolor, linestyle='--', alpha=0.2)
    #     ax1[row, col].plot(xexo, muscleacts_exo[muscle], label='exotendon', color=ecolor, linestyle='--', alpha=0.2)
    #     # not plot the averages for all the natural and all the exo
    #     ax1[row,col].plot(xnat, np.mean(muscleacts_nat[muscle],1), color=ncolor, linewidth=2, label='natural_avg')
    #     ax1[row,col].plot(xexo, np.mean(muscleacts_exo[muscle],1), color=ecolor, linewidth=2, label='exotendon_avg')    
    #     # formatting
    #     ax1[row, col].set_xlabel('% Gait cycle', fontsize=8)
    #     ax1[row, col].set_ylabel('Activation', fontsize=8)
    #     # split the string to get the name for the plot
    #     musclename = muscle.split('_r')[0][10:]
    #     ax1[row, col].set_title(musclename, fontsize=8)
    
    # handles, labels = ax1[0, 0].get_legend_handles_labels()
    # # fig1.legend(handles, labels, loc='upper right')
    # fig1.tight_layout()
    # plt.savefig(os.path.join(analyzedir, 'muscleactivations_' + whichleg + '.png'))


    # # create a figure for the joint moments for natural and exotendon
    # fig2, ax2 = plt.subplots(3, 8, figsize=(20, 8), dpi=500)
    # joints = list(moments_nat.keys())
    # idjoints = list(IDmoments_nat.keys())
    # for i, joint in enumerate(joints):
    #     row = i // 8
    #     col = i % 8
    #     xnat = np.linspace(0, 100, len(moments_nat[joint]))
    #     xexo = np.linspace(0, 100, len(moments_exo[joint]))
        
    #     ax2[row, col].plot(xnat, moments_nat[joint], label='natural', color=ncolor, linestyle='--', alpha=0.2)
    #     ax2[row, col].plot(xexo, moments_exo[joint], label='exotendon', color=ecolor, linestyle='--', alpha=0.2)
    #     # now the means of all the moments
    #     ax2[row, col].plot(xnat, np.mean(moments_nat[joint],1), color=ncolor, linewidth=2, label='natural_avg')
    #     ax2[row, col].plot(xexo, np.mean(moments_exo[joint],1), color=ecolor, linewidth=2, label='exotendon_avg')
    #     # now want the ID moments for comparison
    #     if joint in IDmoments_nat:
    #         # ax2[row, col].plot(xnat, IDmoments_nat[joint], label='ID natural', color='black', linestyle='--', alpha=0.8)
    #         # ax2[row, col].plot(xexo, IDmoments_exo[joint], label='ID exotendon', color='grey', linestyle='--', alpha=0.8)
    #         ax2[row, col].plot(xnat, np.mean(IDmoments_nat[joint],1), label='ID natural', color='black', linestyle='--', alpha=0.8)
    #         ax2[row, col].plot(xexo, np.mean(IDmoments_exo[joint],1), label='ID exotendon', color='grey', linestyle='--', alpha=0.8)
    #     elif 'pelvis_tx' in joint:
    #         ax2[row, col].plot(xnat, np.mean(IDmoments_nat['pelvis_tx_force'],1), label='ID natural', color='black', linestyle='--', alpha=0.8)
    #         ax2[row, col].plot(xexo, np.mean(IDmoments_exo['pelvis_tx_force'],1), label='ID exotendon', color='grey', linestyle='--', alpha=0.8)
    #     elif 'pelvis_ty' in joint:
    #         ax2[row, col].plot(xnat, np.mean(IDmoments_nat['pelvis_ty_force'],1), label='ID natural', color='black', linestyle='--', alpha=0.8)
    #         ax2[row, col].plot(xexo, np.mean(IDmoments_exo['pelvis_ty_force'],1), label='ID exotendon', color='grey', linestyle='--', alpha=0.8)
    #     elif 'pelvis_tz' in joint:
    #         ax2[row, col].plot(xnat, np.mean(IDmoments_nat['pelvis_tz_force'],1), label='ID natural', color='black', linestyle='--', alpha=0.8)
    #         ax2[row, col].plot(xexo, np.mean(IDmoments_exo['pelvis_tz_force'],1), label='ID exotendon', color='grey', linestyle='--', alpha=0.8)
    #     # formatting
    #     ax2[row, col].set_xlabel('% Gait cycle', fontsize=8)
    #     ax2[row, col].set_ylabel('Moment (Nm/kg)', fontsize=8)
    #     ax2[row, col].set_title(joint, fontsize=8)
    # # add a legend to the last subplot
    # handles, labels = ax2[0, 0].get_legend_handles_labels()
    # ax2[-1,-1].axis('off')  # Hide the last subplot
    # ax2[-1,-1].legend(handles, labels, loc='center', fontsize=8)
    # # fig2.legend(handles, labels, loc='upper right')
    # # ax2[0, 0].legend(loc='upper right', fontsize=8)
    # fig2.tight_layout()
    # plt.savefig(analyzedir + '\\jointmoments_' + whichleg + '.png')


    # # now create a figure for the muscle passive forces for natural and exotendon
    # fig3, ax3 = plt.subplots(5, 8, figsize=(20, 12), dpi=500)
    # muscles = list(passiveforces_nat.keys())
    # for i, muscle in enumerate(muscles):
    #     row = i // 8
    #     col = i % 8
    #     xnat = np.linspace(0, 100, len(passiveforces_nat[muscle]))
    #     xexo = np.linspace(0, 100, len(passiveforces_exo[muscle]))
        
    #     ax3[row, col].plot(xnat, passiveforces_nat[muscle], label='natural', color=ncolor, linestyle='--', alpha=0.2)
    #     ax3[row, col].plot(xexo, passiveforces_exo[muscle], label='exotendon', color=ecolor, linestyle='--', alpha=0.2)
    #     # now the averages
    #     ax3[row, col].plot(xnat, np.mean(passiveforces_nat[muscle],1), color=ncolor, linewidth=2, label='natural_avg')
    #     ax3[row, col].plot(xexo, np.mean(passiveforces_exo[muscle],1), color=ecolor, linewidth=2, label='exotendon_avg')
    #     # formatting
    #     ax3[row, col].set_xlabel('% Gait cycle', fontsize=8)
    #     ax3[row, col].set_ylabel('Passive Force (N)', fontsize=8)
    #     # split the strings so that the names are readable
    #     musclename = muscle.split('_r')[0][10:]
    #     ax3[row, col].set_title(musclename, fontsize=8)

    # handles, labels = ax3[0, 0].get_legend_handles_labels()
    # # fig3.legend(handles, labels, loc='upper right')
    # fig3.tight_layout()
    # plt.savefig(analyzedir + '\\passiveforces_' + whichleg + '.png')


    # # now the figure but for the muscle active forces in activeforces_nat and activeforces_exo
    # fig4, ax4 = plt.subplots(5, 8, figsize=(20, 12), dpi=500)
    # muscles = list(activeforces_nat.keys())
    # for i, muscle in enumerate(muscles):
    #     row = i // 8
    #     col = i % 8
    #     xnat = np.linspace(0, 100, len(activeforces_nat[muscle]))
    #     xexo = np.linspace(0, 100, len(activeforces_exo[muscle]))
        
    #     ax4[row, col].plot(xnat, activeforces_nat[muscle], label='natural', color=ncolor, linestyle='--', alpha=0.2)
    #     ax4[row, col].plot(xexo, activeforces_exo[muscle], label='exotendon', color=ecolor, linestyle='--', alpha=0.2)
    #     # now the averages
    #     ax4[row, col].plot(xnat, np.mean(activeforces_nat[muscle],1), color=ncolor, linewidth=2, label='natural_avg')
    #     ax4[row, col].plot(xexo, np.mean(activeforces_exo[muscle],1), color=ecolor, linewidth=2, label='exotendon_avg')
    #     # formatting
    #     ax4[row, col].set_xlabel('% Gait cycle', fontsize=8)
    #     ax4[row, col].set_ylabel('Active Force (N)', fontsize=8)
    #     # split the string to get the name for the plot
    #     musclename = muscle.split('_r')[0][10:]
    #     ax4[row, col].set_title(musclename, fontsize=8)

    # handles, labels = ax4[0, 0].get_legend_handles_labels()
    # # fig4.legend(handles, labels, loc='upper right')
    # fig4.tight_layout()
    # plt.savefig(analyzedir + '\\activeforces_' + whichleg + '.png')


    # # now the total forces
    # fig5, ax5 = plt.subplots(5, 8, figsize=(20, 12), dpi=500)
    # muscles = list(totalforces_nat.keys())
    # for i, muscle in enumerate(muscles):
    #     row = i // 8
    #     col = i % 8
    #     xnat = np.linspace(0, 100, len(totalforces_nat[muscle]))
    #     xexo = np.linspace(0, 100, len(totalforces_exo[muscle]))
        
    #     ax5[row, col].plot(xnat, totalforces_nat[muscle], label='natural', color=ncolor, linestyle='--', alpha=0.2)
    #     ax5[row, col].plot(xexo, totalforces_exo[muscle], label='exotendon', color=ecolor, linestyle='--', alpha=0.2)
    #     # now the averages
    #     ax5[row, col].plot(xnat, np.mean(totalforces_nat[muscle],1), color=ncolor, linewidth=2, label='natural_avg')
    #     ax5[row, col].plot(xexo, np.mean(totalforces_exo[muscle],1), color=ecolor, linewidth=2, label='exotendon_avg')
    #     # formatting
    #     ax5[row, col].set_xlabel('% Gait cycle', fontsize=8)
    #     ax5[row, col].set_ylabel('Total Force (N)', fontsize=8)
    #     # split the string to get the name for the plot
    #     musclename = muscle.split('_r')[0][10:]
    #     ax5[row, col].set_title(musclename, fontsize=8)

    # handles, labels = ax5[0, 0].get_legend_handles_labels()
    # # fig5.legend(handles, labels, loc='upper right')
    # fig5.tight_layout()
    # plt.savefig(analyzedir + '\\totalforces_' + whichleg + '.png')

    # plt.show()


    ################################################################################
    # now the knee contact forces 
    ################################################################################
    print('\n\nNow moving onto the knee contact forces\n')

    # tryin to implement everything below into a plotting function... that way we can just call it for each of the components... ideally.
    datay = {
        'ninterseg_combine': ninterseg_combine,
        'einterseg_combine': einterseg_combine,
        'ntfl_combine': ntfl_combine,
        'etfl_combine': etfl_combine,
        'ngas_combine': ngas_combine,
        'egas_combine': egas_combine,
        'nhams_combine': nhams_combine,
        'ehams_combine': ehams_combine,
        'nquads_combine': nquads_combine,
        'equads_combine': equads_combine,
        'nall_combine': nall_combine,
        'eall_combine': eall_combine,
    }
    plotKneeContactForce('Vertical', analyzedir, welksubjects, ncolor, ecolor, n_timespercent101, e_timespercent101, datay)
    
    datax = {
        'ninterseg_combine': ninterseg_combine_x,
        'einterseg_combine': einterseg_combine_x,
        'ntfl_combine': ntfl_combine_x,
        'etfl_combine': etfl_combine_x,
        'ngas_combine': ngas_combine_x,
        'egas_combine': egas_combine_x,
        'nhams_combine': nhams_combine_x,
        'ehams_combine': ehams_combine_x,
        'nquads_combine': nquads_combine_x,
        'equads_combine': equads_combine_x,
        'nall_combine': nall_combine_x,
        'eall_combine': eall_combine_x,
    }
    # plotKneeContactForce('Anterior-Posterior', analyzedir, welksubjects, ncolor, ecolor, n_timespercent101, e_timespercent101, datax)

    dataz = {
        'ninterseg_combine': ninterseg_combine_z,
        'einterseg_combine': einterseg_combine_z,
        'ntfl_combine': ntfl_combine_z,
        'etfl_combine': etfl_combine_z,
        'ngas_combine': ngas_combine_z,
        'egas_combine': egas_combine_z,
        'nhams_combine': nhams_combine_z,
        'ehams_combine': ehams_combine_z,
        'nquads_combine': nquads_combine_z,
        'equads_combine': equads_combine_z,
        'nall_combine': nall_combine_z,
        'eall_combine': eall_combine_z,
    }
    # plotKneeContactForce('Medial-Lateral', analyzedir, welksubjects, ncolor, ecolor, n_timespercent101, e_timespercent101, dataz)

    # here is where I am going to try and combine the x and z into a shear force single value. 
    datashear = {}
    for each in datax.keys(): 
        datashear[each] = np.sqrt((datax[each]**2) + (dataz[each]**2))
    plotKneeContactForce('Shear', analyzedir, welksubjects, ncolor, ecolor, n_timespercent101, e_timespercent101, datashear)

    # data for the resultant
    dataresultant = {}
    for each in datax.keys():
        dataresultant[each] = np.sqrt((datax[each]**2) + (datay[each]**2) + (dataz[each]**2))
    # plotKneeContactForce('Resultant', analyzedir, welksubjects, ncolor, ecolor, n_timespercent101, e_timespercent101, dataresultant)

    print('assuming we have working plotting functions, this should be the end of the script.')
    sys.exit()


    ###########################################################################
    # plotting for activations, moments, etc. muscle insights for natural and exotendon

    if indresults: 
        # tweak the results directory to print out in an individual folder for subject. 
        analyzedir = os.path.join(analyzedir, welksubjects[0])
        print(analyzedir)
    
    # create a figure for the muscle activations for natural and exotendon
    fig1, ax1 = plt.subplots(5, 8, figsize=(20, 12), dpi=500)
    muscles = list(muscleacts_nat.keys())
    for i, muscle in enumerate(muscles):
        row = i // 8
        col = i % 8
        xnat = np.linspace(0, 100, len(muscleacts_nat[muscle]))
        xexo = np.linspace(0, 100, len(muscleacts_exo[muscle]))
        
        ax1[row, col].plot(xnat, muscleacts_nat[muscle], label='natural', color=ncolor, linestyle='--', alpha=0.2)
        ax1[row, col].plot(xexo, muscleacts_exo[muscle], label='exotendon', color=ecolor, linestyle='--', alpha=0.2)
        # not plot the averages for all the natural and all the exo
        ax1[row,col].plot(xnat, np.mean(muscleacts_nat[muscle],1), color=ncolor, linewidth=2, label='natural_avg')
        ax1[row,col].plot(xexo, np.mean(muscleacts_exo[muscle],1), color=ecolor, linewidth=2, label='exotendon_avg')    
        # formatting
        ax1[row, col].set_xlabel('% Gait cycle', fontsize=8)
        ax1[row, col].set_ylabel('Activation', fontsize=8)
        # split the string to get the name for the plot
        musclename = muscle.split('_r')[0][10:]
        ax1[row, col].set_title(musclename, fontsize=8)
    
    handles, labels = ax1[0, 0].get_legend_handles_labels()
    # fig1.legend(handles, labels, loc='upper right')
    fig1.tight_layout()
    plt.savefig(os.path.join(analyzedir, 'muscleactivations_' + whichleg + '.png'))


    # create a figure for the joint moments for natural and exotendon
    fig2, ax2 = plt.subplots(3, 8, figsize=(20, 8), dpi=500)
    joints = list(moments_nat.keys())
    for i, joint in enumerate(joints):
        row = i // 8
        col = i % 8
        xnat = np.linspace(0, 100, len(moments_nat[joint]))
        xexo = np.linspace(0, 100, len(moments_exo[joint]))
        
        ax2[row, col].plot(xnat, moments_nat[joint], label='natural', color=ncolor, linestyle='--', alpha=0.2)
        ax2[row, col].plot(xexo, moments_exo[joint], label='exotendon', color=ecolor, linestyle='--', alpha=0.2)
        # now the means of all the moments
        ax2[row, col].plot(xnat, np.mean(moments_nat[joint],1), color=ncolor, linewidth=2, label='natural_avg')
        ax2[row, col].plot(xexo, np.mean(moments_exo[joint],1), color=ecolor, linewidth=2, label='exotendon_avg')
        # formatting
        ax2[row, col].set_xlabel('% Gait cycle', fontsize=8)
        ax2[row, col].set_ylabel('Moment (Nm/kg)', fontsize=8)
        ax2[row, col].set_title(joint, fontsize=8)

    handles, labels = ax2[0, 0].get_legend_handles_labels()
    # fig2.legend(handles, labels, loc='upper right')
    fig2.tight_layout()
    plt.savefig(analyzedir + '\\jointmoments_' + whichleg + '.png')


    # now create a figure for the muscle passive forces for natural and exotendon
    fig3, ax3 = plt.subplots(5, 8, figsize=(20, 12), dpi=500)
    muscles = list(passiveforces_nat.keys())
    for i, muscle in enumerate(muscles):
        row = i // 8
        col = i % 8
        xnat = np.linspace(0, 100, len(passiveforces_nat[muscle]))
        xexo = np.linspace(0, 100, len(passiveforces_exo[muscle]))
        
        ax3[row, col].plot(xnat, passiveforces_nat[muscle], label='natural', color=ncolor, linestyle='--', alpha=0.2)
        ax3[row, col].plot(xexo, passiveforces_exo[muscle], label='exotendon', color=ecolor, linestyle='--', alpha=0.2)
        # now the averages
        ax3[row, col].plot(xnat, np.mean(passiveforces_nat[muscle],1), color=ncolor, linewidth=2, label='natural_avg')
        ax3[row, col].plot(xexo, np.mean(passiveforces_exo[muscle],1), color=ecolor, linewidth=2, label='exotendon_avg')
        # formatting
        ax3[row, col].set_xlabel('% Gait cycle', fontsize=8)
        ax3[row, col].set_ylabel('Passive Force (N)', fontsize=8)
        # split the strings so that the names are readable
        musclename = muscle.split('_r')[0][10:]
        ax3[row, col].set_title(musclename, fontsize=8)

    handles, labels = ax3[0, 0].get_legend_handles_labels()
    # fig3.legend(handles, labels, loc='upper right')
    fig3.tight_layout()
    plt.savefig(analyzedir + '\\passiveforces_' + whichleg + '.png')


    # now the figure but for the muscle active forces in activeforces_nat and activeforces_exo
    fig4, ax4 = plt.subplots(5, 8, figsize=(20, 12), dpi=500)
    muscles = list(activeforces_nat.keys())
    for i, muscle in enumerate(muscles):
        row = i // 8
        col = i % 8
        xnat = np.linspace(0, 100, len(activeforces_nat[muscle]))
        xexo = np.linspace(0, 100, len(activeforces_exo[muscle]))
        
        ax4[row, col].plot(xnat, activeforces_nat[muscle], label='natural', color=ncolor, linestyle='--', alpha=0.2)
        ax4[row, col].plot(xexo, activeforces_exo[muscle], label='exotendon', color=ecolor, linestyle='--', alpha=0.2)
        # now the averages
        ax4[row, col].plot(xnat, np.mean(activeforces_nat[muscle],1), color=ncolor, linewidth=2, label='natural_avg')
        ax4[row, col].plot(xexo, np.mean(activeforces_exo[muscle],1), color=ecolor, linewidth=2, label='exotendon_avg')
        # formatting
        ax4[row, col].set_xlabel('% Gait cycle', fontsize=8)
        ax4[row, col].set_ylabel('Active Force (N)', fontsize=8)
        # split the string to get the name for the plot
        musclename = muscle.split('_r')[0][10:]
        ax4[row, col].set_title(musclename, fontsize=8)

    handles, labels = ax4[0, 0].get_legend_handles_labels()
    # fig4.legend(handles, labels, loc='upper right')
    fig4.tight_layout()
    plt.savefig(analyzedir + '\\activeforces_' + whichleg + '.png')


    # now the total forces
    fig5, ax5 = plt.subplots(5, 8, figsize=(20, 12), dpi=500)
    muscles = list(totalforces_nat.keys())
    for i, muscle in enumerate(muscles):
        row = i // 8
        col = i % 8
        xnat = np.linspace(0, 100, len(totalforces_nat[muscle]))
        xexo = np.linspace(0, 100, len(totalforces_exo[muscle]))
        
        ax5[row, col].plot(xnat, totalforces_nat[muscle], label='natural', color=ncolor, linestyle='--', alpha=0.2)
        ax5[row, col].plot(xexo, totalforces_exo[muscle], label='exotendon', color=ecolor, linestyle='--', alpha=0.2)
        # now the averages
        ax5[row, col].plot(xnat, np.mean(totalforces_nat[muscle],1), color=ncolor, linewidth=2, label='natural_avg')
        ax5[row, col].plot(xexo, np.mean(totalforces_exo[muscle],1), color=ecolor, linewidth=2, label='exotendon_avg')
        # formatting
        ax5[row, col].set_xlabel('% Gait cycle', fontsize=8)
        ax5[row, col].set_ylabel('Total Force (N)', fontsize=8)
        # split the string to get the name for the plot
        musclename = muscle.split('_r')[0][10:]
        ax5[row, col].set_title(musclename, fontsize=8)

    handles, labels = ax5[0, 0].get_legend_handles_labels()
    # fig5.legend(handles, labels, loc='upper right')
    fig5.tight_layout()
    plt.savefig(analyzedir + '\\totalforces_' + whichleg + '.png')

    # plt.show()


    ###########################################################################
    # plotting for the joint contacts - natural and exotendon


    if len(welksubjects) == 1:
        ncolors = ['#fee0b6', '#fdae6b', '#fd8d3c', '#e66101']
        ecolors = ['#c7eae5', '#80cdc1', '#35978f', '#01665e']
    

    ###########################################################################
    # figure: segmenting all the muscles between exo and nat 
    ## really nice figure for seeing what is going on, but likely not going to 
    ## be in the paper...
    fig9, ax9 = plt.subplots(2, 4, figsize=(14, 6))# , dpi=300)
    ax9 = ax9.flatten()
    # intersegmental forces average
    for i, curve in enumerate(ninterseg_combine):
         ax9[0].plot(n_timespercent101, curve, label='natural'+str(i), color=ncolors[i] if len(welksubjects) == 1 else ncolor)
    for i, curve in enumerate(einterseg_combine):
        ax9[0].plot(e_timespercent101, curve, label='exotendon'+str(i), color=ecolors[i] if len(welksubjects) == 1 else ecolor)
    ax9[0].set_xlabel('% Gait cycle')
    ax9[0].set_ylabel('Force (BW)')
    ax9[0].set_title('intersegmental')
    # ax9[0].legend()
    
    # tfl forces
    for i, curve in enumerate(ntfl_combine):
        ax9[1].plot(n_timespercent101, curve, label='natural'+str(i), color=ncolors[i] if len(welksubjects) == 1 else ncolor)
    for i, curve in enumerate(etfl_combine):
        ax9[1].plot(e_timespercent101, curve, label='exotendon'+str(i), color=ecolors[i] if len(welksubjects) == 1 else ecolor)
    # ax9[1].plot(n_timespercent101, ntfl_combine, label='natural', color=ncolor)
    # ax9[1].plot(e_timespercent101, etfl_combine, label='exotendon', color=ecolor)
    ax9[1].set_xlabel('% Gait cycle')
    # ax9[1].set_ylabel('Force (BW)')
    # ax9[1].legend()
    ax9[1].set_title('tfl')

    # gastroc forces
    for i, curve in enumerate(ngas_combine):
        ax9[2].plot(n_timespercent101, curve, label='natural'+str(i), color=ncolors[i] if len(welksubjects) == 1 else ncolor)
    for i, curve in enumerate(egas_combine):
        ax9[2].plot(e_timespercent101, curve, label='exotendon'+str(i), color=ecolors[i] if len(welksubjects) == 1 else ecolor)
    # ax9[2].plot(n_timespercent101, ngas_combine, label='natural', color=ncolor)
    # ax9[2].plot(e_timespercent101, egas_combine, label='exotendon', color=ecolor)
    ax9[2].set_xlabel('% Gait cycle')
    # ax9[2].set_ylabel('Force (BW)')
    # ax9[2].legend()
    ax9[2].set_title('gastroc')
    
    # hamstring forces
    for i, curve in enumerate(nhams_combine):
        ax9[3].plot(n_timespercent101, curve, label='natural'+str(i), color=ncolors[i] if len(welksubjects) == 1 else ncolor)
    for i, curve in enumerate(ehams_combine):
        ax9[3].plot(e_timespercent101, curve, label='exotendon'+str(i), color=ecolors[i] if len(welksubjects) == 1 else ecolor)
    # ax9[3].plot(n_timespercent101, nhams_combine, label='natural', color=ncolor)
    # ax9[3].plot(e_timespercent101, ehams_combine, label='exotendon', color=ecolor)
    ax9[3].set_xlabel('% Gait cycle')
    # ax9[3].set_ylabel('Force (BW)')
    # ax9[3].legend()
    ax9[3].set_title('hamstrings')

    # quads forces
    for i, curve in enumerate(nquads_combine):
        ax9[4].plot(n_timespercent101, curve, label='natural'+str(i), color=ncolors[i] if len(welksubjects) == 1 else ncolor)
    for i, curve in enumerate(equads_combine):
        ax9[4].plot(e_timespercent101, curve, label='exotendon'+str(i), color=ecolors[i] if len(welksubjects) == 1 else ecolor)
    # ax9[4].plot(n_timespercent101, nquads_combine, label='natural', color=ncolor)
    # ax9[4].plot(e_timespercent101, equads_combine, label='exotendon', color=ecolor)
    ax9[4].set_xlabel('% Gait cycle')
    # ax9[4].set_ylabel('Force (BW)')
    # ax9[4].legend()
    ax9[4].set_title('quadriceps')

    # # reserve forces
    # ax9[5].plot(n_timespercent101, nreserve_combine, label='natural', color=ncolor)
    # ax9[5].plot(e_timespercent101, ereserve_combine, label='exotendon', color=ecolor)
    # ax9[5].set_xlabel('% Gait cycle')
    # # ax9[5].set_ylabel('Force (BW)')
    # # ax9[5].legend()
    # ax9[5].set_title('reserves')
    
    # # added all forces
    # ax9[5].plot(n_timespercent101, nquads_combine+ nhams_combine+ ngas_combine+ ntfl_combine+ ninterseg_combine+ nreserve_combine, label='natural', color=ncolor)
    # ax9[5].plot(e_timespercent101, equads_combine+ ehams_combine+ egas_combine+ etfl_combine+ einterseg_combine+ ereserve_combine, label='exotendon', color=ecolor)
    # ax9[5].set_xlabel('% Gait cycle')
    # # ax9[6].set_ylabel('Force (BW)')
    # # ax9[6].legend()
    # ax9[5].set_title('Total Vertical Contact')

    # all forces from whole analysis
    for i, curve in enumerate(nall_combine):
        ax9[5].plot(n_timespercent101, curve, label='natural'+str(i), color=ncolors[i] if len(welksubjects) == 1 else ncolor)
    for i, curve in enumerate(eall_combine):
        ax9[5].plot(e_timespercent101, curve, label='exotendon'+str(i), color=ecolors[i] if len(welksubjects) == 1 else ecolor)
    # ax9[5].plot(n_timespercent101, nall_combine, label='natural', color=ncolor)
    # ax9[5].plot(e_timespercent101, eall_combine, label='exotendon', color=ecolor)
    ax9[5].set_xlabel('% Gait cycle')
    # ax9[5].set_ylabel('Force (BW)')
    ax9[5].set_title('Total Vertical Contact')
    # ax9[5].legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10)
    # Hide the last subplot and use it to display the legend   
    ax9[6].axis('off')
    handles, labels = ax9[5].get_legend_handles_labels()
    if len(welksubjects) == 1:
        ax9[6].legend(handles, labels, loc='center', fontsize=14)

    fig9.tight_layout()
    plt.savefig(analyzedir + '\\contact1_' + whichleg + '.png')

    ###########################################################################
    # figure: segmenting all the muscles between exo and nat 
    ## really nice figure for seeing what is going on, but likely not going to 
    ## be in the paper...
    fig10, ax10 = plt.subplots(1,7, figsize=(18,3), dpi=500)
    fontz = 16
    font_properties = {'fontsize': 16, 'fontfamily': 'serif', 'fontname': 'Times New Roman'}
    # tick_font_properties = {'fontfamily': 'serif', 'fontname': 'Times New Roman'}

    # all forces from whole analysis
    ax10[0].fill_between(n_timespercent101, np.mean(nall_combine, 0) - np.std(nall_combine, 0), np.mean(nall_combine, 0) + np.std(nall_combine, 0), color=ncolor, alpha=0.2)
    ax10[0].fill_between(e_timespercent101, np.mean(eall_combine, 0) - np.std(eall_combine, 0), np.mean(eall_combine, 0) + np.std(eall_combine, 0), color=ecolor, alpha=0.2)
    ax10[0].plot(n_timespercent101, np.mean(nall_combine,0), label='natural', color=ncolor)
    ax10[0].plot(e_timespercent101, np.mean(eall_combine,0), label='exotendon', color=ecolor)
    ax10[0].set_xlabel('% Gait cycle', **font_properties)
    ax10[0].set_ylabel('Force (BW)', **font_properties)
    # ax10[0].legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10)
    ax10[0].set_title('Total\nVertical Knee Contact', **font_properties)
    ax10[0].tick_params(axis='both', which='major', labelsize=fontz)

    # quads forces
    ax10[1].fill_between(n_timespercent101, np.mean(nquads_combine, 0) - np.std(nquads_combine, 0), np.mean(nquads_combine, 0) + np.std(nquads_combine, 0), color=ncolor, alpha=0.2)#, label='Natural St. Dev.')
    ax10[1].fill_between(e_timespercent101, np.mean(equads_combine, 0) - np.std(equads_combine, 0), np.mean(equads_combine, 0) + np.std(equads_combine, 0), color=ecolor, alpha=0.2)#, label='Exotendon St. Dev.')
    ax10[1].plot(n_timespercent101, np.mean(nquads_combine, 0), label='Natural (Mean \u00B1 Std.)', color=ncolor)
    ax10[1].plot(e_timespercent101, np.mean(equads_combine, 0), label='Exotendon (Mean \u00B1 Std.)', color=ecolor)
    ax10[1].set_xlabel('% Gait cycle', **font_properties)
    ax10[1].tick_params(axis='both', which='major', labelsize=fontz)
    # ax10[1].set_ylabel('Force (BW)')
    # ax10[1].legend()
    ax10[1].set_title('Contribution of\nQuadriceps', **font_properties)

    # gastroc forces
    ax10[2].fill_between(n_timespercent101, np.mean(ngas_combine, 0) - np.std(ngas_combine, 0), np.mean(ngas_combine, 0) + np.std(ngas_combine, 0), color=ncolor, alpha=0.2)
    ax10[2].fill_between(e_timespercent101, np.mean(egas_combine, 0) - np.std(egas_combine, 0), np.mean(egas_combine, 0) + np.std(egas_combine, 0), color=ecolor, alpha=0.2)
    ax10[2].plot(n_timespercent101, np.mean(ngas_combine, 0), label='natural', color=ncolor)
    ax10[2].plot(e_timespercent101, np.mean(egas_combine, 0), label='exotendon', color=ecolor)
    ax10[2].set_xlabel('% Gait cycle', **font_properties)
    ax10[2].tick_params(axis='both', which='major', labelsize=fontz)
    # ax10[2].set_ylabel('Force (BW)')
    # ax10[2].legend()
    ax10[2].set_title('Contribution of\nGastrocnemius', **font_properties)

    # hamstring forces
    ax10[3].fill_between(n_timespercent101, np.mean(nhams_combine, 0) - np.std(nhams_combine, 0), np.mean(nhams_combine, 0) + np.std(nhams_combine, 0), color=ncolor, alpha=0.2)
    ax10[3].fill_between(e_timespercent101, np.mean(ehams_combine, 0) - np.std(ehams_combine, 0), np.mean(ehams_combine, 0) + np.std(ehams_combine, 0), color=ecolor, alpha=0.2)
    ax10[3].plot(n_timespercent101, np.mean(nhams_combine, 0), label='natural', color=ncolor)
    ax10[3].plot(e_timespercent101, np.mean(ehams_combine, 0), label='exotendon', color=ecolor)
    ax10[3].set_xlabel('% Gait cycle', **font_properties)
    ax10[3].tick_params(axis='both', which='major', labelsize=fontz)
    # ax10[3].set_ylabel('Force (BW)')
    # ax10[3].legend()
    ax10[3].set_title('Contribution of\nHamstrings', **font_properties)

    # tfl forces
    ax10[4].fill_between(n_timespercent101, np.mean(ntfl_combine, 0) - np.std(ntfl_combine, 0), np.mean(ntfl_combine, 0) + np.std(ntfl_combine, 0), color=ncolor, alpha=0.2)
    ax10[4].fill_between(e_timespercent101, np.mean(etfl_combine, 0) - np.std(etfl_combine, 0), np.mean(etfl_combine, 0) + np.std(etfl_combine, 0), color=ecolor, alpha=0.2)
    ax10[4].plot(n_timespercent101, np.mean(ntfl_combine, 0), label='natural', color=ncolor)
    ax10[4].plot(e_timespercent101, np.mean(etfl_combine, 0), label='exotendon', color=ecolor)
    ax10[4].set_xlabel('% Gait cycle', **font_properties)
    ax10[4].tick_params(axis='both', which='major', labelsize=fontz)
    # ax10[4].set_ylabel('Force (BW)')
    # ax10[4].legend()
    ax10[4].set_title('Contribution of\nTensor Fascia Latae', **font_properties)

    # intersegmental forces average
    ax10[5].fill_between(n_timespercent101, np.mean(ninterseg_combine, 0) - np.std(ninterseg_combine, 0), np.mean(ninterseg_combine, 0) + np.std(ninterseg_combine, 0), color=ncolor, alpha=0.2)
    ax10[5].fill_between(e_timespercent101, np.mean(einterseg_combine, 0) - np.std(einterseg_combine, 0), np.mean(einterseg_combine, 0) + np.std(einterseg_combine, 0), color=ecolor, alpha=0.2)
    ax10[5].plot(n_timespercent101, np.mean(ninterseg_combine, 0), label='natural', color=ncolor)
    ax10[5].plot(e_timespercent101, np.mean(einterseg_combine, 0), label='exotendon', color=ecolor)
    ax10[5].set_xlabel('% Gait cycle', **font_properties)
    # ax10[5].set_ylabel('Force (BW)', **font_properties)
    ax10[5].set_title('Contribution of\nIntersegmental Forces', **font_properties)
    ax10[5].tick_params(axis='both', which='major', labelsize=fontz)    

    # # reserve forces
    # ax10[5].plot(n_timespercent101, np.mean(nreserve_combine, 0), label='natural', color=ncolor)
    # ax10[5].plot(e_timespercent101, np.mean(ereserve_combine, 0), label='exotendon', color=ecolor)
    # ax10[5].set_xlabel('% Gait cycle', **font_properties)
    # # ax10[5].set_ylabel('Force (BW)')
    # # ax10[5].legend()
    # ax10[5].set_title('reserves')

    # # added all forces
    # ax10[5].plot(n_timespercent101, np.mean(nquads_combine,0) + np.mean(nhams_combine,0) + np.mean(ngas_combine,0) + np.mean(ntfl_combine,0) + np.mean(ninterseg_combine,0) + np.mean(nreserve_combine,0), label='natural', color=ncolor)
    # ax10[5].plot(e_timespercent101, np.mean(equads_combine,0) + np.mean(ehams_combine,0) + np.mean(egas_combine,0) + np.mean(etfl_combine,0) + np.mean(einterseg_combine,0) + np.mean(ereserve_combine,0), label='exotendon', color=ecolor)
    # ax10[5].set_xlabel('% Gait cycle', **font_properties)
    # # ax10[6].set_ylabel('Force (BW)')
    # # ax10[6].legend()
    # ax10[5].set_title('Total Vertical Contact')

    # Hide the last subplot and use it to display the legend
    ax10[6].axis('off')
    handles, labels = ax10[1].get_legend_handles_labels()
    # ax10[6].legend(handles, labels, loc='center', fontsize=fontz)
    fig10.tight_layout(pad=2.0, w_pad=0.5, h_pad=1.0)
    plt.savefig(analyzedir + '\\contact2_' + whichleg + '.png')

    ###########################################################################
    ### figure: differences between conditions for each and all
    #### this is a simplified look at the same as above. We can see nice stuff, 
    #### but likely not going to be in the paper. 
    fig11, ax11 = plt.subplots(1,7, figsize=(14,3))#, dpi=300)
    # intersegmental forces average
    ax11[0].plot(n_timespercent101, np.mean(einterseg_combine, 0) - np.mean(ninterseg_combine, 0))
    ax11[0].set_xlabel('% Gait cycle')
    ax11[0].set_ylabel('Force (BW)')
    ax11[0].set_title('intersegmental')
    # ax11[0].legend()
    # tfl forces
    ax11[1].plot(n_timespercent101, np.mean(etfl_combine, 0) - np.mean(ntfl_combine, 0))
    ax11[1].set_xlabel('% Gait cycle')
    # ax11[1].set_ylabel('Force (BW)')
    # ax11[1].legend()
    ax11[1].set_title('tfl')

    # gastroc forces
    ax11[2].plot(n_timespercent101, np.mean(egas_combine, 0) - np.mean(ngas_combine, 0))
    ax11[2].set_xlabel('% Gait cycle')
    # ax11[2].set_ylabel('Force (BW)')
    # ax11[2].legend()
    ax11[2].set_title('gastroc')

    # hamstring forces
    ax11[3].plot(n_timespercent101, np.mean(ehams_combine, 0) - np.mean(nhams_combine, 0))
    ax11[3].set_xlabel('% Gait cycle')
    # ax11[3].set_ylabel('Force (BW)')
    # ax11[3].legend()
    ax11[3].set_title('hamstrings')

    # quads forces
    ax11[4].plot(n_timespercent101, np.mean(equads_combine, 0) - np.mean(nquads_combine, 0))
    ax11[4].set_xlabel('% Gait cycle')
    # ax11[4].set_ylabel('Force (BW)')
    # ax11[4].legend()
    ax11[4].set_title('quadriceps')

    # # added all forces
    # ax11[5].plot(n_timespercent101, (np.mean(equads_combine,0) + np.mean(ehams_combine,0) + np.mean(egas_combine,0) + np.mean(etfl_combine,0) + np.mean(einterseg_combine,0)) - (np.mean(nquads_combine,0) + np.mean(nhams_combine,0) + np.mean(ngas_combine,0) + np.mean(ntfl_combine,0) + np.mean(ninterseg_combine,0)))
    # ax11[5].set_xlabel('% Gait cycle')
    # # ax11[5].set_ylabel('Force (BW)')
    # # ax11[5].legend()
    # ax11[5].set_title('Total Vertical Contact')

    # all forces from whole analysis
    ax11[5].plot(n_timespercent101, np.mean(eall_combine,0) - np.mean(nall_combine,0), label='Exotendon diff from Natural')
    ax11[5].set_xlabel('% Gait cycle')
    # ax11[5].set_ylabel('Force (BW)')
    # ax11[5].legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10)
    ax11[5].set_title('Total Vertical Contact')

    # Hide the last subplot and use it to display the legend
    ax11[6].axis('off')
    handles, labels = ax11[5].get_legend_handles_labels()
    ax11[6].legend(handles, labels, loc='center', fontsize=10)
    fig11.tight_layout()
    plt.savefig(analyzedir + '\\contact3_' + whichleg + '.png')
    
    
    ###########################################################################
    ### figure: changes in force segmented together on plot
    #### Possible paper figure for R3. 
    colors4 = ['#648FFF','#785EF0','#DC267F','#FE6100', '#FFB000']
    fig12, ax12 = plt.subplots(1, 2, figsize=(12,4.55), dpi=500)
    ax12[0].axhline(0, color='black', linewidth=1)
    # intersegmental forces average
    ax12[0].plot(n_timespercent101, np.mean(einterseg_combine, 0) - np.mean(ninterseg_combine, 0), label='Intersegmental', color='black', linewidth=3)
    # tfl forces
    ax12[0].plot(n_timespercent101, np.mean(etfl_combine, 0) - np.mean(ntfl_combine, 0), label='Tensor Fasciae Latae', color=colors4[0], linewidth=3)
    # gastroc forces
    ax12[0].plot(n_timespercent101, np.mean(egas_combine, 0) - np.mean(ngas_combine, 0), label='Gastrocnemius', color=colors4[1], linewidth=3)
    # hamstring forces
    ax12[0].plot(n_timespercent101, np.mean(ehams_combine, 0) - np.mean(nhams_combine, 0), label='Hamstrings', color=colors4[2], linewidth=3)
    # quads forces
    ax12[0].plot(n_timespercent101, np.mean(equads_combine, 0) - np.mean(nquads_combine, 0), label='Quadriceps', color=colors4[3], linewidth=3)
    # reserve forces
    # ax12[0].plot(n_timespercent101, np.mean(ereserve_combine, 0) - np.mean(nreserve_combine, 0), label='reserves', color=colors4[4])
    # added all forces
    # ax12[0].plot(n_timespercent101, (np.mean(equads_combine,0) + np.mean(ehams_combine,0) + np.mean(egas_combine,0) + np.mean(etfl_combine,0) + np.mean(einterseg_combine,0) + np.mean(ereserve_combine,0)) - (np.mean(nquads_combine,0) + np.mean(nhams_combine,0) + np.mean(ngas_combine,0) + np.mean(ntfl_combine,0) + np.mean(ninterseg_combine,0) + np.mean(nreserve_combine,0)), label='Total')
    # all forces from whole analysis
    ax12[0].plot(n_timespercent101, np.mean(eall_combine,0) - np.mean(nall_combine,0), label='Total Vertical Contact', linestyle='dashed', color='black', linewidth=3)
    ax12[0].set_xlabel('% Gait cycle', fontsize=16)
    ax12[0].set_ylabel('Vertical knee contact difference (BW)', fontsize=16)
    ax12[0].set_title('Exotendon change in contact force', fontsize=16)
    ax12[0].tick_params(axis='both', which='major', labelsize=16)
    # hide the second subplots and use it for the legend
    ax12[1].axis('off')
    handles, labels = ax12[0].get_legend_handles_labels()
    ax12[1].legend(handles, labels, loc='center', fontsize=16)
    fig12.tight_layout()
    fig12.savefig(analyzedir + '\\contact4_' + whichleg + '.png')

    # TODO: figure out why the difference in total and all added together. 
    # okay so not in how I am adding/averaging. has to be something in how the analysis is done between them.... am I missing something??
    
    fig13, ax13 = plt.subplots(1,3, figsize=(14,5))#, dpi=300)
    # intersegmental forces average - natural
    ax13[0].plot(n_timespercent101, np.mean(ninterseg_combine, 0), label='natural_interseg')
    ax13[0].plot(n_timespercent101, np.mean(ninterseg_combine + ntfl_combine, 0), label='natural_interseg + tfl')
    ax13[0].plot(n_timespercent101, np.mean(ninterseg_combine+ntfl_combine+nhams_combine, 0), label='natural_interseg+tfl+hams')
    ax13[0].plot(n_timespercent101, np.mean(ninterseg_combine+ntfl_combine+nhams_combine+ngas_combine, 0), label='natural_interseg+tfl+hams+gas')
    ax13[0].plot(n_timespercent101, np.mean(ninterseg_combine+ntfl_combine+nhams_combine+ngas_combine+nquads_combine, 0), label='natural_interseg+tfl+hams+gas+quads')
    # ax13[0].plot(n_timespercent101, np.mean(ninterseg_combine+ntfl_combine+nhams_combine+ngas_combine+nquads_combine+nreserve_combine, 0), label='natural_interseg+tfl+hams+gas+quads+reserve')
    ax13[0].plot(n_timespercent101, np.mean(nall_combine, 0), label='nat_all', linestyle='dotted')
    ax13[0].set_xlabel('% Gait cycle')
    ax13[0].set_ylabel('Force (BW)')
    ax13[0].set_title('natural')
    # ax13[0].legend()
    # intersegmental forces average - exotendon
    ax13[1].plot(e_timespercent101, np.mean(einterseg_combine, 0), label='exo_interseg')
    ax13[1].plot(e_timespercent101, np.mean(einterseg_combine+etfl_combine, 0), label='exo_interseg + tfl')
    ax13[1].plot(e_timespercent101, np.mean(einterseg_combine+etfl_combine+ehams_combine, 0), label='exo_interseg+tfl+hams')
    ax13[1].plot(e_timespercent101, np.mean(einterseg_combine+etfl_combine+ehams_combine+egas_combine, 0), label='exo_interseg+tfl+hams+gas')
    ax13[1].plot(e_timespercent101, np.mean(einterseg_combine+etfl_combine+ehams_combine+egas_combine+equads_combine, 0), label='exo_interseg+tfl+hams+gas+quads')
    # ax13[1].plot(e_timespercent101, np.mean(einterseg_combine+etfl_combine+ehams_combine+egas_combine+equads_combine+ereserve_combine, 0), label='exo_interseg+tfl+hams+gas+quads+reserve')
    ax13[1].plot(e_timespercent101, np.mean(eall_combine, 0), label='exo_all', linestyle='dotted')
    ax13[1].set_xlabel('% Gait cycle')
    ax13[1].set_ylabel('Force (BW)')
    ax13[1].set_title('exotendon')
    # ax13[1].legend()

    # Hide the last subplot and use it to display the legend
    ax13[2].axis('off')
    handles, labels = ax13[1].get_legend_handles_labels()
    ax13[2].legend(handles, labels, loc='center', fontsize=14)
    fig13.tight_layout()
    plt.savefig(analyzedir + '\\contact5_' + whichleg + '.png')
    


    ###########################################################################
    ### Figure: total pop average for right leg between nat and exo
    #### Figure in the paper R1... 
    figcon6, axcon6 = plt.subplots(1, 2, figsize=(12,4.55), dpi=500)
    axcon6[0].fill_between(n_timespercent101, np.mean(nall_combine,0)-np.std(nall_combine,0), np.mean(nall_combine,0)+np.std(nall_combine,0), color=ncolorlight, alpha=0.75)
    axcon6[0].fill_between(e_timespercent101, np.mean(eall_combine,0)-np.std(eall_combine,0), np.mean(eall_combine,0)+np.std(eall_combine,0), color=ecolorlight, alpha=0.75)
    axcon6[0].plot(n_timespercent101, np.mean(nall_combine,0), color=ncolor, label='Natural (Mean \u00B1 Std.)', linewidth=3)
    axcon6[0].plot(e_timespercent101, np.mean(eall_combine,0), color=ecolor, label='Exotendon (Mean \u00B1 Std.)', linewidth=3)
    axcon6[0].set_xlabel('% Gait cycle', fontsize=16)
    axcon6[0].set_ylabel('Vertical knee contact force (BW)', fontsize=16)
    axcon6[0].set_title('Total Vertical Contact force', fontsize=16)
    axcon6[0].tick_params(axis='both', which='major', labelsize=16)
    # hide the second subplots and use it for the legend
    axcon6[1].axis('off')
    handles, labels = axcon6[0].get_legend_handles_labels()
    axcon6[1].legend(handles, labels, loc='center', fontsize=16)
    figcon6.tight_layout()
    figcon6.savefig(analyzedir + '\\contact6_' + whichleg + '.png')


    ###########################################################################
    # output the data in an excel sheet
    # Create a dictionary to hold all the data
    data_dict = {
        'ninterseg_combine': ninterseg_combine,
        'einterseg_combine': einterseg_combine,
        'nquads_combine': nquads_combine,
        'equads_combine': equads_combine,
        'nhams_combine': nhams_combine,
        'ehams_combine': ehams_combine,
        'ngas_combine': ngas_combine,
        'egas_combine': egas_combine,
        'ntfl_combine': ntfl_combine,
        'etfl_combine': etfl_combine,
        'nall_combine': nall_combine,
        'eall_combine': eall_combine,
        'nreserve_combine': nreserve_combine,
        'ereserve_combine': ereserve_combine,
        'nnone_combine': nnone_combine,
        'enone_combine': enone_combine
    }

        # 'muscleacts_nat': muscleacts_nat,
        # 'muscleacts_exo': muscleacts_exo,
        # 'moments_nat': moments_nat,
        # 'moments_exo': moments_exo,
        # 'IDmoments_nat': IDmoments_nat,
        # 'IDmoments_exo': IDmoments_exo,
        # 'activeforces_nat': activeforces_nat,
        # 'activeforces_exo': activeforces_exo,
        # 'passiveforces_nat': passiveforces_nat,
        # 'passiveforces_exo': passiveforces_exo,
        # 'totalforces_nat': totalforces_nat,
        # 'totalforces_exo': totalforces_exo

    os.chdir(analyzedir)
    # Create a Pandas Excel writer using XlsxWriter as the engine
    with pd.ExcelWriter('analysis_results'+'_'+whichleg+'.xlsx', engine='xlsxwriter') as writer:
        for key, value in data_dict.items():
            print(key)
            if isinstance(value, dict):
                print('found a dict')
                pdb.set_trace()
                for sub_key, sub_value in value.items():
                    df = pd.DataFrame(sub_value)
                    df.to_excel(writer, sheet_name=f'{key}_{sub_key}')
            else:
                df = pd.DataFrame(value)
                df.to_excel(writer, sheet_name=key)


    pdb.set_trace()
    plt.show()
    print('\n\nbeyond this is the breakdown paper figures.')
    # sys.exit()

    # ###########################################################################
    # # stats for the data
    # # here is the stats for R1 - differences in the peak vertical JCF
    # mean_nall_combine = np.mean(nall_combine,0)
    # std_nall_combine = np.std(nall_combine,0)
    # mean_eall_combine = np.mean(eall_combine,0)
    # std_eall_combine = np.std(eall_combine,0)

    # # get the peaks first and then do the stats on them...
    # peaks_nall_combine = np.max(nall_combine,1)
    # idx_peaks_nall_combine = nall_combine.argmax(1)
    # peaks_eall_combine = np.max(eall_combine,1)
    # idx_peaks_eall_combine = eall_combine.argmax(1)

    # mean_peaks_nall_combine = np.mean(peaks_nall_combine)
    # mean_peaks_eall_combine = np.mean(peaks_eall_combine)
    # std_peaks_nall_combine = np.std(peaks_nall_combine)
    # std_peaks_eall_combine = np.std(peaks_eall_combine)

    # # differences in peaks
    # diff_peaks_all_combine = peaks_nall_combine - peaks_eall_combine
    # mean_diff_peaks_all_combine = np.mean(peaks_nall_combine - peaks_eall_combine)
    # std_diff_peaks_all_combine = np.std(peaks_nall_combine - peaks_eall_combine)

    # # shapiro test for normal on the peaks differences
    # res = scipy.stats.shapiro(diff_peaks_all_combine)
    
    # # t - test for peaks differences is in the excel sheet...
    
    # # percent change in peak
    # perc_diff_peaks_all_combine = (peaks_eall_combine - peaks_nall_combine) / peaks_nall_combine * 100
    # mean_perc_diff_peaks_all_combine = np.mean(perc_diff_peaks_all_combine)
    # std_perc_diff_peaks_all_combine = np.std(perc_diff_peaks_all_combine)
    

    # ###########################################################################
    # more polished figures

    ### Figure: total pop average for right leg for nat - segmented shaded
    plt.figure(dpi=300)
    plt.plot(n_timespercent101, np.mean(ninterseg_combine,0), label='intersegmental', color=ecolor)
    plt.fill_between(n_timespercent101, np.mean(ninterseg_combine,0), color=ecolorlight)
    plt.plot(n_timespercent101, np.mean(ntfl_combine,0) + np.mean(ninterseg_combine,0), label='inter + tfl', color=ncolor3)
    plt.fill_between(n_timespercent101, np.mean(ninterseg_combine,0), np.mean(ntfl_combine,0) + np.mean(ninterseg_combine,0), color=ncolor2)

    plt.plot(n_timespercent101, np.mean(ngas_combine,0) + np.mean(ntfl_combine,0) + np.mean(ninterseg_combine,0), label='inter + tfl + gas', color=ncolor5)
    plt.fill_between(n_timespercent101, np.mean(ntfl_combine,0) + np.mean(ninterseg_combine,0), np.mean(ngas_combine,0) + np.mean(ntfl_combine,0) + np.mean(ninterseg_combine,0), color=ncolor4)

    plt.plot(n_timespercent101, np.mean(nhams_combine,0) + np.mean(ngas_combine,0) + np.mean(ntfl_combine,0) + np.mean(ninterseg_combine,0), label='inter + tfl + gas + hams', color=ncolor7)
    plt.fill_between(n_timespercent101, np.mean(ngas_combine,0) + np.mean(ntfl_combine,0) + np.mean(ninterseg_combine,0), np.mean(nhams_combine,0) + np.mean(ngas_combine,0) + np.mean(ntfl_combine,0) + np.mean(ninterseg_combine,0), color=ncolor6)

    plt.plot(n_timespercent101, np.mean(nquads_combine,0) + np.mean(nhams_combine,0) + np.mean(ngas_combine,0) + np.mean(ntfl_combine,0) + np.mean(ninterseg_combine,0), label='inter + tfl + gas + hams + quads', color=ncolor9)
    plt.fill_between(n_timespercent101, np.mean(nhams_combine,0) + np.mean(ngas_combine,0) + np.mean(ntfl_combine,0) + np.mean(ninterseg_combine,0), np.mean(nquads_combine,0) + np.mean(nhams_combine,0) + np.mean(ngas_combine,0) + np.mean(ntfl_combine,0) + np.mean(ninterseg_combine,0), color=ncolor8)

    plt.legend()
    plt.ylabel('Vertical knee contact force (BW)', fontsize=16)
    plt.xlabel('% Gait cycle', fontsize=16)
    plt.title('Vertical knee contact force - Natural', fontsize=16)    
    plt.tick_params(axis='both', which='major', labelsize=16)
    plt.tight_layout()
    plt.savefig(analyzedir + '\\contact7_Natural_' + whichleg + '.png')
    plt.show()
    
    
    ###########################################################################
    ### Figure: total pop average for right leg for exo - segmented shaded
    plt.figure(dpi=300)
    plt.plot(e_timespercent101, np.mean(einterseg_combine,0), label='intersegmental', color=ecolor)
    plt.fill_between(e_timespercent101, np.mean(einterseg_combine,0), color=ecolorlight)

    plt.plot(e_timespercent101, np.mean(etfl_combine,0) + np.mean(einterseg_combine,0), label='inter + tfl', color=ncolor3)
    plt.fill_between(e_timespercent101, np.mean(einterseg_combine,0), np.mean(etfl_combine,0) + np.mean(einterseg_combine,0), color=ncolor2)

    plt.plot(e_timespercent101, np.mean(egas_combine,0) + np.mean(etfl_combine,0) + np.mean(einterseg_combine,0), label='inter + tfl + gas', color=ncolor5)
    plt.fill_between(e_timespercent101, np.mean(etfl_combine,0) + np.mean(einterseg_combine,0), np.mean(egas_combine,0) + np.mean(etfl_combine,0) + np.mean(einterseg_combine,0), color=ncolor4)

    plt.plot(e_timespercent101, np.mean(ehams_combine,0) + np.mean(egas_combine,0) + np.mean(etfl_combine,0) + np.mean(einterseg_combine,0), label='inter + tfl + gas + hams', color=ncolor7)
    plt.fill_between(e_timespercent101, np.mean(egas_combine,0) + np.mean(etfl_combine,0) + np.mean(einterseg_combine,0), np.mean(ehams_combine,0) + np.mean(egas_combine,0) + np.mean(etfl_combine,0) + np.mean(einterseg_combine,0),color=ncolor6)

    plt.plot(e_timespercent101, np.mean(equads_combine,0) + np.mean(ehams_combine,0) + np.mean(egas_combine,0) + np.mean(etfl_combine,0) + np.mean(einterseg_combine,0), label='inter + tfl + gas + hams + quads', color=ncolor9)
    plt.fill_between(e_timespercent101, np.mean(ehams_combine,0) + np.mean(egas_combine,0) + np.mean(etfl_combine,0) + np.mean(einterseg_combine,0), np.mean(equads_combine,0) + np.mean(ehams_combine,0) + np.mean(egas_combine,0) + np.mean(etfl_combine,0) + np.mean(einterseg_combine,0), color=ncolor8)

    plt.legend()
    plt.ylabel('Vertical knee contact force (BW)', fontsize=16)
    plt.xlabel('% Gait cycle', fontsize=16)
    plt.title('Vertical knee contact force - Exotendon', fontsize=16)
    plt.tick_params(axis='both', which='major', labelsize=16)
    plt.tight_layout()
    plt.savefig(analyzedir + '\\contact8_exo_' + whichleg + '.png')
    plt.show()

    
    ###########################################################################
    # stats for the individual muscle segmented contributions.... check out the excel sheet
    peaks_quads_natural = np.max(nquads_combine,1)
    peaks_gas_natural = np.max(ngas_combine,1)
    peaks_hams_natural = np.max(nhams_combine,1)
    peaks_tfl_natural = np.max(ntfl_combine,1)
    peaks_inter_natural = np.max(ninterseg_combine,1)
    peaks_reserve_natural = np.max(nreserve_combine,1)

    peaks_quads_exo = np.max(equads_combine,1)
    peaks_gas_exo = np.max(egas_combine,1)
    peaks_hams_exo = np.max(ehams_combine,1)
    peaks_tfl_exo = np.max(etfl_combine,1)
    peaks_inter_exo = np.max(einterseg_combine,1)
    peaks_reserve_exo = np.max(ereserve_combine,1)

    
    
    ###########################################################################
    # I don't think this is needed anymore, but keeping for now...

    # ### fig: subplots of nat and exo for each breakdown
    # fig3, ax3 = plt.subplots(1,6)
    # # full forces
    # # ax3[0].plot(n_timespercent101, n_avg_r.mean(0), label='n_full', color='orange')
    # # ax3[0].plot(e_timespercent101, e_avg_r.mean(0), label='e_full', color='purple')
    # # intersegmental forces
    # ax3[1].plot(n_timespercent101, ninterseg_combine.mean(0), label='n_interseg', color='orange')
    # # ax3[1].plot(e_timespercent101, einterseg_combine.mean(0), label='e_interseg', color='purple')
    # ax3[1].set_xlabel('% gait cycle')
    # ax3[1].set_ylabel('BW - intersegmental')
    # # # quads forces
    # # ax3[2].plot(n_timespercent101, nquads_combine.mean(0), label='n_quads', color='orange')
    # # ax3[2].plot(e_timespercent101, equads_combine.mean(0), label='e_quads', color='purple')
    # # ax3[2].set_xlabel('% gait cycle')
    # # ax3[2].set_ylabel('BW - quads')
    # # # hams forces
    # # ax3[3].plot(n_timespercent101, nhams_combine.mean(0), label='n_hams', color='orange')
    # # ax3[3].plot(e_timespercent101, ehams_combine.mean(0), label='e_hams', color='purple')
    # # ax3[3].set_xlabel('% gait cycle')
    # # ax3[3].set_ylabel('BW - hams')
    # # # gastroc forces
    # # ax3[4].plot(n_timespercent101, ngas_combine.mean(0), label='n_gas', color='orange')
    # # ax3[4].plot(e_timespercent101, egas_combine.mean(0), label='e_gas', color='purple')
    # # ax3[4].set_xlabel('% gait cycle')
    # # ax3[4].set_ylabel('BW - gastroc')
    # # # tfl forces
    # # ax3[5].plot(n_timespercent101, ntfl_combine.mean(0), label='n_tfl', color='orange')
    # # ax3[5].plot(e_timespercent101, etfl_combine.mean(0), label='e_tfl', color='purple')
    # # ax3[5].set_xlabel('% gait cycle')
    # # ax3[5].set_ylabel('BW - tfl')

    # plt.legend()
    # plt.show()
    
    pdb.set_trace()
    sys.exit()


    # everything below here is probably old and not needed anymore... keeping for now. 
    '''
    ### fig: straight unedited exo (should include interseg for each) all
    plt.figure()
    # full
    # plt.plot(e_timespercent101, np.mean(n_avg_r, 0), label='n_full')
    # plt.plot(n_timespercent101, np.mean(e_avg_r, 0), label='e_full')
    # quads
    # plt.plot(e_timespercent101, np.mean(nquads_combine, 0), label='n_quads')
    plt.plot(n_timespercent101, np.mean(equads_combine, 0), label='e_quads')
    # gastroc
    # plt.plot(e_timespercent101, np.mean(ngas_combine, 0), label='n_gas')
    plt.plot(n_timespercent101, np.mean(egas_combine, 0), label='e_gas')
    # hams
    # plt.plot(e_timespercent101, np.mean(nhams_combine, 0), label='n_hams')
    plt.plot(n_timespercent101, np.mean(ehams_combine, 0), label='e_hams')
    # tFL
    # plt.plot(e_timespercent101, np.mean(ntfl_combine, 0), label='n_tfl')
    plt.plot(n_timespercent101, np.mean(etfl_combine, 0), label='e_tfl')
    # intersegmental 
    # plt.plot(e_timespercent101, np.mean(einterseg_combine, 0), label='n_intersegmental')
    plt.plot(n_timespercent101, np.mean(ninterseg_combine, 0), label='e_intersegmental')
    plt.legend()
    
    
    ### fig: straight unedited natural (should include interseg for each) all
    plt.figure()
    # full
    # plt.plot(e_timespercent101, np.mean(n_avg_r, 0), label='n_full')
    # plt.plot(n_timespercent101, np.mean(e_avg_r, 0), label='e_full')
    # quads
    plt.plot(e_timespercent101, np.mean(nquads_combine, 0), label='n_quads')
    # plt.plot(n_timespercent101, np.mean(equads_combine, 0), label='e_quads')
    # gastroc
    plt.plot(e_timespercent101, np.mean(ngas_combine, 0), label='n_gas')
    # plt.plot(n_timespercent101, np.mean(egas_combine, 0), label='e_gas')
    # hams
    plt.plot(e_timespercent101, np.mean(nhams_combine, 0), label='n_hams')
    # plt.plot(n_timespercent101, np.mean(ehams_combine, 0), label='e_hams')
    # tFL
    plt.plot(e_timespercent101, np.mean(ntfl_combine, 0), label='n_tfl')
    # plt.plot(n_timespercent101, np.mean(etfl_combine, 0), label='e_tfl')
    # intersegmental 
    plt.plot(e_timespercent101, np.mean(einterseg_combine, 0), label='n_intersegmental')
    # plt.plot(n_timespercent101, np.mean(ninterseg_combine, 0), label='e_intersegmental')
    plt.legend()
    
    
    ### fig: straight unedited (should include interseg for each) all
    plt.figure()
    # full
    # plt.plot(e_timespercent101, np.mean(n_avg_r, 0), label='n_full')
    # plt.plot(n_timespercent101, np.mean(e_avg_r, 0), label='e_full')
    # quads
    plt.plot(e_timespercent101, np.mean(nquads_combine, 0), label='n_quads')
    plt.plot(n_timespercent101, np.mean(equads_combine, 0), label='e_quads')
    # gastroc
    plt.plot(e_timespercent101, np.mean(ngas_combine, 0), label='n_gas')
    plt.plot(n_timespercent101, np.mean(egas_combine, 0), label='e_gas')
    # hams
    plt.plot(e_timespercent101, np.mean(nhams_combine, 0), label='n_hams')
    plt.plot(n_timespercent101, np.mean(ehams_combine, 0), label='e_hams')
    # tFL
    plt.plot(e_timespercent101, np.mean(ntfl_combine, 0), label='n_tfl')
    plt.plot(n_timespercent101, np.mean(etfl_combine, 0), label='e_tfl')
    # intersegmental 
    plt.plot(e_timespercent101, np.mean(einterseg_combine, 0), label='n_intersegmental')
    plt.plot(n_timespercent101, np.mean(ninterseg_combine, 0), label='e_intersegmental')
    plt.legend()
    '''

    

    # going to take out the intersegmental from each and see what happens. 
    '''
    e_avg_r
    
    einterseg_combine
    equads_combine
    egas_combine
    ehams_combine
    etfl_combine
    '''
    
    equadsub = equads_combine - einterseg_combine
    egassub = egas_combine - einterseg_combine
    ehamssub = ehams_combine - einterseg_combine
    etflsub = etfl_combine - einterseg_combine

    nquadsub = nquads_combine - ninterseg_combine
    ngassub = ngas_combine - ninterseg_combine
    nhamssub = nhams_combine - ninterseg_combine
    ntflsub = ntfl_combine - ninterseg_combine
    
    '''
    ### fig: natural building on each other (have to add since they are individualized now)
    plt.figure()
    # full
    # plt.plot(e_timespercent101, np.mean(n_avg_r, 0), label='n_full')
    # plt.plot(n_timespercent101, np.mean(e_avg_r, 0), label='e_full')
    # quads
    plt.plot(e_timespercent101, 
             nquadsub.mean(0) + ngassub.mean(0) + nhamssub.mean(0) + ntflsub.mean(0) + ninterseg_combine.mean(0), 
             label='n_quads')
    # plt.plot(n_timespercent101, 
    #          equadsub.mean(0) + egassub.mean(0) + ehamssub.mean(0) + etflsub.mean(0) + einterseg_combine.mean(0), 
    #          label='e_quads')
    # gastroc
    plt.plot(e_timespercent101, 
             ngassub.mean(0) + nhamssub.mean(0) + ntflsub.mean(0) + ninterseg_combine.mean(0), 
             label='n_gas')
    # plt.plot(n_timespercent101, 
    #          egassub.mean(0) + ehamssub.mean(0) + etflsub.mean(0) + einterseg_combine.mean(0), 
    #          label='e_gas')
    # hams
    plt.plot(e_timespercent101, nhamssub.mean(0) + ntflsub.mean(0) + ninterseg_combine.mean(0),
             label='n_hams')
    # plt.plot(n_timespercent101, ehamssub.mean(0) + etflsub.mean(0) + einterseg_combine.mean(0),
    #          label='e_hams')
    # tFL
    plt.plot(e_timespercent101, ntflsub.mean(0) + ninterseg_combine.mean(0), label='n_tfl')
    # plt.plot(n_timespercent101, etflsub.mean(0) + einterseg_combine.mean(0), label='e_tfl')
    # intersegmental 
    plt.plot(e_timespercent101, einterseg_combine.mean(0), label='n_intersegmental')
    # plt.plot(n_timespercent101, ninterseg_combine.mean(0), label='e_intersegmental')
    plt.legend()

    

    ### fig: natural building on each other (have to add since they are individualized now)
    plt.figure()
    # full
    # plt.plot(e_timespercent101, np.mean(n_avg_r, 0), label='n_full')
    # plt.plot(n_timespercent101, np.mean(e_avg_r, 0), label='e_full')
    # quads
    # plt.plot(e_timespercent101, 
    #          nquadsub.mean(0) + ngassub.mean(0) + nhamssub.mean(0) + ntflsub.mean(0) + ninterseg_combine.mean(0), 
    #          label='n_quads')
    plt.plot(n_timespercent101, 
              equadsub.mean(0) + egassub.mean(0) + ehamssub.mean(0) + etflsub.mean(0) + einterseg_combine.mean(0), 
              label='e_quads')
    # gastroc
    # plt.plot(e_timespercent101, 
    #          ngassub.mean(0) + nhamssub.mean(0) + ntflsub.mean(0) + ninterseg_combine.mean(0), 
    #          label='n_gas')
    plt.plot(n_timespercent101, 
              egassub.mean(0) + ehamssub.mean(0) + etflsub.mean(0) + einterseg_combine.mean(0), 
              label='e_gas')
    # hams
    # plt.plot(e_timespercent101, nhamssub.mean(0) + ntflsub.mean(0) + ninterseg_combine.mean(0),
    #          label='n_hams')
    plt.plot(n_timespercent101, ehamssub.mean(0) + etflsub.mean(0) + einterseg_combine.mean(0),
              label='e_hams')
    # tFL
    # plt.plot(e_timespercent101, ntflsub.mean(0) + ninterseg_combine.mean(0), label='n_tfl')
    plt.plot(n_timespercent101, etflsub.mean(0) + einterseg_combine.mean(0), label='e_tfl')
    # intersegmental 
    # plt.plot(e_timespercent101, einterseg_combine.mean(0), label='n_intersegmental')
    plt.plot(n_timespercent101, ninterseg_combine.mean(0), label='e_intersegmental')
    plt.legend()
    '''

    
    ### fig: both building on each other (have to add since they are individualized now)
    plt.figure(figsize=(11,8.5), dpi=300)
    # full
    # plt.plot(e_timespercent101, np.mean(n_avg_r, 0), label='n_full')
    # plt.plot(n_timespercent101, np.mean(e_avg_r, 0), label='e_full')
    # quads
    plt.plot(e_timespercent101, 
                nquadsub.mean(0) + ngassub.mean(0) + nhamssub.mean(0) + ntflsub.mean(0) + ninterseg_combine.mean(0), 
                label='n_quads', color='#5e3c99')
    plt.plot(n_timespercent101, 
                equadsub.mean(0) + egassub.mean(0) + ehamssub.mean(0) + etflsub.mean(0) + einterseg_combine.mean(0), 
                label='e_quads', linestyle='dashed', color='#5e3c99')
    # gastroc
    plt.plot(e_timespercent101, 
                ngassub.mean(0) + nhamssub.mean(0) + ntflsub.mean(0) + ninterseg_combine.mean(0), 
                label='n_gas', color='#b2abd2')
    plt.plot(n_timespercent101, 
                egassub.mean(0) + ehamssub.mean(0) + etflsub.mean(0) + einterseg_combine.mean(0), 
                label='e_gas', linestyle='dashed', color='#b2abd2')
    # hams
    plt.plot(e_timespercent101, nhamssub.mean(0) + ntflsub.mean(0) + ninterseg_combine.mean(0),
                label='n_hams', color='grey')
    plt.plot(n_timespercent101, ehamssub.mean(0) + etflsub.mean(0) + einterseg_combine.mean(0),
                label='e_hams', linestyle='dashed', color='grey')
    # tFL
    plt.plot(e_timespercent101, ntflsub.mean(0) + ninterseg_combine.mean(0), label='n_tfl', color='#fdb863')
    plt.plot(n_timespercent101, etflsub.mean(0) + einterseg_combine.mean(0), label='e_tfl', linestyle='dashed', color='#fdb863')
    # intersegmental 
    plt.plot(e_timespercent101, einterseg_combine.mean(0), label='n_intersegmental', color='#e66101')
    plt.plot(n_timespercent101, ninterseg_combine.mean(0), label='e_intersegmental', linestyle='dashed', color='#e66101')
    plt.legend(fontsize=20)
    plt.ylabel('BW', fontsize=24)
    plt.xlabel('% gait cycle', fontsize=24)
    plt.xticks(fontsize=24)
    plt.yticks(fontsize=24)

    plt.show()
    
    # think about cumulative loading or some metric- literature # ross miller? Scott U will be good resource
    pdb.set_trace()

   
    ###################################################################
    # NEW PLOTTING SECTION FOR FIGURES IN THE PAPER
    ###################################################################
    
    # 1) tweak the contributions figure for just a breakdown not comparing changes
    ### fig: both building on each other (have to add since they are individualized now)
    plt.figure(figsize=(11,8.5), dpi=300)
    # full
    # plt.plot(e_timespercent101, np.mean(n_avg_r, 0), label='n_full')
    # plt.plot(n_timespercent101, np.mean(e_avg_r, 0), label='e_full')
    # quads
    plt.plot(e_timespercent101, nquadsub.mean(0) + ngassub.mean(0) + nhamssub.mean(0) + ntflsub.mean(0) + ninterseg_combine.mean(0), label='n_quads', color='#d94701')
    # plt.plot(n_timespercent101, 
    #           equadsub.mean(0) + egassub.mean(0) + ehamssub.mean(0) + etflsub.mean(0) + einterseg_combine.mean(0), 
    #           label='e_quads', linestyle='dashed', color='#a63603')
    # gastroc
    plt.plot(e_timespercent101, ngassub.mean(0) + nhamssub.mean(0) + ntflsub.mean(0) + ninterseg_combine.mean(0), label='n_gas', color='#e6550d')
    # plt.plot(n_timespercent101, 
    #           egassub.mean(0) + ehamssub.mean(0) + etflsub.mean(0) + einterseg_combine.mean(0), 
    #           label='e_gas', linestyle='dashed', color='#b2abd2')
    # hams
    plt.plot(e_timespercent101, nhamssub.mean(0) + ntflsub.mean(0) + ninterseg_combine.mean(0),label='n_hams', color='#fd8d3c')
    # plt.plot(n_timespercent101, ehamssub.mean(0) + etflsub.mean(0) + einterseg_combine.mean(0),
    #           label='e_hams', linestyle='dashed', color='grey')
    # tFL
    plt.plot(e_timespercent101, ntflsub.mean(0) + ninterseg_combine.mean(0), label='n_tfl', color='#fdbe85')
    # plt.plot(n_timespercent101, etflsub.mean(0) + einterseg_combine.mean(0), label='e_tfl', linestyle='dashed', color='#fdb863')
    # intersegmental 
    plt.plot(n_timespercent101, ninterseg_combine.mean(0), label='n_intersegmental', color='#5e3c99')
    # plt.plot(e_timespercent101, einterseg_combine.mean(0), label='n_intersegmental', linestyle='dashed', color='#e66101')

    plt.legend(fontsize=20)
    plt.ylabel('BW', fontsize=24)
    plt.xlabel('% gait cycle', fontsize=24)
    plt.xticks(fontsize=24)
    plt.yticks(fontsize=24)

    plt.show()
    
    
    
    # 2) individual subplots for each muscle group to look at changes. 
    ### fig: subplots of nat and exo for each breakdown
    fig500, ax500 = plt.subplots(1,6, figsize=(20,3), dpi=300)
    # full forces
    ax500[0].plot(n_timespercent101, n_avg_r.mean(0), label='n_full', color='orange')
    ax500[0].plot(e_timespercent101, e_avg_r.mean(0), label='e_full', color='purple')
    ax500[0].set_xlabel('% gait cycle')
    ax500[0].set_ylabel('BW - Total')
    # intersegmental forces
    ax500[1].plot(n_timespercent101, ninterseg_combine.mean(0), label='natural', color='orange')
    ax500[1].plot(e_timespercent101, einterseg_combine.mean(0), label='exotendon', color='purple')
    ax500[1].set_xlabel('% gait cycle')
    ax500[1].set_ylabel('BW - intersegmental')
    ax500[1].legend()
    # quads forces
    ax500[2].plot(n_timespercent101, nquadsub.mean(0), label='n_quads', color='orange')
    ax500[2].plot(e_timespercent101, equadsub.mean(0), label='e_quads', color='purple')
    ax500[2].set_xlabel('% gait cycle')
    ax500[2].set_ylabel('BW - quads')
    # hams forces
    ax500[3].plot(n_timespercent101, nhamssub.mean(0), label='n_hams', color='orange')
    ax500[3].plot(e_timespercent101, ehamssub.mean(0), label='e_hams', color='purple')
    ax500[3].set_xlabel('% gait cycle')
    ax500[3].set_ylabel('BW - hams')
    # gastroc forces
    ax500[4].plot(n_timespercent101, ngassub.mean(0), label='n_gas', color='orange')
    ax500[4].plot(e_timespercent101, egassub.mean(0), label='e_gas', color='purple')
    ax500[4].set_xlabel('% gait cycle')
    ax500[4].set_ylabel('BW - gastroc')
    # tfl forces
    ax500[5].plot(n_timespercent101, ntflsub.mean(0), label='n_tfl', color='orange')
    ax500[5].plot(e_timespercent101, etflsub.mean(0), label='e_tfl', color='purple')
    ax500[5].set_xlabel('% gait cycle')
    ax500[5].set_ylabel('BW - tfl')

    # plt.legend('natural', 'exotendon')
    # plt.legend()
    plt.tight_layout()
    plt.show()
    
    
    ###########################################################################
    # figure that does just the changes, not the totals between each group
    ### fig: subplots of nat and exo for each breakdown
    fig501, ax501 = plt.subplots(1,6, figsize=(20,3), dpi=300)
    # full forces
    ax501[0].plot(n_timespercent101, e_avg_r.mean(0) - n_avg_r.mean(0), label='n_full', color='orange')
    # ax501[0].plot(e_timespercent101, e_avg_r.mean(0), label='e_full', color='purple')
    ax501[0].set_xlabel('% gait cycle')
    ax501[0].set_ylabel('BW - Total')
    # intersegmental forces
    ax501[1].plot(n_timespercent101, einterseg_combine.mean(0) - ninterseg_combine.mean(0), label='natural', color='orange')
    # ax501[1].plot(e_timespercent101, einterseg_combine.mean(0), label='exotendon', color='purple')
    ax501[1].set_xlabel('% gait cycle')
    ax501[1].set_ylabel('BW - intersegmental')
    ax501[1].legend()
    # quads forces
    ax501[2].plot(n_timespercent101, equadsub.mean(0) - nquadsub.mean(0), label='n_quads', color='orange')
    # ax501[2].plot(e_timespercent101, equadsub.mean(0), label='e_quads', color='purple')
    ax501[2].set_xlabel('% gait cycle')
    ax501[2].set_ylabel('BW - quads')
    # hams forces
    ax501[3].plot(n_timespercent101, ehamssub.mean(0) - nhamssub.mean(0), label='n_hams', color='orange')
    # ax501[3].plot(e_timespercent101, ehamssub.mean(0), label='e_hams', color='purple')
    ax501[3].set_xlabel('% gait cycle')
    ax501[3].set_ylabel('BW - hams')
    # gastroc forces
    ax501[4].plot(n_timespercent101, egassub.mean(0) - ngassub.mean(0), label='n_gas', color='orange')
    # ax501[4].plot(e_timespercent101, egassub.mean(0), label='e_gas', color='purple')
    ax501[4].set_xlabel('% gait cycle')
    ax501[4].set_ylabel('BW - gastroc')
    # tfl forces
    ax501[5].plot(n_timespercent101, etflsub.mean(0) - ntflsub.mean(0), label='n_tfl', color='orange')
    # ax501[5].plot(e_timespercent101, etflsub.mean(0), label='e_tfl', color='purple')
    ax501[5].set_xlabel('% gait cycle')
    ax501[5].set_ylabel('BW - tfl')


    # plt.legend('natural', 'exotendon')
    # plt.legend()
    plt.tight_layout()
    plt.show()
    
    
    ###########################################################################
    # figure that does just the changes, all on one plot
    ### fig: 
    plt.figure(figsize=(11,8), dpi=300)
    plt.axhline(color='black')
    # full forces
    plt.plot(n_timespercent101, e_avg_r.mean(0) - n_avg_r.mean(0), label='Total force', color='black')
    # quads forces
    plt.plot(n_timespercent101, equadsub.mean(0) - nquadsub.mean(0), label='quadriceps', color='#542788')
    # gastroc forces
    plt.plot(n_timespercent101, egassub.mean(0) - ngassub.mean(0), label='gastrocnemius', color='#998ec3')
    # hams forces
    plt.plot(n_timespercent101, ehamssub.mean(0) - nhamssub.mean(0), label='hamstrings', color='#d8daeb')
    # tfl forces
    plt.plot(n_timespercent101, etflsub.mean(0) - ntflsub.mean(0), label='TFL', color='#f1a340')
    # intersegmental forces
    plt.plot(n_timespercent101, einterseg_combine.mean(0) - ninterseg_combine.mean(0), label='intersegmental', color='#b35806')
    plt.xlabel('% gait cycle', fontsize=24)
    plt.xticks(fontsize=24)
    plt.ylabel('BW', fontsize=24)
    plt.yticks(fontsize=24)
    plt.title('Change in Vertical KCF', fontsize=24)

    # plt.legend('natural', 'exotendon')
    plt.legend(fontsize=20)
    plt.tight_layout()
    plt.show()
    
    ###########################################################################
    # some stats on the individual muscle group segmented out
    # now get the peaks for each of the segmented muscle groups contributions
    npeaksinter = np.max(ninterseg_combine, 1)
    npeaksquad = np.max(nquadsub, 1)
    npeaksgas = np.max(ngassub, 1)
    npeakshams = np.max(nhamssub, 1)
    npeakstfl = np.max(ntflsub, 1)
    epeaksinter = np.max(einterseg_combine, 1)
    epeaksquad = np.max(equadsub, 1)
    epeaksgas = np.max(egassub, 1)
    epeakshams = np.max(ehamssub, 1)
    epeakstfl = np.max(etflsub, 1)
    # others?? - at that point show the figures you have and see what others think    
                

                
                
                
                
                
                
                
                
                
                