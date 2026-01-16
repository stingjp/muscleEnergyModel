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

ncolor1 = '#fff5eb'
ncolor2 = '#fee6ce'
ncolor3 = '#fdd0a2'
ncolor4 = '#fdae6b'
ncolor5 = '#fd8d3c'
ncolor6 = '#f16913'
ncolor7 = '#d94801'
ncolor8 = '#a63603'
ncolor9 = '#7f2704'
ncolor10 = '#662506'
ncolor11 = '#4d1c00'


# plt.rcParams['font.family'] = 'Times New Roman'
# exotendoncolor = '#f1a340'
ecolor = '#5e3c99'
ecolormid = '#b2abd2'
ecolormiddark = '#807dba'
ecolorlight = '#d8daeb'

# function from Nick Bianco - not used in script, but used as reference for moco
def calc_hip_reaction_force(self, root_dir, solution):
    modelProcessor = osim.ModelProcessor()
    model = modelProcessor.process()
    jr = osim.analyzeSpatialVec(model, solution,
                                ['.*hip.*reaction_on_parent.*'])
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
def computeHipContact(trimmodel, initTime, finalTime, trialdir, tag, whichleg, runtool):
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
    tiby = trimjra.getDependentColumn('hip_r_on_femur_r_in_femur_r_fy').to_numpy() # hip_r_on_femur_r_in_femur_r_fy
    # pdb.set_trace()
    # import matplotlib.pyplot as plt
    # plt.figure()
    # plt.plot(np.array(trimjra.getIndependentColumn()), tiby)
    return tiby
# this is the same script ^ but for the prescribed mocoinverse solution
def computeHipContactRedo(trimmodel, initTime, finalTime, trialdir, tag, whichleg, runtool):
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
        tibx = trimjra.getDependentColumn('hip_r_on_femur_r_in_femur_r_fx').to_numpy()
        tiby = trimjra.getDependentColumn('hip_r_on_femur_r_in_femur_r_fy').to_numpy()
        tibz = trimjra.getDependentColumn('hip_r_on_femur_r_in_femur_r_fz').to_numpy()
    elif whichleg == 'left':
        tibx = trimjra.getDependentColumn('hip_l_on_femur_l_in_femur_l_fx').to_numpy()
        tiby = trimjra.getDependentColumn('hip_l_on_femur_l_in_femur_l_fy').to_numpy()
        tibz = trimjra.getDependentColumn('hip_l_on_femur_l_in_femur_l_fz').to_numpy()
    elif whichleg == 'both':
        tibxr = trimjra.getDependentColumn('hip_r_on_femur_r_in_femur_r_fx').to_numpy()
        tibxl = trimjra.getDependentColumn('hip_l_on_femur_l_in_femur_l_fx').to_numpy()
        tibyr = trimjra.getDependentColumn('hip_r_on_femur_r_in_femur_r_fy').to_numpy()
        tibyl = trimjra.getDependentColumn('hip_l_on_femur_l_in_femur_l_fy').to_numpy()
        tibzr = trimjra.getDependentColumn('hip_r_on_femur_r_in_femur_r_fz').to_numpy()
        tibzl = trimjra.getDependentColumn('hip_l_on_femur_l_in_femur_l_fz').to_numpy()
        # tiby = np.array([tibyr, tibyl])

        # # do a little peak finding, then get the index, and orient the index to match the peak in the right
        # # then add the two together.
        # # Find the peak value and index for tibyl
        # peak_tibyl = np.min(tibyl)
        # peak_index_tibyl = np.argmin(tibyl)

        # # Find the peak value and index for tibyr
        # peak_tibyr = np.min(tibyr)
        # peak_index_tibyr = np.argmin(tibyr)

        # # Rotate tibyl to match the peak index of tibyr
        # tibyl_rotated = np.roll(tibyl, peak_index_tibyr - peak_index_tibyl)
        # tibxl_rotated = np.roll(tibxl, peak_index_tibyr - peak_index_tibyl)
        # tibzl_rotated = np.roll(tibzl, peak_index_tibyr - peak_index_tibyl)
        
        ### use the vertical knee as the best signal to align the two sides.
        vertright = trimjra.getDependentColumn('walker_knee_r_on_tibia_r_in_tibia_r_fy').to_numpy()
        vertleft = trimjra.getDependentColumn('walker_knee_l_on_tibia_l_in_tibia_l_fy').to_numpy()

        peak_vertleft = np.min(vertleft)
        peak_index_vertleft = np.argmin(vertleft)
        peak_vertright = np.min(vertright)
        peak_index_vertright = np.argmin(vertright)
        # use these for rotations on tib data
        tibyl_rotated = np.roll(tibyl, peak_index_vertright - peak_index_vertleft)
        tibxl_rotated = np.roll(tibxl, peak_index_vertright - peak_index_vertleft)
        tibzl_rotated = np.roll(tibzl, peak_index_vertright - peak_index_vertleft)

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
def computeHipContactPrescribe(trimmodel, initTime, finalTime, trialdir, tag):
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
    tiby = trimjra.getDependentColumn('hip_r_on_femur_r_in_femur_r_fy').to_numpy()
    # pdb.set_trace()
    # import matplotlib.pyplot as plt
    # plt.figure()
    # plt.plot(np.array(trimjra.getIndependentColumn()), tiby)
    return tiby
# method for computing the individual muscle contributions to knee contact force
def getHipContactributions(trialdir, musclesWanted_split, tag, whichleg, runtool):
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
        jray = computeHipContact(trimmodel, initTime, finalTime, trialdir, tag, whichleg, runtool)
    
    
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
        jray = computeHipContact(trimmodel, initTime, finalTime, trialdir, tag, whichleg, runtool)
    
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
        jray = computeHipContact(trimmodel, initTime, finalTime, trialdir, tag, whichleg, runtool)
    
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
        jray = computeHipContact(trimmodel, initTime, finalTime, trialdir, tag, whichleg, runtool)
        
    
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
        jray = computeHipContact(trimmodel, initTime, finalTime, trialdir, tag, whichleg, runtool)
    
    return jray
# method for computing the individual muscle contributions to knee contact force
def getHipContactributionsRedo(trialdir, musclesWanted_split, tag, whichleg, runtool):
    if not runtool:
        print('\n NOT RUNNING TOOL - JUST COLLECTING OLD RESULTS')
        if musclesWanted_split == []:
            # load the model
            trimmodel = osim.Model('trimmingmodel2_redo_' + tag + '.osim')
            # timings - these don't matter since we are not rerunning, just need filler for the function
            initTime = 0.0
            finalTime = 1.0
            jrax, jray, jraz = computeHipContactRedo(trimmodel, initTime, finalTime, trialdir, tag, whichleg, runtool)
        elif 'all' in musclesWanted_split:
            # load the model
            trimmodel = osim.Model('trimmingmodel2_redo_' + tag + '.osim')
            # timings - these don't matter since we are not rerunning, just need filler for the function
            initTime = 0.0
            finalTime = 1.0
            jrax, jray, jraz = computeHipContactRedo(trimmodel, initTime, finalTime, trialdir, tag, whichleg, runtool)
        elif 'reserve' in musclesWanted_split:
            # load the model
            trimmodel = osim.Model('trimmingmodel2_redo_' + tag + '.osim')
            # timings - these don't matter since we are not rerunning, just need filler for the function
            initTime = 0.0
            finalTime = 1.0
            jrax, jray, jraz = computeHipContactRedo(trimmodel, initTime, finalTime, trialdir, tag, whichleg, runtool)
        elif 'none' in musclesWanted_split:
            # load the model
            trimmodel = osim.Model('trimmingmodel2_redo_' + tag + '.osim')
            # timings - these don't matter since we are not rerunning, just need filler for the function
            initTime = 0.0
            finalTime = 1.0
            jrax, jray, jraz = computeHipContactRedo(trimmodel, initTime, finalTime, trialdir, tag, whichleg, runtool)
        else:
            # load the model
            trimmodel = osim.Model('trimmingmodel2_redo_' + tag + '.osim')
            # timings - these don't matter since we are not rerunning, just need filler for the function
            initTime = 0.0
            finalTime = 1.0
            jrax, jray, jraz = computeHipContactRedo(trimmodel, initTime, finalTime, trialdir, tag, whichleg, runtool)

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
            jrax, jray, jraz = computeHipContactRedo(trimmodel, initTime, finalTime, trialdir, tag, whichleg, runtool)
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
            jrax, jray, jraz = computeHipContactRedo(trimmodel, initTime, finalTime, trialdir, tag, whichleg, runtool)
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
            jrax, jray, jraz = computeHipContactRedo(trimmodel, initTime, finalTime, trialdir, tag, whichleg, runtool)  
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
            jrax, jray, jraz = computeHipContactRedo(trimmodel, initTime, finalTime, trialdir, tag, whichleg, runtool)
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
            jrax, jray, jraz = computeHipContactRedo(trimmodel, initTime, finalTime, trialdir, tag, whichleg, runtool)

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
def getHipContactributionsPrescribe(trialdir, musclesWanted_split, tag):
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
        jray = computeHipContactPrescribe(trimmodel, initTime, finalTime, trialdir, tag)
    
    
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
        jray = computeHipContactPrescribe(trimmodel, initTime, finalTime, trialdir, tag)
    
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
        jray = computeHipContactPrescribe(trimmodel, initTime, finalTime, trialdir, tag)
    
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
        jray = computeHipContactPrescribe(trimmodel, initTime, finalTime, trialdir, tag)
        
    
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
        jray = computeHipContactPrescribe(trimmodel, initTime, finalTime, trialdir, tag)
    
    return jray
# plotting method that works for each component of the hip contact force
def plotHipContactForce(tagcomponent, analyzedir, welksubjects, ncolor, ecolor, n_timespercent101, e_timespercent101, datacombined):
    # start by separating the data...
    ninterseg_combine = datacombined['ninterseg_combine']
    einterseg_combine = datacombined['einterseg_combine']
    nflexors_combine = datacombined['nflexors_combine']
    eflexors_combine = datacombined['eflexors_combine']
    nadductors_combine = datacombined['nadductors_combine']
    eadductors_combine = datacombined['eadductors_combine']
    nabductors_combine = datacombined['nabductors_combine']
    eabductors_combine = datacombined['eabductors_combine']
    nextensors_combine = datacombined['nextensors_combine']
    eextensors_combine = datacombined['eextensors_combine']
    nhams_combine = datacombined['nhams_combine']
    ehams_combine = datacombined['ehams_combine']
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

    # flexors forces
    for i, curve in enumerate(nflexors_combine):
        ax9[1].plot(n_timespercent101, curve, label='natural'+str(i), color=ncolors[i] if len(welksubjects) == 1 else ncolor)
    for i, curve in enumerate(eflexors_combine):
        ax9[1].plot(e_timespercent101, curve, label='exotendon'+str(i), color=ecolors[i] if len(welksubjects) == 1 else ecolor)
    # ax9[1].plot(n_timespercent101, nflexors_combine, label='natural', color=ncolor)
    # ax9[1].plot(e_timespercent101, eflexors_combine, label='exotendon', color=ecolor)
    ax9[1].set_xlabel('% Gait cycle')
    # ax9[1].set_ylabel('Force (BW)')
    # ax9[1].legend()
    ax9[1].set_title('Hip flexors')

    # adductor forces
    for i, curve in enumerate(nadductors_combine):
        ax9[2].plot(n_timespercent101, curve, label='natural'+str(i), color=ncolors[i] if len(welksubjects) == 1 else ncolor)
    for i, curve in enumerate(eadductors_combine):
        ax9[2].plot(e_timespercent101, curve, label='exotendon'+str(i), color=ecolors[i] if len(welksubjects) == 1 else ecolor)
    # ax9[2].plot(n_timespercent101, ngas_combine, label='natural', color=ncolor)
    # ax9[2].plot(e_timespercent101, egas_combine, label='exotendon', color=ecolor)
    ax9[2].set_xlabel('% Gait cycle')
    # ax9[2].set_ylabel('Force (BW)')
    # ax9[2].legend()
    ax9[2].set_title('Adductors')

    # abductors forces
    for i, curve in enumerate(nabductors_combine):
        ax9[3].plot(n_timespercent101, curve, label='natural'+str(i), color=ncolors[i] if len(welksubjects) == 1 else ncolor)
    for i, curve in enumerate(eabductors_combine):
        ax9[3].plot(e_timespercent101, curve, label='exotendon'+str(i), color=ecolors[i] if len(welksubjects) == 1 else ecolor)
    # ax9[3].plot(n_timespercent101, nhams_combine, label='natural', color=ncolor)
    # ax9[3].plot(e_timespercent101, ehams_combine, label='exotendon', color=ecolor)
    ax9[3].set_xlabel('% Gait cycle')
    # ax9[3].set_ylabel('Force (BW)')
    # ax9[3].legend()
    ax9[3].set_title('Abductors')

    # extensor forces
    for i, curve in enumerate(nextensors_combine):
        ax9[4].plot(n_timespercent101, curve, label='natural'+str(i), color=ncolors[i] if len(welksubjects) == 1 else ncolor)
    for i, curve in enumerate(eextensors_combine):
        ax9[4].plot(e_timespercent101, curve, label='exotendon'+str(i), color=ecolors[i] if len(welksubjects) == 1 else ecolor)
    # ax9[4].plot(n_timespercent101, nextensors_combine, label='natural', color=ncolor)
    # ax9[4].plot(e_timespercent101, eextensors_combine, label='exotendon', color=ecolor)
    ax9[4].set_xlabel('% Gait cycle')
    # ax9[4].set_ylabel('Force (BW)')
    # ax9[4].legend()
    ax9[4].set_title('Extensors')

    # hamstring forces
    for i, curve in enumerate(nhams_combine):
        ax9[5].plot(n_timespercent101, curve, label='natural'+str(i), color=ncolors[i] if len(welksubjects) == 1 else ncolor)
    for i, curve in enumerate(ehams_combine):
        ax9[5].plot(e_timespercent101, curve, label='exotendon'+str(i), color=ecolors[i] if len(welksubjects) == 1 else ecolor)
    # ax9[5].plot(n_timespercent101, nhams_combine, label='natural', color=ncolor)
    # ax9[5].plot(e_timespercent101, ehams_combine, label='exotendon', color=ecolor)
    ax9[5].set_xlabel('% Gait cycle')
    # ax9[5].set_ylabel('Force (BW)')
    # ax9[5].legend()
    ax9[5].set_title('Hamstrings')

    # all forces from whole analysis
    for i, curve in enumerate(nall_combine):
        ax9[6].plot(n_timespercent101, curve, label='natural'+str(i), color=ncolors[i] if len(welksubjects) == 1 else ncolor)
    for i, curve in enumerate(eall_combine):
        ax9[6].plot(e_timespercent101, curve, label='exotendon'+str(i), color=ecolors[i] if len(welksubjects) == 1 else ecolor)
    # ax9[6].plot(n_timespercent101, nall_combine, label='natural', color=ncolor)
    # ax9[6].plot(e_timespercent101, eall_combine, label='exotendon', color=ecolor)
    ax9[6].set_xlabel('% Gait cycle')
    # ax9[6].set_ylabel('Force (BW)')
    ax9[6].set_title('Total ' + tagcomponent + ' contact')
    # ax9[6].legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10)
    # Hide the last subplot and use it to display the legend   
    ax9[7].axis('off')
    handles, labels = ax9[6].get_legend_handles_labels()
    if len(welksubjects) == 1:
        ax9[7].legend(handles, labels, loc='center', fontsize=14)

    fig9.tight_layout()
    plt.savefig(analyzedir + '\\hip_contact1_' + whichleg + '_' + tagcomponent + '.png')

    ###########################################################################
    # figure: segmenting all the muscles between exo and nat 
    ## really nice figure for seeing what is going on, but likely not going to 
    ## be in the paper...
    fig10, ax10 = plt.subplots(1,8, figsize=(18,3), dpi=500)
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
    ax10[0].set_title('Total\n' + tagcomponent + ' Hip Contact', **font_properties)
    ax10[0].tick_params(axis='both', which='major', labelsize=fontz)
    
    # flexor forces
    ax10[1].fill_between(n_timespercent101, np.mean(nflexors_combine, 0) - np.std(nflexors_combine, 0), np.mean(nflexors_combine, 0) + np.std(nflexors_combine, 0), color=ncolor, alpha=0.2)#, label='Natural St. Dev.')
    ax10[1].fill_between(e_timespercent101, np.mean(eflexors_combine, 0) - np.std(eflexors_combine, 0), np.mean(eflexors_combine, 0) + np.std(eflexors_combine, 0), color=ecolor, alpha=0.2)#, label='Exotendon St. Dev.')
    ax10[1].plot(n_timespercent101, np.mean(nflexors_combine, 0), label='Natural (Mean \u00B1 Std.)', color=ncolor)
    ax10[1].plot(e_timespercent101, np.mean(eflexors_combine, 0), label='Exotendon (Mean \u00B1 Std.)', color=ecolor)
    ax10[1].set_xlabel('% Gait cycle', **font_properties)
    ax10[1].tick_params(axis='both', which='major', labelsize=fontz)
    # ax10[1].set_ylabel('Force (BW)')
    # ax10[1].legend()
    ax10[1].set_title('Contribution of\nHip Flexors', **font_properties)

    # adductor forces
    ax10[2].fill_between(n_timespercent101, np.mean(nadductors_combine, 0) - np.std(nadductors_combine, 0), np.mean(nadductors_combine, 0) + np.std(nadductors_combine, 0), color=ncolor, alpha=0.2)
    ax10[2].fill_between(e_timespercent101, np.mean(eadductors_combine, 0) - np.std(eadductors_combine, 0), np.mean(eadductors_combine, 0) + np.std(eadductors_combine, 0), color=ecolor, alpha=0.2)
    ax10[2].plot(n_timespercent101, np.mean(nadductors_combine, 0), label='natural', color=ncolor)
    ax10[2].plot(e_timespercent101, np.mean(eadductors_combine, 0), label='exotendon', color=ecolor)
    ax10[2].set_xlabel('% Gait cycle', **font_properties)
    ax10[2].tick_params(axis='both', which='major', labelsize=fontz)
    # ax10[2].set_ylabel('Force (BW)')
    # ax10[2].legend()
    ax10[2].set_title('Contribution of\nHip Adductors', **font_properties)

    # abductors forces
    ax10[3].fill_between(n_timespercent101, np.mean(nabductors_combine, 0) - np.std(nabductors_combine, 0), np.mean(nabductors_combine, 0) + np.std(nabductors_combine, 0), color=ncolor, alpha=0.2)
    ax10[3].fill_between(e_timespercent101, np.mean(eabductors_combine, 0) - np.std(eabductors_combine, 0), np.mean(eabductors_combine, 0) + np.std(eabductors_combine, 0), color=ecolor, alpha=0.2)
    ax10[3].plot(n_timespercent101, np.mean(nabductors_combine, 0), label='natural', color=ncolor)
    ax10[3].plot(e_timespercent101, np.mean(eabductors_combine, 0), label='exotendon', color=ecolor)
    ax10[3].set_xlabel('% Gait cycle', **font_properties)
    ax10[3].tick_params(axis='both', which='major', labelsize=fontz)
    # ax10[3].set_ylabel('Force (BW)')
    # ax10[3].legend()
    ax10[3].set_title('Contribution of\nHip abductors', **font_properties)

    # tfl forces
    ax10[4].fill_between(n_timespercent101, np.mean(nextensors_combine, 0) - np.std(nextensors_combine, 0), np.mean(nextensors_combine, 0) + np.std(nextensors_combine, 0), color=ncolor, alpha=0.2)
    ax10[4].fill_between(e_timespercent101, np.mean(eextensors_combine, 0) - np.std(eextensors_combine, 0), np.mean(eextensors_combine, 0) + np.std(eextensors_combine, 0), color=ecolor, alpha=0.2)
    ax10[4].plot(n_timespercent101, np.mean(nextensors_combine, 0), label='natural', color=ncolor)
    ax10[4].plot(e_timespercent101, np.mean(eextensors_combine, 0), label='exotendon', color=ecolor)
    ax10[4].set_xlabel('% Gait cycle', **font_properties)
    ax10[4].tick_params(axis='both', which='major', labelsize=fontz)
    # ax10[4].set_ylabel('Force (BW)')
    # ax10[4].legend()
    ax10[4].set_title('Contribution of\nHip extensors', **font_properties)

    # hamstring forces
    ax10[5].fill_between(n_timespercent101, np.mean(nhams_combine, 0) - np.std(nhams_combine, 0), np.mean(nhams_combine, 0) + np.std(nhams_combine, 0), color=ncolor, alpha=0.2)
    ax10[5].fill_between(e_timespercent101, np.mean(ehams_combine, 0) - np.std(ehams_combine, 0), np.mean(ehams_combine, 0) + np.std(ehams_combine, 0), color=ecolor, alpha=0.2)
    ax10[5].plot(n_timespercent101, np.mean(nhams_combine, 0), label='natural', color=ncolor)
    ax10[5].plot(e_timespercent101, np.mean(ehams_combine, 0), label='exotendon', color=ecolor)
    ax10[5].set_xlabel('% Gait cycle', **font_properties)
    # ax10[5].set_ylabel('Force (BW)', **font_properties)
    ax10[5].set_title('Contribution of\nHamstrings', **font_properties)
    ax10[5].tick_params(axis='both', which='major', labelsize=fontz)

    # intersegmental forces average
    ax10[6].fill_between(n_timespercent101, np.mean(ninterseg_combine, 0) - np.std(ninterseg_combine, 0), np.mean(ninterseg_combine, 0) + np.std(ninterseg_combine, 0), color=ncolor, alpha=0.2)
    ax10[6].fill_between(e_timespercent101, np.mean(einterseg_combine, 0) - np.std(einterseg_combine, 0), np.mean(einterseg_combine, 0) + np.std(einterseg_combine, 0), color=ecolor, alpha=0.2)
    ax10[6].plot(n_timespercent101, np.mean(ninterseg_combine, 0), label='natural', color=ncolor)
    ax10[6].plot(e_timespercent101, np.mean(einterseg_combine, 0), label='exotendon', color=ecolor)
    ax10[6].set_xlabel('% Gait cycle', **font_properties)
    # ax10[6].set_ylabel('Force (BW)', **font_properties)
    ax10[6].set_title('Contribution of\nIntersegmental Forces', **font_properties)
    ax10[6].tick_params(axis='both', which='major', labelsize=fontz)

    # Hide the last subplot and use it to display the legend
    ax10[7].axis('off')
    handles, labels = ax10[1].get_legend_handles_labels()
    # ax10[6].legend(handles, labels, loc='center', fontsize=fontz)
    fig10.tight_layout(pad=2.0, w_pad=0.5, h_pad=1.0)
    plt.savefig(analyzedir + '\\hip_contact2_' + whichleg + '_' + tagcomponent + '.png')

    ###########################################################################
    ### figure: differences between conditions for each and all
    #### this is a simplified look at the same as above. We can see nice stuff, 
    #### but likely not going to be in the paper. 
    fig11, ax11 = plt.subplots(1,8, figsize=(14,3))#, dpi=300)
    # intersegmental forces average
    ax11[0].plot(n_timespercent101, np.mean(einterseg_combine, 0) - np.mean(ninterseg_combine, 0))
    ax11[0].set_xlabel('% Gait cycle')
    ax11[0].set_ylabel('Force (BW)')
    ax11[0].set_title('intersegmental')
    # ax11[0].legend()

    # flexor forces
    ax11[1].plot(n_timespercent101, np.mean(eflexors_combine, 0) - np.mean(nflexors_combine, 0))
    ax11[1].set_xlabel('% Gait cycle')
    # ax11[1].set_ylabel('Force (BW)')
    # ax11[1].legend()
    ax11[1].set_title('flexors')

    # adductor forces
    ax11[2].plot(n_timespercent101, np.mean(eadductors_combine, 0) - np.mean(nadductors_combine, 0))
    ax11[2].set_xlabel('% Gait cycle')
    # ax11[2].set_ylabel('Force (BW)')
    # ax11[2].legend()
    ax11[2].set_title('adductors')

    # abductors forces
    ax11[3].plot(n_timespercent101, np.mean(eabductors_combine, 0) - np.mean(nabductors_combine, 0))
    ax11[3].set_xlabel('% Gait cycle')
    # ax11[3].set_ylabel('Force (BW)')
    # ax11[3].legend()
    ax11[3].set_title('abductors')

    # extensor forces
    ax11[4].plot(n_timespercent101, np.mean(eextensors_combine, 0) - np.mean(nextensors_combine, 0))
    ax11[4].set_xlabel('% Gait cycle')
    # ax11[4].set_ylabel('Force (BW)')
    # ax11[4].legend()
    ax11[4].set_title('extensors')

    # hamstring forces
    ax11[5].plot(n_timespercent101, np.mean(ehams_combine, 0) - np.mean(nhams_combine, 0))
    ax11[5].set_xlabel('% Gait cycle')
    # ax11[5].set_ylabel('Force (BW)')
    # ax11[5].legend()
    ax11[5].set_title('hamstrings')

    # all forces from whole analysis
    ax11[6].plot(n_timespercent101, np.mean(eall_combine,0) - np.mean(nall_combine,0), label='Exotendon diff from Natural')
    ax11[6].set_xlabel('% Gait cycle')
    # ax11[6].set_ylabel('Force (BW)')
    # ax11[6].legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10)
    ax11[6].set_title('Total ' + tagcomponent + ' contact')

    # Hide the last subplot and use it to display the legend
    ax11[7].axis('off')
    handles, labels = ax11[6].get_legend_handles_labels()
    ax11[7].legend(handles, labels, loc='center', fontsize=10)
    fig11.tight_layout()
    plt.savefig(analyzedir + '\\hip_contact3_' + whichleg + '_' + tagcomponent + '.png')


    ###########################################################################
    ### figure: changes in force segmented together on plot
    #### Possible paper figure for R3. 
    colors4 = ['#648FFF', '#785EF0', '#DC267F', '#FE6100', '#FFB000', '#009E73']
    fig12, ax12 = plt.subplots(1, 2, figsize=(12,4.55), dpi=500)
    ax12[0].axhline(0, color='black', linewidth=1)
    # intersegmental forces average
    ax12[0].plot(n_timespercent101, np.mean(einterseg_combine, 0) - np.mean(ninterseg_combine, 0), label='Intersegmental', color='black', linewidth=3)
    # flexor forces
    ax12[0].plot(n_timespercent101, np.mean(eflexors_combine, 0) - np.mean(nflexors_combine, 0), label='Hip flexors', color=colors4[0], linewidth=3)
    # adductors forces
    ax12[0].plot(n_timespercent101, np.mean(eadductors_combine, 0) - np.mean(nadductors_combine, 0), label='Hip adductors', color=colors4[1], linewidth=3)
    # abductor forces
    ax12[0].plot(n_timespercent101, np.mean(eabductors_combine, 0) - np.mean(nabductors_combine, 0), label='Hip abductors', color=colors4[2], linewidth=3)
    # extensor forces
    ax12[0].plot(n_timespercent101, np.mean(eextensors_combine, 0) - np.mean(nextensors_combine, 0), label='Hip extensors', color=colors4[3], linewidth=3)
    # hamstring forces
    ax12[0].plot(n_timespercent101, np.mean(ehams_combine, 0) - np.mean(nhams_combine, 0), label='Hamstrings', color=colors4[4], linewidth=3)
    # reserve forces
    # ax12[0].plot(n_timespercent101, np.mean(ereserve_combine, 0) - np.mean(nreserve_combine, 0), label='reserves', color=colors4[4])
    # added all forces
    # ax12[0].plot(n_timespercent101, (np.mean(equads_combine,0) + np.mean(ehams_combine,0) + np.mean(egas_combine,0) + np.mean(etfl_combine,0) + np.mean(einterseg_combine,0) + np.mean(ereserve_combine,0)) - (np.mean(nquads_combine,0) + np.mean(nhams_combine,0) + np.mean(ngas_combine,0) + np.mean(ntfl_combine,0) + np.mean(ninterseg_combine,0) + np.mean(nreserve_combine,0)), label='Total')
    # all forces from whole analysis
    ax12[0].plot(n_timespercent101, np.mean(eall_combine,0) - np.mean(nall_combine,0), label='Total ' + tagcomponent + ' contact', linestyle='dashed', color='black', linewidth=3)
    ax12[0].set_xlabel('% Gait cycle', fontsize=16)
    ax12[0].set_ylabel('Vertical hip contact difference (BW)', fontsize=16)
    ax12[0].set_title('Exotendon change in contact force', fontsize=16)
    ax12[0].tick_params(axis='both', which='major', labelsize=16)
    # hide the second subplots and use it for the legend
    ax12[1].axis('off')
    handles, labels = ax12[0].get_legend_handles_labels()
    ax12[1].legend(handles, labels, loc='center', fontsize=16)
    fig12.tight_layout()
    fig12.savefig(analyzedir + '\\hip_contact4_' + whichleg + '_' + tagcomponent + '.png')

    # TODO: figure out why the difference in total and all added together. 
    # okay so not in how I am adding/averaging. has to be something in how the analysis is done between them.... am I missing something??
    
    fig13, ax13 = plt.subplots(1,3, figsize=(14,5))#, dpi=300)
    # intersegmental forces average - natural
    ax13[0].plot(n_timespercent101, np.mean(ninterseg_combine, 0), label='natural_interseg')
    ax13[0].plot(n_timespercent101, np.mean(ninterseg_combine + nflexors_combine, 0), label='natural_interseg + flexors')
    ax13[0].plot(n_timespercent101, np.mean(ninterseg_combine+nflexors_combine+nadductors_combine, 0), label='natural_interseg+flexors+adductors')
    ax13[0].plot(n_timespercent101, np.mean(ninterseg_combine+nflexors_combine+nadductors_combine+nabductors_combine, 0), label='natural_interseg+flexors+adductors+abductors')
    ax13[0].plot(n_timespercent101, np.mean(ninterseg_combine+nflexors_combine+nadductors_combine+nabductors_combine+nextensors_combine, 0), label='natural_interseg+flexors+adductors+abductors+extensors')
    ax13[0].plot(n_timespercent101, np.mean(ninterseg_combine+nflexors_combine+nadductors_combine+nabductors_combine+nextensors_combine+nhams_combine, 0), label='natural_interseg+flexors+adductors+abductors+extensors+hams')
    ax13[0].plot(n_timespercent101, np.mean(nall_combine, 0), label='nat_all', linestyle='dotted')
    ax13[0].set_xlabel('% Gait cycle')
    ax13[0].set_ylabel('Force (BW)')
    ax13[0].set_title('natural')
    # ax13[0].legend()
    # intersegmental forces average - exotendon
    ax13[1].plot(e_timespercent101, np.mean(einterseg_combine, 0), label='exo_interseg')
    ax13[1].plot(e_timespercent101, np.mean(einterseg_combine+eflexors_combine, 0), label='exo_interseg + flexors')
    ax13[1].plot(e_timespercent101, np.mean(einterseg_combine+eflexors_combine+eadductors_combine, 0), label='exo_interseg+flexors+adductors')
    ax13[1].plot(e_timespercent101, np.mean(einterseg_combine+eflexors_combine+eadductors_combine+eabductors_combine, 0), label='exo_interseg+flexors+adductors+abductors')
    ax13[1].plot(e_timespercent101, np.mean(einterseg_combine+eflexors_combine+eadductors_combine+eabductors_combine+eextensors_combine, 0), label='exo_interseg+flexors+adductors+abductors+extensors')
    ax13[1].plot(e_timespercent101, np.mean(einterseg_combine+eflexors_combine+eadductors_combine+eabductors_combine+eextensors_combine+ehams_combine, 0), label='exo_interseg+flexors+adductors+abductors+extensors+hams')
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
    plt.savefig(analyzedir + '\\hip_contact5_' + whichleg + '_' + tagcomponent + '.png')

    ###########################################################################
    ### Figure: total pop average for right leg between nat and exo
    #### Figure in the paper R1... 
    figcon6, axcon6 = plt.subplots(1, 2, figsize=(12,4.55), dpi=500)
    axcon6[0].fill_between(n_timespercent101, np.mean(nall_combine,0)-np.std(nall_combine,0), np.mean(nall_combine,0)+np.std(nall_combine,0), color=ncolorlight, alpha=0.75)
    axcon6[0].fill_between(e_timespercent101, np.mean(eall_combine,0)-np.std(eall_combine,0), np.mean(eall_combine,0)+np.std(eall_combine,0), color=ecolorlight, alpha=0.75)
    axcon6[0].plot(n_timespercent101, np.mean(nall_combine,0), color=ncolor, label='Natural (Mean \u00B1 Std.)', linewidth=3)
    axcon6[0].plot(e_timespercent101, np.mean(eall_combine,0), color=ecolor, label='Exotendon (Mean \u00B1 Std.)', linewidth=3)
    axcon6[0].set_xlabel('% Gait cycle', fontsize=16)
    axcon6[0].set_ylabel(tagcomponent + ' hip contact force (BW)', fontsize=16)
    axcon6[0].set_title('Total ' + tagcomponent + ' hip contact force', fontsize=16)
    axcon6[0].tick_params(axis='both', which='major', labelsize=16)
    # axcon6[0].xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    axcon6[0].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    if tagcomponent == 'Shear':
        axcon6[0].set_ylim([0, 4.5])
    if tagcomponent == 'Vertical':
        axcon6[0].set_ylim([0, 11])
    # hide the second subplots and use it for the legend
    axcon6[1].axis('off')
    handles, labels = axcon6[0].get_legend_handles_labels()
    axcon6[1].legend(handles, labels, loc='center', fontsize=16)
    figcon6.tight_layout()
    figcon6.savefig(analyzedir + '\\hip_contact6_' + whichleg + '_' + tagcomponent + '.png')


    ###########################################################################
    # output the data in an excel sheet
    # Create a dictionary to hold all the data
    data_dict = {
        'ninterseg_combine': ninterseg_combine,
        'einterseg_combine': einterseg_combine,
        'nflexors_combine': nflexors_combine,
        'eflexors_combine': eflexors_combine,
        'nhams_combine': nhams_combine,
        'ehams_combine': ehams_combine,
        'nadductors_combine': nadductors_combine,
        'eadductors_combine': eadductors_combine,
        'nabductors_combine': nabductors_combine,
        'eabductors_combine': eabductors_combine,
        'nextensors_combine': nextensors_combine,
        'eextensors_combine': eextensors_combine,
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
    with pd.ExcelWriter('analysis_results'+'_'+whichleg+'_hip_'+tagcomponent+'.xlsx', engine='xlsxwriter') as writer:
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
    plt.plot(n_timespercent101, np.mean(nflexors_combine,0) + np.mean(ninterseg_combine,0), label='inter + flexors', color=ncolor3)
    plt.fill_between(n_timespercent101, np.mean(ninterseg_combine,0), np.mean(nflexors_combine,0) + np.mean(ninterseg_combine,0), color=ncolor2)

    plt.plot(n_timespercent101, np.mean(nadductors_combine,0) + np.mean(nflexors_combine,0) + np.mean(ninterseg_combine,0), label='inter + flexors + adductors', color=ncolor5)
    plt.fill_between(n_timespercent101, np.mean(nflexors_combine,0) + np.mean(ninterseg_combine,0), np.mean(nadductors_combine,0) + np.mean(nflexors_combine,0) + np.mean(ninterseg_combine,0), color=ncolor4)

    plt.plot(n_timespercent101, np.mean(nhams_combine,0) + np.mean(nadductors_combine,0) + np.mean(nflexors_combine,0) + np.mean(ninterseg_combine,0), label='inter + flexors + adductors + hams', color=ncolor7)
    plt.fill_between(n_timespercent101, np.mean(nadductors_combine,0) + np.mean(nflexors_combine,0) + np.mean(ninterseg_combine,0), np.mean(nhams_combine,0) + np.mean(nadductors_combine,0) + np.mean(nflexors_combine,0) + np.mean(ninterseg_combine,0), color=ncolor6)

    plt.plot(n_timespercent101, np.mean(nabductors_combine,0) + np.mean(nhams_combine,0) + np.mean(nadductors_combine,0) + np.mean(nflexors_combine,0) + np.mean(ninterseg_combine,0), label='inter + flexors + adductors + hams + abductors', color=ncolor9)
    plt.fill_between(n_timespercent101, np.mean(nhams_combine,0) + np.mean(nadductors_combine,0) + np.mean(nflexors_combine,0) + np.mean(ninterseg_combine,0), np.mean(nabductors_combine,0) + np.mean(nhams_combine,0) + np.mean(nadductors_combine,0) + np.mean(nflexors_combine,0) + np.mean(ninterseg_combine,0), color=ncolor8)

    plt.plot(n_timespercent101, np.mean(nextensors_combine,0) + np.mean(nabductors_combine,0) + np.mean(nhams_combine,0) + np.mean(nadductors_combine,0) + np.mean(nflexors_combine,0) + np.mean(ninterseg_combine,0), label='inter + flexors + adductors + hams + abductors + extensors', color=ncolor11)
    plt.fill_between(n_timespercent101, np.mean(nabductors_combine,0) + np.mean(nhams_combine,0) + np.mean(nadductors_combine,0) + np.mean(nflexors_combine,0) + np.mean(ninterseg_combine,0), np.mean(nextensors_combine,0) + np.mean(nabductors_combine,0) + np.mean(nhams_combine,0) + np.mean(nadductors_combine,0) + np.mean(nflexors_combine,0) + np.mean(ninterseg_combine,0), color=ncolor10)

    # plt.legend()
    plt.ylabel(tagcomponent + ' hip contact force (BW)', fontsize=16)
    plt.xlabel('% Gait cycle', fontsize=16)
    plt.title(tagcomponent + ' hip contact force - Natural', fontsize=16)    
    plt.tick_params(axis='both', which='major', labelsize=16)
    plt.tight_layout()
    plt.savefig(analyzedir + '\\hip_contact7_Natural_' + whichleg + '_' + tagcomponent + '_nolegend.png')
    plt.show()
    
    
    ###########################################################################
    ### Figure: total pop average for right leg for exo - segmented shaded
    plt.figure(dpi=300)
    plt.plot(e_timespercent101, np.mean(einterseg_combine,0), label='intersegmental', color=ecolor)
    plt.fill_between(e_timespercent101, np.mean(einterseg_combine,0), color=ecolorlight)

    plt.plot(e_timespercent101, np.mean(eflexors_combine,0) + np.mean(einterseg_combine,0), label='inter + flexors', color=ncolor3)
    plt.fill_between(e_timespercent101, np.mean(einterseg_combine,0), np.mean(eflexors_combine,0) + np.mean(einterseg_combine,0), color=ncolor2)

    plt.plot(e_timespercent101, np.mean(eadductors_combine,0) + np.mean(eflexors_combine,0) + np.mean(einterseg_combine,0), label='inter + flexors + adductors', color=ncolor5)
    plt.fill_between(e_timespercent101, np.mean(eflexors_combine,0) + np.mean(einterseg_combine,0), np.mean(eadductors_combine,0) + np.mean(eflexors_combine,0) + np.mean(einterseg_combine,0), color=ncolor4)

    plt.plot(e_timespercent101, np.mean(ehams_combine,0) + np.mean(eadductors_combine,0) + np.mean(eflexors_combine,0) + np.mean(einterseg_combine,0), label='inter + flexors + adductors + hams', color=ncolor7)
    plt.fill_between(e_timespercent101, np.mean(eadductors_combine,0) + np.mean(eflexors_combine,0) + np.mean(einterseg_combine,0), np.mean(ehams_combine,0) + np.mean(eadductors_combine,0) + np.mean(eflexors_combine,0) + np.mean(einterseg_combine,0),color=ncolor6)

    plt.plot(e_timespercent101, np.mean(eabductors_combine,0) + np.mean(ehams_combine,0) + np.mean(eadductors_combine,0) + np.mean(eflexors_combine,0) + np.mean(einterseg_combine,0), label='inter + flexors + adductors + hams + abductors', color=ncolor9)
    plt.fill_between(e_timespercent101, np.mean(ehams_combine,0) + np.mean(eadductors_combine,0) + np.mean(eflexors_combine,0) + np.mean(einterseg_combine,0), np.mean(eabductors_combine,0) + np.mean(ehams_combine,0) + np.mean(eadductors_combine,0) + np.mean(eflexors_combine,0) + np.mean(einterseg_combine,0), color=ncolor8)

    plt.plot(e_timespercent101, np.mean(eextensors_combine,0) + np.mean(eabductors_combine,0) + np.mean(ehams_combine,0) + np.mean(eadductors_combine,0) + np.mean(eflexors_combine,0) + np.mean(einterseg_combine,0), label='inter + flexors + adductors + hams + abductors + extensors', color=ncolor11)
    plt.fill_between(e_timespercent101, np.mean(eabductors_combine,0) + np.mean(ehams_combine,0) + np.mean(eadductors_combine,0) + np.mean(eflexors_combine,0) + np.mean(einterseg_combine,0), np.mean(eextensors_combine,0) + np.mean(eabductors_combine,0) + np.mean(ehams_combine,0) + np.mean(eadductors_combine,0) + np.mean(eflexors_combine,0) + np.mean(einterseg_combine,0), color=ncolor10)

    # plt.legend()
    plt.ylabel(tagcomponent + ' hip contact force (BW)', fontsize=16)
    plt.xlabel('% Gait cycle', fontsize=16)
    plt.title(tagcomponent + ' hip contact force - Exotendon', fontsize=16)
    plt.tick_params(axis='both', which='major', labelsize=16)
    plt.tight_layout()
    plt.savefig(analyzedir + '\\hip_contact8_exo_' + whichleg + '_' + tagcomponent + '_nolegend.png')
    plt.show()


    ###########################################################################
    # stats for the individual muscle segmented contributions.... check out the excel sheet
    peaks_flexors_natural = np.max(nflexors_combine,1)
    peaks_adductors_natural = np.max(nadductors_combine,1)
    peaks_hams_natural = np.max(nhams_combine,1)
    peaks_abductors_natural = np.max(nabductors_combine,1)
    peaks_extensors_natural = np.max(nextensors_combine,1)
    peaks_inter_natural = np.max(ninterseg_combine,1)
    peaks_reserve_natural = np.max(nreserve_combine,1)

    peaks_flexors_exo = np.max(eflexors_combine,1)
    peaks_adductors_exo = np.max(eadductors_combine,1)
    peaks_hams_exo = np.max(ehams_combine,1)
    peaks_abductors_exo = np.max(eabductors_combine,1)
    peaks_extensors_exo = np.max(eextensors_combine,1)
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
    # previous (combo) results directory

    # resultsdir = 'C:\\Users\\jonstingel\\code\\musclemodel\\testresults\\results\\';
    # analyzedir = 'C:\\Users\\jonstingel\\code\\musclemodel\\testresults\\analysis\\';

    # resultsdir = 'C:\\Users\\jonstingel\\code\\musclemodel\\testresults - Copy\\results\\';
    # analyzedir = 'C:\\Users\\jonstingel\\code\\musclemodel\\testresults - Copy\\analysis\\';

    # resultsdir = 'C:\\Users\\jonstingel\\code\\musclemodel\\testresults - Copy - Copy\\results\\';
    # analyzedir = 'C:\\Users\\jonstingel\\code\\musclemodel\\testresults - Copy - Copy\\analysis\\';

    resultsdir = 'C:\\Users\\jonstingel\\code\\musclemodel\\testresults - Copy - Copy - Copy\\results\\';
    analyzedir = 'C:\\Users\\jonstingel\\code\\musclemodel\\testresults - Copy - Copy - Copy\\analysis\\';

    welkexoconditions = ['welkexo']
    welknaturalconditions = ['welknatural']
    welksubjects = ['welk002','welk010','welk003','welk005','welk008','welk009','welk013'];
    thingstoplot = ['contactForces']
    trials = ['trial01','trial02','trial03','trial04']
    whichleg = 'both'
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
    '''
    # loop the subjects
    for subj in range(len(welksubjects)):
        subject = welksubjects[subj]
        subjdir = os.path.join(resultsdir, subject)

        # create a structure for individual subject stuff
        welknaturalstruct = {}
        welkexostruct = {}

        # # loop through conditions
        # for cond in range(len(welkexoconditions)):
        #     condition = welkexoconditions[cond]
        #     condir = os.path.join(subjdir, condition)

        #     # loop the trials
        #     for tr in range(len(trials)):
        #         trial = trials[tr]
        #         trialdir = os.path.join(condir, trial)

        #         ### now what do we want to do at each of the trials

        #         # need a model
        #         # model = osim.Model(os.path.join(trialdir, 'post_simple_model_all_the_probes_muscletrack.osim'))
        #         model = osim.Model(os.path.join(trialdir, 'simple_model_all_the_probes_adjusted.osim'))
        #         # weld the mtp real quick
        #         modelProcessor = osim.ModelProcessor(model)
        #         weldem = osim.StdVectorString()
        #         weldem.append('mtp_r'); weldem.append('mtp_l')
        #         modelProcessor.append(osim.ModOpReplaceJointsWithWelds(weldem))
        #         model = modelProcessor.process()
        #         model.initSystem()
                
        #         # need the solution
        #         solution = osim.MocoTrajectory(os.path.join(trialdir, 'muscle_statetrack_GRFprescribe_solution_100con.sto'))
        #         statesTable = osim.TimeSeriesTable(os.path.join(trialdir, 'muscletrack_states_100con.sto'))
        #         controlsTable = osim.TimeSeriesTable(os.path.join(trialdir, 'muscletrack_controls_100con.sto'))
        #         fiberLengths = osim.TimeSeriesTable(os.path.join(trialdir, 'analyzemuscles_100conmuscletrack_MuscleAnalysis_FiberLength.sto'))
                
        #         # get a table with just the jointset states
        #         coordsTable = statesTable
        #         statelabels = coordsTable.getColumnLabels()
        #         for coord in range(len(statelabels)):
        #             coordinate = statelabels[coord]
        #             if 'jointset/' not in coordinate:
        #                 coordsTable.removeColumn(coordinate)
                
        #         # get the times - first and last
        #         times = statesTable.getIndependentColumn()
        #         initTime = times[0]
        #         finalTime = times[-1]
                    
        #         # get a table with the combined states and tendon lengths?
        #         combStates = statesTable
        #         comblabels = combStates.getColumnLabels()
        #         # loop the states table and remove the tendon force columns
        #         for col in range(len(comblabels)):
        #             columnname = comblabels[col]
        #             if 'tendon' in columnname:
        #                 combStates.removeColumn(columnname)
                        
        #         comblabels2 = combStates.getColumnLabels()
                
        #         # get the names of all the fiberlengths
        #         fiberlabels = fiberLengths.getColumnLabels()
        #         # loop through them all, grab, and drop into states table, 
        #         for fib in range(len(fiberlabels)):
        #             fiber = fiberlabels[fib]
        #             fiberlength = fiberLengths.getDependentColumn(fiber)
        #             # try to add it to the states table
        #             combStates.appendColumn('/forceset/'+fiber+'/fiber_length', fiberlength)
                
        #         # store the model masses
        #         modelmass = get_model_total_mass(trialdir, 'simple_model_all_the_probes_adjusted.osim')
        #         exostruct_combine[subject] = modelmass                

        #         ######################################################
        #         ## trying the state step approach from HPLers
        #         os.chdir(trialdir)
        #         pdb.set_trace()
                
        #         # load in the base model
        #         jramodelProcessor = osim.ModelProcessor('simple_model_all_the_probes_adjusted.osim')
        #         weldem = osim.StdVectorString();
        #         weldem.append('mtp_r');
        #         weldem.append('mtp_l');
        #         jramodelProcessor.append(osim.ModOpReplaceJointsWithWelds(weldem))
        #         stepmodel = jramodelProcessor.process()                
                
        #         # initialize state
        #         state = stepmodel.initSystem()
        #         # get info all of the force producers
        #         forceSet = stepmodel.getForceSet()
        #         nForces = forceSet.getSize()
        #         actuators = forceSet.updActuators()
        #         nActuators = actuators.getSize()
        #         muscles = stepmodel.updMuscles()
        #         nMuscles = muscles.getSize()
                
        #         # Get names of states
        #         stateNames = stepmodel.getStateVariableNames()
        #         nStates = stateNames.getSize()
                
        #         # good coordinates file. load into arrays
        #         coordValNames = [] # 23
        #         coordSpeedNames = [] # 23
        #         activationNames = [] # 80
        #         fiberlengthNames = []
        #         normTenForceNames = [] # 80
        #         for i in range(len(statelabels)):
        #             if '/value' in statelabels[i]:
        #                 coordValNames.append(statelabels[i])
        #             if '/speed' in statelabels[i]:
        #                 coordSpeedNames.append(statelabels[i])
        #             if '/activation' in statelabels[i]:
        #                 activationNames.append(statelabels[i])
        #             if '/normalized_tendon_force' in statelabels[i]:
        #                 normTenForceNames.append(statelabels[i])
                
        #         # maybe do the same thing for the controls?
        #         controlNames = controlsTable.getColumnLabels()
        #         reserveNames = []
        #         muscleControlNames = []
        #         for i in range(len(controlNames)):
        #             if 'reserve' in controlNames[i]:
        #                 reserveNames.append(controlNames[i])
        #             else:
        #                 muscleControlNames.append(controlNames[i])

        #         # get moments file, load into arrays??
                
        #         # ID prescribed coordinate values??
                
        #         # Joint reaction setup
        #         jointRxn = osim.JointReaction()
        #         jointRxn.setName('jrxnAnalysis100')
        #         wherestr = osim.ArrayStr(); wherestr.append('child')
        #         jointRxn.setInFrame(wherestr) ;
        #         jointRxn.setOnBody(wherestr) ;
                
        #         # jointRxn.setJointNames(jointNames) ;
        #         jointRxn.setStartTime(initTime);
        #         jointRxn.setFinalTime(finalTime)
        #         stepmodel.addAnalysis(jointRxn) ;
        #         jointRxn.setModel(stepmodel) ;
        #         jointRxn.setResultsDir(trialdir)
        #         jointRxn.printToXML('JrxnSetup.xml') ;
        #         time.sleep(0.5)

        #         pdb.set_trace()                
        #         ############################
        #         # another one
        #         jramodel = osim.Model(os.path.join(trialdir, 'post_simple_model_all_the_probes_muscletrack.osim'))

        #         jraoutputs = osim.StdVectorString()
        #         jraoutputs.append('.*walker_knee.*\\|reaction_on_child')
        #         os.chdir(trialdir)
        #         jras = osim.analyzeMocoTrajectorySpatialVec(jramodel, solution, jraoutputs)
        #         # jr = osim.analyzeSpatialVec(model, combStates, controlsTable, ['.*walker_knee.*reaction_on_parent.*'])
        #         time.sleep(0.5)
                
        #         ## now to actually load in the data and do something with it. 
        #         jrastab = jras.flatten(['_mx', '_my', '_mz', '_fx', '_fy', '_fz'])
        #         jrasrx = jrastab.getDependentColumn('/jointset/walker_knee_r|reaction_on_child_fx').to_numpy()
        #         jrasry = jrastab.getDependentColumn('/jointset/walker_knee_r|reaction_on_child_fy').to_numpy()
        #         jrasrz = jrastab.getDependentColumn('/jointset/walker_knee_r|reaction_on_child_fz').to_numpy()
        #         jraslx = jrastab.getDependentColumn('/jointset/walker_knee_l|reaction_on_child_fx').to_numpy()
        #         jrasly = jrastab.getDependentColumn('/jointset/walker_knee_l|reaction_on_child_fy').to_numpy()
        #         jraslz = jrastab.getDependentColumn('/jointset/walker_knee_l|reaction_on_child_fz').to_numpy()
                
        #         import matplotlib.pyplot as plt
        #         plt.figure()
        #         plt.plot(jrasry)
                
        #         ## now transform to the tibia frame from ground
        #         # get the joints
        #         transformmodel = jramodel
        #         transformmodel.finalizeConnections()
        #         joints = transformmodel.getJointSet()
        #         jointr = joints.get('walker_knee_r')
        #         jointl = joints.get('walker_knee_l')
        #         # get child bodies
        #         tibr = jointr.getChildFrame()
        #         tibr_name = tibr.getAbsolutePathString()
        #         tibl = jointl.getChildFrame()
        #         tibl_name = tibl.getAbsolutePathString()
        #         # also get the ground frame
        #         ground = transformmodel.getGround()
        #         # grab the time vector
        #         jrastime = jrastab.getIndependentColumn()
                
        #         # get the base state structure for the model
        #         state = transformmodel.initSystem()
                
        #         newjrasr = np.zeros((len(jrastime), 3))
        #         newjrasl = np.zeros((len(jrastime), 3))
        #         diffr = np.zeros((len(jrastime), 2))
        #         diffl = np.zeros((len(jrastime), 2))
        #         # loop through each time step, 
        #         coords = model.getCoordinateSet()
        #         # pdb.set_trace()
        #         for tim in range(len(solution.getTime())):
        #             # then each coordinate to pose the model and get the transform
        #             for cor in range(coords.getSize()):
        #                 cord = coords.get(cor)
        #                 statecol = statesTable.getDependentColumn(cord.getAbsolutePathString() + '/value').to_numpy()
        #                 cord.setValue(state, statecol[tim])
                        
        #             # realize model position from poses
        #             transformmodel.realizePosition(state)
        #             # then compute the transform for the forces
        #             newjrasr[tim,:] = ground.expressVectorInAnotherFrame(state, osim.Vec3(jrasrx[tim], jrasry[tim], jrasrz[tim]), tibr).to_numpy()
        #             newjrasl[tim,:] = ground.expressVectorInAnotherFrame(state, osim.Vec3(jraslx[tim], jrasly[tim], jraslz[tim]), tibl).to_numpy()
        #             diffr[tim,0] = np.linalg.norm(ground.expressVectorInAnotherFrame(state, osim.Vec3(jrasrx[tim], jrasry[tim], jrasrz[tim]), tibr).to_numpy())
        #             diffr[tim,1] = np.linalg.norm([jrasrx[tim], jrasry[tim], jrasrz[tim]])
        #             diffl[tim,0] = np.linalg.norm(ground.expressVectorInAnotherFrame(state, osim.Vec3(jraslx[tim], jrasly[tim], jraslz[tim]), tibl).to_numpy()) 
        #             diffl[tim,1] = np.linalg.norm([jraslx[tim], jrasly[tim], jraslz[tim]])
                
        #         # pdb.set_trace()
        #         # #############
        #         # # results
        #         # ### TODO clean up based on analysis in the end
                
        #         # # want to focus on just the knees first
        #         # kneeReactions = jras
        #         # jracolumns = jras.getColumnLabels()
        #         # for cols in jracolumns:
        #         #     if 'walker' not in cols: # and '_f' not in cols:
        #         #         kneeReactions.removeColumn(cols)
        #         #     else:
        #         #         print(cols)
        #         # kneeLabels = kneeReactions.getColumnLabels()
        #         # print('\n\nOkay got the knees...\n')
                
        #         # for each in kneeLabels:
        #         #     if '_f' not in each:
        #         #         kneeReactions.removeColumn(each)
        #         #     else:
        #         #         print(each)
        #         # kneeLabels = kneeReactions.getColumnLabels()
                
        #         # # now have a table with just the knee reaction loads
        #         # # create a column that is the net force
        #         # fxr = kneeReactions.getDependentColumn('walker_knee_r_on_tibia_r_in_tibia_r_fx').to_numpy()
        #         # fyr = kneeReactions.getDependentColumn('walker_knee_r_on_tibia_r_in_tibia_r_fy').to_numpy()
        #         # fzr = kneeReactions.getDependentColumn('walker_knee_r_on_tibia_r_in_tibia_r_fz').to_numpy()
        #         # fxl = kneeReactions.getDependentColumn('walker_knee_l_on_tibia_l_in_tibia_l_fx').to_numpy()
        #         # fyl = kneeReactions.getDependentColumn('walker_knee_l_on_tibia_l_in_tibia_l_fy').to_numpy()
        #         # fzl = kneeReactions.getDependentColumn('walker_knee_l_on_tibia_l_in_tibia_l_fz').to_numpy()
                
        #         # norm_r = np.sqrt(fxr**2 + fyr**2 + fzr**2)
        #         # norm_l = np.sqrt(fxl**2 + fyl**2 + fzl**2)
                
        #         # pdb.set_trace()
        #         # # welkexostruct[trial] = (times, norm_r, norm_l)
        #         # welkexostruct[trial] = (times, np.abs(fyr), np.abs(fyl))
        #         welkexostruct[trial] = (times, np.abs(newjrasr[:,1]), np.abs(newjrasl[:,1]))
        #         print(trialdir)

        # now have to loop through the natural side
        # loop through conditions natural
        for cond in range(len(welknaturalconditions)):
            condition = welknaturalconditions[cond]
            condir = os.path.join(subjdir, condition)

            # loop the trials
            for tr in range(len(trials)):
                trial = trials[tr]
                trialdir = os.path.join(condir, trial)
                
                os.chdir(trialdir)
                pdb.set_trace()
                
                ### now what do we want to do at each of the trials

                # # need a model
                # # model = osim.Model(os.path.join(trialdir, 'post_simple_model_all_the_probes_muscletrack.osim'))
                # model = osim.Model(os.path.join(trialdir, 'simple_model_all_the_probes_adjusted.osim'))
                # # weld the mtp real quick
                # modelProcessor = osim.ModelProcessor(model)
                # weldem = osim.StdVectorString()
                # weldem.append('mtp_r'); weldem.append('mtp_l')
                # modelProcessor.append(osim.ModOpReplaceJointsWithWelds(weldem))
                # modelProcessor.append(ModOpAddReserves(1.0));
                # model = modelProcessor.process()
                # model.printToXML('jra_testingModel.osim')
                # model.initSystem()
                
                # need the solution
                solution = osim.MocoTrajectory(os.path.join(trialdir, 'muscle_statetrack_GRFprescribe_solution_100con.sto'))
                statesTable = osim.TimeSeriesTable(os.path.join(trialdir, 'muscletrack_states_100con.sto'))
                controlsTable = osim.TimeSeriesTable(os.path.join(trialdir, 'muscletrack_controls_100con.sto'))
                fiberLengths = osim.TimeSeriesTable(os.path.join(trialdir, 'analyzemuscles_100conmuscletrack_MuscleAnalysis_FiberLength.sto'))
                
                # # get a table with just the jointset states
                # coordsTable = statesTable
                # statelabels = coordsTable.getColumnLabels()
                # for coord in range(len(statelabels)):
                #     coordinate = statelabels[coord]
                #     if 'jointset/' not in coordinate:
                #         coordsTable.removeColumn(coordinate)
                
                # get the times - first and last
                times = statesTable.getIndependentColumn()
                initTime = times[0]
                finalTime = times[-1]
                
                    
                # # get a table with the combined states and tendon lengths?
                # combStates = statesTable
                # comblabels = combStates.getColumnLabels()
                # # loop the states table and remove the tendon force columns
                # for col in range(len(comblabels)):
                #     columnname = comblabels[col]
                #     if 'tendon' in columnname:
                #         combStates.removeColumn(columnname)
                        
                # comblabels2 = combStates.getColumnLabels()
                
                
                # # get the names of all the fiberlengths
                # fiberlabels = fiberLengths.getColumnLabels()
                # # loop through them all, grab, and drop into states table, 
                # for fib in range(len(fiberlabels)):
                #     fiber = fiberlabels[fib]
                #     fiberlength = fiberLengths.getDependentColumn(fiber)
                #     # try to add it to the states table
                #     combStates.appendColumn('/forceset/'+fiber+'/fiber_length', fiberlength)
                
                # # store the model masses
                # modelmass = get_model_total_mass(trialdir, 'simple_model_all_the_probes_adjusted.osim')
                # naturalstruct_combine[subject] = modelmass


                ## create the model - add welds and then force to feet. 
                modelProcessor = osim.ModelProcessor('simple_model_all_the_probes_adjusted.osim')
                weldem = osim.StdVectorString()
                weldem.append('mtp_r'); weldem.append('mtp_l')
                modelProcessor.append(osim.ModOpReplaceJointsWithWelds(weldem))
                modelProcessor.append(osim.ModOpReplaceMusclesWithDeGrooteFregly2016());
                modelProcessor.append(osim.ModOpScaleActiveFiberForceCurveWidthDGF(1.5));
                modelProcessor.append(osim.ModOpFiberDampingDGF(0.01));
                # modelProcessor.append(ModOpTendonComplianceDynamicsModeDGF('implicit'));
                
                modelProcessor.append(osim.ModOpAddReserves(1.0))
                
                
                model = modelProcessor.process()
                # add in the external loads rather than specifying.
                force_expressed_in_body = ['ground', 'ground']
                point_identifier = ['rCOP_', 'lCOP_']
                point_expressed_in_body = ['ground', 'ground']
                force_identifier = ['rF_', 'lF_']
                applied_to_body = ['calcn_r', 'calcn_l']
                grfNames = ['GRF_r', 'GRF_l']
                dataSource = osim.Storage('ground_reaction.mot')
                for i in range(len(grfNames)):
                    newForce = osim.ExternalForce()
                    newForce.setName(grfNames[i])
                    newForce.set_applied_to_body(applied_to_body[i])
                    newForce.set_force_expressed_in_body(force_expressed_in_body[i])
                    newForce.set_force_identifier(force_identifier[i])
                    newForce.set_point_expressed_in_body(point_expressed_in_body[i]) ;
                    newForce.set_point_identifier(point_identifier[i]) ;
                    newForce.setDataSource(dataSource) ;
                    model.addForce(newForce) ;
                    
                
                
                model.initSystem()
                model.printToXML('jratestingmodel.osim')

                modelstates = model.getStateVariableNames();
                
                
                # testsolution = solution
                # statesTest = statesTable
                # statescols = statesTest.getColumnLabels()
                # fibertest = fiberLengths
                # fibercols = fibertest.getColumnLabels()

                # newfibercols = []
                # for i in range(len(fibercols)):
                #     print(fibercols[i])
                #     newfibercols.append('/forceset/' + fibercols[i] + '/fiber_length')
                # fibertest.setColumnLabels(newfibercols)
                # testsolution.insertStatesTrajectory(fibertest)
                # osim.STOFileAdapter.write(testsolution.exportToStatesTable(), 'testfibsolution.sto')

                
                # # grab a version of the table to manipulate
                # statestest = testsolution.exportToStatesTable()
                
                # loop the states from the model to 




                # # try the analysis
                jr_tool = osim.AnalyzeTool()
                jr_tool.setName('jr_analysis_100con')
                # jr_tool.setModelFilename(os.path.join(trialdir, 'post_simple_model_all_the_probes_muscletrack.osim'))
                
                                
                
                # jr_tool.setStatesFileName(os.path.join(trialdir, 'muscletrack_states_100con.sto'))
                statesStorage = osim.Storage('muscletrack_states_100con.sto')
                jr_tool.setStatesStorage(statesStorage)
                # jr_tool.setStatesFileName('testfibsolution.sto')

                # figure out if need the loads or not
                # jr_tool.setExternalLoadsFileName('grf_walk.xml')
                jr_tool.updControllerSet().cloneAndAppend(osim.PrescribedController(os.path.join(trialdir, 'muscletrack_controls_100con.sto')))
                jra = osim.JointReaction()
                jra.setName('jra')
                wherestr = osim.ArrayStr(); wherestr.append('child')
                jra.setInFrame(wherestr)
                
                jr_tool.updAnalysisSet().cloneAndAppend(jra)
                jr_tool.setInitialTime(initTime)
                jr_tool.setFinalTime(finalTime)
                jr_tool.setResultsDir(trialdir)
                
                
                # jr_tool.setModelFilename('jratestingmodel.osim')
                model.addAnalysis(jra)
                jr_tool.setModel(model)
                
                
                
                ## uncomment to rerun the analysis
                # jr_tool.printToXML(os.path.join(trialdir, 'jr_setup.xml'))
                # time.sleep(0.5)
                # jr_tool = osim.AnalyzeTool(os.path.join(trialdir, 'jr_setup.xml'))
                jr_tool.run()
                time.sleep(0.5)

                testtable = osim.TimeSeriesTable('jr_analysis_100con_jra_ReactionLoads.sto')
                testtib = testtable.getDependentColumn('walker_knee_r_on_tibia_r_in_tibia_r_fy').to_numpy()

                import matplotlib.pyplot as plt
                plt.figure()
                plt.plot(testtib)


                ### anmother one
                # Load your MocoTrajectory file
                trajectory_file = "muscle_statetrack_grfprescribe_solution_100con.sto"
                trajectory = osim.MocoTrajectory(trajectory_file)
                
                # Load your OpenSim model (with the knee joint defined)
                # model = opensim.Model("your_model.osim")
                
                # Get the knee joint
                knee_joint = model.getJointSet().get("walker_knee_r")  # Replace with your knee joint name
                
                # Create a storage to store the knee joint reaction forces
                forces_storage = opensim.Storage()
                
                # Get the time values from the trajectory
                time_column = trajectory.getTimeMat()
                
                # Loop through each time point in the trajectory
                for time in time_column:
                    # Get the states at the current time point
                    state = trajectory.getStatesTrajectory().get(time)
                
                    # Set the model state to the current state
                    model.setState(state)
                
                    # Calculate the knee joint reaction forces
                    knee_reaction_forces = knee_joint.getReactionForce(state)
                
                    # Append the knee reaction forces to the storage
                    forces_storage.append(time, knee_reaction_forces)
                
                # Save the reaction forces to a file (optional)
                forces_storage.printToFile("knee_reaction_forces.sto")

                
                



                # # jr = osim.analyzeSpatialVec(model, combStates, controlsTable, ['.*walker_knee.*reaction_on_parent.*'])

                # ## now to actually load in the data and do something with it. 
                # jras = osim.TimeSeriesTable(os.path.join(trialdir, 'jr_analysis_100con_jra_ReactionLoads.sto'))
                
                ############################
                # another one
                jramodel = osim.Model(os.path.join(trialdir, 'post_simple_model_all_the_probes_muscletrack.osim'))

                jraoutputs = osim.StdVectorString()
                jraoutputs.append('.*walker_knee.*\\|reaction_on_child')
                os.chdir(trialdir)
                jras = osim.analyzeMocoTrajectorySpatialVec(jramodel, solution, jraoutputs)
                # jr = osim.analyzeSpatialVec(model, combStates, controlsTable, ['.*walker_knee.*reaction_on_parent.*'])
                time.sleep(0.5)
                
                
                
                ## now to actually load in the data and do something with it. 
                jrastab = jras.flatten(['_mx', '_my', '_mz', '_fx', '_fy', '_fz'])
                jrasrx = jrastab.getDependentColumn('/jointset/walker_knee_r|reaction_on_child_fx').to_numpy()
                jrasry = jrastab.getDependentColumn('/jointset/walker_knee_r|reaction_on_child_fy').to_numpy()
                jrasrz = jrastab.getDependentColumn('/jointset/walker_knee_r|reaction_on_child_fz').to_numpy()
                jraslx = jrastab.getDependentColumn('/jointset/walker_knee_l|reaction_on_child_fx').to_numpy()
                jrasly = jrastab.getDependentColumn('/jointset/walker_knee_l|reaction_on_child_fy').to_numpy()
                jraslz = jrastab.getDependentColumn('/jointset/walker_knee_l|reaction_on_child_fz').to_numpy()
                
                ## now transform to the tibia frame from ground
                # get the joints
                transformmodel = jramodel
                transformmodel.finalizeConnections()
                joints = transformmodel.getJointSet()
                jointr = joints.get('walker_knee_r')
                jointl = joints.get('walker_knee_l')
                # get child bodies
                tibr = jointr.getChildFrame()
                tibr_name = tibr.getAbsolutePathString()
                tibl = jointl.getChildFrame()
                tibl_name = tibl.getAbsolutePathString()
                # also get the ground frame
                ground = transformmodel.getGround()
                # grab the time vector
                jrastime = jrastab.getIndependentColumn()
                
                # get the base state structure for the model
                state = transformmodel.initSystem()
                
                newjrasr = np.zeros((len(jrastime), 3))
                newjrasl = np.zeros((len(jrastime), 3))
                diffr = np.zeros((len(jrastime), 2))
                diffl = np.zeros((len(jrastime), 2))
                # loop through each time step, 
                coords = model.getCoordinateSet()
                # pdb.set_trace()
                for tim in range(len(solution.getTime())):
                    # then each coordinate to pose the model and get the transform
                    for cor in range(coords.getSize()):
                        cord = coords.get(cor)
                        statecol = statesTable.getDependentColumn(cord.getAbsolutePathString() + '/value').to_numpy()
                        cord.setValue(state, statecol[tim])
                        
                    # realize model position from poses
                    transformmodel.realizePosition(state)
                    # then compute the transform for the forces
                    newjrasr[tim,:] = ground.expressVectorInAnotherFrame(state, osim.Vec3(jrasrx[tim], jrasry[tim], jrasrz[tim]), tibr).to_numpy()
                    newjrasl[tim,:] = ground.expressVectorInAnotherFrame(state, osim.Vec3(jraslx[tim], jrasly[tim], jraslz[tim]), tibl).to_numpy()
                    diffr[tim,0] = np.linalg.norm(ground.expressVectorInAnotherFrame(state, osim.Vec3(jrasrx[tim], jrasry[tim], jrasrz[tim]), tibr).to_numpy())
                    diffr[tim,1] = np.linalg.norm([jrasrx[tim], jrasry[tim], jrasrz[tim]])
                    diffl[tim,0] = np.linalg.norm(ground.expressVectorInAnotherFrame(state, osim.Vec3(jraslx[tim], jrasly[tim], jraslz[tim]), tibl).to_numpy()) 
                    diffl[tim,1] = np.linalg.norm([jraslx[tim], jrasly[tim], jraslz[tim]])
                
                # pdb.set_trace()
                # #############
                # # results
                # ### TODO clean up based on analysis in the end
                
                # # want to focus on just the knees first
                # kneeReactions = jras
                # jracolumns = jras.getColumnLabels()
                # for cols in jracolumns:
                #     if 'walker' not in cols: # and '_f' not in cols:
                #         kneeReactions.removeColumn(cols)
                #     else:
                #         print(cols)
                # kneeLabels = kneeReactions.getColumnLabels()
                # print('\n\nOkay got the knees...\n')
                
                # for each in kneeLabels:
                #     if '_f' not in each:
                #         kneeReactions.removeColumn(each)
                #     else:
                #         print(each)
                # kneeLabels = kneeReactions.getColumnLabels()
                
                # # now have a table with just the knee reaction loads
                # # create a column that is the net force
                # fxr = kneeReactions.getDependentColumn('walker_knee_r_on_tibia_r_in_tibia_r_fx').to_numpy()
                # fyr = kneeReactions.getDependentColumn('walker_knee_r_on_tibia_r_in_tibia_r_fy').to_numpy()
                # fzr = kneeReactions.getDependentColumn('walker_knee_r_on_tibia_r_in_tibia_r_fz').to_numpy()
                # fxl = kneeReactions.getDependentColumn('walker_knee_l_on_tibia_l_in_tibia_l_fx').to_numpy()
                # fyl = kneeReactions.getDependentColumn('walker_knee_l_on_tibia_l_in_tibia_l_fy').to_numpy()
                # fzl = kneeReactions.getDependentColumn('walker_knee_l_on_tibia_l_in_tibia_l_fz').to_numpy()
                
                # norm_r = np.sqrt(fxr**2 + fyr**2 + fzr**2)
                # norm_l = np.sqrt(fxl**2 + fyl**2 + fzl**2)
                
                # pdb.set_trace()

                # welknaturalstruct[trial] = (times, norm_r, norm_l)
                # welknaturalstruct[trial] = (times, np.abs(fyr), np.abs(fyl))
                welknaturalstruct[trial] = (times, np.abs(newjrasr[:,1]), np.abs(newjrasl[:,1]))
                print(trialdir)

                
        # now add this subject struct to the larger combined one
        welkexostruct_combine[subject] = welkexostruct
        welknaturalstruct_combine[subject] = welknaturalstruct
            
    ###########################################################################
    ## analyze the joint reaction forces at the knee
    ###########################################################################
    import matplotlib.pyplot as plt
    fig1, ax1 = plt.subplots(1,7)
    count = 0
    # start with only the right legs
    n_holdr = np.zeros((4, 101))
    e_holdr = np.zeros((4, 101))
    n_holdl = np.zeros((4, 101))
    e_holdl = np.zeros((4, 101))
    
    n_avg_r = np.zeros(((len(trials)*len(welksubjects)), 101))
    n_avg_l = np.zeros(((len(trials)*len(welksubjects)), 101))
    e_avg_r = np.zeros(((len(trials)*len(welksubjects)), 101))
    e_avg_l = np.zeros(((len(trials)*len(welksubjects)), 101))
    pop = 0
    
    ## percent gait cycle for each subject
    for key in welkexostruct_combine.keys():
        # loop subjects first
        print(key)
        
        for tri in welkexostruct_combine[key].keys():
            print(tri)
            # loop through each trial and convert them all to an interpolated 
            # vector that is percent gait cycle
            e_temptime = np.array(welkexostruct_combine[key][tri][0])
            e_tempnormr = np.array(welkexostruct_combine[key][tri][1])
            e_tempnorml = np.array(welkexostruct_combine[key][tri][2])
            
            # now interpolate them all 
            e_timespercent = (e_temptime - e_temptime[0]) / (e_temptime[-1] - e_temptime[0])*100
            e_timespercent101 = np.arange(0,101,1)
            e_normr_interp = np.interp(e_timespercent101, e_timespercent, e_tempnormr)
            e_norml_interp = np.interp(e_timespercent101, e_timespercent, e_tempnorml)
            
            # rewrite struct with percent values and vectors
            welkexostruct_combine[key][tri] = (e_timespercent101, e_normr_interp, e_norml_interp)
            
            # these will be averaged over trials later
            e_holdr[int(tri[-1])-1,:] = e_normr_interp
            e_holdl[int(tri[-1])-1,:] = e_norml_interp
            # the full struct for all the subjects
            e_avg_r[pop, :] = e_normr_interp/exostruct_combine[key]
            e_avg_l[pop, :] = e_norml_interp/exostruct_combine[key]
            
            # and the natural counterparts
            n_temptime = np.array(welknaturalstruct_combine[key][tri][0])
            n_tempnormr = np.array(welknaturalstruct_combine[key][tri][1])
            n_tempnorml = np.array(welknaturalstruct_combine[key][tri][2])
    
            n_timespercent = (n_temptime - n_temptime[0]) / (n_temptime[-1] - n_temptime[0])*100
            n_timespercent101 = np.arange(0,101,1)
            n_normr_interp = np.interp(n_timespercent101, n_timespercent, n_tempnormr)
            n_norml_interp = np.interp(n_timespercent101, n_timespercent, n_tempnorml)
            
            # rewrite struct with percent values and vectors
            welknaturalstruct_combine[key][tri] = (n_timespercent101, n_normr_interp, n_norml_interp)
            
            # these will be averaged over trials later
            n_holdr[int(tri[-1])-1,:] = n_normr_interp
            n_holdl[int(tri[-1])-1,:] = n_norml_interp
            # population one
            n_avg_r[pop, :] = n_normr_interp/naturalstruct_combine[key]
            n_avg_l[pop, :] = n_norml_interp/naturalstruct_combine[key]
            pop += 1

        # do the averaging for the average for each subject
        exostruct_avg[key] = [e_holdr.mean(0), e_holdl.mean(0)]
        naturalstruct_avg[key] = [n_holdr.mean(0), n_holdl.mean(0)]
            
        # plot the natural and the exo cases for each subject, 
        ax1[count].plot(n_timespercent101, n_holdr.mean(0), color='orange', label='natural_r')
        ax1[count].plot(n_timespercent101, n_holdl.mean(0), color='orange', label='natural_l')
        ax1[count].plot(e_timespercent101, e_holdr.mean(0), color='purple', label='exotendon_r')
        ax1[count].plot(e_timespercent101, e_holdl.mean(0), color='purple', label='exotendon_l')

        count += 1
        
    plt.show()
    
    ## take a population average with some body normalization
    plt.figure(figsize=(11,8), dpi=300)
    plt.plot(n_timespercent101, n_avg_r.mean(0), color='#e66101', label='natural_r')
    plt.plot(n_timespercent101, n_avg_l.mean(0), color='#fdb863', label='natural_l')
    plt.plot(e_timespercent101, e_avg_r.mean(0), color='#5e3c99', label='exotendon_r')
    plt.plot(e_timespercent101, e_avg_l.mean(0), color='#b2abd2', label='exotendon_l')
    
    plt.ylabel('BW', fontsize=24)
    plt.yticks(fontsize=24)
    plt.xlabel('% Gait cycle', fontsize=24)
    plt.xticks(fontsize=24)
    plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0, fontsize=20)
    plt.title('Vertical Knee Contact Force', fontsize=24)
    plt.show()
    
    
    #################################################################
    # stats for the total knee contact forces
    
    # first do stats on the peak vertical knee contact force
    npeaksvr = np.max(n_avg_r, 1)
    epeaksvr = np.max(e_avg_r, 1)
    # take these values into excel sheet to play around/test normality
    # then do a t-test (see the excel sheet for all the stats)
    '''
    
    ################################################################
    # do more analysis! yay 
    ################################################################
    # think about segmenting between the muscles and their individual contributions
    # figure out how to get each of the muscles and do the analysis for hip contact of just that muscle
    # need to figure out how to do joint reaction and isolate just the muscles
    
    # 1) load the model, load a verson of the model, set passives to zero
        # then set all controls but one we want to zero
        # then all activation states to zero that we don't want
    # data for y: vertical forces
    einterseg_combine = np.zeros((len(welksubjects)*len(trials), 101))
    ninterseg_combine = np.zeros((len(welksubjects)*len(trials), 101))
    eflexors_combine    = np.zeros((len(welksubjects)*len(trials), 101))
    nflexors_combine    = np.zeros((len(welksubjects)*len(trials), 101))
    eadductors_combine      = np.zeros((len(welksubjects)*len(trials), 101))
    nadductors_combine      = np.zeros((len(welksubjects)*len(trials), 101))
    eabductors_combine      = np.zeros((len(welksubjects)*len(trials), 101))
    nabductors_combine      = np.zeros((len(welksubjects)*len(trials), 101))
    eextensors_combine      = np.zeros((len(welksubjects)*len(trials), 101))
    nextensors_combine      = np.zeros((len(welksubjects)*len(trials), 101))
    ehams_combine     = np.zeros((len(welksubjects)*len(trials), 101))
    nhams_combine     = np.zeros((len(welksubjects)*len(trials), 101))
    eall_combine      = np.zeros((len(welksubjects)*len(trials), 101))
    nall_combine      = np.zeros((len(welksubjects)*len(trials), 101))
    ereserve_combine  = np.zeros((len(welksubjects)*len(trials), 101))
    nreserve_combine  = np.zeros((len(welksubjects)*len(trials), 101))
    enone_combine     = np.zeros((len(welksubjects)*len(trials), 101))
    nnone_combine     = np.zeros((len(welksubjects)*len(trials), 101))
    
    # data for x: anterior-posterior forces
    einterseg_combine_x = np.zeros((len(welksubjects)*len(trials), 101))
    ninterseg_combine_x = np.zeros((len(welksubjects)*len(trials), 101))
    eflexors_combine_x    = np.zeros((len(welksubjects)*len(trials), 101))
    nflexors_combine_x    = np.zeros((len(welksubjects)*len(trials), 101))
    eadductors_combine_x      = np.zeros((len(welksubjects)*len(trials), 101))
    nadductors_combine_x      = np.zeros((len(welksubjects)*len(trials), 101))
    eabductors_combine_x      = np.zeros((len(welksubjects)*len(trials), 101))
    nabductors_combine_x      = np.zeros((len(welksubjects)*len(trials), 101))
    eextensors_combine_x      = np.zeros((len(welksubjects)*len(trials), 101))
    nextensors_combine_x      = np.zeros((len(welksubjects)*len(trials), 101))
    ehams_combine_x     = np.zeros((len(welksubjects)*len(trials), 101))
    nhams_combine_x     = np.zeros((len(welksubjects)*len(trials), 101))
    eall_combine_x      = np.zeros((len(welksubjects)*len(trials), 101))
    nall_combine_x      = np.zeros((len(welksubjects)*len(trials), 101))
    ereserve_combine_x  = np.zeros((len(welksubjects)*len(trials), 101))
    nreserve_combine_x  = np.zeros((len(welksubjects)*len(trials), 101))
    enone_combine_x     = np.zeros((len(welksubjects)*len(trials), 101))
    nnone_combine_x     = np.zeros((len(welksubjects)*len(trials), 101))

    # data for z: medial-lateral forces
    einterseg_combine_z = np.zeros((len(welksubjects)*len(trials), 101))
    ninterseg_combine_z = np.zeros((len(welksubjects)*len(trials), 101))
    eflexors_combine_z    = np.zeros((len(welksubjects)*len(trials), 101))
    nflexors_combine_z    = np.zeros((len(welksubjects)*len(trials), 101))
    eadductors_combine_z      = np.zeros((len(welksubjects)*len(trials), 101))
    nadductors_combine_z      = np.zeros((len(welksubjects)*len(trials), 101))
    eabductors_combine_z      = np.zeros((len(welksubjects)*len(trials), 101))
    nabductors_combine_z      = np.zeros((len(welksubjects)*len(trials), 101))
    eextensors_combine_z      = np.zeros((len(welksubjects)*len(trials), 101))
    nextensors_combine_z      = np.zeros((len(welksubjects)*len(trials), 101))
    ehams_combine_z     = np.zeros((len(welksubjects)*len(trials), 101))
    nhams_combine_z     = np.zeros((len(welksubjects)*len(trials), 101))
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
        musclesWanted['hipFlexors'] = ['recfem_r', 'psoas_r', 'iliacus_r']
        musclesWanted['hipAdductors']   = ['addbrev_r','addlong_r','addmagDist_r','addmagIsch_r','addmagMid_r','addmagProx_r']
        musclesWanted['hipAbductors']   = ['glmed1_r','glmed2_r','glmed3_r','glmin1_r','glmin2_r','glmin3_r','piri_r','tfl_r']
        musclesWanted['hipExtensors']   = ['glmax1_r','glmax2_r','glmax3_r']
        musclesWanted['hams']  = ['bflh_r', 'grac_r', 'sart_r', 'semimem_r', 'semiten_r']
        musclesWanted['all'] = ['all']
        musclesWanted['reserve'] = ['reserve']
        musclesWanted['none'] = ['none']
    elif whichleg == 'left':
        musclesWanted['inter'] = []
        musclesWanted['hipFlexors'] = ['recfem_l', 'psoas_l', 'iliacus_l']
        musclesWanted['hipAdductors']   = ['addbrev_l','addlong_l','addmagDist_l','addmagIsch_l','addmagMid_l','addmagProx_l']
        musclesWanted['hipAbductors']   = ['glmed1_l','glmed2_l','glmed3_l','glmin1_l','glmin2_l','glmin3_l','piri_l','tfl_l']
        musclesWanted['hipExtensors']   = ['glmax1_l','glmax2_l','glmax3_l']
        musclesWanted['hams']  = ['bflh_l', 'grac_l', 'sart_l', 'semimem_l', 'semiten_l']
        musclesWanted['all'] = ['all']
        musclesWanted['reserve'] = ['reserve']
        musclesWanted['none'] = ['none']
    elif whichleg == 'both':
        musclesWanted['inter'] = []
        musclesWanted['hipFlexors'] = ['recfem_r', 'psoas_r', 'iliacus_r', 'recfem_l', 'psoas_l', 'iliacus_l']
        musclesWanted['hipAdductors']   = ['addbrev_r','addlong_r','addmagDist_r','addmagIsch_r','addmagMid_r','addmagProx_r','addbrev_l','addlong_l','addmagDist_l','addmagIsch_l','addmagMid_l','addmagProx_l']
        musclesWanted['hipAbductors']   = ['glmed1_r','glmed2_r','glmed3_r','glmin1_r','glmin2_r','glmin3_r','piri_r','tfl_r','glmed1_l','glmed2_l','glmed3_l','glmin1_l','glmin2_l','glmin3_l','piri_l','tfl_l']
        musclesWanted['hipExtensors']   = ['glmax1_r','glmax2_r','glmax3_r','glmax1_l','glmax2_l','glmax3_l']
        musclesWanted['hams']  = ['bflh_r', 'grac_r', 'sart_r', 'semimem_r', 'semiten_r','bflh_l', 'grac_l', 'sart_l', 'semimem_l', 'semiten_l']
        musclesWanted['all'] = ['all']
        musclesWanted['reserve'] = ['reserve']
        musclesWanted['none'] = ['none']
    # TODO
    
    # plot the stuff
    # figure out subtracting the intersegmental from each of the others?
    # pdb.set_trace()

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
                # write a method that gives the hip contact forces, 
                # pass in model, and the muscles that you want
                # pdb.set_trace()
                os.chdir(trialdir)                
                

                try:
                    if oldnotredo:
                        ## okay now going to focus on the figures that I actually wanted 
                        jrasrflexors = getHipContactributions(trialdir, musclesWanted['hipFlexors'], 'hipFlexors', whichleg, runtool)
                        jrasradductors = getHipContactributions(trialdir, musclesWanted['hipAdductors'], 'hipAdductors', whichleg, runtool)
                        jrasrabductors = getHipContactributions(trialdir, musclesWanted['hipAbductors'], 'hipAbductors', whichleg, runtool)
                        jrasrextensors = getHipContactributions(trialdir, musclesWanted['hipExtensors'], 'hipExtensors', whichleg, runtool)
                        jrasrhams = getHipContactributions(trialdir, musclesWanted['hams'], 'hams', whichleg, runtool)
                        jrasrinter = getHipContactributions(trialdir, musclesWanted['inter'], 'inter', whichleg, runtool)
                        jrasrall = getHipContactributions(trialdir, musclesWanted['all'], 'all', whichleg, runtool)
                        jrasrinterreserve = getHipContactributions(trialdir, musclesWanted['reserve'], 'reserve', whichleg, runtool)
                        jrasrnone = getHipContactributions(trialdir, musclesWanted['none'], 'none', whichleg, runtool)
                    else:
                        ## okay now going to focus on the figures that I actually wanted 
                        jrasrflexorsx, jrasrflexors, jrasrflexorsz = getHipContactributionsRedo(trialdir, musclesWanted['hipFlexors'], 'hipFlexors', whichleg, runtool)
                        jrasradductorsx, jrasradductors, jrasradductorsz = getHipContactributionsRedo(trialdir, musclesWanted['hipAdductors'], 'hipAdductors', whichleg, runtool)
                        jrasrabductorsx, jrasrabductors, jrasrabductorsz = getHipContactributionsRedo(trialdir, musclesWanted['hipAbductors'], 'hipAbductors', whichleg, runtool)
                        jrasrextensorsx, jrasrextensors, jrasrextensorsz = getHipContactributionsRedo(trialdir, musclesWanted['hipExtensors'], 'hipExtensors', whichleg, runtool)
                        jrasrhamsx, jrasrhams, jrasrhamsz = getHipContactributionsRedo(trialdir, musclesWanted['hams'], 'hams', whichleg, runtool)
                        jrasrinterx, jrasrinter, jrasrinterz = getHipContactributionsRedo(trialdir, musclesWanted['inter'], 'inter', whichleg, runtool)
                        jrasrallx, jrasrall, jrasrallz = getHipContactributionsRedo(trialdir, musclesWanted['all'], 'all', whichleg, runtool)
                        jrasrinterreservex, jrasrinterreserve, jrasrinterreservez = getHipContactributionsRedo(trialdir, musclesWanted['reserve'], 'reserve', whichleg, runtool)
                        jrasrnonex, jrasrnone, jrasrnonez = getHipContactributionsRedo(trialdir, musclesWanted['none'], 'none', whichleg, runtool)
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

                jrasrflexorsonly = jrasrflexors - jrasrnone
                jrasradductorsonly = jrasradductors - jrasrnone
                jrasrabductorsonly = jrasrabductors - jrasrnone                
                jrasrextensorsonly = jrasrextensors - jrasrnone                
                jrasrhamsonly = jrasrhams - jrasrnone

                # get percentages
                e_timespercent101 = np.arange(0,101,1)
                e_times = np.arange(0,len(jrasrflexorsonly),1)
                e_timesinterp = np.linspace(0,len(e_times), 103)

                # get something in BW and interp to 100% gait cycle points. 
                jrasrflexorsonly101 = -1*(np.interp(e_timesinterp, e_times, jrasrflexorsonly)) / (modelmass * 9.81)
                jrasradductorsonly101 = -1*(np.interp(e_timesinterp, e_times, jrasradductorsonly)) / (modelmass * 9.81)
                jrasrabductorsonly101 = -1*(np.interp(e_timesinterp, e_times, jrasrabductorsonly)) / (modelmass * 9.81)
                jrasrextensorsonly101 = -1*(np.interp(e_timesinterp, e_times, jrasrextensorsonly)) / (modelmass * 9.81)
                jrasrhamsonly101 = -1*(np.interp(e_timesinterp, e_times, jrasrhamsonly)) / (modelmass * 9.81)
                jrasrinteronly101 = -1*(np.interp(e_timesinterp, e_times, jrasrinteronly)) / (modelmass * 9.81)
                jrasrallonly101 = -1*(np.interp(e_timesinterp, e_times, jrasrall)) / (modelmass * 9.81)
                jrasrreserveonly101 = -1*(np.interp(e_timesinterp, e_times, jrasrreserveonly)) / (modelmass * 9.81)
                jrasrnoneonly101 = -1*(np.interp(e_timesinterp, e_times, jrasrnone)) / (modelmass * 9.81)
                
                # data for y: vertical forces
                einterseg_combine[spot, :] = jrasrinteronly101[:-2]
                eflexors_combine[spot, :] = jrasrflexorsonly101[:-2]
                eadductors_combine[spot,:] = jrasradductorsonly101[:-2]
                eabductors_combine[spot,:] = jrasrabductorsonly101[:-2]
                eextensors_combine[spot,:] = jrasrextensorsonly101[:-2]
                ehams_combine[spot,:] = jrasrhamsonly101[:-2]
                eall_combine[spot,:] = jrasrallonly101[:-2]
                ereserve_combine[spot,:] = jrasrreserveonly101[:-2]
                enone_combine[spot,:] = jrasrnoneonly101[:-2]
                

                ### X: anterior-posterior forces
                # subtract out the inter segmental
                jrasrinteronlyx = jrasrinterreservex
                jrasrreserveonlyx = jrasrinterx - jrasrinterreservex
                jrasrnoneonlyx = jrasrnonex

                jrasrflexorsonlyx = jrasrflexorsx - jrasrnonex
                jrasradductorsonlyx = jrasradductorsx - jrasrnonex
                jrasrabductorsonlyx = jrasrabductorsx - jrasrnonex                
                jrasrextensorsonlyx = jrasrextensorsx - jrasrnonex                
                jrasrhamsonlyx = jrasrhamsx - jrasrnonex

                # get percentages
                e_timespercent101 = np.arange(0,101,1)
                e_timesx = np.arange(0,len(jrasrflexorsonlyx),1)
                e_timesinterpx = np.linspace(0,len(e_timesx), 103)

                # get something in BW and interp to 100% gait cycle points. 
                jrasrflexorsonly101x = -1*(np.interp(e_timesinterpx, e_timesx, jrasrflexorsonlyx)) / (modelmass * 9.81)
                jrasradductorsonly101x = -1*(np.interp(e_timesinterpx, e_timesx, jrasradductorsonlyx)) / (modelmass * 9.81)
                jrasrabductorsonly101x = -1*(np.interp(e_timesinterpx, e_timesx, jrasrabductorsonlyx)) / (modelmass * 9.81)
                jrasrextensorsonly101x = -1*(np.interp(e_timesinterpx, e_timesx, jrasrextensorsonlyx)) / (modelmass * 9.81)
                jrasrhamsonly101x = -1*(np.interp(e_timesinterpx, e_timesx, jrasrhamsonlyx)) / (modelmass * 9.81)
                jrasrinteronly101x = -1*(np.interp(e_timesinterpx, e_timesx, jrasrinteronlyx)) / (modelmass * 9.81)
                jrasrallonly101x = -1*(np.interp(e_timesinterpx, e_timesx, jrasrallx)) / (modelmass * 9.81)
                jrasrreserveonly101x = -1*(np.interp(e_timesinterpx, e_timesx, jrasrreserveonlyx)) / (modelmass * 9.81)
                jrasrnoneonly101x = -1*(np.interp(e_timesinterpx, e_timesx, jrasrnonex)) / (modelmass * 9.81)

                # data for x: anterior-posterior forces
                einterseg_combine_x[spot, :] = jrasrinteronly101x[:-2]
                eflexors_combine_x[spot, :] = jrasrflexorsonly101x[:-2]
                eadductors_combine_x[spot,:] = jrasradductorsonly101x[:-2]
                eabductors_combine_x[spot,:] = jrasrabductorsonly101x[:-2]
                eextensors_combine_x[spot,:] = jrasrextensorsonly101x[:-2]
                ehams_combine_x[spot,:] = jrasrhamsonly101x[:-2]
                eall_combine_x[spot,:] = jrasrallonly101x[:-2]
                ereserve_combine_x[spot,:] = jrasrreserveonly101x[:-2]
                enone_combine_x[spot,:] = jrasrnoneonly101x[:-2]
                

                ### Z: medial-lateral forces
                # subtract out the inter segmental
                jrasrinteronlyz = jrasrinterreservez
                jrasrreserveonlyz = jrasrinterz - jrasrinterreservez
                jrasrnoneonlyz = jrasrnonez

                jrasrflexorsonlyz = jrasrflexorsz - jrasrnonez
                jrasradductorsonlyz = jrasradductorsz - jrasrnonez
                jrasrabductorsonlyz = jrasrabductorsz - jrasrnonez                
                jrasrextensorsonlyz = jrasrextensorsz - jrasrnonez                
                jrasrhamsonlyz = jrasrhamsz - jrasrnonez

                # get percentages
                e_timespercent101 = np.arange(0,101,1)
                e_timesz = np.arange(0,len(jrasrflexorsonlyz),1)
                e_timesinterpz = np.linspace(0,len(e_timesz), 103)

                # get something in BW and interp to 100% gait cycle points. 
                jrasrflexorsonly101z = -1*(np.interp(e_timesinterpz, e_timesz, jrasrflexorsonlyz)) / (modelmass * 9.81)
                jrasradductorsonly101z = -1*(np.interp(e_timesinterpz, e_timesz, jrasradductorsonlyz)) / (modelmass * 9.81)
                jrasrabductorsonly101z = -1*(np.interp(e_timesinterpz, e_timesz, jrasrabductorsonlyz)) / (modelmass * 9.81)
                jrasrextensorsonly101z = -1*(np.interp(e_timesinterpz, e_timesz, jrasrextensorsonlyz)) / (modelmass * 9.81)
                jrasrhamsonly101z = -1*(np.interp(e_timesinterpz, e_timesz, jrasrhamsonlyz)) / (modelmass * 9.81)
                jrasrinteronly101z = -1*(np.interp(e_timesinterpz, e_timesz, jrasrinteronlyz)) / (modelmass * 9.81)
                jrasrallonly101z = -1*(np.interp(e_timesinterpz, e_timesz, jrasrallz)) / (modelmass * 9.81)
                jrasrreserveonly101z = -1*(np.interp(e_timesinterpz, e_timesz, jrasrreserveonlyz)) / (modelmass * 9.81)
                jrasrnoneonly101z = -1*(np.interp(e_timesinterpz, e_timesz, jrasrnonez)) / (modelmass * 9.81)
                
                # data for z: medial-lateral forces
                einterseg_combine_z[spot, :] = jrasrinteronly101z[:-2]
                eflexors_combine_z[spot, :] = jrasrflexorsonly101z[:-2]
                eadductors_combine_z[spot,:] = jrasradductorsonly101z[:-2]
                eabductors_combine_z[spot,:] = jrasrabductorsonly101z[:-2]
                eextensors_combine_z[spot,:] = jrasrextensorsonly101z[:-2]
                ehams_combine_z[spot,:] = jrasrhamsonly101z[:-2]
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
                # write a method that gives the hip contact forces, 
                # pass in model, and the muscles that you want
                os.chdir(trialdir)                
                # pdb.set_trace()
                
                
                try:
                    if oldnotredo:
                        ## okay now going to focus on the figures that I actually wanted 
                        jrasrflexors = getHipContactributions(trialdir, musclesWanted['hipFlexors'], 'hipFlexors', whichleg, runtool)
                        jrasradductors = getHipContactributions(trialdir, musclesWanted['hipAdductors'], 'hipAdductors', whichleg, runtool)
                        jrasrabductors = getHipContactributions(trialdir, musclesWanted['hipAbductors'], 'hipAbductors', whichleg, runtool)
                        jrasrextensors = getHipContactributions(trialdir, musclesWanted['hipExtensors'], 'hipExtensors', whichleg, runtool)
                        jrasrhams = getHipContactributions(trialdir, musclesWanted['hams'], 'hams', whichleg, runtool)
                        jrasrinter = getHipContactributions(trialdir, musclesWanted['inter'], 'inter', whichleg, runtool)
                        jrasrall = getHipContactributions(trialdir, musclesWanted['all'], 'all', whichleg, runtool)
                        jrasrinterreserve = getHipContactributions(trialdir, musclesWanted['reserve'], 'reserve', whichleg, runtool)
                        jrasrnone = getHipContactributions(trialdir, musclesWanted['none'], 'none', whichleg, runtool)
                    else:
                        ## okay now going to focus on the figures that I actually wanted 
                        jrasrflexorsx, jrasrflexors, jrasrflexorsz = getHipContactributionsRedo(trialdir, musclesWanted['hipFlexors'], 'hipFlexors', whichleg, runtool)
                        jrasradductorsx, jrasradductors, jrasradductorsz = getHipContactributionsRedo(trialdir, musclesWanted['hipAdductors'], 'hipAdductors', whichleg, runtool)
                        jrasrabductorsx, jrasrabductors, jrasrabductorsz = getHipContactributionsRedo(trialdir, musclesWanted['hipAbductors'], 'hipAbductors', whichleg, runtool)
                        jrasrextensorsx, jrasrextensors, jrasrextensorsz = getHipContactributionsRedo(trialdir, musclesWanted['hipExtensors'], 'hipExtensors', whichleg, runtool)
                        jrasrhamsx, jrasrhams, jrasrhamsz = getHipContactributionsRedo(trialdir, musclesWanted['hams'], 'hams', whichleg, runtool)
                        jrasrinterx, jrasrinter, jrasrinterz = getHipContactributionsRedo(trialdir, musclesWanted['inter'], 'inter', whichleg, runtool)
                        jrasrallx, jrasrall, jrasrallz = getHipContactributionsRedo(trialdir, musclesWanted['all'], 'all', whichleg, runtool)
                        jrasrinterreservex, jrasrinterreserve, jrasrinterreservez = getHipContactributionsRedo(trialdir, musclesWanted['reserve'], 'reserve', whichleg, runtool)
                        jrasrnonex, jrasrnone, jrasrnonez = getHipContactributionsRedo(trialdir, musclesWanted['none'], 'none', whichleg, runtool)
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
                jrasrflexorsonly = jrasrflexors - jrasrnone
                jrasradductorsonly = jrasradductors - jrasrnone
                jrasrabductorsonly = jrasrabductors - jrasrnone                
                jrasrextensorsonly = jrasrextensors - jrasrnone                
                jrasrhamsonly = jrasrhams - jrasrnone

                # get percentages
                n_timespercent101 = np.arange(0,101,1)
                n_times = np.arange(0,len(jrasrflexorsonly),1)
                n_timesinterp = np.linspace(0,len(n_times), 103)

                # get something in BW and interp to 100% gait cycle points. 
                jrasrflexorsonly101 = -1*(np.interp(n_timesinterp, n_times, jrasrflexorsonly)) / (modelmass * 9.81)
                jrasradductorsonly101 = -1*(np.interp(n_timesinterp, n_times, jrasradductorsonly)) / (modelmass * 9.81)
                jrasrabductorsonly101 = -1*(np.interp(n_timesinterp, n_times, jrasrabductorsonly)) / (modelmass * 9.81)
                jrasrextensorsonly101 = -1*(np.interp(n_timesinterp, n_times, jrasrextensorsonly)) / (modelmass * 9.81)
                jrasrhamsonly101 = -1*(np.interp(n_timesinterp, n_times, jrasrhamsonly)) / (modelmass * 9.81)
                jrasrinteronly101 = -1*(np.interp(n_timesinterp, n_times, jrasrinteronly)) / (modelmass * 9.81)
                jrasrallonly101 = -1*(np.interp(n_timesinterp, n_times, jrasrall)) / (modelmass * 9.81)
                jrasrreserveonly101 = -1*(np.interp(n_timesinterp, n_times, jrasrreserveonly)) / (modelmass * 9.81)
                jrasrnoneonly101 = -1*(np.interp(n_timesinterp, n_times, jrasrnone)) / (modelmass * 9.81)
                
                # natural combine into the big structure
                ninterseg_combine[spot, :] = jrasrinteronly101[:-2]
                nflexors_combine[spot, :] = jrasrflexorsonly101[:-2]
                nadductors_combine[spot,:] = jrasradductorsonly101[:-2]
                nabductors_combine[spot,:] = jrasrabductorsonly101[:-2]
                nextensors_combine[spot,:] = jrasrextensorsonly101[:-2]
                nhams_combine[spot,:] = jrasrhamsonly101[:-2]
                nall_combine[spot,:] = jrasrallonly101[:-2]
                nreserve_combine[spot,:] = jrasrreserveonly101[:-2]
                nnone_combine[spot,:] = jrasrnoneonly101[:-2]
                

                ### X: anterior-posterior forces
                # subtract out the inter segmental
                jrasrinteronlyx = jrasrinterreservex
                jrasrreserveonlyx = jrasrinterx - jrasrinterreservex
                jrasrnoneonlyx = jrasrnonex

                jrasrflexorsonlyx = jrasrflexorsx - jrasrnonex
                jrasradductorsonlyx = jrasradductorsx - jrasrnonex
                jrasrabductorsonlyx = jrasrabductorsx - jrasrnonex                
                jrasrextensorsonlyx = jrasrextensorsx - jrasrnonex                
                jrasrhamsonlyx = jrasrhamsx - jrasrnonex

                # get percentages
                n_timespercent101 = np.arange(0,101,1)
                n_timesx = np.arange(0,len(jrasrflexorsonlyx),1)
                n_timesinterpx = np.linspace(0,len(n_timesx), 103)

                # get something in BW and interp to 100% gait cycle points. 
                jrasrflexorsonly101x = -1*(np.interp(n_timesinterpx, n_timesx, jrasrflexorsonlyx)) / (modelmass * 9.81)
                jrasradductorsonly101x = -1*(np.interp(n_timesinterpx, n_timesx, jrasradductorsonlyx)) / (modelmass * 9.81)
                jrasrabductorsonly101x = -1*(np.interp(n_timesinterpx, n_timesx, jrasrabductorsonlyx)) / (modelmass * 9.81)
                jrasrextensorsonly101x = -1*(np.interp(n_timesinterpx, n_timesx, jrasrextensorsonlyx)) / (modelmass * 9.81)
                jrasrhamsonly101x = -1*(np.interp(n_timesinterpx, n_timesx, jrasrhamsonlyx)) / (modelmass * 9.81)
                jrasrinteronly101x = -1*(np.interp(n_timesinterpx, n_timesx, jrasrinteronlyx)) / (modelmass * 9.81)
                jrasrallonly101x = -1*(np.interp(n_timesinterpx, n_timesx, jrasrallx)) / (modelmass * 9.81)
                jrasrreserveonly101x = -1*(np.interp(n_timesinterpx, n_timesx, jrasrreserveonlyx)) / (modelmass * 9.81)
                jrasrnoneonly101x = -1*(np.interp(n_timesinterpx, n_timesx, jrasrnonex)) / (modelmass * 9.81)
                
                # natural combine into the big structure
                ninterseg_combine_x[spot, :] = jrasrinteronly101x[:-2]
                nflexors_combine_x[spot, :] = jrasrflexorsonly101x[:-2]
                nadductors_combine_x[spot,:] = jrasradductorsonly101x[:-2]
                nabductors_combine_x[spot,:] = jrasrabductorsonly101x[:-2]
                nextensors_combine_x[spot,:] = jrasrextensorsonly101x[:-2]
                nhams_combine_x[spot,:] = jrasrhamsonly101x[:-2]
                nall_combine_x[spot,:] = jrasrallonly101x[:-2]
                nreserve_combine_x[spot,:] = jrasrreserveonly101x[:-2]
                nnone_combine_x[spot,:] = jrasrnoneonly101x[:-2]
                

                ### Z: medial-lateral forces
                # subtract out the inter segmental
                jrasrinteronlyz = jrasrinterreservez
                jrasrreserveonlyz = jrasrinterz - jrasrinterreservez
                jrasrnoneonlyz = jrasrnonez

                jrasrflexorsonlyz = jrasrflexorsz - jrasrnonez
                jrasradductorsonlyz = jrasradductorsz - jrasrnonez
                jrasrabductorsonlyz = jrasrabductorsz - jrasrnonez                
                jrasrextensorsonlyz = jrasrextensorsz - jrasrnonez                
                jrasrhamsonlyz = jrasrhamsz - jrasrnonez

                # get percentages
                n_timespercent101 = np.arange(0,101,1)
                n_timesz = np.arange(0,len(jrasrflexorsonlyz),1)
                n_timesinterpz = np.linspace(0,len(n_timesz), 103)

                # get something in BW and interp to 100% gait cycle points. 
                jrasrflexorsonly101z = -1*(np.interp(n_timesinterpz, n_timesz, jrasrflexorsonlyz)) / (modelmass * 9.81)
                jrasradductorsonly101z = -1*(np.interp(n_timesinterpz, n_timesz, jrasradductorsonlyz)) / (modelmass * 9.81)
                jrasrabductorsonly101z = -1*(np.interp(n_timesinterpz, n_timesz, jrasrabductorsonlyz)) / (modelmass * 9.81)
                jrasrextensorsonly101z = -1*(np.interp(n_timesinterpz, n_timesz, jrasrextensorsonlyz)) / (modelmass * 9.81)
                jrasrhamsonly101z = -1*(np.interp(n_timesinterpz, n_timesz, jrasrhamsonlyz)) / (modelmass * 9.81)
                jrasrinteronly101z = -1*(np.interp(n_timesinterpz, n_timesz, jrasrinteronlyz)) / (modelmass * 9.81)
                jrasrallonly101z = -1*(np.interp(n_timesinterpz, n_timesz, jrasrallz)) / (modelmass * 9.81)
                jrasrreserveonly101z = -1*(np.interp(n_timesinterpz, n_timesz, jrasrreserveonlyz)) / (modelmass * 9.81)
                jrasrnoneonly101z = -1*(np.interp(n_timesinterpz, n_timesz, jrasrnonez)) / (modelmass * 9.81)
                
                # natural combine into the big structure
                ninterseg_combine_z[spot, :] = jrasrinteronly101z[:-2]
                nflexors_combine_z[spot, :] = jrasrflexorsonly101z[:-2]
                nadductors_combine_z[spot,:] = jrasradductorsonly101z[:-2]
                nabductors_combine_z[spot,:] = jrasrabductorsonly101z[:-2]
                nextensors_combine_z[spot,:] = jrasrextensorsonly101z[:-2]
                nhams_combine_z[spot,:] = jrasrhamsonly101z[:-2]
                nall_combine_z[spot,:] = jrasrallonly101z[:-2]
                nreserve_combine_z[spot,:] = jrasrreserveonly101z[:-2]
                nnone_combine_z[spot,:] = jrasrnoneonly101z[:-2]
                ## increase the spot - count of trials
                spot += 1



    ########################################################################################
    # plotting for activations, moments, etc. muscle insights for natural and exotendon
    ########################################################################################

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
    plt.savefig(os.path.join(analyzedir, 'hip_muscleactivations_' + whichleg + '.png'))


    # create a figure for the joint moments for natural and exotendon
    fig2, ax2 = plt.subplots(3, 8, figsize=(20, 8), dpi=500)
    joints = list(moments_nat.keys())
    idjoints = list(IDmoments_nat.keys())
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
        # now want the ID moments for comparison
        if joint in IDmoments_nat:
            # ax2[row, col].plot(xnat, IDmoments_nat[joint], label='ID natural', color='black', linestyle='--', alpha=0.8)
            # ax2[row, col].plot(xexo, IDmoments_exo[joint], label='ID exotendon', color='grey', linestyle='--', alpha=0.8)
            ax2[row, col].plot(xnat, np.mean(IDmoments_nat[joint],1), label='ID natural', color='black', linestyle='--', alpha=0.8)
            ax2[row, col].plot(xexo, np.mean(IDmoments_exo[joint],1), label='ID exotendon', color='grey', linestyle='--', alpha=0.8)
        elif 'pelvis_tx' in joint:
            ax2[row, col].plot(xnat, np.mean(IDmoments_nat['pelvis_tx_force'],1), label='ID natural', color='black', linestyle='--', alpha=0.8)
            ax2[row, col].plot(xexo, np.mean(IDmoments_exo['pelvis_tx_force'],1), label='ID exotendon', color='grey', linestyle='--', alpha=0.8)
        elif 'pelvis_ty' in joint:
            ax2[row, col].plot(xnat, np.mean(IDmoments_nat['pelvis_ty_force'],1), label='ID natural', color='black', linestyle='--', alpha=0.8)
            ax2[row, col].plot(xexo, np.mean(IDmoments_exo['pelvis_ty_force'],1), label='ID exotendon', color='grey', linestyle='--', alpha=0.8)
        elif 'pelvis_tz' in joint:
            ax2[row, col].plot(xnat, np.mean(IDmoments_nat['pelvis_tz_force'],1), label='ID natural', color='black', linestyle='--', alpha=0.8)
            ax2[row, col].plot(xexo, np.mean(IDmoments_exo['pelvis_tz_force'],1), label='ID exotendon', color='grey', linestyle='--', alpha=0.8)
        # formatting
        ax2[row, col].set_xlabel('% Gait cycle', fontsize=8)
        ax2[row, col].set_ylabel('Moment (Nm/kg)', fontsize=8)
        ax2[row, col].set_title(joint, fontsize=8)
    # add a legend to the last subplot
    handles, labels = ax2[0, 0].get_legend_handles_labels()
    ax2[-1,-1].axis('off')  # Hide the last subplot
    ax2[-1,-1].legend(handles, labels, loc='center', fontsize=8)
    # fig2.legend(handles, labels, loc='upper right')
    # ax2[0, 0].legend(loc='upper right', fontsize=8)
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
    plt.savefig(analyzedir + '\\hip_passiveforces_' + whichleg + '.png')


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
    plt.savefig(analyzedir + '\\hip_activeforces_' + whichleg + '.png')


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
    plt.savefig(analyzedir + '\\hip_totalforces_' + whichleg + '.png')

    plt.show()


    ################################################################################
    # now the knee contact forces 
    ################################################################################
    print('\n\nNow moving onto the hip contact forces\n')

    # tryin to implement everything below into a plotting function... that way we can just call it for each of the components... ideally.
    datay = {
        'ninterseg_combine': ninterseg_combine,
        'einterseg_combine': einterseg_combine,
        'nflexors_combine': nflexors_combine,
        'eflexors_combine': eflexors_combine,
        'nadductors_combine': nadductors_combine,
        'eadductors_combine': eadductors_combine,
        'nabductors_combine': nabductors_combine,
        'eabductors_combine': eabductors_combine,
        'nextensors_combine': nextensors_combine,
        'eextensors_combine': eextensors_combine,
        'nhams_combine': nhams_combine,
        'ehams_combine': ehams_combine,
        'nall_combine': nall_combine,
        'eall_combine': eall_combine,
    }
    plotHipContactForce('Vertical', analyzedir, welksubjects, ncolor, ecolor, n_timespercent101, e_timespercent101, datay)
    
    datax = {
        'ninterseg_combine': ninterseg_combine_x,
        'einterseg_combine': einterseg_combine_x,
        'nflexors_combine': nflexors_combine_x,
        'eflexors_combine': eflexors_combine_x,
        'nadductors_combine': nadductors_combine_x,
        'eadductors_combine': eadductors_combine_x,
        'nabductors_combine': nabductors_combine_x,
        'eabductors_combine': eabductors_combine_x,
        'nextensors_combine': nextensors_combine_x,
        'eextensors_combine': eextensors_combine_x,
        'nhams_combine': nhams_combine_x,
        'ehams_combine': ehams_combine_x,
        'nall_combine': nall_combine_x,
        'eall_combine': eall_combine_x,
    }
    # plotHipContactForce('Anterior-Posterior', analyzedir, welksubjects, ncolor, ecolor, n_timespercent101, e_timespercent101, datax)

    dataz = {
        'ninterseg_combine': ninterseg_combine_z,
        'einterseg_combine': einterseg_combine_z,
        'nflexors_combine': nflexors_combine_z,
        'eflexors_combine': eflexors_combine_z,
        'nadductors_combine': nadductors_combine_z,
        'eadductors_combine': eadductors_combine_z,
        'nabductors_combine': nabductors_combine_z,
        'eabductors_combine': eabductors_combine_z,
        'nextensors_combine': nextensors_combine_z,
        'eextensors_combine': eextensors_combine_z,
        'nhams_combine': nhams_combine_z,
        'ehams_combine': ehams_combine_z,
        'nall_combine': nall_combine_z,
        'eall_combine': eall_combine_z,
    }
    # plotHipContactForce('Medial-Lateral', analyzedir, welksubjects, ncolor, ecolor, n_timespercent101, e_timespercent101, dataz)

    # here is where I am going to try and combine the x and z into a shear force single value. 
    datashear = {}
    for each in datax.keys(): 
        datashear[each] = np.sqrt((datax[each]**2) + (dataz[each]**2))
    plotHipContactForce('Shear', analyzedir, welksubjects, ncolor, ecolor, n_timespercent101, e_timespercent101, datashear)

    # data for the resultant
    dataresultant = {}
    for each in datax.keys():
        dataresultant[each] = np.sqrt((datax[each]**2) + (datay[each]**2) + (dataz[each]**2))
    # plotHipContactForce('Resultant', analyzedir, welksubjects, ncolor, ecolor, n_timespercent101, e_timespercent101, dataresultant)

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
    plt.savefig(os.path.join(analyzedir, 'muscleactivations_hip_' + whichleg + '.png'))


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
    plt.savefig(analyzedir + '\\jointmoments_hip_' + whichleg + '.png')


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
    plt.savefig(analyzedir + '\\passiveforces_hip_' + whichleg + '.png')


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
    plt.savefig(analyzedir + '\\activeforces_hip_' + whichleg + '.png')


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
    plt.savefig(analyzedir + '\\totalforces_hip_' + whichleg + '.png')

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
    for i, curve in enumerate(nextensors_combine):
        ax9[1].plot(n_timespercent101, curve, label='natural'+str(i), color=ncolors[i] if len(welksubjects) == 1 else ncolor)
    for i, curve in enumerate(eextensors_combine):
        ax9[1].plot(e_timespercent101, curve, label='exotendon'+str(i), color=ecolors[i] if len(welksubjects) == 1 else ecolor)
    # ax9[1].plot(n_timespercent101, ntfl_combine, label='natural', color=ncolor)
    # ax9[1].plot(e_timespercent101, etfl_combine, label='exotendon', color=ecolor)
    ax9[1].set_xlabel('% Gait cycle')
    # ax9[1].set_ylabel('Force (BW)')
    # ax9[1].legend()
    ax9[1].set_title('tfl')

    # gastroc forces
    for i, curve in enumerate(nadductors_combine):
        ax9[2].plot(n_timespercent101, curve, label='natural'+str(i), color=ncolors[i] if len(welksubjects) == 1 else ncolor)
    for i, curve in enumerate(eadductors_combine):
        ax9[2].plot(e_timespercent101, curve, label='exotendon'+str(i), color=ecolors[i] if len(welksubjects) == 1 else ecolor)
    # ax9[2].plot(n_timespercent101, ngas_combine, label='natural', color=ncolor)
    # ax9[2].plot(e_timespercent101, egas_combine, label='exotendon', color=ecolor)
    ax9[2].set_xlabel('% Gait cycle')
    # ax9[2].set_ylabel('Force (BW)')
    # ax9[2].legend()
    ax9[2].set_title('gastroc')
    
    # hamstring forces
    for i, curve in enumerate(nabductors_combine):
        ax9[3].plot(n_timespercent101, curve, label='natural'+str(i), color=ncolors[i] if len(welksubjects) == 1 else ncolor)
    for i, curve in enumerate(eabductors_combine):
        ax9[3].plot(e_timespercent101, curve, label='exotendon'+str(i), color=ecolors[i] if len(welksubjects) == 1 else ecolor)
    # ax9[3].plot(n_timespercent101, nhams_combine, label='natural', color=ncolor)
    # ax9[3].plot(e_timespercent101, ehams_combine, label='exotendon', color=ecolor)
    ax9[3].set_xlabel('% Gait cycle')
    # ax9[3].set_ylabel('Force (BW)')
    # ax9[3].legend()
    ax9[3].set_title('hamstrings')

    # quads forces
    for i, curve in enumerate(nhams_combine):
        ax9[4].plot(n_timespercent101, curve, label='natural'+str(i), color=ncolors[i] if len(welksubjects) == 1 else ncolor)
    for i, curve in enumerate(ehams_combine):
        ax9[4].plot(e_timespercent101, curve, label='exotendon'+str(i), color=ecolors[i] if len(welksubjects) == 1 else ecolor)
    # ax9[4].plot(n_timespercent101, nquads_combine, label='natural', color=ncolor)
    # ax9[4].plot(e_timespercent101, equads_combine, label='exotendon', color=ecolor)
    ax9[4].set_xlabel('% Gait cycle')
    # ax9[4].set_ylabel('Force (BW)')
    # ax9[4].legend()
    ax9[4].set_title('quadriceps')

    # flexors forces
    for i, curve in enumerate(nflexors_combine):
        ax9[4].plot(n_timespercent101, curve, label='natural'+str(i), color=ncolors[i] if len(welksubjects) == 1 else ncolor)
    for i, curve in enumerate(eflexors_combine):
        ax9[4].plot(e_timespercent101, curve, label='exotendon'+str(i), color=ecolors[i] if len(welksubjects) == 1 else ecolor)
    # ax9[4].plot(n_timespercent101, nquads_combine, label='natural', color=ncolor)
    # ax9[4].plot(e_timespercent101, equads_combine, label='exotendon', color=ecolor)
    ax9[4].set_xlabel('% Gait cycle')
    # ax9[4].set_ylabel('Force (BW)')
    # ax9[4].legend()
    ax9[4].set_title('quadriceps')

    # reserve forces
    for i, curve in enumerate(nreserve_combine):
        ax9[5].plot(n_timespercent101, curve, label='natural'+str(i), color=ncolors[i] if len(welksubjects) == 1 else ncolor)
    for i, curve in enumerate(ereserve_combine):
        ax9[5].plot(e_timespercent101, curve, label='exotendon'+str(i), color=ecolors[i] if len(welksubjects) == 1 else ecolor)
    # ax9[5].plot(n_timespercent101, nreserve_combine, label='natural', color=ncolor)
    # ax9[5].plot(e_timespercent101, ereserve_combine, label='exotendon', color=ecolor)
    ax9[5].set_xlabel('% Gait cycle')
    # ax9[5].set_ylabel('Force (BW)')
    # ax9[5].legend()
    ax9[5].set_title('reserves')
    
    # # added all forces
    # ax9[5].plot(n_timespercent101, nquads_combine+ nhams_combine+ ngas_combine+ ntfl_combine+ ninterseg_combine+ nreserve_combine, label='natural', color=ncolor)
    # ax9[5].plot(e_timespercent101, equads_combine+ ehams_combine+ egas_combine+ etfl_combine+ einterseg_combine+ ereserve_combine, label='exotendon', color=ecolor)
    # ax9[5].set_xlabel('% Gait cycle')
    # # ax9[6].set_ylabel('Force (BW)')
    # # ax9[6].legend()
    # ax9[5].set_title('Total vertical contact')

    # all forces from whole analysis
    for i, curve in enumerate(nall_combine):
        ax9[6].plot(n_timespercent101, curve, label='natural'+str(i), color=ncolors[i] if len(welksubjects) == 1 else ncolor)
    for i, curve in enumerate(eall_combine):
        ax9[6].plot(e_timespercent101, curve, label='exotendon'+str(i), color=ecolors[i] if len(welksubjects) == 1 else ecolor)
    # ax9[6].plot(n_timespercent101, nall_combine, label='natural', color=ncolor)
    # ax9[6].plot(e_timespercent101, eall_combine, label='exotendon', color=ecolor)
    ax9[6].set_xlabel('% Gait cycle')
    # ax9[6].set_ylabel('Force (BW)')
    ax9[6].set_title('Total vertical contact')
    # ax9[6].legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10)
    # Hide the last subplot and use it to display the legend   
    ax9[7].axis('off')
    handles, labels = ax9[5].get_legend_handles_labels()
    if len(welksubjects) == 1:
        ax9[7].legend(handles, labels, loc='center', fontsize=14)

    fig9.tight_layout()
    plt.savefig(analyzedir + '\\contact1_hip_' + whichleg + '.png')


    ###########################################################################
    # figure: segmenting all the muscles between exo and nat
    ## really nice figure for seeing what is going on, but likely not going to
    ## be in the paper...
    fig10, ax10 = plt.subplots(1,8, figsize=(18,3), dpi=500)
    fontz = 16
    font_properties = {'fontsize': 16, 'fontfamily': 'serif', 'fontname': 'Times New Roman'}
    # tick_font_properties = {'fontfamily': 'serif', 'fontname': 'Times New Roman'}
    
    # all forces average
    ax10[0].fill_between(n_timespercent101, np.mean(nall_combine,0) - np.std(nall_combine,0), np.mean(nall_combine,0) + np.std(nall_combine,0), color=ncolor, alpha=0.2)
    ax10[0].fill_between(e_timespercent101, np.mean(eall_combine,0) - np.std(eall_combine,0), np.mean(eall_combine,0) + np.std(eall_combine,0), color=ecolor, alpha=0.2)
    ax10[0].plot(n_timespercent101, np.mean(nall_combine, 0), label='Natural (Mean \u00B1 Std.)', color=ncolor)
    ax10[0].plot(e_timespercent101, np.mean(eall_combine, 0), label='Exotendon (Mean \u00B1 Std.)', color=ecolor)
    ax10[0].set_xlabel('% Gait cycle', **font_properties)
    ax10[0].set_ylabel('Force (BW)', **font_properties)
    # ax10[0].legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10)
    ax10[0].set_title('Total\nVertical Hip Contact', **font_properties)
    ax10[0].tick_params(axis='both', which='both', labelsize=fontz)
    
    # extensors forces
    ax10[1].fill_between(n_timespercent101, np.mean(nextensors_combine,0) - np.std(nextensors_combine,0), np.mean(nextensors_combine,0) + np.std(nextensors_combine,0), color=ncolor, alpha=0.2)
    ax10[1].fill_between(e_timespercent101, np.mean(eextensors_combine,0) - np.std(eextensors_combine,0), np.mean(eextensors_combine,0) + np.std(eextensors_combine,0), color=ecolor, alpha=0.2)
    ax10[1].plot(n_timespercent101, np.mean(nextensors_combine, 0), label='natural', color=ncolor)
    ax10[1].plot(e_timespercent101, np.mean(eextensors_combine, 0), label='exotendon', color=ecolor)
    ax10[1].set_xlabel('% Gait cycle', **font_properties)
    # ax10[1].set_ylabel('Force (BW)')
    # ax10[1].legend()
    ax10[1].tick_params(axis='both', which='major', labelsize=fontz)
    ax10[1].set_title('Contribution of\nHip Extensors', **font_properties)

    # adductor forces
    ax10[2].fill_between(n_timespercent101, np.mean(nadductors_combine,0) - np.std(nadductors_combine,0), np.mean(nadductors_combine,0) + np.std(nadductors_combine,0), color=ncolor, alpha=0.2)
    ax10[2].fill_between(e_timespercent101, np.mean(eadductors_combine,0) - np.std(eadductors_combine,0), np.mean(eadductors_combine,0) + np.std(eadductors_combine,0), color=ecolor, alpha=0.2)
    ax10[2].plot(n_timespercent101, np.mean(nadductors_combine, 0), label='natural', color=ncolor)
    ax10[2].plot(e_timespercent101, np.mean(eadductors_combine, 0), label='exotendon', color=ecolor)
    ax10[2].set_xlabel('% Gait cycle', **font_properties)
    ax10[2].tick_params(axis='both', which='major', labelsize=fontz)
    # ax10[2].set_ylabel('Force (BW)')
    # ax10[2].legend()
    ax10[2].set_title('Contribution of\nHip Adductors', **font_properties)
    
    # abductors forces
    ax10[3].fill_between(n_timespercent101, np.mean(nabductors_combine,0) - np.std(nabductors_combine,0), np.mean(nabductors_combine,0) + np.std(nabductors_combine,0), color=ncolor, alpha=0.2)
    ax10[3].fill_between(e_timespercent101, np.mean(eabductors_combine,0) - np.std(eabductors_combine,0), np.mean(eabductors_combine,0) + np.std(eabductors_combine,0), color=ecolor, alpha=0.2)
    ax10[3].plot(n_timespercent101, np.mean(nabductors_combine, 0), label='natural', color=ncolor)
    ax10[3].plot(e_timespercent101, np.mean(eabductors_combine, 0), label='exotendon', color=ecolor)
    ax10[3].set_xlabel('% Gait cycle', **font_properties)
    ax10[3].tick_params(axis='both', which='major', labelsize=fontz)
    # ax10[2].set_ylabel('Force (BW)')
    # ax10[2].legend()
    ax10[3].set_title('Contribution of\nHip Abductors', **font_properties)
    
    # hamstring forces
    ax10[4].fill_between(n_timespercent101, np.mean(nhams_combine,0) - np.std(nhams_combine,0), np.mean(nhams_combine,0) + np.std(nhams_combine,0), color=ncolor, alpha=0.2)
    ax10[4].fill_between(e_timespercent101, np.mean(ehams_combine,0) - np.std(ehams_combine,0), np.mean(ehams_combine,0) + np.std(ehams_combine,0), color=ecolor, alpha=0.2)
    ax10[4].plot(n_timespercent101, np.mean(nhams_combine, 0), label='natural', color=ncolor)
    ax10[4].plot(e_timespercent101, np.mean(ehams_combine, 0), label='exotendon', color=ecolor)
    ax10[4].set_xlabel('% Gait cycle', **font_properties)
    ax10[4].tick_params(axis='both', which='major', labelsize=fontz)
    # ax10[3].set_ylabel('Force (BW)')
    # ax10[3].legend()
    ax10[4].set_title('Contribution of\nHamstrings', **font_properties)

    # flexors forces
    ax10[5].fill_between(n_timespercent101, np.mean(nflexors_combine,0) - np.std(nflexors_combine,0), np.mean(nflexors_combine,0) + np.std(nflexors_combine,0), color=ncolor, alpha=0.2)
    ax10[5].fill_between(e_timespercent101, np.mean(eflexors_combine,0) - np.std(eflexors_combine,0), np.mean(eflexors_combine,0) + np.std(eflexors_combine,0), color=ecolor, alpha=0.2)
    ax10[5].plot(n_timespercent101, np.mean(nflexors_combine, 0), label='natural', color=ncolor)
    ax10[5].plot(e_timespercent101, np.mean(eflexors_combine, 0), label='exotendon', color=ecolor)
    ax10[5].set_xlabel('% Gait cycle', **font_properties)
    ax10[5].tick_params(axis='both', which='major', labelsize=fontz)
    # ax10[4].set_ylabel('Force (BW)')
    # ax10[4].legend()
    ax10[5].set_title('Contribution of\nHip Flexors', **font_properties)

    # intersegmental forces
    ax10[6].fill_between(n_timespercent101, np.mean(ninterseg_combine,0) - np.std(ninterseg_combine,0), np.mean(ninterseg_combine,0) + np.std(ninterseg_combine,0), color=ncolor, alpha=0.2)
    ax10[6].fill_between(e_timespercent101, np.mean(einterseg_combine,0) - np.std(einterseg_combine,0), np.mean(einterseg_combine,0) + np.std(einterseg_combine,0), color=ecolor, alpha=0.2)
    ax10[6].plot(n_timespercent101, np.mean(ninterseg_combine, 0), label='natural', color=ncolor)
    ax10[6].plot(e_timespercent101, np.mean(einterseg_combine, 0), label='exotendon', color=ecolor)
    ax10[6].set_xlabel('% Gait cycle', **font_properties)
    ax10[6].tick_params(axis='both', which='major', labelsize=fontz)
    # ax10[6].set_ylabel('Force (BW)')
    # ax10[6].legend()
    ax10[6].set_title('Contribution of\nIntersegmental Forces', **font_properties)
    
    # # added all forces
    # ax10[6].plot(n_timespercent101, np.mean(nextensors_combine,0) + np.mean(nhams_combine,0) + np.mean(nadductors_combine,0) + np.mean(nabductors_combine,0) + np.mean(nflexors_combine,0) + np.mean(ninterseg_combine,0) + np.mean(nreserve_combine,0), label='natural', color=ncolor)
    # ax10[6].plot(e_timespercent101, np.mean(eextensors_combine,0) + np.mean(ehams_combine,0) + np.mean(eadductors_combine,0) + np.mean(eabductors_combine,0) + np.mean(eflexors_combine,0) + np.mean(einterseg_combine,0) + np.mean(ereserve_combine,0), label='exotendon', color=ecolor)
    # ax10[6].set_xlabel('% Gait cycle')
    # # ax10[6].set_ylabel('Force (BW)')
    # # ax10[6].legend()
    # ax10[6].set_title('Total vertical contact added')

    # # all forces from whole analysis
    # ax10[7].plot(n_timespercent101, np.mean(nall_combine,0), label='natural', color=ncolor)
    # ax10[7].plot(e_timespercent101, np.mean(eall_combine,0), label='exotendon', color=ecolor)
    # ax10[7].set_xlabel('% Gait cycle')
    # # ax10[7].set_ylabel('Force (BW)')
    # ax10[7].legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10)
    # ax10[7].set_title('Total vertical contact')

    # Hide the last subplot and use it to display the legend
    ax10[7].axis('off')
    handles, labels = ax10[0].get_legend_handles_labels()
    # ax10[7].legend(handles, labels, loc='center', fontsize=14)
    fig10.tight_layout(pad=2.0, w_pad=0.5, h_pad=1.0)
    plt.savefig(analyzedir + '\\contact2_hip_' + whichleg + '.png')
    
    ###########################################################################
    ### figure: differences between conditions for each and all
    #### this is a simplified look at the same as above. We can see nice stuff, 
    #### but likely not going to be in the paper. 
        
    fig11, ax11 = plt.subplots(1,8, figsize=(18,3), dpi=300)
    # intersegmental forces average
    ax11[0].plot(n_timespercent101, np.mean(einterseg_combine, 0) - np.mean(ninterseg_combine, 0))
    ax11[0].set_xlabel('% Gait cycle')
    ax11[0].set_ylabel('Force (BW)')
    ax11[0].set_title('intersegmental')
    # ax11[0].legend()
    
    # extensor forces
    ax11[1].plot(n_timespercent101, np.mean(eextensors_combine, 0) - np.mean(nextensors_combine, 0))
    ax11[1].set_xlabel('% Gait cycle')
    # ax11[1].set_ylabel('Force (BW)')
    # ax11[1].legend()
    ax11[1].set_title('hip extensors')

    # adductor forces
    ax11[2].plot(n_timespercent101, np.mean(eadductors_combine, 0) - np.mean(nadductors_combine, 0))
    ax11[2].set_xlabel('% Gait cycle')
    # ax11[2].set_ylabel('Force (BW)')
    # ax11[2].legend()
    ax11[2].set_title('hip adductors')
    
    # abductor forces
    ax11[3].plot(n_timespercent101, np.mean(eabductors_combine, 0) - np.mean(nabductors_combine, 0))
    ax11[3].set_xlabel('% Gait cycle')
    # ax11[2].set_ylabel('Force (BW)')
    # ax11[2].legend()
    ax11[3].set_title('hip abductors')
    
    # hamstring forces
    ax11[4].plot(n_timespercent101, np.mean(ehams_combine, 0) - np.mean(nhams_combine, 0))
    ax11[4].set_xlabel('% Gait cycle')
    # ax11[3].set_ylabel('Force (BW)')
    # ax11[3].legend()
    ax11[4].set_title('hamstrings')

    # flexor forces
    ax11[5].plot(n_timespercent101, np.mean(eflexors_combine, 0) - np.mean(nflexors_combine, 0))
    ax11[5].set_xlabel('% Gait cycle')
    # ax11[4].set_ylabel('Force (BW)')
    # ax11[4].legend()
    ax11[5].set_title('hip flexors')

    # # added all forces
    # ax11[6].plot(n_timespercent101, (np.mean(eextensors_combine,0) + np.mean(eadductors_combine,0) + np.mean(eabductors_combine,0) + np.mean(ehams_combine,0) + np.mean(eflexors_combine,0) + np.mean(einterseg_combine,0)) - (np.mean(nextensors_combine,0) + np.mean(nhams_combine,0) + np.mean(nadductors_combine,0) + np.mean(nabductors_combine,0) + np.mean(nflexors_combine,0) + np.mean(ninterseg_combine,0)))
    # ax11[6].set_xlabel('% Gait cycle')
    # # ax11[5].set_ylabel('Force (BW)')
    # # ax11[5].legend()
    # ax11[6].set_title('Total vertical contact added')

    # all forces from whole analysis
    ax11[6].plot(n_timespercent101, np.mean(eall_combine,0) - np.mean(nall_combine,0))
    ax11[6].set_xlabel('% Gait cycle')
    # ax11[6].set_ylabel('Force (BW)')
    ax11[6].legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10)
    ax11[6].set_title('Total vertical contact')

    ax11[7].axis('off')
    handles, labels = ax11[6].get_legend_handles_labels()
    ax11[6].legend(handles, labels, loc='center', fontsize=10)
    fig11.tight_layout()
    plt.savefig(analyzedir + '\\contact3_hip_' + whichleg + '.png')
    
    
    ###########################################################################
    ### figure: changes in force segmented together on plot
    #### Possible paper figure for R3. 
    colors4 = ['#648FFF','#785EF0','#DC267F','#FE6100', '#FFB000']
    fig12, ax12 = plt.subplots(1, 2, figsize=(12,4.55), dpi=500)
    ax12[0].axhline(0, color='black', linewidth=1)
    # intersegmental forces average
    ax12[0].plot(n_timespercent101, np.mean(einterseg_combine, 0) - np.mean(ninterseg_combine, 0), label='Intersegmental', color='black', linewidth=3)
    # extensor forces
    ax12[0].plot(n_timespercent101, np.mean(eextensors_combine, 0) - np.mean(nextensors_combine, 0), label='Hip extensors', color=colors4[0], linewidth=3)
    # adductor forces
    ax12[0].plot(n_timespercent101, np.mean(eadductors_combine, 0) - np.mean(nadductors_combine, 0), label='Hip adductors', color=colors4[1], linewidth=3)
    # abductor forces
    ax12[0].plot(n_timespercent101, np.mean(eabductors_combine, 0) - np.mean(nabductors_combine, 0), label='Hip abductors', color=colors4[2], linewidth=3)
    # hamstring forces
    ax12[0].plot(n_timespercent101, np.mean(ehams_combine, 0) - np.mean(nhams_combine, 0), label='Hamstrings', color=colors4[3], linewidth=3)
    # flexor forces
    ax12[0].plot(n_timespercent101, np.mean(eflexors_combine, 0) - np.mean(nflexors_combine, 0), label='Hip flexors', color=colors4[4], linewidth=3)
    # reserve forces
    # ax12[0].plot(n_timespercent101, np.mean(ereserve_combine, 0) - np.mean(nreserve_combine, 0), label='reserves')
    # added all forces
    # ax12[0].plot(n_timespercent101, (np.mean(eextensors_combine,0) + np.mean(ehams_combine,0) + np.mean(eadductors_combine,0) + np.mean(eabductors_combine,0) + np.mean(eflexors_combine,0) + np.mean(einterseg_combine,0) + np.mean(ereserve_combine,0)) - (np.mean(nextensors_combine,0) + np.mean(nhams_combine,0) + np.mean(nadductors_combine,0) + np.mean(nabductors_combine,0) + np.mean(nflexors_combine,0) + np.mean(ninterseg_combine,0) + np.mean(nreserve_combine,0)), label='Total')
    # all forces from whole analysis
    ax12[0].plot(n_timespercent101, np.mean(eall_combine,0) - np.mean(nall_combine,0), label='Total vertical contact', linestyle='dashed', color='black', linewidth=3)
    ax12[0].set_xlabel('% Gait cycle', fontsize=16)
    ax12[0].set_ylabel('Vertical hip contact difference (BW)', fontsize=16)
    ax12[0].set_title('Exotendon change in contact force', fontsize=16)
    ax12[0].tick_params(axis='both', which='major', labelsize=16)
    # hide the other subplot and use it for the legend
    ax12[1].axis('off')
    handles, labels = ax12[0].get_legend_handles_labels()
    ax12[1].legend(handles, labels, loc='center', fontsize=16)
    fig12.tight_layout()
    fig12.savefig(analyzedir + '\\contact4_hip_' + whichleg + '.png')

    ###########################################################################
    # TODO: figure out why the difference in total and all added together. 
    # okay so not in how I am adding/averaging. has to be something in how the analysis is done between them.... am I missing something??
        
    fig13, ax13 = plt.subplots(1,3, figsize=(14,5))#, dpi=300)
    # intersegmental forces average - natural
    ax13[0].plot(n_timespercent101, np.mean(ninterseg_combine, 0), label='natural_interseg')
    ax13[0].plot(n_timespercent101, np.mean(ninterseg_combine + nextensors_combine, 0), label='natural_interseg + extensors')
    ax13[0].plot(n_timespercent101, np.mean(ninterseg_combine+nextensors_combine+nhams_combine, 0), label='natural_interseg+extensors+hams')
    ax13[0].plot(n_timespercent101, np.mean(ninterseg_combine+nextensors_combine+nhams_combine+nadductors_combine, 0), label='natural_interseg+extensors+hams+adductors')
    ax13[0].plot(n_timespercent101, np.mean(ninterseg_combine+nextensors_combine+nhams_combine+nadductors_combine+nabductors_combine, 0), label='natural_interseg+extensors+hams+adductors+abductors')
    ax13[0].plot(n_timespercent101, np.mean(ninterseg_combine+nextensors_combine+nhams_combine+nadductors_combine+nabductors_combine+nflexors_combine, 0), label='natural_interseg+extensors+hams+adductors+abductors+flexors')
    ax13[0].plot(n_timespercent101, np.mean(ninterseg_combine+nextensors_combine+nhams_combine+nadductors_combine+nabductors_combine+nflexors_combine+nreserve_combine, 0), label='natural_interseg+extensors+hams+adductors+abductors+flexors+reserve')
    ax13[0].plot(n_timespercent101, np.mean(nall_combine, 0), label='nat_all', linestyle='dotted')
    ax13[0].set_xlabel('% Gait cycle')
    ax13[0].set_ylabel('Force (BW)')
    ax13[0].set_title('natural')
    # ax13[0].legend()
    # intersegmental forces average - exotendon
    ax13[1].plot(e_timespercent101, np.mean(einterseg_combine, 0), label='exo_interseg')
    ax13[1].plot(e_timespercent101, np.mean(einterseg_combine+eextensors_combine, 0), label='exo_interseg + extensors')
    ax13[1].plot(e_timespercent101, np.mean(einterseg_combine+eextensors_combine+ehams_combine, 0), label='exo_interseg+extensors+hams')
    ax13[1].plot(e_timespercent101, np.mean(einterseg_combine+eextensors_combine+ehams_combine+eadductors_combine, 0), label='exo_interseg+extensors+hams+adductors')
    ax13[1].plot(e_timespercent101, np.mean(einterseg_combine+eextensors_combine+ehams_combine+eadductors_combine+eabductors_combine, 0), label='exo_interseg+extensors+hams+adductors+abductors')
    ax13[1].plot(e_timespercent101, np.mean(einterseg_combine+eextensors_combine+ehams_combine+eadductors_combine+eabductors_combine+eflexors_combine, 0), label='exo_interseg+extensors+hams+adductors+abductors+flexors')
    ax13[1].plot(e_timespercent101, np.mean(einterseg_combine+eextensors_combine+ehams_combine+eadductors_combine+eabductors_combine+eflexors_combine+ereserve_combine, 0), label='exo_interseg+extensors+hams+adductors+abductors+flexors+reserve')
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
    plt.savefig(analyzedir + '\\contact5_hip_' + whichleg + '.png')
    
    
    # ###########################################################################
    # # figures like ^ but individuals
    # fig131, ax131 = plt.subplots(1,2, figsize=(10,8), dpi=300)
    # # intersegmental forces average - natural welk002
    # ax131[0].plot(n_timespercent101, np.mean(ninterseg_combine[0:4,:], 0), label='natural_interseg')
    # ax131[0].plot(n_timespercent101, np.mean(ninterseg_combine[0:4,:] + nextensors_combine[0:4,:], 0), label='natural_interseg + extensors')
    # ax131[0].plot(n_timespercent101, np.mean(ninterseg_combine[0:4,:]+nextensors_combine[0:4,:]+nhams_combine[0:4,:], 0), label='natural_interseg+extensors+hams')
    # ax131[0].plot(n_timespercent101, np.mean(ninterseg_combine[0:4,:]+nextensors_combine[0:4,:]+nhams_combine[0:4,:]+nadductors_combine[0:4,:], 0), label='natural_interseg+extensors+hams+adductors')
    # ax131[0].plot(n_timespercent101, np.mean(ninterseg_combine[0:4,:]+nextensors_combine[0:4,:]+nhams_combine[0:4,:]+nadductors_combine[0:4,:]+nabductors_combine[0:4,:], 0), label='natural_interseg+extensors+hams+adductors+abductors')
    # ax131[0].plot(n_timespercent101, np.mean(ninterseg_combine[0:4,:]+nextensors_combine[0:4,:]+nhams_combine[0:4,:]+nadductors_combine[0:4,:]+nabductors_combine[0:4,:]+nflexors_combine[0:4,:], 0), label='natural_interseg+extensors+hams+adductors+abductors+flexors')
    # ax131[0].plot(n_timespercent101, np.mean(ninterseg_combine[0:4,:]+nextensors_combine[0:4,:]+nhams_combine[0:4,:]+nadductors_combine[0:4,:]+nabductors_combine[0:4,:]+nflexors_combine[0:4,:]+nreserve_combine[0:4,:], 0), label='natural_interseg+extensors+hams+adductors+abductors+flexors+reserve')
    # ax131[0].plot(n_timespercent101, np.mean(nall_combine[0:4,:], 0), label='nat_all', linestyle='dotted')
    # ax131[0].set_xlabel('% Gait cycle')
    # ax131[0].set_ylabel('Force (BW)')
    # ax131[0].set_title('welk002 natural')
    # ax131[0].legend()
    # # intersegmental forces average - exotendon
    # ax131[1].plot(e_timespercent101, np.mean(einterseg_combine[0:4,:], 0), label='exo_interseg')
    # ax131[1].plot(e_timespercent101, np.mean(einterseg_combine[0:4,:]+eextensors_combine[0:4,:], 0), label='exo_interseg + extensors')
    # ax131[1].plot(e_timespercent101, np.mean(einterseg_combine[0:4,:]+eextensors_combine[0:4,:]+ehams_combine[0:4,:], 0), label='exo_interseg+extensors+hams')
    # ax131[1].plot(e_timespercent101, np.mean(einterseg_combine[0:4,:]+eextensors_combine[0:4,:]+ehams_combine[0:4,:]+eadductors_combine[0:4,:], 0), label='exo_interseg+extensors+hams+adductors')
    # ax131[1].plot(e_timespercent101, np.mean(einterseg_combine[0:4,:]+eextensors_combine[0:4,:]+ehams_combine[0:4,:]+eadductors_combine[0:4,:]+eabductors_combine[0:4,:], 0), label='exo_interseg+extensors+hams+adductors+abductors')
    # ax131[1].plot(e_timespercent101, np.mean(einterseg_combine[0:4,:]+eextensors_combine[0:4,:]+ehams_combine[0:4,:]+eadductors_combine[0:4,:]+eabductors_combine[0:4,:]+eflexors_combine[0:4,:], 0), label='exo_interseg+extensors+hams+adductors+abductors+flexors')
    # ax131[1].plot(e_timespercent101, np.mean(einterseg_combine[0:4,:]+eextensors_combine[0:4,:]+ehams_combine[0:4,:]+eadductors_combine[0:4,:]+eabductors_combine[0:4,:]+eflexors_combine[0:4,:]+ereserve_combine[0:4,:], 0), label='exo_interseg+extensors+hams+adductors+abductors+flexors+reserve')
    # ax131[1].plot(e_timespercent101, np.mean(eall_combine[0:4,:], 0), label='exo_all', linestyle='dotted')
    # ax131[1].set_xlabel('% Gait cycle')
    # ax131[1].set_ylabel('Force (BW)')
    # ax131[1].set_title('welk002 exotendon')
    # ax131[1].legend()
    # fig131.tight_layout()
    # # welk003
    # # figures like ^ but individuals
    # fig132, ax132 = plt.subplots(1,2, figsize=(10,8), dpi=300)
    # # intersegmental forces average - natural welk002
    # ax132[0].plot(n_timespercent101, np.mean(ninterseg_combine[4:8,:], 0), label='natural_interseg')
    # ax132[0].plot(n_timespercent101, np.mean(ninterseg_combine[4:8,:] + nextensors_combine[4:8,:], 0), label='natural_interseg + extensors')
    # ax132[0].plot(n_timespercent101, np.mean(ninterseg_combine[4:8,:]+nextensors_combine[4:8,:]+nhams_combine[4:8,:], 0), label='natural_interseg+extensors+hams')
    # ax132[0].plot(n_timespercent101, np.mean(ninterseg_combine[4:8,:]+nextensors_combine[4:8,:]+nhams_combine[4:8,:]+nadductors_combine[4:8,:], 0), label='natural_interseg+extensors+hams+adductors')
    # ax132[0].plot(n_timespercent101, np.mean(ninterseg_combine[4:8,:]+nextensors_combine[4:8,:]+nhams_combine[4:8,:]+nadductors_combine[4:8,:]+nabductors_combine[4:8,:], 0), label='natural_interseg+extensors+hams+adductors+abductors')
    # ax132[0].plot(n_timespercent101, np.mean(ninterseg_combine[4:8,:]+nextensors_combine[4:8,:]+nhams_combine[4:8,:]+nadductors_combine[4:8,:]+nabductors_combine[4:8,:]+nflexors_combine[4:8,:], 0), label='natural_interseg+extensors+hams+adductors+abductors+flexors')
    # ax132[0].plot(n_timespercent101, np.mean(ninterseg_combine[4:8,:]+nextensors_combine[4:8,:]+nhams_combine[4:8,:]+nadductors_combine[4:8,:]+nabductors_combine[4:8,:]+nflexors_combine[4:8,:]+nreserve_combine[4:8,:], 0), label='natural_interseg+extensors+hams+adductors+abductors+flexors+reserve')
    # ax132[0].plot(n_timespercent101, np.mean(nall_combine[4:8,:], 0), label='nat_all', linestyle='dotted')
    # ax132[0].set_xlabel('% Gait cycle')
    # ax132[0].set_ylabel('Force (BW)')
    # ax132[0].set_title('welk003 natural')
    # ax132[0].legend()
    # # intersegmental forces average - exotendon
    # ax132[1].plot(e_timespercent101, np.mean(einterseg_combine[4:8,:], 0), label='exo_interseg')
    # ax132[1].plot(e_timespercent101, np.mean(einterseg_combine[4:8,:]+eextensors_combine[4:8,:], 0), label='exo_interseg + extensors')
    # ax132[1].plot(e_timespercent101, np.mean(einterseg_combine[4:8,:]+eextensors_combine[4:8,:]+ehams_combine[4:8,:], 0), label='exo_interseg+extensors+hams')
    # ax132[1].plot(e_timespercent101, np.mean(einterseg_combine[4:8,:]+eextensors_combine[4:8,:]+ehams_combine[4:8,:]+eadductors_combine[4:8,:], 0), label='exo_interseg+extensors+hams+adductors')
    # ax132[1].plot(e_timespercent101, np.mean(einterseg_combine[4:8,:]+eextensors_combine[4:8,:]+ehams_combine[4:8,:]+eadductors_combine[4:8,:]+eabductors_combine[4:8,:], 0), label='exo_interseg+extensors+hams+adductors+abductors')
    # ax132[1].plot(e_timespercent101, np.mean(einterseg_combine[4:8,:]+eextensors_combine[4:8,:]+ehams_combine[4:8,:]+eadductors_combine[4:8,:]+eabductors_combine[4:8,:]+eflexors_combine[4:8,:], 0), label='exo_interseg+extensors+hams+adductors+abductors+flexors')
    # ax132[1].plot(e_timespercent101, np.mean(einterseg_combine[4:8,:]+eextensors_combine[4:8,:]+ehams_combine[4:8,:]+eadductors_combine[4:8,:]+eabductors_combine[4:8,:]+eflexors_combine[4:8,:]+ereserve_combine[4:8,:], 0), label='exo_interseg+extensors+hams+adductors+abductors+flexors+reserve')
    # ax132[1].plot(e_timespercent101, np.mean(eall_combine[4:8,:], 0), label='exo_all', linestyle='dotted')
    # ax132[1].set_xlabel('% Gait cycle')
    # ax132[1].set_ylabel('Force (BW)')
    # ax132[1].set_title('welk003 exotendon')
    # ax132[1].legend()
    # fig132.tight_layout()
    # # welk005
    
    # # figures like ^ but individuals
    # fig133, ax133 = plt.subplots(1,2, figsize=(10,8), dpi=300)
    # # intersegmental forces average - natural welk002
    # ax133[0].plot(n_timespercent101, np.mean(ninterseg_combine[8:12,:], 0), label='natural_interseg')
    # ax133[0].plot(n_timespercent101, np.mean(ninterseg_combine[8:12,:] + nextensors_combine[8:12,:], 0), label='natural_interseg + extensors')
    # ax133[0].plot(n_timespercent101, np.mean(ninterseg_combine[8:12,:]+nextensors_combine[8:12,:]+nhams_combine[8:12,:], 0), label='natural_interseg+extensors+hams')
    # ax133[0].plot(n_timespercent101, np.mean(ninterseg_combine[8:12,:]+nextensors_combine[8:12,:]+nhams_combine[8:12,:]+nadductors_combine[8:12,:], 0), label='natural_interseg+extensors+hams+adductors')
    # ax133[0].plot(n_timespercent101, np.mean(ninterseg_combine[8:12,:]+nextensors_combine[8:12,:]+nhams_combine[8:12,:]+nadductors_combine[8:12,:]+nabductors_combine[8:12,:], 0), label='natural_interseg+extensors+hams+adductors+abductors')
    # ax133[0].plot(n_timespercent101, np.mean(ninterseg_combine[8:12,:]+nextensors_combine[8:12,:]+nhams_combine[8:12,:]+nadductors_combine[8:12,:]+nabductors_combine[8:12,:]+nflexors_combine[8:12,:], 0), label='natural_interseg+extensors+hams+adductors+abductors+flexors')
    # ax133[0].plot(n_timespercent101, np.mean(ninterseg_combine[8:12,:]+nextensors_combine[8:12,:]+nhams_combine[8:12,:]+nadductors_combine[8:12,:]+nabductors_combine[8:12,:]+nflexors_combine[8:12,:]+nreserve_combine[8:12,:], 0), label='natural_interseg+extensors+hams+adductors+abductors+flexors+reserve')
    # ax133[0].plot(n_timespercent101, np.mean(nall_combine[8:12,:], 0), label='nat_all', linestyle='dotted')
    # ax133[0].set_xlabel('% Gait cycle')
    # ax133[0].set_ylabel('Force (BW)')
    # ax133[0].set_title('welk005 natural')
    # ax133[0].legend()
    # # intersegmental forces average - exotendon
    # ax133[1].plot(e_timespercent101, np.mean(einterseg_combine[8:12,:], 0), label='exo_interseg')
    # ax133[1].plot(e_timespercent101, np.mean(einterseg_combine[8:12,:]+eextensors_combine[8:12,:], 0), label='exo_interseg + extensors')
    # ax133[1].plot(e_timespercent101, np.mean(einterseg_combine[8:12,:]+eextensors_combine[8:12,:]+ehams_combine[8:12,:], 0), label='exo_interseg+extensors+hams')
    # ax133[1].plot(e_timespercent101, np.mean(einterseg_combine[8:12,:]+eextensors_combine[8:12,:]+ehams_combine[8:12,:]+eadductors_combine[8:12,:], 0), label='exo_interseg+extensors+hams+adductors')
    # ax133[1].plot(e_timespercent101, np.mean(einterseg_combine[8:12,:]+eextensors_combine[8:12,:]+ehams_combine[8:12,:]+eadductors_combine[8:12,:]+eabductors_combine[8:12,:], 0), label='exo_interseg+extensors+hams+adductors+abductors')
    # ax133[1].plot(e_timespercent101, np.mean(einterseg_combine[8:12,:]+eextensors_combine[8:12,:]+ehams_combine[8:12,:]+eadductors_combine[8:12,:]+eabductors_combine[8:12,:]+eflexors_combine[8:12,:], 0), label='exo_interseg+extensors+hams+adductors+abductors+flexors')
    # ax133[1].plot(e_timespercent101, np.mean(einterseg_combine[8:12,:]+eextensors_combine[8:12,:]+ehams_combine[8:12,:]+eadductors_combine[8:12,:]+eabductors_combine[8:12,:]+eflexors_combine[8:12,:]+ereserve_combine[8:12,:], 0), label='exo_interseg+extensors+hams+adductors+abductors+flexors+reserve')
    # ax133[1].plot(e_timespercent101, np.mean(eall_combine[8:12,:], 0), label='exo_all', linestyle='dotted')
    # ax133[1].set_xlabel('% Gait cycle')
    # ax133[1].set_ylabel('Force (BW)')
    # ax133[1].set_title('welk005 exotendon')
    # ax133[1].legend()
    # fig133.tight_layout()
    # # welk008
    # # figures like ^ but individuals
    # fig134, ax134 = plt.subplots(1,2, figsize=(10,8), dpi=300)
    # # intersegmental forces average - natural welk002
    # ax134[0].plot(n_timespercent101, np.mean(ninterseg_combine[12:16,:], 0), label='natural_interseg')
    # ax134[0].plot(n_timespercent101, np.mean(ninterseg_combine[12:16,:] + nextensors_combine[12:16,:], 0), label='natural_interseg + extensors')
    # ax134[0].plot(n_timespercent101, np.mean(ninterseg_combine[12:16,:]+nextensors_combine[12:16,:]+nhams_combine[12:16,:], 0), label='natural_interseg+extensors+hams')
    # ax134[0].plot(n_timespercent101, np.mean(ninterseg_combine[12:16,:]+nextensors_combine[12:16,:]+nhams_combine[12:16,:]+nadductors_combine[12:16,:], 0), label='natural_interseg+extensors+hams+adductors')
    # ax134[0].plot(n_timespercent101, np.mean(ninterseg_combine[12:16,:]+nextensors_combine[12:16,:]+nhams_combine[12:16,:]+nadductors_combine[12:16,:]+nabductors_combine[12:16,:], 0), label='natural_interseg+extensors+hams+adductors+abductors')
    # ax134[0].plot(n_timespercent101, np.mean(ninterseg_combine[12:16,:]+nextensors_combine[12:16,:]+nhams_combine[12:16,:]+nadductors_combine[12:16,:]+nabductors_combine[12:16,:]+nflexors_combine[12:16,:], 0), label='natural_interseg+extensors+hams+adductors+abductors+flexors')
    # ax134[0].plot(n_timespercent101, np.mean(ninterseg_combine[12:16,:]+nextensors_combine[12:16,:]+nhams_combine[12:16,:]+nadductors_combine[12:16,:]+nabductors_combine[12:16,:]+nflexors_combine[12:16,:]+nreserve_combine[12:16,:], 0), label='natural_interseg+extensors+hams+adductors+abductors+flexors+reserve')
    # ax134[0].plot(n_timespercent101, np.mean(nall_combine[12:16,:], 0), label='nat_all', linestyle='dotted')
    # ax134[0].set_xlabel('% Gait cycle')
    # ax134[0].set_ylabel('Force (BW)')
    # ax134[0].set_title('welk008 natural')
    # ax134[0].legend()
    # # intersegmental forces average - exotendon
    # ax134[1].plot(e_timespercent101, np.mean(einterseg_combine[12:16,:], 0), label='exo_interseg')
    # ax134[1].plot(e_timespercent101, np.mean(einterseg_combine[12:16,:]+eextensors_combine[12:16,:], 0), label='exo_interseg + extensors')
    # ax134[1].plot(e_timespercent101, np.mean(einterseg_combine[12:16,:]+eextensors_combine[12:16,:]+ehams_combine[12:16,:], 0), label='exo_interseg+extensors+hams')
    # ax134[1].plot(e_timespercent101, np.mean(einterseg_combine[12:16,:]+eextensors_combine[12:16,:]+ehams_combine[12:16,:]+eadductors_combine[12:16,:], 0), label='exo_interseg+extensors+hams+adductors')
    # ax134[1].plot(e_timespercent101, np.mean(einterseg_combine[12:16,:]+eextensors_combine[12:16,:]+ehams_combine[12:16,:]+eadductors_combine[12:16,:]+eabductors_combine[12:16,:], 0), label='exo_interseg+extensors+hams+adductors+abductors')
    # ax134[1].plot(e_timespercent101, np.mean(einterseg_combine[12:16,:]+eextensors_combine[12:16,:]+ehams_combine[12:16,:]+eadductors_combine[12:16,:]+eabductors_combine[12:16,:]+eflexors_combine[12:16,:], 0), label='exo_interseg+extensors+hams+adductors+abductors+flexors')
    # ax134[1].plot(e_timespercent101, np.mean(einterseg_combine[12:16,:]+eextensors_combine[12:16,:]+ehams_combine[12:16,:]+eadductors_combine[12:16,:]+eabductors_combine[12:16,:]+eflexors_combine[12:16,:]+ereserve_combine[12:16,:], 0), label='exo_interseg+extensors+hams+adductors+abductors+flexors+reserve')
    # ax134[1].plot(e_timespercent101, np.mean(eall_combine[12:16,:], 0), label='exo_all', linestyle='dotted')
    # ax134[1].set_xlabel('% Gait cycle')
    # ax134[1].set_ylabel('Force (BW)')
    # ax134[1].set_title('welk008 exotendon')
    # ax134[1].legend()
    # fig134.tight_layout()
    # # welk009
    # # figures like ^ but individuals
    # fig135, ax135 = plt.subplots(1,2, figsize=(10,8), dpi=300)
    # # intersegmental forces average - natural welk002
    # ax135[0].plot(n_timespercent101, np.mean(ninterseg_combine[16:20,:], 0), label='natural_interseg')
    # ax135[0].plot(n_timespercent101, np.mean(ninterseg_combine[16:20,:] + nextensors_combine[16:20,:], 0), label='natural_interseg + extensors')
    # ax135[0].plot(n_timespercent101, np.mean(ninterseg_combine[16:20,:]+nextensors_combine[16:20,:]+nhams_combine[16:20,:], 0), label='natural_interseg+extensors+hams')
    # ax135[0].plot(n_timespercent101, np.mean(ninterseg_combine[16:20,:]+nextensors_combine[16:20,:]+nhams_combine[16:20,:]+nadductors_combine[16:20,:], 0), label='natural_interseg+extensors+hams+adductors')
    # ax135[0].plot(n_timespercent101, np.mean(ninterseg_combine[16:20,:]+nextensors_combine[16:20,:]+nhams_combine[16:20,:]+nadductors_combine[16:20,:]+nabductors_combine[16:20,:], 0), label='natural_interseg+extensors+hams+adductors+abductors')
    # ax135[0].plot(n_timespercent101, np.mean(ninterseg_combine[16:20,:]+nextensors_combine[16:20,:]+nhams_combine[16:20,:]+nadductors_combine[16:20,:]+nabductors_combine[16:20,:]+nflexors_combine[16:20,:], 0), label='natural_interseg+extensors+hams+adductors+abductors+flexors')
    # ax135[0].plot(n_timespercent101, np.mean(ninterseg_combine[16:20,:]+nextensors_combine[16:20,:]+nhams_combine[16:20,:]+nadductors_combine[16:20,:]+nabductors_combine[16:20,:]+nflexors_combine[16:20,:]+nreserve_combine[16:20,:], 0), label='natural_interseg+extensors+hams+adductors+abductors+flexors+reserve')
    # ax135[0].plot(n_timespercent101, np.mean(nall_combine[16:20,:], 0), label='nat_all', linestyle='dotted')
    # ax135[0].set_xlabel('% Gait cycle')
    # ax135[0].set_ylabel('Force (BW)')
    # ax135[0].set_title('welk009 natural')
    # ax135[0].legend()
    # # intersegmental forces average - exotendon
    # ax135[1].plot(e_timespercent101, np.mean(einterseg_combine[16:20,:], 0), label='exo_interseg')
    # ax135[1].plot(e_timespercent101, np.mean(einterseg_combine[16:20,:]+eextensors_combine[16:20,:], 0), label='exo_interseg + extensors')
    # ax135[1].plot(e_timespercent101, np.mean(einterseg_combine[16:20,:]+eextensors_combine[16:20,:]+ehams_combine[16:20,:], 0), label='exo_interseg+extensors+hams')
    # ax135[1].plot(e_timespercent101, np.mean(einterseg_combine[16:20,:]+eextensors_combine[16:20,:]+ehams_combine[16:20,:]+eadductors_combine[16:20,:], 0), label='exo_interseg+extensors+hams+adductors')
    # ax135[1].plot(e_timespercent101, np.mean(einterseg_combine[16:20,:]+eextensors_combine[16:20,:]+ehams_combine[16:20,:]+eadductors_combine[16:20,:]+eabductors_combine[16:20,:], 0), label='exo_interseg+extensors+hams+adductors+abductors')
    # ax135[1].plot(e_timespercent101, np.mean(einterseg_combine[16:20,:]+eextensors_combine[16:20,:]+ehams_combine[16:20,:]+eadductors_combine[16:20,:]+eabductors_combine[16:20,:]+eflexors_combine[16:20,:], 0), label='exo_interseg+extensors+hams+adductors+abductors+flexors')
    # ax135[1].plot(e_timespercent101, np.mean(einterseg_combine[16:20,:]+eextensors_combine[16:20,:]+ehams_combine[16:20,:]+eadductors_combine[16:20,:]+eabductors_combine[16:20,:]+eflexors_combine[16:20,:]+ereserve_combine[16:20,:], 0), label='exo_interseg+extensors+hams+adductors+abductors+flexors+reserve')
    # ax135[1].plot(e_timespercent101, np.mean(eall_combine[16:20,:], 0), label='exo_all', linestyle='dotted')
    # ax135[1].set_xlabel('% Gait cycle')
    # ax135[1].set_ylabel('Force (BW)')
    # ax135[1].set_title('welk009 exotendon')
    # ax135[1].legend()
    # fig135.tight_layout()
    # # welk010
    # # figures like ^ but individuals
    # fig136, ax136 = plt.subplots(1,2, figsize=(10,8), dpi=300)
    # # intersegmental forces average - natural welk002
    # ax136[0].plot(n_timespercent101, np.mean(ninterseg_combine[20:24,:], 0), label='natural_interseg')
    # ax136[0].plot(n_timespercent101, np.mean(ninterseg_combine[20:24,:] + nextensors_combine[20:24,:], 0), label='natural_interseg + extensors')
    # ax136[0].plot(n_timespercent101, np.mean(ninterseg_combine[20:24,:]+nextensors_combine[20:24,:]+nhams_combine[20:24,:], 0), label='natural_interseg+extensors+hams')
    # ax136[0].plot(n_timespercent101, np.mean(ninterseg_combine[20:24,:]+nextensors_combine[20:24,:]+nhams_combine[20:24,:]+nadductors_combine[20:24,:], 0), label='natural_interseg+extensors+hams+adductors')
    # ax136[0].plot(n_timespercent101, np.mean(ninterseg_combine[20:24,:]+nextensors_combine[20:24,:]+nhams_combine[20:24,:]+nadductors_combine[20:24,:]+nabductors_combine[20:24,:], 0), label='natural_interseg+extensors+hams+adductors+abductors')
    # ax136[0].plot(n_timespercent101, np.mean(ninterseg_combine[20:24,:]+nextensors_combine[20:24,:]+nhams_combine[20:24,:]+nadductors_combine[20:24,:]+nabductors_combine[20:24,:]+nflexors_combine[20:24,:], 0), label='natural_interseg+extensors+hams+adductors+abductors+flexors')
    # ax136[0].plot(n_timespercent101, np.mean(ninterseg_combine[20:24,:]+nextensors_combine[20:24,:]+nhams_combine[20:24,:]+nadductors_combine[20:24,:]+nabductors_combine[20:24,:]+nflexors_combine[20:24,:]+nreserve_combine[20:24,:], 0), label='natural_interseg+extensors+hams+adductors+abductors+flexors+reserve')
    # ax136[0].plot(n_timespercent101, np.mean(nall_combine[20:24,:], 0), label='nat_all', linestyle='dotted')
    # ax136[0].set_xlabel('% Gait cycle')
    # ax136[0].set_ylabel('Force (BW)')
    # ax136[0].set_title('welk010 natural')
    # ax136[0].legend()
    # # intersegmental forces average - exotendon
    # ax136[1].plot(e_timespercent101, np.mean(einterseg_combine[20:24,:], 0), label='exo_interseg')
    # ax136[1].plot(e_timespercent101, np.mean(einterseg_combine[20:24,:]+eextensors_combine[20:24,:], 0), label='exo_interseg + extensors')
    # ax136[1].plot(e_timespercent101, np.mean(einterseg_combine[20:24,:]+eextensors_combine[20:24,:]+ehams_combine[20:24,:], 0), label='exo_interseg+extensors+hams')
    # ax136[1].plot(e_timespercent101, np.mean(einterseg_combine[20:24,:]+eextensors_combine[20:24,:]+ehams_combine[20:24,:]+eadductors_combine[20:24,:], 0), label='exo_interseg+extensors+hams+adductors')
    # ax136[1].plot(e_timespercent101, np.mean(einterseg_combine[20:24,:]+eextensors_combine[20:24,:]+ehams_combine[20:24,:]+eadductors_combine[20:24,:]+eabductors_combine[20:24,:], 0), label='exo_interseg+extensors+hams+adductors+abductors')
    # ax136[1].plot(e_timespercent101, np.mean(einterseg_combine[20:24,:]+eextensors_combine[20:24,:]+ehams_combine[20:24,:]+eadductors_combine[20:24,:]+eabductors_combine[20:24,:]+eflexors_combine[20:24,:], 0), label='exo_interseg+extensors+hams+adductors+abductors+flexors')
    # ax136[1].plot(e_timespercent101, np.mean(einterseg_combine[20:24,:]+eextensors_combine[20:24,:]+ehams_combine[20:24,:]+eadductors_combine[20:24,:]+eabductors_combine[20:24,:]+eflexors_combine[20:24,:]+ereserve_combine[20:24,:], 0), label='exo_interseg+extensors+hams+adductors+abductors+flexors+reserve')
    # ax136[1].plot(e_timespercent101, np.mean(eall_combine[20:24,:], 0), label='exo_all', linestyle='dotted')
    # ax136[1].set_xlabel('% Gait cycle')
    # ax136[1].set_ylabel('Force (BW)')
    # ax136[1].set_title('welk010 exotendon')
    # ax136[1].legend()
    # fig136.tight_layout()
    # # welk013
    # # figures like ^ but individuals
    # fig137, ax137 = plt.subplots(1,2, figsize=(10,8), dpi=300)
    # # intersegmental forces average - natural welk002
    # ax137[0].plot(n_timespercent101, np.mean(ninterseg_combine[24:28,:], 0), label='natural_interseg')
    # ax137[0].plot(n_timespercent101, np.mean(ninterseg_combine[24:28,:] + nextensors_combine[24:28,:], 0), label='natural_interseg + extensors')
    # ax137[0].plot(n_timespercent101, np.mean(ninterseg_combine[24:28,:]+nextensors_combine[24:28,:]+nhams_combine[24:28,:], 0), label='natural_interseg+extensors+hams')
    # ax137[0].plot(n_timespercent101, np.mean(ninterseg_combine[24:28,:]+nextensors_combine[24:28,:]+nhams_combine[24:28,:]+nadductors_combine[24:28,:], 0), label='natural_interseg+extensors+hams+adductors')
    # ax137[0].plot(n_timespercent101, np.mean(ninterseg_combine[24:28,:]+nextensors_combine[24:28,:]+nhams_combine[24:28,:]+nadductors_combine[24:28,:]+nabductors_combine[24:28,:], 0), label='natural_interseg+extensors+hams+adductors+abductors')
    # ax137[0].plot(n_timespercent101, np.mean(ninterseg_combine[24:28,:]+nextensors_combine[24:28,:]+nhams_combine[24:28,:]+nadductors_combine[24:28,:]+nabductors_combine[24:28,:]+nflexors_combine[24:28,:], 0), label='natural_interseg+extensors+hams+adductors+abductors+flexors')
    # ax137[0].plot(n_timespercent101, np.mean(ninterseg_combine[24:28,:]+nextensors_combine[24:28,:]+nhams_combine[24:28,:]+nadductors_combine[24:28,:]+nabductors_combine[24:28,:]+nflexors_combine[24:28,:]+nreserve_combine[24:28,:], 0), label='natural_interseg+extensors+hams+adductors+abductors+flexors+reserve')
    # ax137[0].plot(n_timespercent101, np.mean(nall_combine[24:28,:], 0), label='nat_all', linestyle='dotted')
    # ax137[0].set_xlabel('% Gait cycle')
    # ax137[0].set_ylabel('Force (BW)')
    # ax137[0].set_title('welk013 natural')
    # ax137[0].legend()
    # # intersegmental forces average - exotendon
    # ax137[1].plot(e_timespercent101, np.mean(einterseg_combine[24:28,:], 0), label='exo_interseg')
    # ax137[1].plot(e_timespercent101, np.mean(einterseg_combine[24:28,:]+eextensors_combine[24:28,:], 0), label='exo_interseg + extensors')
    # ax137[1].plot(e_timespercent101, np.mean(einterseg_combine[24:28,:]+eextensors_combine[24:28,:]+ehams_combine[24:28,:], 0), label='exo_interseg+extensors+hams')
    # ax137[1].plot(e_timespercent101, np.mean(einterseg_combine[24:28,:]+eextensors_combine[24:28,:]+ehams_combine[24:28,:]+eadductors_combine[24:28,:], 0), label='exo_interseg+extensors+hams+adductors')
    # ax137[1].plot(e_timespercent101, np.mean(einterseg_combine[24:28,:]+eextensors_combine[24:28,:]+ehams_combine[24:28,:]+eadductors_combine[24:28,:]+eabductors_combine[24:28,:], 0), label='exo_interseg+extensors+hams+adductors+abductors')
    # ax137[1].plot(e_timespercent101, np.mean(einterseg_combine[24:28,:]+eextensors_combine[24:28,:]+ehams_combine[24:28,:]+eadductors_combine[24:28,:]+eabductors_combine[24:28,:]+eflexors_combine[24:28,:], 0), label='exo_interseg+extensors+hams+adductors+abductors+flexors')
    # ax137[1].plot(e_timespercent101, np.mean(einterseg_combine[24:28,:]+eextensors_combine[24:28,:]+ehams_combine[24:28,:]+eadductors_combine[24:28,:]+eabductors_combine[24:28,:]+eflexors_combine[24:28,:]+ereserve_combine[24:28,:], 0), label='exo_interseg+extensors+hams+adductors+abductors+flexors+reserve')
    # ax137[1].plot(e_timespercent101, np.mean(eall_combine[24:28,:], 0), label='exo_all', linestyle='dotted')
    # ax137[1].set_xlabel('% Gait cycle')
    # ax137[1].set_ylabel('Force (BW)')
    # ax137[1].set_title('welk013 exotendon')
    # ax137[1].legend()
    # fig137.tight_layout()
    


    ###########################################################################
    ### Figure: total pop average for right leg between nat and exo
    #### Figure in the paper R1... 
    figcon6, axcon6 = plt.subplots(1, 2, figsize=(12,4.55), dpi=500)
    axcon6[0].fill_between(n_timespercent101, np.mean(nall_combine,0)-np.std(nall_combine,0), np.mean(nall_combine,0)+np.std(nall_combine,0), color=ncolorlight, alpha=0.75)
    axcon6[0].fill_between(e_timespercent101, np.mean(eall_combine,0)-np.std(eall_combine,0), np.mean(eall_combine,0)+np.std(eall_combine,0), color=ecolorlight, alpha=0.75)
    axcon6[0].plot(n_timespercent101, np.mean(nall_combine,0), color=ncolor, label='Natural (Mean \u00B1 Std.)', linewidth=3)
    axcon6[0].plot(e_timespercent101, np.mean(eall_combine,0), color=ecolor, label='Exotendon (Mean \u00B1 Std.)', linewidth=3)
    axcon6[0].set_xlabel('% Gait cycle', fontsize=16)
    axcon6[0].set_ylabel('Vertical hip contact force (BW)', fontsize=16)
    axcon6[0].set_title('Total vertical contact force', fontsize=16)
    axcon6[0].tick_params(axis='both', which='major', labelsize=16)
    # hide the second subplot axis and use it for legend
    axcon6[1].axis('off')
    handles, labels = axcon6[0].get_legend_handles_labels()
    axcon6[1].legend(handles, labels, loc='center', fontsize=16)
    figcon6.tight_layout()
    figcon6.savefig(analyzedir + '\\contact6_hip_' + whichleg + '.png')


    ###########################################################################
    # output the data in an excel sheet
    # Create a dictionary to hold the data
    data_dict = {
        'ninterseg_combine': ninterseg_combine,
        'einterseg_combine': einterseg_combine,
        'nextensors_combine': nextensors_combine,
        'eextensors_combine': eextensors_combine,
        'nadductors_combine': nadductors_combine,
        'eadductors_combine': eadductors_combine,
        'nabductors_combine': nabductors_combine,
        'eabductors_combine': eabductors_combine,
        'nhams_combine': nhams_combine,
        'ehams_combine': ehams_combine,
        'nflexors_combine': nflexors_combine,
        'eflexors_combine': eflexors_combine,
        'nreserve_combine': nreserve_combine,
        'ereserve_combine': ereserve_combine,
        'nall_combine': nall_combine,
        'eall_combine': eall_combine,
        'nnone_combine': nnone_combine,
        'enone_combine': enone_combine,
    }

    # 

    os.chdir(analyzedir)
    # Create a Pandas Excel writer using XlsxWriter as the engine.
    with pd.ExcelWriter('analysis_results_hip'+'_'+whichleg+'.xlsx', engine='xlsxwriter') as writer:
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
    plt.fill_between(n_timespercent101, 
                     np.mean(ninterseg_combine,0), color=ecolorlight)
    
    plt.plot(n_timespercent101, np.mean(nextensors_combine,0) + np.mean(ninterseg_combine,0), 
             label='inter + extensors', color=ncolor2)
    plt.fill_between(n_timespercent101, np.mean(ninterseg_combine,0), 
                     np.mean(nextensors_combine,0) + np.mean(ninterseg_combine,0), color=ncolor1)
    
    plt.plot(n_timespercent101, np.mean(nadductors_combine,0) + np.mean(nextensors_combine,0) + np.mean(ninterseg_combine,0), 
             label='inter + extensors + adductors', color=ncolor4)
    plt.fill_between(n_timespercent101, np.mean(nextensors_combine,0) + np.mean(ninterseg_combine,0), 
                     np.mean(nadductors_combine,0) + np.mean(nextensors_combine,0) + np.mean(ninterseg_combine,0), color=ncolor3)
    
    plt.plot(n_timespercent101, np.mean(nabductors_combine,0) + np.mean(nadductors_combine,0) + np.mean(nextensors_combine,0) + np.mean(ninterseg_combine,0), 
             label='inter + extensors + adductors + abductors', color=ncolor6)
    plt.fill_between(n_timespercent101, np.mean(nextensors_combine,0) + np.mean(ninterseg_combine,0) + np.mean(nadductors_combine,0), 
                     np.mean(nabductors_combine,0) + np.mean(nadductors_combine,0) + np.mean(nextensors_combine,0) + np.mean(ninterseg_combine,0), color=ncolor5)

    plt.plot(n_timespercent101, np.mean(nhams_combine,0) + np.mean(nabductors_combine,0) + np.mean(nadductors_combine,0) + np.mean(nextensors_combine,0) + np.mean(ninterseg_combine,0), 
             label='inter + extensors + adductors + abductors + hams', color=ncolor8)
    plt.fill_between(n_timespercent101, np.mean(nextensors_combine,0) + np.mean(ninterseg_combine,0) + np.mean(nadductors_combine,0) + np.mean(nabductors_combine,0), 
                     np.mean(nhams_combine,0) + np.mean(nabductors_combine,0) + np.mean(nadductors_combine,0) + np.mean(nextensors_combine,0) + np.mean(ninterseg_combine,0), color=ncolor7)
    
    plt.plot(n_timespercent101, np.mean(nflexors_combine,0) + np.mean(nabductors_combine,0) + np.mean(nhams_combine,0) + np.mean(nadductors_combine,0) + np.mean(nextensors_combine,0) + np.mean(ninterseg_combine,0), 
             label='inter + extensors + adductors + abductors + hams + flexors', color=ncolor10)
    plt.fill_between(n_timespercent101, np.mean(nextensors_combine,0) + np.mean(ninterseg_combine,0) + np.mean(nadductors_combine,0) + np.mean(nabductors_combine,0) + np.mean(nhams_combine,0), 
                     np.mean(nflexors_combine,0) + np.mean(nhams_combine,0) + np.mean(nabductors_combine,0) + np.mean(nadductors_combine,0) + np.mean(nextensors_combine,0) + np.mean(ninterseg_combine,0), color=ncolor9)
    
    # plt.legend()
    plt.ylabel('Vertical hip contact force (BW)', fontsize=16)
    plt.xlabel('% Gait cycle', fontsize=16)
    plt.title('Vertical hip contact force - Natural', fontsize=16)
    plt.tick_params(axis='both', which='major', labelsize=16)
    plt.tight_layout()
    plt.savefig(analyzedir + '\\contact7_hip_natural_' + whichleg + '.png')
    plt.show()
    
    ###########################################################################
    ### Figure: total pop average for right leg for exo - segmented shaded
    plt.figure(dpi=300)
    plt.plot(e_timespercent101, np.mean(einterseg_combine,0), label='intersegmental', color=ecolor)
    plt.fill_between(e_timespercent101, 
                     np.mean(einterseg_combine,0), color=ecolorlight)
    
    plt.plot(e_timespercent101, np.mean(eextensors_combine,0) + np.mean(einterseg_combine,0), 
             label='inter + extensors', color=ncolor2)
    plt.fill_between(e_timespercent101, np.mean(einterseg_combine,0), 
                     np.mean(eextensors_combine,0) + np.mean(einterseg_combine,0), color=ncolor1)
    
    plt.plot(e_timespercent101, np.mean(eadductors_combine,0) + np.mean(eextensors_combine,0) + np.mean(einterseg_combine,0), 
             label='inter + extensors + adductors', color=ncolor4)
    plt.fill_between(e_timespercent101, np.mean(eextensors_combine,0) + np.mean(einterseg_combine,0), 
                     np.mean(eadductors_combine,0) + np.mean(eextensors_combine,0) + np.mean(einterseg_combine,0), color=ncolor3)
    
    plt.plot(e_timespercent101, np.mean(eabductors_combine,0) + np.mean(eadductors_combine,0) + np.mean(eextensors_combine,0) + np.mean(einterseg_combine,0), 
             label='inter + extensors + adductors + abductors', color=ncolor6)
    plt.fill_between(e_timespercent101, np.mean(eextensors_combine,0) + np.mean(einterseg_combine,0) + np.mean(eadductors_combine,0), 
                     np.mean(eabductors_combine,0) + np.mean(eadductors_combine,0) + np.mean(eextensors_combine,0) + np.mean(einterseg_combine,0), color=ncolor5)

    plt.plot(e_timespercent101, np.mean(ehams_combine,0) + np.mean(eabductors_combine,0) + np.mean(eadductors_combine,0) + np.mean(eextensors_combine,0) + np.mean(einterseg_combine,0), 
             label='inter + extensors + adductors + abductors + hams', color=ncolor8)
    plt.fill_between(e_timespercent101, np.mean(eextensors_combine,0) + np.mean(einterseg_combine,0) + np.mean(eadductors_combine,0) + np.mean(eabductors_combine,0), 
                     np.mean(ehams_combine,0) + np.mean(eabductors_combine,0) + np.mean(eadductors_combine,0) + np.mean(eextensors_combine,0) + np.mean(einterseg_combine,0), color=ncolor7)
    
    plt.plot(e_timespercent101, np.mean(eflexors_combine,0) + np.mean(eabductors_combine,0) + np.mean(ehams_combine,0) + np.mean(eadductors_combine,0) + np.mean(eextensors_combine,0) + np.mean(einterseg_combine,0), 
             label='inter + extensors + adductors + abductors + hams + flexors', color=ncolor10)
    plt.fill_between(e_timespercent101, np.mean(eextensors_combine,0) + np.mean(einterseg_combine,0) + np.mean(eadductors_combine,0) + np.mean(eabductors_combine,0) + np.mean(ehams_combine,0), 
                     np.mean(eflexors_combine,0) + np.mean(ehams_combine,0) + np.mean(eabductors_combine,0) + np.mean(eadductors_combine,0) + np.mean(eextensors_combine,0) + np.mean(einterseg_combine,0), color=ncolor9)
    
    plt.legend()
    plt.ylabel('Vertical hip contact force (BW)', fontsize=16)
    plt.xlabel('% Gait cycle', fontsize=16)
    plt.title('Vertical hip contact force - exotendon', fontsize=16)
    plt.tick_params(axis='both', which='major', labelsize=16)
    plt.tight_layout()
    plt.savefig(analyzedir + '\\contact8_hip_exo_' + whichleg + '.png')
    plt.show()

    
    
    ###########################################################################
    # stats for the individual muscle segmented contributions.... check out the excel sheet
    peaks_extensors_natural = np.max(nextensors_combine,1)
    peaks_adductors_natural = np.max(nadductors_combine,1)
    peaks_hams_natural = np.max(nhams_combine,1)
    peaks_abductors_natural = np.max(nabductors_combine,1)
    peaks_flexors_natural = np.max(nflexors_combine,1)
    peaks_inter_natural = np.max(ninterseg_combine,1)
    peaks_reserve_natural = np.max(nreserve_combine,1)

    peaks_extensors_exo = np.max(eextensors_combine,1)
    peaks_adductors_exo = np.max(eadductors_combine,1)
    peaks_hams_exo = np.max(ehams_combine,1)
    peaks_abductors_exo = np.max(eabductors_combine,1)
    peaks_flexors_exo = np.max(eflexors_combine,1)
    peaks_inter_exo = np.max(einterseg_combine,1)
    peaks_reserve_exo = np.max(ereserve_combine,1)

    # copy what you need and go to excel... 
    
    ###########################################################################
    
    pdb.set_trace()
    sys.exit()
    

    ###########################################################################
    ### fig: subplots of nat and exo for each breakdown
    fig3, ax3 = plt.subplots(1,6)
    # full forces
    # ax3[0].plot(n_timespercent101, n_avg_r.mean(0), label='n_full', color='orange')
    # ax3[0].plot(e_timespercent101, e_avg_r.mean(0), label='e_full', color='purple')
    # intersegmental forces
    ax3[1].plot(n_timespercent101, ninterseg_combine.mean(0), label='n_interseg', color='orange')
    # ax3[1].plot(e_timespercent101, einterseg_combine.mean(0), label='e_interseg', color='purple')
    ax3[1].set_xlabel('% gait cycle')
    ax3[1].set_ylabel('BW - intersegmental')
    # # quads forces
    # ax3[2].plot(n_timespercent101, nquads_combine.mean(0), label='n_quads', color='orange')
    # ax3[2].plot(e_timespercent101, equads_combine.mean(0), label='e_quads', color='purple')
    # ax3[2].set_xlabel('% gait cycle')
    # ax3[2].set_ylabel('BW - quads')
    # # hams forces
    # ax3[3].plot(n_timespercent101, nhams_combine.mean(0), label='n_hams', color='orange')
    # ax3[3].plot(e_timespercent101, ehams_combine.mean(0), label='e_hams', color='purple')
    # ax3[3].set_xlabel('% gait cycle')
    # ax3[3].set_ylabel('BW - hams')
    # # gastroc forces
    # ax3[4].plot(n_timespercent101, ngas_combine.mean(0), label='n_gas', color='orange')
    # ax3[4].plot(e_timespercent101, egas_combine.mean(0), label='e_gas', color='purple')
    # ax3[4].set_xlabel('% gait cycle')
    # ax3[4].set_ylabel('BW - gastroc')
    # # tfl forces
    # ax3[5].plot(n_timespercent101, ntfl_combine.mean(0), label='n_tfl', color='orange')
    # ax3[5].plot(e_timespercent101, etfl_combine.mean(0), label='e_tfl', color='purple')
    # ax3[5].set_xlabel('% gait cycle')
    # ax3[5].set_ylabel('BW - tfl')
    
    plt.legend()
    plt.show()
    
    pdb.set_trace()
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

    '''
    pulling numbers specifically from the force reporter for muscles to gut-check 
    the magnitudes that they had, comparing them to their contributions
    This is not exact matching, but rather should be within the ball park to get 
    a good idea of if the analysis is properly applying all the forces etc. 
    
    CASE STUDY FOR WELK002, WELKEXO, TRIAL01
    
    #### FORCE REPORTER FROM ANALYSIS TOOL ####
    TFL:
        TFL: ========== 340
    Gastroc: ========== 1463
        medial = 766
        lateral = 697
    Quads: ============ 6117
        vaslat = 3985
        vasmed = 1154
        vasint = 697
        recfem = 281
    hamstrings: ======= 1521
        bflh = 254
        bfsh = 144
        gracilis = 64
        sartorius = 67
        semimem = 899
        semiten = 93
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    #### CONTRIBUTION FORCES FROM JRA ####
    TFL ===== 335
    GASTROC = 1420
    QUADS === 5824
    HAMS ==== 1180
    
    ######################################################################
    # SEEMS LIKE IT IS RELATIVELY GOOD FIT. 
    ######################################################################
    
    
    CASE STUDY FOR WELK013, WELKNATURAL TRIAL04 
    #### FORCE REPORTER FROM ANALYSIS TOOL ####
    TFL:
        TFL: ========== 
    Gastroc: ========== 
        medial = 
        lateral = 
    Quads: ============ 
        vaslat = 
        vasmed = 
        vasint = 
        recfem = 
    hamstrings: ======= 
        bflh = 
        bfsh = 
        gracilis = 
        sartorius = 
        semimem = 
        semiten = 
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    #### CONTRIBUTION FORCES FROM JRA #### DOING THE NORM
    TFL ===== 
    GASTROC = 
    QUADS === 
    HAMS ==== 
    
    '''
    
    
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
    plt.plot(e_timespercent101, 
             nquadsub.mean(0) + ngassub.mean(0) + nhamssub.mean(0) + ntflsub.mean(0) + ninterseg_combine.mean(0), 
             label='n_quads', color='#d94701')
    # plt.plot(n_timespercent101, 
    #           equadsub.mean(0) + egassub.mean(0) + ehamssub.mean(0) + etflsub.mean(0) + einterseg_combine.mean(0), 
    #           label='e_quads', linestyle='dashed', color='#a63603')
    # gastroc
    plt.plot(e_timespercent101, 
             ngassub.mean(0) + nhamssub.mean(0) + ntflsub.mean(0) + ninterseg_combine.mean(0), 
             label='n_gas', color='#e6550d')
    # plt.plot(n_timespercent101, 
    #           egassub.mean(0) + ehamssub.mean(0) + etflsub.mean(0) + einterseg_combine.mean(0), 
    #           label='e_gas', linestyle='dashed', color='#b2abd2')
    # hams
    plt.plot(e_timespercent101, nhamssub.mean(0) + ntflsub.mean(0) + ninterseg_combine.mean(0),
             label='n_hams', color='#fd8d3c')
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
                
           
                
                
                
                
                
                
                
                
                
                
                