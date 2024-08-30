import os
import opensim as osim
import numpy as np
import pdb
import time
import matplotlib.pyplot as plt
import scipy


# naturalcolor = '#fdb863'
ncolor = '#e66101'
ncolorlight = '#fee0b6'

ncolor1 = '#fff5eb'
ncolor2 = '#fee6ce'
ncolor3 = '#fdd0a2'
ncolor4 = '#fdae6b'
ncolor5 = '#fd8d3c'
ncolor6 = '#f16913'
ncolor7 = '#d94801'
ncolor8 = '#a63603'
ncolor9 = '#7f2704'



# exotendoncolor = '#f1a340'
ecolor = '#5e3c99'
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
    return modelmass*grav
# function that actually runs the analysis to compute knee contact force, and
# transforms it to the tibia frame
def computeKneeContact(trimmodel, initTime, finalTime, trialdir, tag):
    # '''
    # intersegmental forces - method 2
    # try the analysis
    jr_tool = osim.AnalyzeTool()
    jr_tool.setName('jr_analysis_2Dpred')
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
    jr_tool.run()
    time.sleep(0.5)
    # '''
    # figure out how to do an intersegmental with proper value of forces showing up.
    trimjra = osim.TimeSeriesTable('jr_analysis_2Dpred_jra_' + tag + '_ReactionLoads.sto')
    trimjralabels = trimjra.getColumnLabels()
    tiby = trimjra.getDependentColumn('walker_knee_r_on_tibia_r_in_tibia_r_fy').to_numpy()
    # pdb.set_trace()
    # import matplotlib.pyplot as plt
    # plt.figure()
    # plt.plot(np.array(trimjra.getIndependentColumn()), tiby)
    return tiby

# method for computing the individual muscle contributions to knee contact force
def getKneeContactributions(trialdir, musclesWanted_split, tag):
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
        statesStorage = osim.Storage('muscletrack_states_2Dpred.sto')
        statesTable = osim.TimeSeriesTable('muscletrack_states_2Dpred.sto')
        stateslabels = statesTable.getColumnLabels()
        statesTableTrim = osim.TimeSeriesTable('muscletrack_states_2Dpred.sto')
        statesTableTrim2 = osim.TimeSeriesTable('muscletrack_states_2Dpred.sto')
        
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
        controlsTable = osim.TimeSeriesTable('muscletrack_controls_2Dpred.sto')
        controlslabels = controlsTable.getColumnLabels()
        controlsTableTrim = osim.TimeSeriesTable('muscletrack_controls_2Dpred.sto')
        
        # get a version of the controls that is trimmed down. 
        for con in controlslabels:
            if 'reserve' not in con and 'lumbar' not in con:
                controlsTableTrim.removeColumn(con)
        osim.STOFileAdapter.write(controlsTableTrim, 'trimmingControls_' + tag + '.sto')
        
        # get a version of the model with no muscles in it
        # muscmodel = osim.Model('post_simple_model_all_the_probes_muscletrack.osim')
        trimmodelprocessor = osim.ModelProcessor('mk12_rv1_dgf.osim')
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
        jray = computeKneeContact(trimmodel, initTime, finalTime, trialdir, tag)
    
    
    elif 'all' in musclesWanted_split:
        statesStorage = osim.Storage('muscletrack_states_2Dpred.sto')
        statesTable = osim.TimeSeriesTable('muscletrack_states_2Dpred.sto')
        stateslabels = statesTable.getColumnLabels()
        statesTableTrim = osim.TimeSeriesTable('muscletrack_states_2Dpred.sto')
        statesTableTrim2 = osim.TimeSeriesTable('muscletrack_states_2Dpred.sto')
        # in this case, we want all the muscles in the model
        osim.STOFileAdapter.write(statesTableTrim, 'trimmingStates_' + tag + '.sto')
        
        for stat in stateslabels:
            if 'speed' in stat:
                statesTableTrim2.removeColumn(stat)
        osim.STOFileAdapter.write(statesTableTrim2, 'trimmingStates2_' + tag + '.sto')


        # get a version of the controls that matches. 
        controlsTable = osim.TimeSeriesTable('muscletrack_controls_2Dpred.sto')
        controlslabels = controlsTable.getColumnLabels()
        controlsTableTrim = osim.TimeSeriesTable('muscletrack_controls_2Dpred.sto')
        
        # in this case we want all the controls, not getting rid of any muscles
        osim.STOFileAdapter.write(controlsTableTrim, 'trimmingControls_' + tag + '.sto')
        
        
        # get a version of the model with no muscles in it
        # muscmodel = osim.Model('post_simple_model_all_the_probes_muscletrack.osim')
        trimmodel = osim.Model('mk12_rv1_dgf.osim')
        
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
        jray = computeKneeContact(trimmodel, initTime, finalTime, trialdir, tag)
    
    elif 'reserve' in musclesWanted_split:
        statesStorage = osim.Storage('muscletrack_states_2Dpred.sto')
        statesTable = osim.TimeSeriesTable('muscletrack_states_2Dpred.sto')
        stateslabels = statesTable.getColumnLabels()
        statesTableTrim = osim.TimeSeriesTable('muscletrack_states_2Dpred.sto')
        statesTableTrim2 = osim.TimeSeriesTable('muscletrack_states_2Dpred.sto')
        
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
        controlsTable = osim.TimeSeriesTable('muscletrack_controls_2Dpred.sto')
        controlslabels = controlsTable.getColumnLabels()
        controlsTableTrim = osim.TimeSeriesTable('muscletrack_controls_2Dpred.sto')
        
        # get a version of the controls that is trimmed down. 
        for con in controlslabels:
            if 'lumbar' not in con:
                controlsTableTrim.removeColumn(con)
        osim.STOFileAdapter.write(controlsTableTrim, 'trimmingControls_' + tag + '.sto')
        
        # get a version of the model with no muscles in it (or reserves?)
        # muscmodel = osim.Model('post_simple_model_all_the_probes_muscletrack.osim')
        trimmodelprocessor = osim.ModelProcessor('mk12_rv1_dgf.osim')
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
        jray = computeKneeContact(trimmodel, initTime, finalTime, trialdir, tag)
    
    elif 'none' in musclesWanted_split:
        statesStorage = osim.Storage('muscletrack_states_2Dpred.sto')
        statesTable = osim.TimeSeriesTable('muscletrack_states_2Dpred.sto')
        stateslabels = statesTable.getColumnLabels()
        statesTableTrim = osim.TimeSeriesTable('muscletrack_states_2Dpred.sto')
        statesTableTrim2 = osim.TimeSeriesTable('muscletrack_states_2Dpred.sto')
        
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
        controlsTable = osim.TimeSeriesTable('muscletrack_controls_2Dpred.sto')
        controlslabels = controlsTable.getColumnLabels()
        controlsTableTrim = osim.TimeSeriesTable('muscletrack_controls_2Dpred.sto')
        
        # get a version of the controls that is trimmed down. 
        for con in controlslabels:
            controlsTableTrim.removeColumn(con)
        osim.STOFileAdapter.write(controlsTableTrim, 'trimmingControls_' + tag + '.sto')
        
        # get a version of the model with no muscles in it (or reserves?)
        # muscmodel = osim.Model('post_simple_model_all_the_probes_muscletrack.osim')
        trimmodelprocessor = osim.ModelProcessor('mk12_rv1_dgf.osim')
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
        jray = computeKneeContact(trimmodel, initTime, finalTime, trialdir, tag)
        
    
    else:
        statesStorage = osim.Storage('muscletrack_states_2Dpred.sto')
        statesTable = osim.TimeSeriesTable('muscletrack_states_2Dpred.sto')
        stateslabels = statesTable.getColumnLabels()
        statesTableTrim = osim.TimeSeriesTable('muscletrack_states_2Dpred.sto')
        statesTableTrim2 = osim.TimeSeriesTable('muscletrack_states_2Dpred.sto')
        
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
        controlsTable = osim.TimeSeriesTable('muscletrack_controls_2Dpred.sto')
        controlslabels = controlsTable.getColumnLabels()
        controlsTableTrim = osim.TimeSeriesTable('muscletrack_controls_2Dpred.sto')
        
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
        trimmodel = osim.Model('mk12_rv1_dgf.osim')
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
        jray = computeKneeContact(trimmodel, initTime, finalTime, trialdir, tag)
    
    return jray
# method for computing the individual muscle contributions to knee contact force
def getKneeContactributions001(trialdir, musclesWanted_split, tag):
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
    if musclesWanted_split == []:
        statesStorage = osim.Storage('muscletrack_states_100con001_2.sto')
        statesTable = osim.TimeSeriesTable('muscletrack_states_100con001_2.sto')
        stateslabels = statesTable.getColumnLabels()
        statesTableTrim = osim.TimeSeriesTable('muscletrack_states_100con001_2.sto')
        statesTableTrim2 = osim.TimeSeriesTable('muscletrack_states_100con001_2.sto')
        
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
        controlsTable = osim.TimeSeriesTable('muscletrack_controls_100con001_2.sto')
        controlslabels = controlsTable.getColumnLabels()
        controlsTableTrim = osim.TimeSeriesTable('muscletrack_controls_100con001_2.sto')
        
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
        jray = computeKneeContact(trimmodel, initTime, finalTime, trialdir, tag)

        
    else:
        statesStorage = osim.Storage('muscletrack_states_100con001_2.sto')
        statesTable = osim.TimeSeriesTable('muscletrack_states_100con001_2.sto')
        stateslabels = statesTable.getColumnLabels()
        statesTableTrim = osim.TimeSeriesTable('muscletrack_states_100con001_2.sto')
        statesTableTrim2 = osim.TimeSeriesTable('muscletrack_states_100con001_2.sto')
        
        # musclesWanted_split = ['1'2'3'4']
        # get a trimmed set of states that is just the joint angles speeds, and whatever muscles you want. 
        for stat in stateslabels:
            if 'forceset' in stat:
                getrid = True
                for want in musclesWanted_split:
                    if want in stat:
                        # want this one
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
        controlsTable = osim.TimeSeriesTable('muscletrack_controls_100con001_2.sto')
        controlslabels = controlsTable.getColumnLabels()
        controlsTableTrim = osim.TimeSeriesTable('muscletrack_controls_100con001_2.sto')
        
        # get a version of the controls that is trimmed down. 
        for con in controlslabels:
            if 'reserve' not in con and 'lumbar' not in con:
                getrid = True
                for want in musclesWanted_split:
                    if want in con:
                        # want this
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
            if 'reserve' not in fo.getName() and 'lumbar' not in fo.getName():
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
        jray = computeKneeContact(trimmodel, initTime, finalTime, trialdir, tag)
    
    return jray
# method for computing the individual muscle contributions to knee contact force
def getKneeContactributions0001(trialdir, musclesWanted_split, tag):
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
    if musclesWanted_split == []:
        statesStorage = osim.Storage('muscletrack_states_100con0001_2.sto')
        statesTable = osim.TimeSeriesTable('muscletrack_states_100con0001_2.sto')
        stateslabels = statesTable.getColumnLabels()
        statesTableTrim = osim.TimeSeriesTable('muscletrack_states_100con0001_2.sto')
        statesTableTrim2 = osim.TimeSeriesTable('muscletrack_states_100con0001_2.sto')
        
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
        controlsTable = osim.TimeSeriesTable('muscletrack_controls_100con0001_2.sto')
        controlslabels = controlsTable.getColumnLabels()
        controlsTableTrim = osim.TimeSeriesTable('muscletrack_controls_100con0001_2.sto')
        
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
        jray = computeKneeContact(trimmodel, initTime, finalTime, trialdir, tag)

        
    else:
        statesStorage = osim.Storage('muscletrack_states_100con0001_2.sto')
        statesTable = osim.TimeSeriesTable('muscletrack_states_100con0001_2.sto')
        stateslabels = statesTable.getColumnLabels()
        statesTableTrim = osim.TimeSeriesTable('muscletrack_states_100con0001_2.sto')
        statesTableTrim2 = osim.TimeSeriesTable('muscletrack_states_100con0001_2.sto')
        
        # musclesWanted_split = ['1'2'3'4']
        # get a trimmed set of states that is just the joint angles speeds, and whatever muscles you want. 
        for stat in stateslabels:
            if 'forceset' in stat:
                getrid = True
                for want in musclesWanted_split:
                    if want in stat:
                        # want this one
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
        controlsTable = osim.TimeSeriesTable('muscletrack_controls_100con0001_2.sto')
        controlslabels = controlsTable.getColumnLabels()
        controlsTableTrim = osim.TimeSeriesTable('muscletrack_controls_100con0001_2.sto')
        
        # get a version of the controls that is trimmed down. 
        for con in controlslabels:
            if 'reserve' not in con and 'lumbar' not in con:
                getrid = True
                for want in musclesWanted_split:
                    if want in con:
                        # want this
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
            if 'reserve' not in fo.getName() and 'lumbar' not in fo.getName():
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
        jray = computeKneeContact(trimmodel, initTime, finalTime, trialdir, tag)
    
    return jray
def getKneeContactributions01_wt1(trialdir, musclesWanted_split, tag):
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
    if musclesWanted_split == []:
        statesStorage = osim.Storage('muscletrack_states_100con01_wt1.sto')
        statesTable = osim.TimeSeriesTable('muscletrack_states_100con01_wt1.sto')
        stateslabels = statesTable.getColumnLabels()
        statesTableTrim = osim.TimeSeriesTable('muscletrack_states_100con01_wt1.sto')
        statesTableTrim2 = osim.TimeSeriesTable('muscletrack_states_100con01_wt1.sto')
        
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
        controlsTable = osim.TimeSeriesTable('muscletrack_controls_100con01_wt1.sto')
        controlslabels = controlsTable.getColumnLabels()
        controlsTableTrim = osim.TimeSeriesTable('muscletrack_controls_100con01_wt1.sto')
        
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
        jray = computeKneeContact(trimmodel, initTime, finalTime, trialdir, tag)

        
    else:
        statesStorage = osim.Storage('muscletrack_states_100con01_wt1.sto')
        statesTable = osim.TimeSeriesTable('muscletrack_states_100con01_wt1.sto')
        stateslabels = statesTable.getColumnLabels()
        statesTableTrim = osim.TimeSeriesTable('muscletrack_states_100con01_wt1.sto')
        statesTableTrim2 = osim.TimeSeriesTable('muscletrack_states_100con01_wt1.sto')
        
        # musclesWanted_split = ['1'2'3'4']
        # get a trimmed set of states that is just the joint angles speeds, and whatever muscles you want. 
        for stat in stateslabels:
            if 'forceset' in stat:
                getrid = True
                for want in musclesWanted_split:
                    if want in stat:
                        # want this one
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
        controlsTable = osim.TimeSeriesTable('muscletrack_controls_100con01_wt1.sto')
        controlslabels = controlsTable.getColumnLabels()
        controlsTableTrim = osim.TimeSeriesTable('muscletrack_controls_100con01_wt1.sto')
        
        # get a version of the controls that is trimmed down. 
        for con in controlslabels:
            if 'reserve' not in con and 'lumbar' not in con:
                getrid = True
                for want in musclesWanted_split:
                    if want in con:
                        # want this
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
            if 'reserve' not in fo.getName() and 'lumbar' not in fo.getName():
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
        jray = computeKneeContact(trimmodel, initTime, finalTime, trialdir, tag)
    
    return jray
def getKneeContactributions01_wt2(trialdir, musclesWanted_split, tag):
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
    if musclesWanted_split == []:
        statesStorage = osim.Storage('muscletrack_states_100con01_wt2.sto')
        statesTable = osim.TimeSeriesTable('muscletrack_states_100con01_wt2.sto')
        stateslabels = statesTable.getColumnLabels()
        statesTableTrim = osim.TimeSeriesTable('muscletrack_states_100con01_wt2.sto')
        statesTableTrim2 = osim.TimeSeriesTable('muscletrack_states_100con01_wt2.sto')
        
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
        controlsTable = osim.TimeSeriesTable('muscletrack_controls_100con01_wt2.sto')
        controlslabels = controlsTable.getColumnLabels()
        controlsTableTrim = osim.TimeSeriesTable('muscletrack_controls_100con01_wt2.sto')
        
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
        jray = computeKneeContact(trimmodel, initTime, finalTime, trialdir, tag)

        
    else:
        statesStorage = osim.Storage('muscletrack_states_100con01_wt2.sto')
        statesTable = osim.TimeSeriesTable('muscletrack_states_100con01_wt2.sto')
        stateslabels = statesTable.getColumnLabels()
        statesTableTrim = osim.TimeSeriesTable('muscletrack_states_100con01_wt2.sto')
        statesTableTrim2 = osim.TimeSeriesTable('muscletrack_states_100con01_wt2.sto')
        
        # musclesWanted_split = ['1'2'3'4']
        # get a trimmed set of states that is just the joint angles speeds, and whatever muscles you want. 
        for stat in stateslabels:
            if 'forceset' in stat:
                getrid = True
                for want in musclesWanted_split:
                    if want in stat:
                        # want this one
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
        controlsTable = osim.TimeSeriesTable('muscletrack_controls_100con01_wt2.sto')
        controlslabels = controlsTable.getColumnLabels()
        controlsTableTrim = osim.TimeSeriesTable('muscletrack_controls_100con01_wt2.sto')
        
        # get a version of the controls that is trimmed down. 
        for con in controlslabels:
            if 'reserve' not in con and 'lumbar' not in con:
                getrid = True
                for want in musclesWanted_split:
                    if want in con:
                        # want this
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
            if 'reserve' not in fo.getName() and 'lumbar' not in fo.getName():
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
        jray = computeKneeContact(trimmodel, initTime, finalTime, trialdir, tag)
    
    return jray



if __name__ == '__main__':
    # now to define all the setup that he has and is required
    basedir = os.getcwd()
    repodir = 'G:\\Shared drives\\Exotendon\\muscleModel\\muscleEnergyModel';
    resultsdir = os.path.join(repodir, '..\\results');

    welkexoconditions = ['welkexo']
    welknaturalconditions = ['welknatural']
    welksubjects = ['welk005'] # ['welk002','welk003','welk005','welk008','welk009','welk010','welk013'];
    # welksubjects = ['example3DWalking']
    thingstoplot = ['contactForces'];
    trials = ['015_3ActMet']    #'trial01','trial02','trial03','trial04']


    # get some results structures going
    welknaturalstruct_combine = {}
    welkexostruct_combine = {}
    naturalstruct_combine = {}
    exostruct_combine = {}
    naturalstruct_avg = {}
    exostruct_avg = {}


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
    
    
    
    # define the muscles that we want in each of the splits
    musclesWanted = {}
    musclesWanted['inter'] = []
    musclesWanted['quads'] = ['vaslat_r', 'vasmed_r', 'vasint_r', 'recfem_r']
    musclesWanted['hams']  = ['bflh_r', 'bfsh_r', 'grac_r', 'sart_r', 'semimem_r', 'semiten_r']
    musclesWanted['gas']   = ['gaslat_r', 'gasmed_r']
    musclesWanted['tfl']   = ['tfl_r']
    musclesWanted['all'] = ['all']
    musclesWanted['reserve'] = ['reserve']
    musclesWanted['none'] = ['none']
    
    # TODO
    # plot the stuff
    # figure out subtracting the intersegmental from each of the others?
    

    pdb.set_trace()

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
                modelmass = get_model_total_mass(trialdir, 'mk12_rv1_dgf.osim')
                naturalstruct_combine[subject] = modelmass
                
                ### now what do we want to do at each of the trials
                # write a method that gives the knee contact forces, 
                # pass in model, and the muscles that you want
                # pdb.set_trace()
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

                ## okay now going to focus on the figures that I actually wanted 
                jrasrquads = getKneeContactributions(trialdir, musclesWanted['quads'], 'quads')
                jrasrhams = getKneeContactributions(trialdir, musclesWanted['hams'], 'hams')
                jrasrgas = getKneeContactributions(trialdir, musclesWanted['gas'], 'gas')
                jrasrtfl = getKneeContactributions(trialdir, musclesWanted['tfl'], 'tfl')
                jrasrinter = getKneeContactributions(trialdir, musclesWanted['inter'], 'inter')
                jrasrall = getKneeContactributions(trialdir, musclesWanted['all'], 'all')
                jrasrinterreserve = getKneeContactributions(trialdir, musclesWanted['reserve'], 'reserve')
                jrasrnone = getKneeContactributions(trialdir, musclesWanted['none'], 'none')
                
                # pdb.set_trace()
                
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
                jrasrquadsonly101 = -1*(np.interp(e_timesinterp, e_times, jrasrquadsonly)) / modelmass
                jrasrhamsonly101 = -1*(np.interp(e_timesinterp, e_times, jrasrhamsonly)) / modelmass
                jrasrgasonly101 = -1*(np.interp(e_timesinterp, e_times, jrasrgasonly)) / modelmass
                jrasrtflonly101 = -1*(np.interp(e_timesinterp, e_times, jrasrtflonly)) / modelmass
                jrasrinteronly101 = -1*(np.interp(e_timesinterp, e_times, jrasrinteronly)) / modelmass
                jrasrallonly101 = -1*(np.interp(e_timesinterp, e_times, jrasrall)) / modelmass
                jrasrreserveonly101 = -1*(np.interp(e_timesinterp, e_times, jrasrreserveonly)) / modelmass
                jrasrnoneonly101 = -1*(np.interp(e_timesinterp, e_times, jrasrnone)) / modelmass
                
                
                
                einterseg_combine[spot, :] = jrasrinteronly101[:-2]
                equads_combine[spot, :] = jrasrquadsonly101[:-2]
                ehams_combine[spot,:] = jrasrhamsonly101[:-2]
                egas_combine[spot,:] = jrasrgasonly101[:-2]
                etfl_combine[spot,:] = jrasrtflonly101[:-2]
                eall_combine[spot,:] = jrasrallonly101[:-2]
                ereserve_combine[spot,:] = jrasrreserveonly101[:-2]
                enone_combine[spot,:] = jrasrnoneonly101[:-2]
                
                # ## increase the spot - count of trials                
                spot += 1
                
                ## TODO: method for all the muscles ie. don't remove any
    # '''        
    # now the natural conditions
    # now the natural 
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
                modelmass = get_model_total_mass(trialdir, 'mk12_rv1_dgf.osim')
                naturalstruct_combine[subject] = modelmass
                
                ### now what do we want to do at each of the trials
                # write a method that gives the knee contact forces, 
                # pass in model, and the muscles that you want
                os.chdir(trialdir)                
                # pdb.set_trace()
                
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
                pdb.set_trace()
                ## okay now going to focus on the figures that I actually wanted 
                jrasrquads = getKneeContactributions(trialdir, musclesWanted['quads'], 'quads')
                jrasrhams = getKneeContactributions(trialdir, musclesWanted['hams'], 'hams')
                jrasrgas = getKneeContactributions(trialdir, musclesWanted['gas'], 'gas')
                jrasrtfl = getKneeContactributions(trialdir, musclesWanted['tfl'], 'tfl')
                jrasrinter = getKneeContactributions(trialdir, musclesWanted['inter'], 'inter')
                jrasrall = getKneeContactributions(trialdir, musclesWanted['all'], 'all')
                jrasrinterreserve = getKneeContactributions(trialdir, musclesWanted['reserve'], 'reserve')
                jrasrnone = getKneeContactributions(trialdir, musclesWanted['none'], 'none')
                
                # important: interreserve has the reserves removed, where inter includes them still. interreserves is the only one that removes the reserves...                
                
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
                jrasrquadsonly101 = -1*(np.interp(n_timesinterp, n_times, jrasrquadsonly)) / modelmass
                jrasrhamsonly101 = -1*(np.interp(n_timesinterp, n_times, jrasrhamsonly)) / modelmass
                jrasrgasonly101 = -1*(np.interp(n_timesinterp, n_times, jrasrgasonly)) / modelmass
                jrasrtflonly101 = -1*(np.interp(n_timesinterp, n_times, jrasrtflonly)) / modelmass
                jrasrinteronly101 = -1*(np.interp(n_timesinterp, n_times, jrasrinteronly)) / modelmass
                jrasrallonly101 = -1*(np.interp(n_timesinterp, n_times, jrasrall)) / modelmass
                jrasrreserveonly101 = -1*(np.interp(n_timesinterp, n_times, jrasrreserveonly)) / modelmass
                jrasrnoneonly101 = -1*(np.interp(n_timesinterp, n_times, jrasrnoneonly)) / modelmass
                
                # natural combine into the big structure
                ninterseg_combine[spot,:] = jrasrinteronly101[:-2]
                nquads_combine[spot,:] = jrasrquadsonly101[:-2]
                nhams_combine[spot,:] = jrasrhamsonly101[:-2]
                ngas_combine[spot,:] = jrasrgasonly101[:-2]
                ntfl_combine[spot,:] = jrasrtflonly101[:-2]
                nall_combine[spot,:] = jrasrallonly101[:-2]
                nreserve_combine[spot,:] = jrasrreserveonly101[:-2]
                nnone_combine[spot,:] = jrasrnoneonly101[:-2]
                
                # ## increase the spot - count of trials                
                spot += 1


    
    pdb.set_trace()
    
    # TODO: stacking contributions figure
    # plt.figure()
   
    
    ###########################################################################
    # figure: segmenting all the muscles between exo and nat 
    ## really nice figure for seeing what is going on, but likely not going to 
    ## be in the paper... 
    fig10, ax10 = plt.subplots(1,7, figsize=(18,3), dpi=300)
    # intersegmental forces average
    ax10[0].plot(n_timespercent101, np.mean(ninterseg_combine, 0), label='natural', color=ncolor)
    ax10[0].plot(e_timespercent101, np.mean(einterseg_combine, 0), label='exotendon', color=ecolor)
    ax10[0].set_xlabel('% Gait cycle')
    ax10[0].set_ylabel('Force (BW)')
    ax10[0].set_title('intersegmental')
    # ax10[0].legend()
    # tfl forces
    ax10[1].plot(n_timespercent101, np.mean(ntfl_combine, 0), label='natural', color=ncolor)
    ax10[1].plot(e_timespercent101, np.mean(etfl_combine, 0), label='exotendon', color=ecolor)
    ax10[1].set_xlabel('% Gait cycle')
    # ax10[1].set_ylabel('Force (BW)')
    # ax10[1].legend()
    ax10[1].set_title('tfl')

    # gastroc forces
    ax10[2].plot(n_timespercent101, np.mean(ngas_combine, 0), label='natural', color=ncolor)
    ax10[2].plot(e_timespercent101, np.mean(egas_combine, 0), label='exotendon', color=ecolor)
    ax10[2].set_xlabel('% Gait cycle')
    # ax10[2].set_ylabel('Force (BW)')
    # ax10[2].legend()
    ax10[2].set_title('gastroc')
    
    # hamstring forces
    ax10[3].plot(n_timespercent101, np.mean(nhams_combine, 0), label='natural', color=ncolor)
    ax10[3].plot(e_timespercent101, np.mean(ehams_combine, 0), label='exotendon', color=ecolor)
    ax10[3].set_xlabel('% Gait cycle')
    # ax10[3].set_ylabel('Force (BW)')
    # ax10[3].legend()
    ax10[3].set_title('hamstrings')

    # quads forces
    ax10[4].plot(n_timespercent101, np.mean(nquads_combine, 0), label='natural', color=ncolor)
    ax10[4].plot(e_timespercent101, np.mean(equads_combine, 0), label='exotendon', color=ecolor)
    ax10[4].set_xlabel('% Gait cycle')
    # ax10[4].set_ylabel('Force (BW)')
    # ax10[4].legend()
    ax10[4].set_title('quadriceps')

    # # reserve forces
    # ax10[5].plot(n_timespercent101, np.mean(nreserve_combine, 0), label='natural', color=ncolor)
    # ax10[5].plot(e_timespercent101, np.mean(ereserve_combine, 0), label='exotendon', color=ecolor)
    # ax10[5].set_xlabel('% Gait cycle')
    # # ax10[5].set_ylabel('Force (BW)')
    # # ax10[5].legend()
    # ax10[5].set_title('reserves')
    
    # added all forces
    ax10[5].plot(n_timespercent101, np.mean(nquads_combine,0) + np.mean(nhams_combine,0) + np.mean(ngas_combine,0) + np.mean(ntfl_combine,0) + np.mean(ninterseg_combine,0) + np.mean(nreserve_combine,0), label='natural', color=ncolor)
    ax10[5].plot(e_timespercent101, np.mean(equads_combine,0) + np.mean(ehams_combine,0) + np.mean(egas_combine,0) + np.mean(etfl_combine,0) + np.mean(einterseg_combine,0) + np.mean(ereserve_combine,0), label='exotendon', color=ecolor)
    ax10[5].set_xlabel('% Gait cycle')
    # ax10[6].set_ylabel('Force (BW)')
    # ax10[6].legend()
    ax10[5].set_title('Total vertical contact')

    # all forces from whole analysis
    ax10[6].plot(n_timespercent101, np.mean(nall_combine,0), label='natural', color=ncolor)
    ax10[6].plot(e_timespercent101, np.mean(eall_combine,0), label='exotendon', color=ecolor)
    ax10[6].set_xlabel('% Gait cycle')
    # ax10[7].set_ylabel('Force (BW)')
    ax10[6].legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10)
    ax10[6].set_title('Total vertical contact')

    fig10.tight_layout()
    
    ###########################################################################
    ### figure: differences between conditions for each and all
    #### this is a simplified look at the same as above. We can see nice stuff, 
    #### but likely not going to be in the paper. 
    fig11, ax11 = plt.subplots(1,7, figsize=(18,3), dpi=300)
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

    # added all forces
    ax11[5].plot(n_timespercent101, (np.mean(equads_combine,0) + np.mean(ehams_combine,0) + np.mean(egas_combine,0) + np.mean(etfl_combine,0) + np.mean(einterseg_combine,0)) - (np.mean(nquads_combine,0) + np.mean(nhams_combine,0) + np.mean(ngas_combine,0) + np.mean(ntfl_combine,0) + np.mean(ninterseg_combine,0)))
    ax11[5].set_xlabel('% Gait cycle')
    # ax11[5].set_ylabel('Force (BW)')
    # ax11[5].legend()
    ax11[5].set_title('Total vertical contact')

    # all forces from whole analysis
    ax11[6].plot(n_timespercent101, np.mean(eall_combine,0) - np.mean(nall_combine,0))
    ax11[6].set_xlabel('% Gait cycle')
    # ax11[6].set_ylabel('Force (BW)')
    ax11[6].legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10)
    ax11[6].set_title('Total vertical contact')

    fig11.tight_layout()
    
    
    ###########################################################################
    ### figure: changes in force segmented together on plot
    #### Possible paper figure for R3. 
    fig12 = plt.figure(figsize=(7,4), dpi=300)
    # intersegmental forces average
    plt.plot(n_timespercent101, np.mean(einterseg_combine, 0) - np.mean(ninterseg_combine, 0), label='intersegmental')
    # tfl forces
    plt.plot(n_timespercent101, np.mean(etfl_combine, 0) - np.mean(ntfl_combine, 0), label='tfl')
    # gastroc forces
    plt.plot(n_timespercent101, np.mean(egas_combine, 0) - np.mean(ngas_combine, 0), label='gastroc')
    # hamstring forces
    plt.plot(n_timespercent101, np.mean(ehams_combine, 0) - np.mean(nhams_combine, 0), label='hamstrings')
    # quads forces
    plt.plot(n_timespercent101, np.mean(equads_combine, 0) - np.mean(nquads_combine, 0), label='quadriceps')
    # reserve forces
    # plt.plot(n_timespercent101, np.mean(ereserve_combine, 0) - np.mean(nreserve_combine, 0), label='reserves')
    # added all forces
    # plt.plot(n_timespercent101, (np.mean(equads_combine,0) + np.mean(ehams_combine,0) + np.mean(egas_combine,0) + np.mean(etfl_combine,0) + np.mean(einterseg_combine,0) + np.mean(ereserve_combine,0)) - (np.mean(nquads_combine,0) + np.mean(nhams_combine,0) + np.mean(ngas_combine,0) + np.mean(ntfl_combine,0) + np.mean(ninterseg_combine,0) + np.mean(nreserve_combine,0)), label='Total')
    # all forces from whole analysis
    plt.plot(n_timespercent101, np.mean(eall_combine,0) - np.mean(nall_combine,0), label='Total vertical contact', linestyle='dashed')
    plt.xlabel('% Gait cycle')
    plt.ylabel('Force (BW)')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10)
    plt.title('Exotendon change in contact force')
    plt.tight_layout()

    ###########################################################################
    # TODO: figure out why the difference in total and all added together. 
    # okay so not in how I am adding/averaging. has to be something in how the analysis is done between them.... am I missing something??
    
    fig13, ax13 = plt.subplots(1,2, figsize=(10,8), dpi=300)
    # intersegmental forces average - natural
    ax13[0].plot(n_timespercent101, np.mean(ninterseg_combine, 0), label='natural_interseg')
    ax13[0].plot(n_timespercent101, np.mean(ninterseg_combine + ntfl_combine, 0), label='natural_interseg + tfl')
    ax13[0].plot(n_timespercent101, np.mean(ninterseg_combine+ntfl_combine+nhams_combine, 0), label='natural_interseg+tfl+hams')
    ax13[0].plot(n_timespercent101, np.mean(ninterseg_combine+ntfl_combine+nhams_combine+ngas_combine, 0), label='natural_interseg+tfl+hams+gas')
    ax13[0].plot(n_timespercent101, np.mean(ninterseg_combine+ntfl_combine+nhams_combine+ngas_combine+nquads_combine, 0), label='natural_interseg+tfl+hams+gas+quads')
    ax13[0].plot(n_timespercent101, np.mean(ninterseg_combine+ntfl_combine+nhams_combine+ngas_combine+nquads_combine+nreserve_combine, 0), label='natural_interseg+tfl+hams+gas+quads+reserve')
    ax13[0].plot(n_timespercent101, np.mean(nall_combine, 0), label='nat_all', linestyle='dotted')
    ax13[0].set_xlabel('% Gait cycle')
    ax13[0].set_ylabel('Force (BW)')
    ax13[0].set_title('natural')
    ax13[0].legend()
    # intersegmental forces average - exotendon
    ax13[1].plot(e_timespercent101, np.mean(einterseg_combine, 0), label='exo_interseg')
    ax13[1].plot(e_timespercent101, np.mean(einterseg_combine+etfl_combine, 0), label='exo_interseg + tfl')
    ax13[1].plot(e_timespercent101, np.mean(einterseg_combine+etfl_combine+ehams_combine, 0), label='exo_interseg+tfl+hams')
    ax13[1].plot(e_timespercent101, np.mean(einterseg_combine+etfl_combine+ehams_combine+egas_combine, 0), label='exo_interseg+tfl+hams+gas')
    ax13[1].plot(e_timespercent101, np.mean(einterseg_combine+etfl_combine+ehams_combine+egas_combine+equads_combine, 0), label='exo_interseg+tfl+hams+gas+quads')
    ax13[1].plot(e_timespercent101, np.mean(einterseg_combine+etfl_combine+ehams_combine+egas_combine+equads_combine+ereserve_combine, 0), label='exo_interseg+tfl+hams+gas+quads+reserve')
    ax13[1].plot(e_timespercent101, np.mean(eall_combine, 0), label='exo_all', linestyle='dotted')
    ax13[1].set_xlabel('% Gait cycle')
    ax13[1].set_ylabel('Force (BW)')
    ax13[1].set_title('exotendon')
    ax13[1].legend()
    fig13.tight_layout()
    
    
    ###########################################################################
    # figures like ^ but individuals
    fig131, ax131 = plt.subplots(1,2, figsize=(10,8), dpi=300)
    # intersegmental forces average - natural welk002
    ax131[0].plot(n_timespercent101, np.mean(ninterseg_combine[0:4,:], 0), label='natural_interseg')
    ax131[0].plot(n_timespercent101, np.mean(ninterseg_combine[0:4,:] + ntfl_combine[0:4,:], 0), label='natural_interseg + tfl')
    ax131[0].plot(n_timespercent101, np.mean(ninterseg_combine[0:4,:]+ntfl_combine[0:4,:]+nhams_combine[0:4,:], 0), label='natural_interseg+tfl+hams')
    ax131[0].plot(n_timespercent101, np.mean(ninterseg_combine[0:4,:]+ntfl_combine[0:4,:]+nhams_combine[0:4,:]+ngas_combine[0:4,:], 0), label='natural_interseg+tfl+hams+gas')
    ax131[0].plot(n_timespercent101, np.mean(ninterseg_combine[0:4,:]+ntfl_combine[0:4,:]+nhams_combine[0:4,:]+ngas_combine[0:4,:]+nquads_combine[0:4,:], 0), label='natural_interseg+tfl+hams+gas+quads')
    ax131[0].plot(n_timespercent101, np.mean(ninterseg_combine[0:4,:]+ntfl_combine[0:4,:]+nhams_combine[0:4,:]+ngas_combine[0:4,:]+nquads_combine[0:4,:]+nreserve_combine[0:4,:], 0), label='natural_interseg+tfl+hams+gas+quads+reserve')
    ax131[0].plot(n_timespercent101, np.mean(nall_combine[0:4,:], 0), label='nat_all', linestyle='dotted')
    ax131[0].set_xlabel('% Gait cycle')
    ax131[0].set_ylabel('Force (BW)')
    ax131[0].set_title('welk002 natural')
    ax131[0].legend()
    # intersegmental forces average - exotendon
    ax131[1].plot(e_timespercent101, np.mean(einterseg_combine[0:4,:], 0), label='exo_interseg')
    ax131[1].plot(e_timespercent101, np.mean(einterseg_combine[0:4,:]+etfl_combine[0:4,:], 0), label='exo_interseg + tfl')
    ax131[1].plot(e_timespercent101, np.mean(einterseg_combine[0:4,:]+etfl_combine[0:4,:]+ehams_combine[0:4,:], 0), label='exo_interseg+tfl+hams')
    ax131[1].plot(e_timespercent101, np.mean(einterseg_combine[0:4,:]+etfl_combine[0:4,:]+ehams_combine[0:4,:]+egas_combine[0:4,:], 0), label='exo_interseg+tfl+hams+gas')
    ax131[1].plot(e_timespercent101, np.mean(einterseg_combine[0:4,:]+etfl_combine[0:4,:]+ehams_combine[0:4,:]+egas_combine[0:4,:]+equads_combine[0:4,:], 0), label='exo_interseg+tfl+hams+gas+quads')
    ax131[1].plot(e_timespercent101, np.mean(einterseg_combine[0:4,:]+etfl_combine[0:4,:]+ehams_combine[0:4,:]+egas_combine[0:4,:]+equads_combine[0:4,:]+ereserve_combine[0:4,:], 0), label='exo_interseg+tfl+hams+gas+quads+reserve')
    ax131[1].plot(e_timespercent101, np.mean(eall_combine[0:4,:], 0), label='exo_all', linestyle='dotted')
    ax131[1].set_xlabel('% Gait cycle')
    ax131[1].set_ylabel('Force (BW)')
    ax131[1].set_title('welk002 exotendon')
    ax131[1].legend()
    fig131.tight_layout()
    
    # welk003
    fig132, ax132 = plt.subplots(1,2, figsize=(10,8), dpi=300)
    # intersegmental forces average - natural welk003
    ax132[0].plot(n_timespercent101, np.mean(ninterseg_combine[4:8,:], 0), label='natural_interseg')
    ax132[0].plot(n_timespercent101, np.mean(ninterseg_combine[4:8,:] + ntfl_combine[4:8,:], 0), label='natural_interseg + tfl')
    ax132[0].plot(n_timespercent101, np.mean(ninterseg_combine[4:8,:]+ntfl_combine[4:8,:]+nhams_combine[4:8,:], 0), label='natural_interseg+tfl+hams')
    ax132[0].plot(n_timespercent101, np.mean(ninterseg_combine[4:8,:]+ntfl_combine[4:8,:]+nhams_combine[4:8,:]+ngas_combine[4:8,:], 0), label='natural_interseg+tfl+hams+gas')
    ax132[0].plot(n_timespercent101, np.mean(ninterseg_combine[4:8,:]+ntfl_combine[4:8,:]+nhams_combine[4:8,:]+ngas_combine[4:8,:]+nquads_combine[4:8,:], 0), label='natural_interseg+tfl+hams+gas+quads')
    ax132[0].plot(n_timespercent101, np.mean(ninterseg_combine[4:8,:]+ntfl_combine[4:8,:]+nhams_combine[4:8,:]+ngas_combine[4:8,:]+nquads_combine[4:8,:]+nreserve_combine[4:8,:], 0), label='natural_interseg+tfl+hams+gas+quads+reserve')
    ax132[0].plot(n_timespercent101, np.mean(nall_combine[4:8,:], 0), label='nat_all', linestyle='dotted')
    ax132[0].set_xlabel('% Gait cycle')
    ax132[0].set_ylabel('Force (BW)')
    ax132[0].set_title('welk003 natural')
    ax132[0].legend()
    # intersegmental forces average - exotendon
    ax132[1].plot(e_timespercent101, np.mean(einterseg_combine[4:8,:], 0), label='exo_interseg')
    ax132[1].plot(e_timespercent101, np.mean(einterseg_combine[4:8,:]+etfl_combine[4:8,:], 0), label='exo_interseg + tfl')
    ax132[1].plot(e_timespercent101, np.mean(einterseg_combine[4:8,:]+etfl_combine[4:8,:]+ehams_combine[4:8,:], 0), label='exo_interseg+tfl+hams')
    ax132[1].plot(e_timespercent101, np.mean(einterseg_combine[4:8,:]+etfl_combine[4:8,:]+ehams_combine[4:8,:]+egas_combine[4:8,:], 0), label='exo_interseg+tfl+hams+gas')
    ax132[1].plot(e_timespercent101, np.mean(einterseg_combine[4:8,:]+etfl_combine[4:8,:]+ehams_combine[4:8,:]+egas_combine[4:8,:]+equads_combine[4:8,:], 0), label='exo_interseg+tfl+hams+gas+quads')
    ax132[1].plot(e_timespercent101, np.mean(einterseg_combine[4:8,:]+etfl_combine[4:8,:]+ehams_combine[4:8,:]+egas_combine[4:8,:]+equads_combine[4:8,:]+ereserve_combine[4:8,:], 0), label='exo_interseg+tfl+hams+gas+quads+reserve')
    ax132[1].plot(e_timespercent101, np.mean(eall_combine[4:8,:], 0), label='exo_all', linestyle='dotted')
    ax132[1].set_xlabel('% Gait cycle')
    ax132[1].set_ylabel('Force (BW)')
    ax132[1].set_title('welk003 exotendon')
    ax132[1].legend()
    fig132.tight_layout()
    # welk005
    fig133, ax133 = plt.subplots(1,2, figsize=(10,8), dpi=300)
    # intersegmental forces average - natural welk005
    ax133[0].plot(n_timespercent101, np.mean(ninterseg_combine[8:12,:], 0), label='natural_interseg')
    ax133[0].plot(n_timespercent101, np.mean(ninterseg_combine[8:12,:] + ntfl_combine[8:12,:], 0), label='natural_interseg + tfl')
    ax133[0].plot(n_timespercent101, np.mean(ninterseg_combine[8:12,:]+ntfl_combine[8:12,:]+nhams_combine[8:12,:], 0), label='natural_interseg+tfl+hams')
    ax133[0].plot(n_timespercent101, np.mean(ninterseg_combine[8:12,:]+ntfl_combine[8:12,:]+nhams_combine[8:12,:]+ngas_combine[8:12,:], 0), label='natural_interseg+tfl+hams+gas')
    ax133[0].plot(n_timespercent101, np.mean(ninterseg_combine[8:12,:]+ntfl_combine[8:12,:]+nhams_combine[8:12,:]+ngas_combine[8:12,:]+nquads_combine[8:12,:], 0), label='natural_interseg+tfl+hams+gas+quads')
    ax133[0].plot(n_timespercent101, np.mean(ninterseg_combine[8:12,:]+ntfl_combine[8:12,:]+nhams_combine[8:12,:]+ngas_combine[8:12,:]+nquads_combine[8:12,:]+nreserve_combine[8:12,:], 0), label='natural_interseg+tfl+hams+gas+quads+reserve')
    ax133[0].plot(n_timespercent101, np.mean(nall_combine[8:12,:], 0), label='nat_all', linestyle='dotted')
    ax133[0].set_xlabel('% Gait cycle')
    ax133[0].set_ylabel('Force (BW)')
    ax133[0].set_title('welk005 natural')
    ax133[0].legend()
    # intersegmental forces average - exotendon
    ax133[1].plot(e_timespercent101, np.mean(einterseg_combine[8:12,:], 0), label='exo_interseg')
    ax133[1].plot(e_timespercent101, np.mean(einterseg_combine[8:12,:]+etfl_combine[8:12,:], 0), label='exo_interseg + tfl')
    ax133[1].plot(e_timespercent101, np.mean(einterseg_combine[8:12,:]+etfl_combine[8:12,:]+ehams_combine[8:12,:], 0), label='exo_interseg+tfl+hams')
    ax133[1].plot(e_timespercent101, np.mean(einterseg_combine[8:12,:]+etfl_combine[8:12,:]+ehams_combine[8:12,:]+egas_combine[8:12,:], 0), label='exo_interseg+tfl+hams+gas')
    ax133[1].plot(e_timespercent101, np.mean(einterseg_combine[8:12,:]+etfl_combine[8:12,:]+ehams_combine[8:12,:]+egas_combine[8:12,:]+equads_combine[8:12,:], 0), label='exo_interseg+tfl+hams+gas+quads')
    ax133[1].plot(e_timespercent101, np.mean(einterseg_combine[8:12,:]+etfl_combine[8:12,:]+ehams_combine[8:12,:]+egas_combine[8:12,:]+equads_combine[8:12,:]+ereserve_combine[8:12,:], 0), label='exo_interseg+tfl+hams+gas+quads+reserve')
    ax133[1].plot(e_timespercent101, np.mean(eall_combine[8:12,:], 0), label='exo_all', linestyle='dotted')
    ax133[1].set_xlabel('% Gait cycle')
    ax133[1].set_ylabel('Force (BW)')
    ax133[1].set_title('welk005 exotendon')
    ax133[1].legend()
    fig133.tight_layout()
    # welk008
    fig134, ax134 = plt.subplots(1,2, figsize=(10,8), dpi=300)
    # intersegmental forces average - natural welk008
    ax134[0].plot(n_timespercent101, np.mean(ninterseg_combine[12:16,:], 0), label='natural_interseg')
    ax134[0].plot(n_timespercent101, np.mean(ninterseg_combine[12:16,:] + ntfl_combine[12:16,:], 0), label='natural_interseg + tfl')
    ax134[0].plot(n_timespercent101, np.mean(ninterseg_combine[12:16,:]+ntfl_combine[12:16,:]+nhams_combine[12:16,:], 0), label='natural_interseg+tfl+hams')
    ax134[0].plot(n_timespercent101, np.mean(ninterseg_combine[12:16,:]+ntfl_combine[12:16,:]+nhams_combine[12:16,:]+ngas_combine[12:16,:], 0), label='natural_interseg+tfl+hams+gas')
    ax134[0].plot(n_timespercent101, np.mean(ninterseg_combine[12:16,:]+ntfl_combine[12:16,:]+nhams_combine[12:16,:]+ngas_combine[12:16,:]+nquads_combine[12:16,:], 0), label='natural_interseg+tfl+hams+gas+quads')
    ax134[0].plot(n_timespercent101, np.mean(ninterseg_combine[12:16,:]+ntfl_combine[12:16,:]+nhams_combine[12:16,:]+ngas_combine[12:16,:]+nquads_combine[12:16,:]+nreserve_combine[12:16,:], 0), label='natural_interseg+tfl+hams+gas+quads+reserve')
    ax134[0].plot(n_timespercent101, np.mean(nall_combine[12:16,:], 0), label='nat_all', linestyle='dotted')
    ax134[0].set_xlabel('% Gait cycle')
    ax134[0].set_ylabel('Force (BW)')
    ax134[0].set_title('welk008 natural')
    ax134[0].legend()
    # intersegmental forces average - exotendon
    ax134[1].plot(e_timespercent101, np.mean(einterseg_combine[12:16,:], 0), label='exo_interseg')
    ax134[1].plot(e_timespercent101, np.mean(einterseg_combine[12:16,:]+etfl_combine[12:16,:], 0), label='exo_interseg + tfl')
    ax134[1].plot(e_timespercent101, np.mean(einterseg_combine[12:16,:]+etfl_combine[12:16,:]+ehams_combine[12:16,:], 0), label='exo_interseg+tfl+hams')
    ax134[1].plot(e_timespercent101, np.mean(einterseg_combine[12:16,:]+etfl_combine[12:16,:]+ehams_combine[12:16,:]+egas_combine[12:16,:], 0), label='exo_interseg+tfl+hams+gas')
    ax134[1].plot(e_timespercent101, np.mean(einterseg_combine[12:16,:]+etfl_combine[12:16,:]+ehams_combine[12:16,:]+egas_combine[12:16,:]+equads_combine[12:16,:], 0), label='exo_interseg+tfl+hams+gas+quads')
    ax134[1].plot(e_timespercent101, np.mean(einterseg_combine[12:16,:]+etfl_combine[12:16,:]+ehams_combine[12:16,:]+egas_combine[12:16,:]+equads_combine[12:16,:]+ereserve_combine[12:16,:], 0), label='exo_interseg+tfl+hams+gas+quads+reserve')
    ax134[1].plot(e_timespercent101, np.mean(eall_combine[12:16,:], 0), label='exo_all', linestyle='dotted')
    ax134[1].set_xlabel('% Gait cycle')
    ax134[1].set_ylabel('Force (BW)')
    ax134[1].set_title('welk008 exotendon')
    ax134[1].legend()
    fig134.tight_layout()
    # welk009
    fig135, ax135 = plt.subplots(1,2, figsize=(10,8), dpi=300)
    # intersegmental forces average - natural welk009
    ax135[0].plot(n_timespercent101, np.mean(ninterseg_combine[16:20,:], 0), label='natural_interseg')
    ax135[0].plot(n_timespercent101, np.mean(ninterseg_combine[16:20,:] + ntfl_combine[16:20,:], 0), label='natural_interseg + tfl')
    ax135[0].plot(n_timespercent101, np.mean(ninterseg_combine[16:20,:]+ntfl_combine[16:20,:]+nhams_combine[16:20,:], 0), label='natural_interseg+tfl+hams')
    ax135[0].plot(n_timespercent101, np.mean(ninterseg_combine[16:20,:]+ntfl_combine[16:20,:]+nhams_combine[16:20,:]+ngas_combine[16:20,:], 0), label='natural_interseg+tfl+hams+gas')
    ax135[0].plot(n_timespercent101, np.mean(ninterseg_combine[16:20,:]+ntfl_combine[16:20,:]+nhams_combine[16:20,:]+ngas_combine[16:20,:]+nquads_combine[16:20,:], 0), label='natural_interseg+tfl+hams+gas+quads')
    ax135[0].plot(n_timespercent101, np.mean(ninterseg_combine[16:20,:]+ntfl_combine[16:20,:]+nhams_combine[16:20,:]+ngas_combine[16:20,:]+nquads_combine[16:20,:]+nreserve_combine[16:20,:], 0), label='natural_interseg+tfl+hams+gas+quads+reserve')
    ax135[0].plot(n_timespercent101, np.mean(nall_combine[16:20,:], 0), label='nat_all', linestyle='dotted')
    ax135[0].set_xlabel('% Gait cycle')
    ax135[0].set_ylabel('Force (BW)')
    ax135[0].set_title('welk009 natural')
    ax135[0].legend()
    # intersegmental forces average - exotendon
    ax135[1].plot(e_timespercent101, np.mean(einterseg_combine[16:20,:], 0), label='exo_interseg')
    ax135[1].plot(e_timespercent101, np.mean(einterseg_combine[16:20,:]+etfl_combine[16:20,:], 0), label='exo_interseg + tfl')
    ax135[1].plot(e_timespercent101, np.mean(einterseg_combine[16:20,:]+etfl_combine[16:20,:]+ehams_combine[16:20,:], 0), label='exo_interseg+tfl+hams')
    ax135[1].plot(e_timespercent101, np.mean(einterseg_combine[16:20,:]+etfl_combine[16:20,:]+ehams_combine[16:20,:]+egas_combine[16:20,:], 0), label='exo_interseg+tfl+hams+gas')
    ax135[1].plot(e_timespercent101, np.mean(einterseg_combine[16:20,:]+etfl_combine[16:20,:]+ehams_combine[16:20,:]+egas_combine[16:20,:]+equads_combine[16:20,:], 0), label='exo_interseg+tfl+hams+gas+quads')
    ax135[1].plot(e_timespercent101, np.mean(einterseg_combine[16:20,:]+etfl_combine[16:20,:]+ehams_combine[16:20,:]+egas_combine[16:20,:]+equads_combine[16:20,:]+ereserve_combine[16:20,:], 0), label='exo_interseg+tfl+hams+gas+quads+reserve')
    ax135[1].plot(e_timespercent101, np.mean(eall_combine[16:20,:], 0), label='exo_all', linestyle='dotted')
    ax135[1].set_xlabel('% Gait cycle')
    ax135[1].set_ylabel('Force (BW)')
    ax135[1].set_title('welk009 exotendon')
    ax135[1].legend()
    fig135.tight_layout()
    # welk010
    fig136, ax136 = plt.subplots(1,2, figsize=(10,8), dpi=300)
    # intersegmental forces average - natural welk010
    ax136[0].plot(n_timespercent101, np.mean(ninterseg_combine[20:24,:], 0), label='natural_interseg')
    ax136[0].plot(n_timespercent101, np.mean(ninterseg_combine[20:24,:] + ntfl_combine[20:24,:], 0), label='natural_interseg + tfl')
    ax136[0].plot(n_timespercent101, np.mean(ninterseg_combine[20:24,:]+ntfl_combine[20:24,:]+nhams_combine[20:24,:], 0), label='natural_interseg+tfl+hams')
    ax136[0].plot(n_timespercent101, np.mean(ninterseg_combine[20:24,:]+ntfl_combine[20:24,:]+nhams_combine[20:24,:]+ngas_combine[20:24,:], 0), label='natural_interseg+tfl+hams+gas')
    ax136[0].plot(n_timespercent101, np.mean(ninterseg_combine[20:24,:]+ntfl_combine[20:24,:]+nhams_combine[20:24,:]+ngas_combine[20:24,:]+nquads_combine[20:24,:], 0), label='natural_interseg+tfl+hams+gas+quads')
    ax136[0].plot(n_timespercent101, np.mean(ninterseg_combine[20:24,:]+ntfl_combine[20:24,:]+nhams_combine[20:24,:]+ngas_combine[20:24,:]+nquads_combine[20:24,:]+nreserve_combine[20:24,:], 0), label='natural_interseg+tfl+hams+gas+quads+reserve')
    ax136[0].plot(n_timespercent101, np.mean(nall_combine[20:24,:], 0), label='nat_all', linestyle='dotted')
    ax136[0].set_xlabel('% Gait cycle')
    ax136[0].set_ylabel('Force (BW)')
    ax136[0].set_title('welk010 natural')
    ax136[0].legend()
    # intersegmental forces average - exotendon
    ax136[1].plot(e_timespercent101, np.mean(einterseg_combine[20:24,:], 0), label='exo_interseg')
    ax136[1].plot(e_timespercent101, np.mean(einterseg_combine[20:24,:]+etfl_combine[20:24,:], 0), label='exo_interseg + tfl')
    ax136[1].plot(e_timespercent101, np.mean(einterseg_combine[20:24,:]+etfl_combine[20:24,:]+ehams_combine[20:24,:], 0), label='exo_interseg+tfl+hams')
    ax136[1].plot(e_timespercent101, np.mean(einterseg_combine[20:24,:]+etfl_combine[20:24,:]+ehams_combine[20:24,:]+egas_combine[20:24,:], 0), label='exo_interseg+tfl+hams+gas')
    ax136[1].plot(e_timespercent101, np.mean(einterseg_combine[20:24,:]+etfl_combine[20:24,:]+ehams_combine[20:24,:]+egas_combine[20:24,:]+equads_combine[20:24,:], 0), label='exo_interseg+tfl+hams+gas+quads')
    ax136[1].plot(e_timespercent101, np.mean(einterseg_combine[20:24,:]+etfl_combine[20:24,:]+ehams_combine[20:24,:]+egas_combine[20:24,:]+equads_combine[20:24,:]+ereserve_combine[20:24,:], 0), label='exo_interseg+tfl+hams+gas+quads+reserve')
    ax136[1].plot(e_timespercent101, np.mean(eall_combine[20:24,:], 0), label='exo_all', linestyle='dotted')
    ax136[1].set_xlabel('% Gait cycle')
    ax136[1].set_ylabel('Force (BW)')
    ax136[1].set_title('welk010 exotendon')
    ax136[1].legend()
    fig136.tight_layout()
    # welk013
    fig137, ax137 = plt.subplots(1,2, figsize=(10,8), dpi=300)
    # intersegmental forces average - natural welk013
    ax137[0].plot(n_timespercent101, np.mean(ninterseg_combine[24:28,:], 0), label='natural_interseg')
    ax137[0].plot(n_timespercent101, np.mean(ninterseg_combine[24:28,:] + ntfl_combine[24:28,:], 0), label='natural_interseg + tfl')
    ax137[0].plot(n_timespercent101, np.mean(ninterseg_combine[24:28,:]+ntfl_combine[24:28,:]+nhams_combine[24:28,:], 0), label='natural_interseg+tfl+hams')
    ax137[0].plot(n_timespercent101, np.mean(ninterseg_combine[24:28,:]+ntfl_combine[24:28,:]+nhams_combine[24:28,:]+ngas_combine[24:28,:], 0), label='natural_interseg+tfl+hams+gas')
    ax137[0].plot(n_timespercent101, np.mean(ninterseg_combine[24:28,:]+ntfl_combine[24:28,:]+nhams_combine[24:28,:]+ngas_combine[24:28,:]+nquads_combine[24:28,:], 0), label='natural_interseg+tfl+hams+gas+quads')
    ax137[0].plot(n_timespercent101, np.mean(ninterseg_combine[24:28,:]+ntfl_combine[24:28,:]+nhams_combine[24:28,:]+ngas_combine[24:28,:]+nquads_combine[24:28,:]+nreserve_combine[24:28,:], 0), label='natural_interseg+tfl+hams+gas+quads+reserve')
    ax137[0].plot(n_timespercent101, np.mean(nall_combine[24:28,:], 0), label='nat_all', linestyle='dotted')
    ax137[0].set_xlabel('% Gait cycle')
    ax137[0].set_ylabel('Force (BW)')
    ax137[0].set_title('welk013 natural')
    ax137[0].legend()
    # intersegmental forces average - exotendon
    ax137[1].plot(e_timespercent101, np.mean(einterseg_combine[24:28,:], 0), label='exo_interseg')
    ax137[1].plot(e_timespercent101, np.mean(einterseg_combine[24:28,:]+etfl_combine[24:28,:], 0), label='exo_interseg + tfl')
    ax137[1].plot(e_timespercent101, np.mean(einterseg_combine[24:28,:]+etfl_combine[24:28,:]+ehams_combine[24:28,:], 0), label='exo_interseg+tfl+hams')
    ax137[1].plot(e_timespercent101, np.mean(einterseg_combine[24:28,:]+etfl_combine[24:28,:]+ehams_combine[24:28,:]+egas_combine[24:28,:], 0), label='exo_interseg+tfl+hams+gas')
    ax137[1].plot(e_timespercent101, np.mean(einterseg_combine[24:28,:]+etfl_combine[24:28,:]+ehams_combine[24:28,:]+egas_combine[24:28,:]+equads_combine[24:28,:], 0), label='exo_interseg+tfl+hams+gas+quads')
    ax137[1].plot(e_timespercent101, np.mean(einterseg_combine[24:28,:]+etfl_combine[24:28,:]+ehams_combine[24:28,:]+egas_combine[24:28,:]+equads_combine[24:28,:]+ereserve_combine[24:28,:], 0), label='exo_interseg+tfl+hams+gas+quads+reserve')
    ax137[1].plot(e_timespercent101, np.mean(eall_combine[24:28,:], 0), label='exo_all', linestyle='dotted')
    ax137[1].set_xlabel('% Gait cycle')
    ax137[1].set_ylabel('Force (BW)')
    ax137[1].set_title('welk013 exotendon')
    ax137[1].legend()
    fig137.tight_layout()
    
    
    
    pdb.set_trace()
    ###########################################################################
    ### Figure: total pop average for right leg between nat and exo
    #### Figure in the paper R1... 
    plt.figure(dpi=300)
    plt.fill_between(n_timespercent101, np.mean(nall_combine,0)-np.std(nall_combine,0), np.mean(nall_combine,0)+np.std(nall_combine,0), color=ncolorlight)
    plt.fill_between(e_timespercent101, np.mean(eall_combine,0)-np.std(eall_combine,0), np.mean(eall_combine,0)+np.std(eall_combine,0), color=ecolorlight)
    plt.plot(n_timespercent101, np.mean(nall_combine,0), color=ncolor, label='natural')
    plt.plot(e_timespercent101, np.mean(eall_combine,0), color=ecolor, label='exotendon')
    plt.xlabel('% Gait cycle')
    plt.ylabel('Force (BW)')
    plt.title('Total vertical contact force')
    plt.legend(loc='upper right')
    plt.tight_layout()
    
    # here is the stats for R1 - differences in the peak vertical JCF
    mean_nall_combine = np.mean(nall_combine,0)
    std_nall_combine = np.std(nall_combine,0)
    mean_eall_combine = np.mean(eall_combine,0)
    std_eall_combine = np.std(eall_combine,0)
    
    # get the peaks first and then do the stats on them...
    peaks_nall_combine = np.max(nall_combine,1)
    idx_peaks_nall_combine = nall_combine.argmax(1)
    peaks_eall_combine = np.max(eall_combine,1)
    idx_peaks_eall_combine = eall_combine.argmax(1)
    
    mean_peaks_nall_combine = np.mean(peaks_nall_combine)
    mean_peaks_eall_combine = np.mean(peaks_eall_combine)
    std_peaks_nall_combine = np.std(peaks_nall_combine)
    std_peaks_eall_combine = np.std(peaks_eall_combine)
    
    # differences in peaks
    diff_peaks_all_combine = peaks_nall_combine - peaks_eall_combine
    mean_diff_peaks_all_combine = np.mean(peaks_nall_combine - peaks_eall_combine)
    std_diff_peaks_all_combine = np.std(peaks_nall_combine - peaks_eall_combine)
    
    # shapiro test for normal on the peaks differences
    res = scipy.stats.shapiro(diff_peaks_all_combine)
    
    # t - test for peaks differences is in the excel sheet...
    
    # percent change in peak
    perc_diff_peaks_all_combine = (peaks_eall_combine - peaks_nall_combine) / peaks_nall_combine * 100
    mean_perc_diff_peaks_all_combine = np.mean(perc_diff_peaks_all_combine)
    std_perc_diff_peaks_all_combine = np.std(perc_diff_peaks_all_combine)
    
    ###########################################################################
    ### Figure: total pop average for right leg for nat - segmented shaded
    plt.figure(dpi=300)
    plt.plot(n_timespercent101, np.mean(ninterseg_combine,0), label='intersegmental', color=ecolor)
    plt.fill_between(n_timespercent101, 
                     np.mean(ninterseg_combine,0), color=ecolorlight)
    
    plt.plot(n_timespercent101, np.mean(ntfl_combine,0) + np.mean(ninterseg_combine,0), 
             label='inter + tfl', color=ncolor3)
    plt.fill_between(n_timespercent101, np.mean(ninterseg_combine,0), 
                     np.mean(ntfl_combine,0) + np.mean(ninterseg_combine,0), color=ncolor2)
    
    plt.plot(n_timespercent101, np.mean(ngas_combine,0) + np.mean(ntfl_combine,0) + np.mean(ninterseg_combine,0), 
             label='inter + tfl + gas', color=ncolor5)
    plt.fill_between(n_timespercent101, np.mean(ntfl_combine,0) + np.mean(ninterseg_combine,0), 
                     np.mean(ngas_combine,0) + np.mean(ntfl_combine,0) + np.mean(ninterseg_combine,0), color=ncolor4)
    
    plt.plot(n_timespercent101, np.mean(nhams_combine,0) + np.mean(ngas_combine,0) + np.mean(ntfl_combine,0) + np.mean(ninterseg_combine,0), 
             label='inter + tfl + gas + hams', color=ncolor7)
    plt.fill_between(n_timespercent101, np.mean(ngas_combine,0) + np.mean(ntfl_combine,0) + np.mean(ninterseg_combine,0), 
                     np.mean(nhams_combine,0) + np.mean(ngas_combine,0) + np.mean(ntfl_combine,0) + np.mean(ninterseg_combine,0),
                     color=ncolor6)
    
    plt.plot(n_timespercent101, np.mean(nquads_combine,0) + np.mean(nhams_combine,0) + np.mean(ngas_combine,0) + np.mean(ntfl_combine,0) + np.mean(ninterseg_combine,0), 
             label='inter + tfl + gas + hams + quads', color=ncolor9)
    plt.fill_between(n_timespercent101, np.mean(nhams_combine,0) + np.mean(ngas_combine,0) + np.mean(ntfl_combine,0) + np.mean(ninterseg_combine,0), 
                     np.mean(nquads_combine,0) + np.mean(nhams_combine,0) + np.mean(ngas_combine,0) + np.mean(ntfl_combine,0) + np.mean(ninterseg_combine,0), 
                     color=ncolor8)
    
    plt.legend()
    plt.ylabel('Force (BW)')
    plt.xlabel('% Gait cycle')
    plt.title('Vertical knee contact force - Natural')    
    
    
    ###########################################################################
    ### Figure: total pop average for right leg for exo - segmented shaded
    plt.figure(dpi=300)
    plt.plot(e_timespercent101, np.mean(einterseg_combine,0), label='intersegmental', color=ecolor)
    plt.fill_between(e_timespercent101, 
                     np.mean(einterseg_combine,0), color=ecolorlight)
    
    plt.plot(e_timespercent101, np.mean(etfl_combine,0) + np.mean(einterseg_combine,0), 
             label='inter + tfl', color=ncolor3)
    plt.fill_between(e_timespercent101, np.mean(einterseg_combine,0), 
                     np.mean(etfl_combine,0) + np.mean(einterseg_combine,0), color=ncolor2)
    
    plt.plot(e_timespercent101, np.mean(egas_combine,0) + np.mean(etfl_combine,0) + np.mean(einterseg_combine,0), 
             label='inter + tfl + gas', color=ncolor5)
    plt.fill_between(e_timespercent101, np.mean(etfl_combine,0) + np.mean(einterseg_combine,0), 
                     np.mean(egas_combine,0) + np.mean(etfl_combine,0) + np.mean(einterseg_combine,0), color=ncolor4)
    
    plt.plot(e_timespercent101, np.mean(ehams_combine,0) + np.mean(egas_combine,0) + np.mean(etfl_combine,0) + np.mean(einterseg_combine,0), 
             label='inter + tfl + gas + hams', color=ncolor7)
    plt.fill_between(e_timespercent101, np.mean(egas_combine,0) + np.mean(etfl_combine,0) + np.mean(einterseg_combine,0), 
                     np.mean(ehams_combine,0) + np.mean(egas_combine,0) + np.mean(etfl_combine,0) + np.mean(einterseg_combine,0),
                     color=ncolor6)
    
    plt.plot(e_timespercent101, np.mean(equads_combine,0) + np.mean(ehams_combine,0) + np.mean(egas_combine,0) + np.mean(etfl_combine,0) + np.mean(einterseg_combine,0), 
             label='inter + tfl + gas + hams + quads', color=ncolor9)
    plt.fill_between(e_timespercent101, np.mean(ehams_combine,0) + np.mean(egas_combine,0) + np.mean(etfl_combine,0) + np.mean(einterseg_combine,0), 
                     np.mean(equads_combine,0) + np.mean(ehams_combine,0) + np.mean(egas_combine,0) + np.mean(etfl_combine,0) + np.mean(einterseg_combine,0), 
                     color=ncolor8)
    
    plt.legend()
    plt.ylabel('Force (BW)')
    plt.xlabel('% Gait cycle')
    plt.title('Vertical knee contact force - Exotendon')    
    
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
    
    pdb.set_trace()
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
                
           
                
                
                
                
                
                
                
                
                
                
                