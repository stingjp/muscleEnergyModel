'''
Jon Stingel
2024/11/26
Script for batching through all the python versions of simulations. 
'''

import os
from shutil import copyfile
from shutil import copy
from shutil import copytree
os.add_dll_directory('C:/Users/jonstingel/opensim-core-4.5.1-2024-08-23-cf3ef35/bin')
import opensim as osim
import OsimUtilityfunctions as ouf
import time
import pdb
import re
import argparse
from examplePolynomialPathFitter_plotting import (plot_coordinate_samples, plot_path_lengths, plot_moment_arms)
import numpy as np

# start with the blanket analyze subject function to call other simulations. 
def analyzeSubject(subject, condition, trial, whatfailed, trackGRF, halfcycle, fitpaths, wantpaths, jointreact, guessmin, guess100, guessIK, guessPrev):
    # set up the paths
    print('working on Subject-condition-trial...')
    repodir = 'C:\\Users\\jonstingel\\code\\muscleModel\\muscleEnergyModel\\'
    resultsbasedir = os.path.join(repodir,'..\\results\\')
    analysisbasedir = os.path.join(repodir,'..\\analysis\\')
    # get the current working directory
    workingdir = os.getcwd()
    sctdir = os.path.join(resultsbasedir,subject,condition,trial)
    
    
    # get the trial name
    trialname = os.path.basename(workingdir)
    # go up one level
    os.chdir('..')
    # get the condition name
    condname = os.path.basename(os.getcwd())
    # go up one level
    os.chdir('..')
    # get the subject name
    subjectname = os.path.basename(os.getcwd())
    # get the experiment name
    experimentname = subjectname[0:4]
    # go back to the working directory
    os.chdir(workingdir)
    # create a list of issues
    Issues = []
    # muscleStateTrackGRFPrescribe_secondpass(repodir, subjectname, condname, trialname)
    whatfailed = muscleStateTrackGRFPrescribe_thirdpass(repodir, subjectname, condname, trialname, whatfailed, trackGRF, fitpaths, wantpaths, jointreact, guessmin, guess100, guessIK, guessPrev)
    # whatfailed = torqueStateTrackGRFTrack(repodir, subjectname, condname, trialname, whatfailed, trackGRF, halfcycle)
    return whatfailed

# analyze mimic to just do some post processing
def analyzeSubject_post(subject, condition, trial, exoTensions):
    # set up the paths
    print('working on Subject-condition-trial...')
    repodir = 'C:\\Users\\jonstingel\\code\\muscleModel\\muscleEnergyModel\\'
    resultsbasedir = os.path.join(repodir,'..\\results\\')
    analysisbasedir = os.path.join(repodir,'..\\analysis\\')
    # get the current working directory
    workingdir = os.getcwd()
    sctdir = os.path.join(resultsbasedir,subject,condition,trial)
    
    subjdir = os.path.join(repodir, '..\\results\\' + subject)
    
    # get the trial name
    trialname = os.path.basename(workingdir)
    # go up one level
    os.chdir('..')
    # get the condition name
    condname = os.path.basename(os.getcwd())
    # go up one level
    os.chdir('..')
    # get the subject name
    subjectname = os.path.basename(os.getcwd())
    # get the experiment name
    experimentname = subjectname[0:4]
    # go back to the working directory
    os.chdir(workingdir)
    # create a list of issues
    Issues = []
    # # muscleStateTrackGRFPrescribe_secondpass(repodir, subjectname, condname, trialname)
    # # muscleStateTrackGRFPrescribe_thirdpass(repodir, subjectname, condname, trialname)
    # # solution1 = osim.MocoTrajectory('muscle_statetrack_grfprescribe_solution_100con_py.sto')
    # solution2 = osim.MocoTrajectory('muscle_statetrack_grfprescribe_solution_redoarms_py.sto')

    # # ID plotter
    # # ouf.IDplotter(osim.TimeSeriesTable('muscletrack_moments_py.sto'), 'muscletrack', True)
    # ouf.IDplotter(osim.TimeSeriesTable('muscletrack_redo_moments_py.sto'), 'muscletrack_redo', True, [subject, condition, trial])

    # get exotendon tension plots/data
    exoTensions = ouf.getExotendonTension(subjectname, trialname, workingdir, exoTensions)



    ## lets see if we can get the polynomial paths fit and a model outputs
    # modelProcessor = osim.ModelProcessor('simple_model_all_the_probes.osim')
    # modelProcessor = osim.ModelProcessor('13pre_polyfit.osim')
    # newmodel = replaceMusclePaths(modelProcessor, 'results_IK_redoarms.mot', newfit=True, subjdir=subjdir)

# analyze mimic to just do some post processing
def cumulativeLoading(subject, condition, trial):
    # set up the paths
    print('working on Subject-condition-trial...')
    repodir = 'C:\\Users\\jonstingel\\code\\muscleModel\\muscleEnergyModel\\'
    resultsbasedir = os.path.join(repodir,'..\\results\\')
    analysisbasedir = os.path.join(repodir,'..\\analysis\\')
    # get the current working directory
    workingdir = os.getcwd()
    sctdir = os.path.join(resultsbasedir,subject,condition,trial)
    
    subjdir = os.path.join(repodir, '..\\results\\' + subject)
    
    # get the trial name
    trialname = os.path.basename(workingdir)
    # go up one level
    os.chdir('..')
    # get the condition name
    condname = os.path.basename(os.getcwd())
    # go up one level
    os.chdir('..')
    # get the subject name
    subjectname = os.path.basename(os.getcwd())
    # get the experiment name
    experimentname = subjectname[0:4]
    # go back to the working directory
    os.chdir(workingdir)
    # create a list of issues
    Issues = []
    
    # import matplotlib.pyplot as plt
    # import pdb; pdb.set_trace()

    
    # figure out how to access the GRF
    GRFFile = os.path.join(workingdir, 'ground_reaction.mot')
    grfData = osim.TimeSeriesTable(GRFFile)
    # get time
    grfTime = grfData.getIndependentColumn()
    # get vertical
    grfVertical_r = grfData.getDependentColumn('rF_y').to_numpy()
    grfVertical_l = grfData.getDependentColumn('lF_y').to_numpy()

    # figure out how to get just stance parts, and time
    grfY_r = [grfVertical_r[i] for i in range(len(grfVertical_r)//2) if grfVertical_r[i] > 0]
    stance_r = [grfTime[i] for i in range(len(grfVertical_r)//2) if grfVertical_r[i] > 0]
    grfY_l = [grfVertical_l[i] for i in range(len(grfVertical_l)//2,len(grfVertical_l)) if grfVertical_l[i] > 0]
    stance_l = [grfTime[i] for i in range(len(grfVertical_l)//2,len(grfVertical_l)) if grfVertical_l[i] > 0]

    # extra time
    stance_extra = [grfTime[i] for i in range(len(grfVertical_r)//2,len(grfVertical_r)) if grfVertical_r[i] > 0]

    # then compute stance durations... .
    stancetime = np.average([stance_r[-1] - stance_r[0], stance_l[-1] - stance_l[0]])

    # then load in the actual knee force data
    kneeForceFile = os.path.join(workingdir, 'jr_analysis_redo_jra_redo_all_ReactionLoads.sto')
    kneeForceData = osim.TimeSeriesTable(kneeForceFile)
    jraTime = kneeForceData.getIndependentColumn()
    kneeForce_r = kneeForceData.getDependentColumn('walker_knee_r_on_tibia_r_in_tibia_r_fy').to_numpy()
    kneeForce_l = kneeForceData.getDependentColumn('walker_knee_l_on_tibia_l_in_tibia_l_fy').to_numpy()
    
    # # get shortened knee force data based on start and end of stance times. 
    # kneeForce_r_stance = [kneeForce_r[i] for i in range(len(kneeForce_r)) if jraTime[i] >= stance_r[0] and jraTime[i] <= stance_r[-1]]
    # jraShortTime_r = [jraTime[i] for i in range(len(jraTime)) if jraTime[i] >= stance_r[0] and jraTime[i] <= stance_r[-1]]
    # kneeForce_l_stance = [kneeForce_l[i] for i in range(len(kneeForce_l)) if jraTime[i] >= stance_l[0] and jraTime[i] <= stance_l[-1]]
    # jraShortTime_l = [jraTime[i] for i in range(len(jraTime)) if jraTime[i] >= stance_l[0] and jraTime[i] <= stance_l[-1]]

    # # then compute cumulative loading - integral of force, divided by stance time, and divided by stride length
    # if len(jraShortTime_r) > 1 and len(kneeForce_r_stance) > 1:
    #     kneeForce_r_avg = np.trapz(kneeForce_r_stance, jraShortTime_r) #/ (jraShortTime_r[-1] - jraShortTime_r[0])
    # else:
    #     kneeForce_r_integral = 0
    # if len(jraShortTime_l) > 1 and len(kneeForce_l_stance) > 1:
    #     kneeForce_l_avg = np.trapz(kneeForce_l_stance, jraShortTime_l) #/ (jraShortTime_l[-1] - jraShortTime_l[0])
    # else:
    #     kneeForce_l_integral = 0
    

    # FULL STRIDESSSSSSS 
    kneeForce_r_avg = np.trapz(kneeForce_r, jraTime)
    kneeForce_l_avg = np.trapz(kneeForce_l, jraTime)


    # now think about step length and dividing by that as well. 
    avgSpeed = 2.7 # m/s
    # stepTime = np.average([stance_l[0] - stance_r[0], stance_extra[0] - stance_l[0]])
    # stepLength = avgSpeed * stepTime # m
    strideTime = jraTime[-1] - jraTime[0] # seconds
    strideLength = avgSpeed * strideTime # m
    # then return them
    cumulativeLoad = -np.average([kneeForce_r_avg, kneeForce_l_avg]) / (strideLength * strideTime)

    # load model and get mass
    modelFile = os.path.join(workingdir, 'simple_model_all_the_probes.osim')
    model = osim.Model(modelFile)
    s = model.initSystem()
    modelmass = model.getTotalMass(s)
    # import pdb; pdb.set_trace()

    return cumulativeLoad/(modelmass * 9.81) # in multiples of body weight

# muscle driven state tracking simulation - second pass 
def muscleStateTrackGRFPrescribe_secondpass(repodir, subjectname, conditionname, trialname):
    # create the tracking problem
    track = osim.MocoTrack()
    track.setName("muscle_statetrack_grfprescribe")
    # construct a ModelProcessor and add it to the tool.
    modelProcessor = osim.ModelProcessor("simple_model_all_the_probes_adjusted.osim")
    modelProcessor.append(osim.ModOpAddExternalLoads("grf_walk.xml"))
    weldem = osim.StdVectorString()
    weldem.append('mtp_r')
    weldem.append('mtp_l')
    modelProcessor.append(osim.ModOpReplaceJointsWithWelds(weldem))
    modelProcessor.append(osim.ModOpReplaceMusclesWithDeGrooteFregly2016())
    modelProcessor.append(osim.ModOpScaleActiveFiberForceCurveWidthDGF(1.5))
    modelProcessor.append(osim.ModOpAddReserves(1.0))
    basemodel = modelProcessor.process()
    basemodel.printToXML('basemodel_simple_model_all_the_probes.osim')
    basemodel = ouf.probeActivate(basemodel)
    basemodel.initSystem()
    basemuscles = basemodel.updMuscles()
    numBaseMuscles = basemuscles.getSize()
    # for m = 0:numBaseMuscles-1
    #     set tendon compliance on for certain muscles
    #     if lopt > lst want stiff (ignore)
    #     get the muscle
    #     basemusc = basemuscles.get(m)
    #     get lopt
    #     baselopt = basemusc.getOptimalFiberLength()
    #     get lst
    #     baselst = basemusc.getTendonSlackLength()
    #     set compliance if lopt > lst
    #     if baselopt < baselst
    #         basemusc.set_ignore_tendon_compliance(false)

    modelProcessorDC = osim.ModelProcessor(basemodel)
    modelProcessorDC.append(osim.ModOpFiberDampingDGF(0.01))
    # modelProcessorDC.append(osim.ModOpAddReserves(1, 2.5, True))
    modelProcessorDC.append(osim.ModOpTendonComplianceDynamicsModeDGF('implicit'))
    track.setModel(modelProcessorDC)
    # construct a TableProcessor of the coordinate data and pass it to the tracking tool.
    tableProcessor = osim.TableProcessor('coordinates_updated.mot')
    tableProcessor.append(osim.TabOpLowPassFilter(6))
    tableProcessor.append(osim.TabOpUseAbsoluteStateNames())
    track.setStatesReference(tableProcessor)
    # prescribeTable = osim.TableProcessor('muscleprescribe_states.sto')
    tempkintable = osim.TimeSeriesTable('coordinates_updated.mot')
    track.set_states_global_tracking_weight(100)
    track.set_allow_unused_references(True)
    track.set_track_reference_position_derivatives(True)
    # set specific weights for the individual weight set
    coordinateweights = osim.MocoWeightSet()
    coordinateweights.cloneAndAppend(osim.MocoWeight("pelvis_tx", 1e5))
    coordinateweights.cloneAndAppend(osim.MocoWeight("pelvis_ty", 1e7))
    coordinateweights.cloneAndAppend(osim.MocoWeight("pelvis_tz", 1e3))
    coordinateweights.cloneAndAppend(osim.MocoWeight("pelvis_list", 1e6))
    coordinateweights.cloneAndAppend(osim.MocoWeight("pelvis_rotation", 1e6))
    coordinateweights.cloneAndAppend(osim.MocoWeight("pelvis_tilt", 1e6))
    coordinateweights.cloneAndAppend(osim.MocoWeight("hip_rotation_r", 1e-6))
    coordinateweights.cloneAndAppend(osim.MocoWeight("hip_rotation_l", 1e-6))
    coordinateweights.cloneAndAppend(osim.MocoWeight("hip_adduction_r", 1e5))
    coordinateweights.cloneAndAppend(osim.MocoWeight("hip_adduction_l", 1e5))
    coordinateweights.cloneAndAppend(osim.MocoWeight("ankle_angle_r", 1e2))
    coordinateweights.cloneAndAppend(osim.MocoWeight("ankle_angle_l", 1e2))
    coordinateweights.cloneAndAppend(osim.MocoWeight("subtalar_angle_r", 1e-6))
    coordinateweights.cloneAndAppend(osim.MocoWeight("subtalar_angle_l", 1e-6))
    coordinateweights.cloneAndAppend(osim.MocoWeight('lumber_extension', 1e3))
    coordinateweights.cloneAndAppend(osim.MocoWeight('lumber_bending', 1e3))
    coordinateweights.cloneAndAppend(osim.MocoWeight('lumber_rotation', 1e3))
    track.set_states_weight_set(coordinateweights)

    # get the gait cycle timings
    cycles = ouf.subjectGaitTimings()
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

    # set the times and mesh interval, mesh points are computed internally.
    track.set_initial_time(gait_start)
    track.set_final_time(gait_end)
    track.set_mesh_interval(0.03)
    # initialize and set goals
    study = track.initialize()
    problem = study.updProblem()
    # effort goal
    effort = osim.MocoControlGoal.safeDownCast(problem.updGoal('control_effort'))
    effort.setWeight(0.5)
    # initial activation goals
    initactivationgoal = osim.MocoInitialActivationGoal('init_activation')
    initactivationgoal.setWeight(10)
    problem.addGoal(initactivationgoal)
    # put large weight on the pelvis CoordinateActuators, which act as the
    # residual, or 'hand-of-god' forces which we would like to keep small
    model = modelProcessorDC.process()
    model.printToXML('post_simple_model_all_the_probes_muscletrack.osim')
    model.initSystem()
    forceSet = model.getForceSet()
    for i in range(forceSet.getSize()):
        forcePath = forceSet.get(i).getAbsolutePathString()
        if 'pelvis' in forcePath:
            print('need to dial in the pelvis actuators...')
            # effort.setWeightForControl(forcePath, 1000)
            # if 'pelvis_ty' in forcePath:
            #     effort.setWeightForControl(forcePath, 1e8)
            # if 'hip_rotation' in forcePath:
            #     effort.setWeightForControl(forcePath, 1e4)
        elif 'reserve' in forcePath and 'subtalar' in forcePath:
            effort.setWeightForControl(forcePath, 100)
        elif 'reserve' in forcePath and 'hip_rotation' in forcePath:
            effort.setWeightForControl(forcePath, 10)
        # if 'hip_rotation' in forcePath:
        #    effort.setWeightForControl(forcePath, 10)

    # set an initial guess up
    # twosteptraj = osim.MocoTrajectory('muscle_stateprescribe_grfprescribe_solution.sto')
    # twosteptraj = osim.MocoTrajectory('muscle_statetrack_grfprescribe_solution_100con.sto')
    twosteptraj = osim.MocoTrajectory('muscle_statetrack_grfprescribe_solution_100con.sto')    
    steps = twosteptraj.getNumTimes()    
    # solver changes. 
    solver = osim.MocoCasADiSolver.safeDownCast(study.updSolver())
    solver.resetProblem(problem)
    solver.set_optim_convergence_tolerance(1e-2)
    solver.set_optim_constraint_tolerance(1e-4)
    # solver.set_optim_finite_difference_scheme('forward')
    # solver.set_optim_finite_difference_scheme('central')

    guess = solver.createGuess('bounds')
    guess.write('boundsguess.sto')
    # solver.setGuess(guess)
    randomguess = osim.MocoTrajectory('boundsguess.sto')
    randomguess.resampleWithNumTimes(steps)
    # go through and overwrite the states first
    randomstatenames = randomguess.getStateNames()
    # this will cover joint values, speeds, muscle activations, and norm tendon force
    for s in range(len(randomstatenames)):
        statename = randomstatenames[s]
        # temprandom = randomguess.getStateMat(statename);
        temp2step = twosteptraj.getStateMat(statename)
        randomguess.setState(statename,temp2step)
    # go through all the controls - excitations
    randomcontrolnames = randomguess.getControlNames()
    for c in range(len(randomcontrolnames)):
        controlname = randomcontrolnames[c]
        # temprandom = randomguess.getControlMat(controlname)
        temp2step = twosteptraj.getControlMat(controlname)
        randomguess.setControl(controlname, temp2step)
    # go through others??
    # randomparamnames = randomguess.getParameterNames()
    # this is empty in the normal condition
    # multipliers
    randommultnames = randomguess.getMultiplierNames()
    for m in range(len(randommultnames)):
        multname = randommultnames[m]
        # temprandom = randomguess.getMultiplierMat(multname)
        try:
            temp2step = twosteptraj.getMultiplierMat(multname)
            randomguess.setMultiplier(multname, temp2step)
        except:
            print('did not have the multiplier in the 2 step problem solution')
    # now for the implicit derivatives
    randomderivnames = randomguess.getDerivativeNames()
    for d in range(len(randomderivnames)):
        derivname = randomderivnames[d]
        # temprandom = randomguess.getDerivativeMat(derivname)
        temp2step = twosteptraj.getDerivativeMat(derivname)
        randomguess.setDerivative(derivname, temp2step)

    # now set the guess for the solver
    solver.setGuess(randomguess)
    # solve and visualize
    solution = study.solve()
    solution.write('muscle_statetrack_grfprescribe_solution_100con_py.sto')
    osim.STOFileAdapter.write(solution.exportToControlsTable(), 'muscletrack_controls_100con_py.sto')
    osim.STOFileAdapter.write(solution.exportToStatesTable(), 'muscletrack_states_100con_py.sto')

    # solution = osim.MocoTrajectory('muscle_statetrack_grfprescribe_solution_100con_py.sto')

    # do some post analysis
    analyzeStrings_forcePaths = osim.StdVectorString();
    analyzeStrings_forcePaths.append('.*externalloads.*');
    # analyzeStrings_forcePaths.append('/contactHeel_r');
    # analyzeStrings_forcePaths.append('/contactLateralMidfoot_r');
    # analyzeStrings_forcePaths.append('/contactMedialToe_r');
    # analyzeStrings_forcePaths.append('/contactMedialMidfoot_r');
    # analyzeStrings_forcePaths.append('/contactHeel_l');
    # analyzeStrings_forcePaths.append('/contactLateralMidfoot_l');
    # analyzeStrings_forcePaths.append('/contactMedialToe_l');
    # analyzeStrings_forcePaths.append('/contactMedialMidfoot_l');
    analyzeStrings_forcePaths.append('.*contact.*');
    table_jointMoments = study.calcGeneralizedForces(solution, analyzeStrings_forcePaths);
    osim.STOFileAdapter.write(table_jointMoments, 'muscletrack_moments_py.sto');
    ouf.IDplotter(osim.TimeSeriesTable('muscletrack_moments_py.sto'), 'muscletrack', False)

    # need to write the state, control, and other contact files. 


    # study.visualize(solution)

    # post analysis etc. 
    # solution1 = osim.MocoTrajectory('muscle_statetrack_grfprescribe_solution_100con.sto')
    # solution2 = osim.MocoTrajectory('muscle_statetrack_grfprescribe_solution_100con_py.sto')


    # likely just use the matlab infrastructure to do the post analysis.
    # otherwise have to update everything to python - not worth the time likely for post analysis. s
    return

# muscle driven state tracking simulation = third pass (testing)
def muscleStateTrackGRFPrescribe_thirdpass(repodir, subjectname, conditionname, trialname, whatfailed, trackGRF, fitpaths, wantpaths, jointreact, guessmin, guess100, guessIK, guessPrev):
    
    # create the tracking problem
    track = osim.MocoTrack()
    track.setName("muscle_statetrack_grfprescribe")
    # construct a ModelProcessor and add it to the tool.

    weldem = osim.StdVectorString()
    subjdir = os.path.join(repodir, '..\\results\\' + subjectname)

    if not trackGRF:
        print('not tracking GRF')
        time.sleep(1)
        modelProcessor = osim.ModelProcessor("simple_model_all_the_probes.osim")
        weldem.append('subtalar_r')
        weldem.append('subtalar_l')
        weldem.append('radius_hand_r')
        weldem.append('radius_hand_l')
        weldem.append('mtp_r')
        weldem.append('mtp_l')
        modelProcessor.append(osim.ModOpReplaceJointsWithWelds(weldem))
        if wantpaths:
            modelProcessor = osim.ModelProcessor(replaceMusclePaths(modelProcessor, 'results_IK_redoarms.mot', newfit=fitpaths, subjdir=subjdir))
        modelProcessor.append(osim.ModOpAddExternalLoads("grf_walk.xml"))
        
    else: 
        print('tracking GRF')
        time.sleep(1)
        modelProcessor = osim.ModelProcessor("simple_model_all_the_probes_spheres.osim")

        weldem.append('subtalar_r')
        weldem.append('subtalar_l')
        weldem.append('radius_hand_r')
        weldem.append('radius_hand_l')
        weldem.append('mtp_r')
        weldem.append('mtp_l')
        modelProcessor.append(osim.ModOpReplaceJointsWithWelds(weldem))
        if wantpaths:
            modelProcessor = osim.ModelProcessor(replaceMusclePaths(modelProcessor, 'results_IK_redoarms.mot', newfit=fitpaths, subjdir=subjdir))
        modelProcessor.append(osim.ModOpAddExternalLoads("grf_walk.xml"))
    
    modelProcessor.append(osim.ModOpReplaceMusclesWithDeGrooteFregly2016())
    modelProcessor.append(osim.ModOpScaleActiveFiberForceCurveWidthDGF(1.5))
    modelProcessor.append(osim.ModOpAddReserves(1.0))
    basemodel = modelProcessor.process()
    basemodel.printToXML('basemodel_simple_model_all_the_probes_redo.osim')
    basemodel = ouf.probeActivate(basemodel)
    basemodel.initSystem()
    basemuscles = basemodel.updMuscles()
    numBaseMuscles = basemuscles.getSize()
    # for m = 0:numBaseMuscles-1
    #     set tendon compliance on for certain muscles
    #     if lopt > lst want stiff (ignore)
    #     get the muscle
    #     basemusc = basemuscles.get(m)
    #     get lopt
    #     baselopt = basemusc.getOptimalFiberLength()
    #     get lst
    #     baselst = basemusc.getTendonSlackLength()
    #     set compliance if lopt > lst
    #     if baselopt < baselst
    #         basemusc.set_ignore_tendon_compliance(false)

    modelProcessorDC = osim.ModelProcessor(basemodel)
    modelProcessorDC.append(osim.ModOpFiberDampingDGF(0.01))
    # modelProcessorDC.append(osim.ModOpAddReserves(1, 2.5, True))
    modelProcessorDC.append(osim.ModOpTendonComplianceDynamicsModeDGF('implicit'))
    track.setModel(modelProcessorDC)

    # construct a TableProcessor of the coordinate data and pass it to the tracking tool.
    tableProcessor = osim.TableProcessor('results_IK_redoarms.mot')
    tableProcessor.append(osim.TabOpLowPassFilter(15))
    tableProcessor.append(osim.TabOpUseAbsoluteStateNames())
    track.setStatesReference(tableProcessor)
    # prescribeTable = osim.TableProcessor('muscleprescribe_states.sto')
    tempkintable = osim.TimeSeriesTable('results_IK_redoarms.mot')
    track.set_states_global_tracking_weight(100) 
    track.set_allow_unused_references(True)
    track.set_track_reference_position_derivatives(True)
    # set specific weights for the individual weight set
    coordinateweights = osim.MocoWeightSet()
    coordinateweights.cloneAndAppend(osim.MocoWeight("pelvis_tx", 1e5))
    coordinateweights.cloneAndAppend(osim.MocoWeight("pelvis_ty", 1e7))
    coordinateweights.cloneAndAppend(osim.MocoWeight("pelvis_tz", 1e3))
    coordinateweights.cloneAndAppend(osim.MocoWeight("pelvis_list", 1e6))
    coordinateweights.cloneAndAppend(osim.MocoWeight("pelvis_rotation", 1e6))
    coordinateweights.cloneAndAppend(osim.MocoWeight("pelvis_tilt", 1e6))
    coordinateweights.cloneAndAppend(osim.MocoWeight("hip_rotation_r", 1e-6))
    coordinateweights.cloneAndAppend(osim.MocoWeight("hip_rotation_l", 1e-6))
    coordinateweights.cloneAndAppend(osim.MocoWeight("hip_adduction_r", 1e8))
    coordinateweights.cloneAndAppend(osim.MocoWeight("hip_adduction_l", 1e8))
    coordinateweights.cloneAndAppend(osim.MocoWeight("ankle_angle_r", 1e2))
    coordinateweights.cloneAndAppend(osim.MocoWeight("ankle_angle_l", 1e2))
    coordinateweights.cloneAndAppend(osim.MocoWeight("subtalar_angle_r", 1e-6))
    coordinateweights.cloneAndAppend(osim.MocoWeight("subtalar_angle_l", 1e-6))
    coordinateweights.cloneAndAppend(osim.MocoWeight('lumber_extension', 1e3))
    coordinateweights.cloneAndAppend(osim.MocoWeight('lumber_bending', 1e3))
    coordinateweights.cloneAndAppend(osim.MocoWeight('lumber_rotation', 1e3))
    # track.set_states_weight_set(coordinateweights)

    # get the gait cycle timings
    cycles = ouf.subjectGaitTimings()
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

    # set the times and mesh interval, mesh points are computed internally.
    track.set_initial_time(gait_start)
    track.set_final_time(gait_end)
    track.set_mesh_interval(0.03)
    # initialize and set goals
    study = track.initialize()
    problem = study.updProblem()
    # effort goal
    effort = osim.MocoControlGoal.safeDownCast(problem.updGoal('control_effort'))
    effort.setWeight(0.5)
    # initial activation goals
    initactivationgoal = osim.MocoInitialActivationGoal('init_activation')
    initactivationgoal.setWeight(10)
    problem.addGoal(initactivationgoal)
    # put large weight on the pelvis CoordinateActuators, which act as the
    # residual, or 'hand-of-god' forces which we would like to keep small
    model = modelProcessorDC.process()
    if wantpaths:
        model.printToXML('post_simple_model_all_the_probes_muscletrack_redoarms_poly.osim')
    else: 
        model.printToXML('post_simple_model_all_the_probes_muscletrack_redo.osim')
    model.initSystem()
    forceSet = model.getForceSet()
    for i in range(forceSet.getSize()):
        forcePath = forceSet.get(i).getAbsolutePathString()
        if 'pelvis' in forcePath:
            print('need to dial in the pelvis actuators...')
            effort.setWeightForControl(forcePath, 1e-5)
            # if 'pelvis_ty' in forcePath:
            #     effort.setWeightForControl(forcePath, 1e8)
            # if 'hip_rotation' in forcePath:
            #     effort.setWeightForControl(forcePath, 1e4)
        # elif 'reserve' in forcePath and 'pelvis' in forcePath:
        #     effort.setWeightForControl(forcePath, 1e-4)
        elif 'reserve' in forcePath and 'subtalar' in forcePath:
            effort.setWeightForControl(forcePath, 100)
        elif 'reserve' in forcePath and 'hip_rotation' in forcePath:
            effort.setWeightForControl(forcePath, 10)
        # if 'hip_rotation' in forcePath:
        #    effort.setWeightForControl(forcePath, 10)
        # elif 'vas' in forcePath:
        #     print(forcePath)
        #     effort.setWeightForControl(forcePath, 2.5)

    # set up the moment tracking goal
    # test a moment tracking goal from the id moments
    # Add a joint moment tracking goal to the problem.
    jointMomentTracking = osim.MocoGeneralizedForceTrackingGoal('joint_moment_tracking', 7.0) # type: ignore
    # low-pass filter the data at 10 Hz. The reference data should use the 
    # same column label format as the output of the Inverse Dynamics Tool.
    jointMomentRef = osim.TableProcessor('./IDactual/inverse_dynamics.sto')
    # jointMomentRef.append(osim.TabOpLowPassFilter(10))
    jointMomentTracking.setReference(jointMomentRef)
    # Set the force paths that will be applied to the model to compute the
    # generalized forces. Usually these are the external loads and actuators 
    # (e.g., muscles) should be excluded, but any model force can be included 
    # or excluded. Gravitational force is applied by default.
    # Regular expression are supported when setting the force paths.
    forcePaths = osim.StdVectorString()
    forcePaths.append('.*externalloads.*')
    forcePaths.append('.*contact.*')
    jointMomentTracking.setForcePaths(forcePaths)
    # Allow unused columns in the reference data.
    jointMomentTracking.setAllowUnusedReferences(True)
    # Normalize the tracking error for each generalized for by the maximum 
    # absolute value in the reference data for that generalized force.
    jointMomentTracking.setNormalizeTrackingError(True)
    # Ignore coordinates that are locked, prescribed, or coupled to other
    # coordinates via CoordinateCouplerConstraints (true by default).
    jointMomentTracking.setIgnoreConstrainedCoordinates(True)
    # Do not track generalized forces associated with pelvis residuals.
    # testing only tracking the knee and the ankle moments, and letting everything else free
    jointMomentTracking.setWeightForGeneralizedForcePattern('.*pelvis.*', 0)
    jointMomentTracking.setWeightForGeneralizedForcePattern('.*mtp.*', 0)
    jointMomentTracking.setWeightForGeneralizedForcePattern('.*subtalar.*', 0)
    jointMomentTracking.setWeightForGeneralizedForcePattern('.*radius_hand.*', 0)
    jointMomentTracking.setWeightForGeneralizedForcePattern('.*knee.*', 200)
    jointMomentTracking.setWeightForGeneralizedForcePattern('.*beta.*', 0)
    jointMomentTracking.setWeightForGeneralizedForcePattern('.*ankle.*', 200)
    jointMomentTracking.setWeightForGeneralizedForcePattern('.*hip.*', 0)
    jointMomentTracking.setWeightForGeneralizedForcePattern('.*lumbar.*', 0)
    jointMomentTracking.setWeightForGeneralizedForcePattern('.*arm.*', 0)
    jointMomentTracking.setWeightForGeneralizedForcePattern('.*elbow.*', 0)
    jointMomentTracking.setWeightForGeneralizedForcePattern('.*pro.*', 0)
    jointMomentTracking.setWeightForGeneralizedForcePattern('.*wrist.*', 0)
    problem.addGoal(jointMomentTracking)

    ### set up a joint reaction goal to minimize knee joint contact... 
    jointReaction_r = osim.MocoJointReactionGoal('joint_reaction_r', 0.05)#0.0075) #0.5 0.05 before
    # jointpath = osim.StdVectorString(); jointpath.append('/jointset/walker_knee_r')
    # loadframe = osim.StdVectorString(); loadframe.append('child')
    # framepaths = osim.StdVectorString(); framepaths.append('/bodyset/tibia_r')
    whichForces_r = osim.StdVectorString(); whichForces_r.append('force-y')
    jointReaction_r.setJointPath('/jointset/walker_knee_r')
    jointReaction_r.setLoadsFrame('child')
    jointReaction_r.setExpressedInFramePath('/bodyset/tibia_r')
    jointReaction_r.setReactionMeasures(whichForces_r)
    ### set up a joint reaction goal to minimize knee joint contact... 
    jointReaction_l = osim.MocoJointReactionGoal('joint_reaction_l', 0.0075) #0.5 0.05 before
    # jointpath = osim.StdVectorString(); jointpath.append('/jointset/walker_knee_r')
    # loadframe = osim.StdVectorString(); loadframe.append('child')
    # framepaths = osim.StdVectorString(); framepaths.append('/bodyset/tibia_r')
    whichForces_l = osim.StdVectorString(); whichForces_l.append('force-y')
    jointReaction_l.setJointPath('/jointset/walker_knee_l')
    jointReaction_l.setLoadsFrame('child')
    jointReaction_l.setExpressedInFramePath('/bodyset/tibia_l')
    jointReaction_l.setReactionMeasures(whichForces_l)
    if jointreact:
        problem.addGoal(jointReaction_r)
        problem.addGoal(jointReaction_l)


    ### grf tracking goal... if we specify that we want it... (not default for this project.)
    GRFTrackingWeight = 1
    # set up the GRF tracking goal
    # if GRFTrackingWeight != 0:
    if trackGRF:
        # % Track the right and left vertical and fore-aft ground reaction forces.
        contactTracking = osim.MocoContactTrackingGoal('contact', GRFTrackingWeight)
        # what data are we tracking, GRF exp, or from tight tracking results


        contactTracking.setExternalLoadsFile('grf_walk.xml')
        # if trackIK:
        #     contactTracking.setExternalLoadsFile('grf_walk_nat_1.xml')
        # else:
        #     # contactTracking.setExternalLoadsFile('grf_walk_nat_1_tight.xml'); # grf_walk - Copy
        #     ## current work around for time changing... not easy way to access the xml and adjust the names and things... 
        #     # contactTracking.setExternalLoadsFile('grf_walk_nat_1_tight_poly_' + str(finalTime*2)[2:] + '.xml')
        #     # contactTracking.setExternalLoadsFile('grf_walk_nat_1_extratight_poly_' + str(finalTime*2)[2:] + '.xml')
        #     contactTracking.setExternalLoadsFile('grf_walk_nat_1_9tight_poly_' + str(finalTime*2)[2:] + '.xml')

        forceNamesRightFoot = osim.StdVectorString();
        forceNamesRightFoot.append('/contactHeel_r');
        # forceNamesRightFoot.append('/forceset/contactLateralRearfoot_r');
        forceNamesRightFoot.append('/contactLateralMidfoot_r');
        # forceNamesRightFoot.append('/contactLateralToe_r');
        forceNamesRightFoot.append('/contactMedialToe_r');
        forceNamesRightFoot.append('/contactMedialMidfoot_r');
        # contactTracking.addContactGroup(forceNamesRightFoot, 'Right_GRF');
        contactTrackingSplitRight = osim.MocoContactTrackingGoalGroup(forceNamesRightFoot, 'Right_GRF');
        contactTrackingSplitRight.append_alternative_frame_paths('/bodyset/toes_r')
        contactTracking.addContactGroup(contactTrackingSplitRight);

        forceNamesLeftFoot = osim.StdVectorString();
        forceNamesLeftFoot.append('/contactHeel_l');
        # forceNamesLeftFoot.append('/forceset/contactLateralRearfoot_l');
        forceNamesLeftFoot.append('/contactLateralMidfoot_l');
        # forceNamesLeftFoot.append('/contactLateralToe_l');
        forceNamesLeftFoot.append('/contactMedialToe_l');
        forceNamesLeftFoot.append('/contactMedialMidfoot_l');
        # contactTracking.addContactGroup(forceNamesLeftFoot, 'Left_GRF');
        contactTrackingSplitLeft = osim.MocoContactTrackingGoalGroup(forceNamesLeftFoot, 'Left_GRF');
        contactTrackingSplitLeft.append_alternative_frame_paths('/bodyset/toes_l')
        contactTracking.addContactGroup(contactTrackingSplitLeft);
        
        contactTracking.setProjection('plane');
        contactTracking.setProjectionVector(osim.Vec3(0, 0, 1));

        # contactTracking.setDivideByDuration(True)
        contactTracking.setDivideByMass(True)
        problem.addGoal(contactTracking);


    # guess100 = True then we want the original 100con as IG

    wantguess = True
    # set an initial guess up
    if wantguess:
        ## the try does the most recent solution, the except does the 100con solution from previous study.
        
        # twosteptraj = osim.MocoTrajectory('muscle_statetrack_grfprescribe_solution_100con.sto')
        if not guess100:
            if guessmin:
                try: 
                    twosteptraj = osim.MocoTrajectory('muscle_statetrack_grfprescribe_solution_redoarms_jointreact_py.sto')
                    twoname = 'muscle_statetrack_grfprescribe_solution_redoarms_jointreact_py.sto'
                except:
                    twosteptraj = osim.MocoTrajectory('muscle_statetrack_grfprescribe_solution_redoarms_py.sto')
                    twoname = 'muscle_statetrack_grfprescribe_solution_redoarms_py.sto'
            
            else:
                if trackGRF:
                    twosteptraj = osim.MocoTrajectory('muscle_statetrack_grftrack_solution_redoarms_poly_py.sto')
                    twoname = 'muscle_statetrack_grftrack_solution_redoarms_poly_py.sto'
                else:    
                    if not wantpaths:
                        twosteptraj = osim.MocoTrajectory('muscle_statetrack_grfprescribe_solution_redoarms_py.sto')
                        twoname = 'muscle_statetrack_grfprescribe_solution_redoarms_py.sto'
                    else: 
                        twosteptraj = osim.MocoTrajectory('muscle_statetrack_grfprescribe_solution_redoarms_poly_py.sto')
                        twoname = 'muscle_statetrack_grfprescribe_solution_redoarms_poly_py.sto'
        else:
            twosteptraj = osim.MocoTrajectory('muscle_statetrack_grfprescribe_solution_100con.sto')
            twoname = 'muscle_statetrack_grfprescribe_solution_100con.sto'
        
        print('twosteptraj: ' + twoname)
        # import pdb; pdb.set_trace()
        # twosteptraj = osim.MocoTrajectory('thirdpass_IG.sto')    
        steps = twosteptraj.getNumTimes()
        # solver changes. 
        solver = osim.MocoCasADiSolver.safeDownCast(study.updSolver())
        solver.resetProblem(problem)
        solver.set_optim_convergence_tolerance(1e-2)
        solver.set_optim_constraint_tolerance(1e-4)
        solver.set_parallel(16)
        # set the number of intervals potentially 
        # solver.set_num_mesh_intervals(24)

        # solver.set_optim_max_iterations(3)
        # solver.set_optim_finite_difference_scheme('forward')
        # solver.set_optim_finite_difference_scheme('central')
        if not trackGRF: 
            solver.set_optim_max_iterations(6000)
        else: 
            solver.set_optim_max_iterations(6000)
            
        guess = solver.createGuess('bounds')
        guess.write('boundsguess.sto')
        # solver.setGuess(guess)
        randomguess = osim.MocoTrajectory('boundsguess.sto')
        if randomguess.getNumTimes() != steps:
            print('resampling the guess')
            try:
                randomguess.resampleWithNumTimes(steps)
            except:
                print('could not resample the guess')
                print(os.getcwd())
                pdb.set_trace()
                return
        if not guessIK:
            # go through and overwrite the states first
            randomstatenames = randomguess.getStateNames()
            # this will cover joint values, speeds, muscle activations, and norm tendon force
            for s in range(len(randomstatenames)):
                statename = randomstatenames[s]
                # temprandom = randomguess.getStateMat(statename);
                try:
                    temp2step = twosteptraj.getStateMat(statename)
                    randomguess.setState(statename,temp2step)
                except:
                    print('did not have the state in the 2 step problem solution - keeping random. ')
                # temp2step = twosteptraj.getStateMat(statename)
                # randomguess.setState(statename,temp2step)
        elif guessIK:
            # go through and overwrite the states first
            randomstatenames = randomguess.getStateNames()
            # this will cover muscle activations, and norm tendon force (all kinematics stuff should be from IK)
            for s in range(len(randomstatenames)):
                statename = randomstatenames[s]
                # temprandom = randomguess.getStateMat(statename);
                if 'value' not in statename and 'speed' not in statename:
                    # this is a state that is not a coordinate, so we can overwrite it with the previous solution
                    try:
                        temp2step = twosteptraj.getStateMat(statename)
                        randomguess.setState(statename,temp2step)
                    except:
                        print('did not have the state in the 2 step problem solution - keeping random. ')
                    # temp2step = twosteptraj.getStateMat(statename)
                    # randomguess.setState(statename,temp2step)
                else: 
                    # this is a coordinate state, so we should not overwrite it with the previous solution
                    # see if we can input the ik data or coordinate values. 
                    # tempkintable
                    try:
                        print(statename)
                        import pdb; pdb.set_trace()
                        iknamefull = statename.split('/')[-2]
                        ikstate = tempkintable.getDependentColumn(iknamefull)
                        # if ikstate is not None:
                        #     # print('setting the state from the IK data')
                        #     randomguess.setState(statename, ikstate)
                        # else:
                        #     print('did not have the state in the IK data - keeping random. ')
                    except: 
                        print('did not have the state in the IK data - keeping random. ')
                        # import pdb; pdb.set_trace()
                        #
                    pass    







        # go through all the controls - excitations
        randomcontrolnames = randomguess.getControlNames()
        for c in range(len(randomcontrolnames)):
            controlname = randomcontrolnames[c]
            # temprandom = randomguess.getControlMat(controlname)
            try:
                temp2step = twosteptraj.getControlMat(controlname)
                randomguess.setControl(controlname, temp2step)
            except:
                print('did not have the control in the 2 step problem solution - keeping random. ')
            # temp2step = twosteptraj.getControlMat(controlname)
            # randomguess.setControl(controlname, temp2step)
        # go through others??
        # randomparamnames = randomguess.getParameterNames()
        # this is empty in the normal condition
        # multipliers
        randommultnames = randomguess.getMultiplierNames()
        for m in range(len(randommultnames)):
            multname = randommultnames[m]
            # temprandom = randomguess.getMultiplierMat(multname)
            try:
                temp2step = twosteptraj.getMultiplierMat(multname)
                randomguess.setMultiplier(multname, temp2step)
            except:
                print('did not have the multiplier in the 2 step problem solution')
        # now for the implicit derivatives
        randomderivnames = randomguess.getDerivativeNames()
        for d in range(len(randomderivnames)):
            derivname = randomderivnames[d]
            # temprandom = randomguess.getDerivativeMat(derivname)
            try:
                temp2step = twosteptraj.getDerivativeMat(derivname)
                randomguess.setDerivative(derivname, temp2step)
            except:
                print('did not have the derivative in the 2 step problem solution')
            # temp2step = twosteptraj.getDerivativeMat(derivname)
            # randomguess.setDerivative(derivname, temp2step)

        # now set the guess for the solver
        solver.setGuess(randomguess)


    # # solve and visualize
    # try:
    #     solution = study.solve()
    #     if wantpaths: 
    #         if trackGRF:
    #             if jointreact: 
    #                 solution.write('muscle_statetrack_grftrack_solution_redoarms_jointreact_poly_py.sto')
    #                 osim.STOFileAdapter.write(solution.exportToControlsTable(), 'muscletrack_redo_controls_jointreact_poly_py.sto')
    #                 osim.STOFileAdapter.write(solution.exportToStatesTable(), 'muscletrack_redo_states_jointreact_poly_py.sto')
    #             else:
    #                 solution.write('muscle_statetrack_grftrack_solution_redoarms_poly_py.sto')
    #                 osim.STOFileAdapter.write(solution.exportToControlsTable(), 'muscletrack_redo_controls_poly_py.sto')
    #                 osim.STOFileAdapter.write(solution.exportToStatesTable(), 'muscletrack_redo_states_poly_py.sto')
    #         else:
    #             if jointreact: 
    #                 solution.write('muscle_statetrack_grfprescribe_solution_redoarms_jointreact_poly_py.sto')
    #                 osim.STOFileAdapter.write(solution.exportToControlsTable(), 'muscletrack_redo_controls_jointreact_poly_py.sto')
    #                 osim.STOFileAdapter.write(solution.exportToStatesTable(), 'muscletrack_redo_states_jointreact_poly_py.sto')
    #             else:
    #                 solution.write('muscle_statetrack_grfprescribe_solution_redoarms_poly_py.sto')
    #                 osim.STOFileAdapter.write(solution.exportToControlsTable(), 'muscletrack_redo_controls_poly_py.sto')
    #                 osim.STOFileAdapter.write(solution.exportToStatesTable(), 'muscletrack_redo_states_poly_py.sto')
    #     else:
    #         if trackGRF:
    #             if jointreact: 
    #                 solution.write('muscle_statetrack_grftrack_solution_redoarms_jointreact_py.sto')
    #                 osim.STOFileAdapter.write(solution.exportToControlsTable(), 'muscletrack_redo_controls_jointreact_py.sto')
    #                 osim.STOFileAdapter.write(solution.exportToStatesTable(), 'muscletrack_redo_states_jointreact_py.sto')
    #             else:
    #                 solution.write('muscle_statetrack_grftrack_solution_redoarms_py.sto')
    #                 osim.STOFileAdapter.write(solution.exportToControlsTable(), 'muscletrack_redo_controls_py.sto')
    #                 osim.STOFileAdapter.write(solution.exportToStatesTable(), 'muscletrack_redo_states_py.sto')

    #         else:
    #             if jointreact: 
    #                 solution.write('muscle_statetrack_grfprescribe_solution_redoarms_jointreact_py.sto')
    #                 osim.STOFileAdapter.write(solution.exportToControlsTable(), 'muscletrack_redo_controls_jointreact_py.sto')
    #                 osim.STOFileAdapter.write(solution.exportToStatesTable(), 'muscletrack_redo_states_jointreact_py.sto')
    #             else:
    #                 solution.write('muscle_statetrack_grfprescribe_solution_redoarms_py.sto')
    #                 osim.STOFileAdapter.write(solution.exportToControlsTable(), 'muscletrack_redo_controls_py.sto')
    #                 osim.STOFileAdapter.write(solution.exportToStatesTable(), 'muscletrack_redo_states_py.sto')
        
    #     print('ran the base')
    # except:
    #     print(os.getcwd())
    #     print('could not solve the problem')
    #     whatfailed[subjectname + '_' + conditionname + '_' + trialname] = os.getcwd()
    #     solution = solution.unseal()
    #     solution.write('muscle_statetrack_grfprescribe_solution_unseal_redoarms_py.sto')
    #     return whatfailed
    
    solution = osim.MocoTrajectory('muscle_statetrack_grfprescribe_solution_redoarms_py.sto')

    # osim.STOFileAdapter.write(solution.exportToControlsTable(), 'muscletrack_redo_controls_py.sto')
    # osim.STOFileAdapter.write(solution.exportToStatesTable(), 'muscletrack_redo_states_py.sto')

    # try:
    #     solution = osim.MocoTrajectory('muscle_statetrack_grfprescribe_solution_redoarms_py.sto')
    # except:
    #     print('could not read the solution')
    #     whatfailed[subjectname + '_' + conditionname + '_' + trialname] = os.getcwd()
    #     return whatfailed

    # create file with exotendon outputs
    analyzeStrings_exotendon = osim.StdVectorString();
    analyzeStrings_exotendon.append('.*HOBL.*');
    table_exotendon = study.analyze(solution, analyzeStrings_exotendon);
    osim.STOFileAdapter.write(table_exotendon, 'muscletrack_redo_exotendon_py.sto');



    # do some post analysis
    analyzeStrings_forcePaths = osim.StdVectorString();
    analyzeStrings_forcePaths.append('.*externalloads.*');
    # analyzeStrings_forcePaths.append('/contactHeel_r');
    # analyzeStrings_forcePaths.append('/contactLateralMidfoot_r');
    # analyzeStrings_forcePaths.append('/contactMedialToe_r');
    # analyzeStrings_forcePaths.append('/contactMedialMidfoot_r');
    # analyzeStrings_forcePaths.append('/contactHeel_l');
    # analyzeStrings_forcePaths.append('/contactLateralMidfoot_l');
    # analyzeStrings_forcePaths.append('/contactMedialToe_l');
    # analyzeStrings_forcePaths.append('/contactMedialMidfoot_l');
    analyzeStrings_forcePaths.append('.*contact.*');
    table_jointMoments = study.calcGeneralizedForces(solution, analyzeStrings_forcePaths);
    osim.STOFileAdapter.write(table_jointMoments, 'muscletrack_redo_moments_py.sto');
    ouf.IDplotter(osim.TimeSeriesTable('muscletrack_redo_moments_py.sto'), 'muscletrack_redo', False, [subjectname, conditionname, trialname])

    # add muscle forces to the post analysis
    analyzeStrings_muscleForces = osim.StdVectorString();
    analyzeStrings_muscleForces.append('.*active_fiber_force');
    analyzeStrings_muscleForces.append('.*passive_fiber_force');
    analyzeStrings_muscleForces.append('.*fiber_force_along_tendon');
    table_muscleForces = study.analyze(solution, analyzeStrings_muscleForces);
    osim.STOFileAdapter.write(table_muscleForces, 'muscletrack_redo_muscleforces_py.sto');

    if trackGRF: 
        # % Extract ground reaction forces
        # % ==============================
        contact_r = osim.StdVectorString();
        contact_l = osim.StdVectorString();
        contact_r.append('/contactHeel_r');
        contact_r.append('/contactLateralMidfoot_r');
        contact_r.append('/contactMedialToe_r');
        contact_r.append('/contactMedialMidfoot_r');
        contact_l.append('/contactHeel_l');
        contact_l.append('/contactLateralMidfoot_l');
        contact_l.append('/contactMedialToe_l');
        contact_l.append('/contactMedialMidfoot_l');

        externalForcesTableFlat = osim.createExternalLoadsTableForGait(model, 
                                        solution,contact_r,contact_l);
        osim.STOFileAdapter.write(externalForcesTableFlat, 
                                    'muscle_statetrack_grftrack_grf_solution_redo_py.sto');
        fullstrideGRF = solution;
        fullstrideGRF.insertStatesTrajectory(externalForcesTableFlat);
        fullstrideGRF.write('muscle_statetrack_grftrack_grf_solution_redo_py.sto'); 
            # study.visualize(solution)

    # post analysis etc. 
    # solution1 = osim.MocoTrajectory('muscle_statetrack_grfprescribe_solution_100con.sto')
    # solution2 = osim.MocoTrajectory('muscle_statetrack_grfprescribe_solution_100con_py.sto')

    # likely just use the matlab infrastructure to do the post analysis.
    # otherwise have to update everything to python - not worth the time likely for post analysis. s
    return whatfailed

###################################################################################################
# sets up moco track problem with a torque driven model, and tracks GRF as well. 
def torqueStateTrackGRFTrack(repodir, subjectname, conditionname, trialname, whatfailed, trackGRF, halfcycle):
    # establish a few weights for the problem
    kinematicsWeight = 40
    GRFTrackingWeight = 1e1
    effortWeight = 1e-2
    momentWeight = 10
    stepsize = 0.03
    trackGRF = True


    # create the tracking problem
    track = osim.MocoTrack()
    track.setName("torque_statetrack_grftrack")
    # construct a ModelProcessor and add it to the tool.

    weldem = osim.StdVectorString()

    if not trackGRF:
        print('not tracking GRF')
        modelProcessor = osim.ModelProcessor("simple_model_all_the_probes.osim")
        # modelProcessor.append(osim.ModOpAddExternalLoads("grf_walk.xml"))
        modelProcessor.append(osim.ModOpAddExternalLoads("grf_walk - Copy.xml"))

        weldem.append('mtp_r')
        weldem.append('mtp_l')
    else: 
        print('tracking GRF')
        modelProcessor = osim.ModelProcessor("simple_model_all_the_probes_spheres.osim")

    weldem.append('subtalar_r')
    weldem.append('subtalar_l')
    weldem.append('radius_hand_r')
    weldem.append('radius_hand_l')
    modelProcessor.append(osim.ModOpReplaceJointsWithWelds(weldem))
    modelProcessor.append(osim.ModOpRemoveMuscles())
    modelProcessor.append(osim.ModOpAddReserves(250))
    # go through and remove specific coordinate actuators. 
    testmodel = modelProcessor.process()
    forceset = testmodel.updForceSet()
    count = 0
    for i in range(forceset.getSize()):
        force = forceset.get(i-count)
        if 'reserve_jointset_mtp' in force.getName():# or 'pelvis' in force.getName():
            forceset.remove(i-count)
            count += 1
        # if 'pelvis' in force.getName():
        #     force = osim.CoordinateActuator.safeDownCast(force)
        #     force.set_optimal_force(100.0)
    modelProcessor = osim.ModelProcessor(testmodel)
    torquemodel = modelProcessor.process()
    torquemodel.printToXML('torquemodel_simple_model_all_the_probes_redo.osim')
    track.setModel(modelProcessor)
    torquemodel.initSystem()


    # get the kinematics data that we are tracking
    tableProcessor = osim.TableProcessor('results_IK_redoarms.mot')
    tableProcessor.append(osim.TabOpLowPassFilter(15))
    tableProcessor.append(osim.TabOpUseAbsoluteStateNames())
    track.setStatesReference(tableProcessor)
    track.set_states_global_tracking_weight(kinematicsWeight)
    track.set_allow_unused_references(True)
    track.set_track_reference_position_derivatives(True)
    
    # set up specific weights for individual coordinates
    coordinateweights = osim.MocoWeightSet()
    coordinateweights.cloneAndAppend(osim.MocoWeight("pelvis_tx", 0.01))
    coordinateweights.cloneAndAppend(osim.MocoWeight("pelvis_ty", 0))
    coordinateweights.cloneAndAppend(osim.MocoWeight("pelvis_tz", 0.01))
    coordinateweights.cloneAndAppend(osim.MocoWeight("pelvis_list", 0.01))
    coordinateweights.cloneAndAppend(osim.MocoWeight("pelvis_rotation", 0.01))
    coordinateweights.cloneAndAppend(osim.MocoWeight("pelvis_tilt", 0.01))
    # coordinateweights.cloneAndAppend(osim.MocoWeight("hip_rotation_r", 0))
    # coordinateweights.cloneAndAppend(osim.MocoWeight("hip_rotation_l", 0))
    coordinateweights.cloneAndAppend(osim.MocoWeight("hip_adduction_r", 2))
    coordinateweights.cloneAndAppend(osim.MocoWeight("hip_adduction_l", 2))
    coordinateweights.cloneAndAppend(osim.MocoWeight("hip_flexion_r", 4))
    coordinateweights.cloneAndAppend(osim.MocoWeight("hip_flexion_l", 4))
    coordinateweights.cloneAndAppend(osim.MocoWeight("knee_angle_r", 5))
    coordinateweights.cloneAndAppend(osim.MocoWeight("knee_angle_l", 5))
    coordinateweights.cloneAndAppend(osim.MocoWeight("ankle_angle_r", 6))
    coordinateweights.cloneAndAppend(osim.MocoWeight("ankle_angle_l", 6))
    track.set_states_weight_set(coordinateweights)

    # get the individual subject-condition-trial timings
    cycles = ouf.subjectGaitTimings()
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

    # set the times and mesh interval, mesh points are computed internally.
    track.set_initial_time(gait_start)
    if halfcycle:
        gait_end = gait_end - (0.5*(gait_end - gait_start)) 
    track.set_final_time(gait_end)
    # track.set_mesh_interval(0.01)
    
    # initialize and set goals
    study = track.initialize()
    problem = study.updProblem()

    # create a symmetry/periodicity goal
    if halfcycle: 
        # % Symmetry (to permit simulating only one step)
        symmetryGoal = osim.MocoPeriodicityGoal('symmetryGoal');
        problem.addGoal(symmetryGoal);
        # % Symmetric coordinate values (except for pelvis_tx) and speeds
        for i in range(torquemodel.getNumStateVariables()):
            currentStateName = str(torquemodel.getStateVariableNames().getitem(i));
            if str.startswith(currentStateName , '/jointset') and 'beta' not in currentStateName:
                # print('\niiii')
                # print(currentStateName)
                # if 'hip_rotation_r' in currentStateName:
                #     pair = osim.MocoPeriodicityGoalPair(currentStateName, re.sub('hip_rotation_r', 'hip_rotation_l', currentStateName));
                #     symmetryGoal.addStatePair(pair);
                # if 'hip_rotation_l' in currentStateName:
                #     pair = osim.MocoPeriodicityGoalPair(currentStateName, re.sub('hip_rotation_l', 'hip_rotation_r', currentStateName));
                #     symmetryGoal.addStatePair(pair);
                # right and left limb coordinates
                if currentStateName.endswith('_r/value') or currentStateName.endswith('_r/speed') and 'hip_rotation' not in currentStateName:
                    pair = osim.MocoPeriodicityGoalPair(currentStateName,re.sub('_r/', '_l/', currentStateName));
                    symmetryGoal.addStatePair(pair);
                    # print('1 - rights and lefts - pair')
                    # print(currentStateName)
                    # print(re.sub('_r', '_l', currentStateName))
                if currentStateName.endswith('_l/value') or currentStateName.endswith('_l/speed') and 'hip_rotation' not in currentStateName:
                    pair = osim.MocoPeriodicityGoalPair(currentStateName,re.sub('_l/', '_r/', currentStateName)); 
                    symmetryGoal.addStatePair(pair);
                    # print('2 - lefts and rights - pair')
                    # print(currentStateName)
                    # print(re.sub('_l', '_r', currentStateName))
                # pelvis tilt
                if currentStateName.endswith('_tilt/value') or currentStateName.endswith('_tilt/speed'):
                    symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(currentStateName));
                    # print('3 - pelvis tilt - pair value and speed')
                    # print(currentStateName)
                # pelvis list
                if currentStateName.endswith('_list/value') or currentStateName.endswith('_list/speed'):
                    if currentStateName.endswith('/value'):
                        symmetryGoal.addNegatedStatePair(osim.MocoPeriodicityGoalPair(currentStateName));
                        # print('4 - pelvis list - negated pair value')
                        # print(currentStateName)
                    if currentStateName.endswith('/speed'):
                        symmetryGoal.addNegatedStatePair(osim.MocoPeriodicityGoalPair(currentStateName));
                        # print('4 - pelvis list - speed pair')
                        # print(currentStateName)
                # pelvis rotation
                if currentStateName.endswith('pelvis_rotation/value') or currentStateName.endswith('pelvis_rotation/speed'):
                    if currentStateName.endswith('/value'):
                        symmetryGoal.addNegatedStatePair(osim.MocoPeriodicityGoalPair(currentStateName));
                        # print('5 - pelvis rotation - negated pair value')
                        # print(currentStateName)
                    if currentStateName.endswith('/speed'):
                        symmetryGoal.addNegatedStatePair(osim.MocoPeriodicityGoalPair(currentStateName));
                        # print('5 - pelvis rotation - pair speed')
                        # print(currentStateName)
                # pelvis ty symmetry
                if currentStateName.endswith('_ty/value') or currentStateName.endswith('_ty/speed'):
                    symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(currentStateName));
                    # print('6 - pelvis ty - pair value and speed')
                    # print(currentStateName)
                # pelvis tx
                if currentStateName.endswith('_tx/speed'): # overground so not value
                    symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(currentStateName));
                    # print('7 - pelvis tx - pair speed')
                    # print(currentStateName)
                # # pelvis Tz
                # if currentStateName.endswith('_tz/value') or currentStateName.endswith('_tz/speed'):
                #     if currentStateName.endswith('value'):
                #         symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(currentStateName));
                #         print('8 - pelvis tz - pair value')
                #         print(currentStateName)
                #     if currentStateName.endswith('speed'):
                #         symmetryGoal.addNegatedStatePair(osim.MocoPeriodicityGoalPair(currentStateName));
                #         print('8 - pelvis tz - negated pair speed')
                #         print(currentStateName) 
                # lumbar extension 
                if currentStateName.endswith('_extension/value') or currentStateName.endswith('_extension/speed'):
                    symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(currentStateName));
                    # print('9 - lumbar extension - pair value and speed')
                    # print(currentStateName)
                # lumbar bending 
                if currentStateName.endswith('_bending/value') or currentStateName.endswith('_bending/speed'):
                    if currentStateName.endswith('/value'):
                        symmetryGoal.addNegatedStatePair(osim.MocoPeriodicityGoalPair(currentStateName));
                        # print('10 - lumbar bending - negated value pair')
                        # print(currentStateName)
                    if currentStateName.endswith('/speed'):
                        symmetryGoal.addNegatedStatePair(osim.MocoPeriodicityGoalPair(currentStateName));
                        # print('10 - lumbar bending - speed pair')
                        # print(currentStateName)
                # lumbar rotation
                if currentStateName.endswith('lumbar_rotation/value') or currentStateName.endswith('lumbar_rotation/speed'):
                    if currentStateName.endswith('/value'):
                        symmetryGoal.addNegatedStatePair(osim.MocoPeriodicityGoalPair(currentStateName));
                        # print('11 - lumbar rotation - negated pair value')
                        # print(currentStateName)
                    if currentStateName.endswith('_rotation/speed'):
                        symmetryGoal.addNegatedStatePair(osim.MocoPeriodicityGoalPair(currentStateName));
                        # print('11 - lumbar rotation - pair speed')
                        # print(currentStateName)

            if 'beta' in currentStateName:
                # print('\niiii')
                # print(currentStateName)
                print('beta not included in symmetry')

        print('\n\naaaaaaaaaaaaaaaaaaaahahaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa')
        # % Symmetric muscle activations
        for i in range(torquemodel.getNumStateVariables()):
            currentStateName = str(torquemodel.getStateVariableNames().getitem(i));
            
            if str.endswith(currentStateName, '/normalized_tendon_force'):
                # print('\nhhhhhh')
                # print(currentStateName)
            # TODO tendon forces symmetry
                if str.endswith(currentStateName, '_r/normalized_tendon_force'):
                    pair = osim.MocoPeriodicityGoalPair(currentStateName, re.sub('_r','_l', currentStateName));
                    symmetryGoal.addStatePair(pair);
                    # print('norm tendon pairs')
                    # print(re.sub('_r', '_l', currentStateName))
                if str.endswith(currentStateName, '_l/normalized_tendon_force'):
                    pair = osim.MocoPeriodicityGoalPair(currentStateName, re.sub('_l','_r', currentStateName));
                    symmetryGoal.addStatePair(pair);
                    # print('norm tendon pairs')
                    # print(re.sub('_l', '_r', currentStateName))

            if str.endswith(currentStateName,'/activation'):
                # print('\naaaa')
                # print(currentStateName)
                # activations squared for this actuator
                ## testing out a higher weight for the activations on the quads to get lower joint reaction forces. 


                if currentStateName.endswith('_r/activation'):
                    pair = osim.MocoPeriodicityGoalPair(currentStateName, re.sub('_r','_l', currentStateName));
                    symmetryGoal.addStatePair(pair);
                    # print('a1 - rights and lefts - pair')
                    # # print(currentStateName)
                    # print(re.sub('_r', '_l', currentStateName))
                if currentStateName.endswith('_l/activation'):
                    pair = osim.MocoPeriodicityGoalPair(currentStateName, re.sub('_l', '_r', currentStateName));
                    symmetryGoal.addStatePair(pair);
                    # print('a2 - lefts and rights - pair')
                    # # print(currentStateName)
                    # print(re.sub('_l', '_r', currentStateName))
                # bending , rotation , tz, list - these are gonna be different 
                # if 'Bend' in currentStateName or 'Rot' in currentStateName:
                #     symmetryGoal.addNegatedStatePair(osim.MocoPeriodicityGoalPair(currentStateName));
                #     print('a3 - lumbar bending, lumbar rotation actuators - negated pair')
                #     print(currentStateName)
                # if 'Ext' in currentStateName:
                #     symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(currentStateName));
                #     print('a3 - lumbar extension actuator -  pair')
                #     print(currentStateName)

        print('\nccccccccccccccccccccccccccccccccccccccccccccccccccc')
        # now get the controls and do symmetry for them as well
        # controlTab = torquemodel.getControlsTable();
        # controlnames = controlTab.getColumnLabels();
        # for i in range(len(controlnames)):
        #     currentcontrol = controlnames[i]
        #     print('/controllerset/' + currentcontrol)
        #     if str.endswith(currentcontrol, '_r'):
        #         pair = Moco
        forceSet = torquemodel.getForceSet();
        for i in range(forceSet.getSize()):
            currentforce = forceSet.get(i).getAbsolutePathString();
            # print('\nccccc')
            # print(currentforce)
            if str.endswith(currentforce, '_r') and 'Passive' not in currentforce and 'mtp' not in currentforce and '_rot' not in currentforce:
                pair = osim.MocoPeriodicityGoalPair(currentforce, re.sub('_r', '_l', currentforce));
                symmetryGoal.addControlPair(pair);
                # print(re.sub('_r', '_l', currentforce))
            if str.endswith(currentforce, '_l') and 'Passive' not in currentforce and 'mtp' not in currentforce and '_rot' not in currentforce:
                pair = osim.MocoPeriodicityGoalPair(currentforce, re.sub('_l', '_r', currentforce));
                symmetryGoal.addControlPair(pair);
                # print(re.sub('_l', '_r', currentforce))
    
    # we set up periodicity for the full gait cycle instead of symmetry between the limbs. 
    else:
        print('full cycle simulation')
        # Constrain the states and controls to be periodic.
        periodicityGoal = osim.MocoPeriodicityGoal("periodicity")
        for i in range(torquemodel.getNumStateVariables()):
            currentStateName = str(torquemodel.getStateVariableNames().getitem(i))
            if 'pelvis_tx/value' not in currentStateName and 'Passive' not in currentStateName and 'mtp' not in currentStateName and '_rot' not in currentStateName:
                periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair(currentStateName))
            
        forceSet = torquemodel.getForceSet()
        for i in range(forceSet.getSize()):
            forcePath = forceSet.get(i).getAbsolutePathString()
            if 'Passive' not in forcePath and 'mtp' not in forcePath and '_rot' not in forcePath:
                periodicityGoal.addControlPair(osim.MocoPeriodicityGoalPair(forcePath))

        problem.addGoal(periodicityGoal)



    # effort goal
    effort = osim.MocoControlGoal.safeDownCast(problem.updGoal('control_effort'))
    effort.setWeight(effortWeight)
    # initial activation goals
    initactivationgoal = osim.MocoInitialActivationGoal('init_activation')
    initactivationgoal.setWeight(10)
    problem.addGoal(initactivationgoal)
    # put large weight on the pelvis CoordinateActuators, which act as the
    # residual, or 'hand-of-god' forces which we would like to keep small
    model = torquemodel
    model.initSystem()
    forceSet = model.getForceSet()
    for i in range(forceSet.getSize()):
        forcePath = forceSet.get(i).getAbsolutePathString()
        if 'pelvis' in forcePath:
            print('need to dial in the pelvis actuators...')
            effort.setWeightForControl(forcePath, 1e2)
        elif 'reserve' in forcePath and 'subtalar' in forcePath:
            effort.setWeightForControl(forcePath, 10)
        elif 'reserve' in forcePath and 'hip_rotation' in forcePath:
            effort.setWeightForControl(forcePath, 10)
    
    # Add a joint moment tracking goal to the problem.
    jointMomentTracking = osim.MocoGeneralizedForceTrackingGoal('joint_moment_tracking', momentWeight) # type: ignore
    # low-pass filter the data at 15 Hz. The reference data should use the
    # same column label format as the output of the Inverse Dynamics Tool.
    jointMomentRef = osim.TableProcessor('./IDactual/inverse_dynamics.sto')
    jointMomentRef.append(osim.TabOpLowPassFilter(15))
    jointMomentTracking.setReference(jointMomentRef)
    # Set the force paths that will be applied to the model to compute the
    # generalized forces. Usually these are the external loads and actuators
    # (e.g., muscles) should be excluded, but any model force can be included
    # or excluded. Gravitational force is applied by default.
    # Regular expression are supported when setting the force paths.
    forcePaths = osim.StdVectorString()
    forcePaths.append('.*externalloads.*')
    forcePaths.append('.*contact.*')
    jointMomentTracking.setForcePaths(forcePaths)
    jointMomentTracking.setAllowUnusedReferences(True)
    # Normalize the tracking error for each generalized for by the maximum
    # absolute value in the reference data for that generalized force.
    jointMomentTracking.setNormalizeTrackingError(True)
    # Ignore coordinates that are locked, prescribed, or coupled to other
    # coordinates via CoordinateCouplerConstraints (true by default).
    jointMomentTracking.setIgnoreConstrainedCoordinates(True)
    # Do not track generalized forces associated with pelvis residuals.
    # testing only tracking the knee and the ankle moments, and letting everything else free
    jointMomentTracking.setWeightForGeneralizedForcePattern('.*pelvis.*', 0)
    jointMomentTracking.setWeightForGeneralizedForcePattern('.*mtp.*', 0)
    jointMomentTracking.setWeightForGeneralizedForcePattern('.*subtalar.*', 0)
    jointMomentTracking.setWeightForGeneralizedForcePattern('.*radius_hand.*', 0)
    jointMomentTracking.setWeightForGeneralizedForcePattern('.*knee.*', 800)
    jointMomentTracking.setWeightForGeneralizedForcePattern('.*beta.*', 0)
    jointMomentTracking.setWeightForGeneralizedForcePattern('.*ankle.*', 200)
    jointMomentTracking.setWeightForGeneralizedForcePattern('.*hip.*', 100)
    jointMomentTracking.setWeightForGeneralizedForcePattern('.*lumbar.*', 0)
    jointMomentTracking.setWeightForGeneralizedForcePattern('.*arm.*', 0)
    jointMomentTracking.setWeightForGeneralizedForcePattern('.*elbow.*', 0)
    jointMomentTracking.setWeightForGeneralizedForcePattern('.*pro.*', 0)
    jointMomentTracking.setWeightForGeneralizedForcePattern('.*wrist.*', 0)
    problem.addGoal(jointMomentTracking)    
    
    # set up the GRF tracking goal
    if trackGRF: 
        # % Track the right and left vertical and fore-aft ground reaction forces.
        contactTracking = osim.MocoContactTrackingGoal('contact', GRFTrackingWeight)
        # what data are we tracking, GRF exp, or from tight tracking results
        contactTracking.setExternalLoadsFile('grf_walk.xml')

        forceNamesRightFoot = osim.StdVectorString();
        forceNamesRightFoot.append('/contactHeel_r');
        forceNamesRightFoot.append('/contactLateralMidfoot_r');
        forceNamesRightFoot.append('/contactMedialToe_r');
        forceNamesRightFoot.append('/contactMedialMidfoot_r');
        contactTrackingSplitRight = osim.MocoContactTrackingGoalGroup(forceNamesRightFoot, 'Right_GRF');
        contactTrackingSplitRight.append_alternative_frame_paths('/bodyset/toes_r')
        contactTracking.addContactGroup(contactTrackingSplitRight);

        forceNamesLeftFoot = osim.StdVectorString();
        forceNamesLeftFoot.append('/contactHeel_l');
        forceNamesLeftFoot.append('/contactLateralMidfoot_l');
        forceNamesLeftFoot.append('/contactMedialToe_l');
        forceNamesLeftFoot.append('/contactMedialMidfoot_l');
        contactTrackingSplitLeft = osim.MocoContactTrackingGoalGroup(forceNamesLeftFoot, 'Left_GRF');
        contactTrackingSplitLeft.append_alternative_frame_paths('/bodyset/toes_l')
        contactTracking.addContactGroup(contactTrackingSplitLeft);
        
        contactTracking.setProjection('plane');
        contactTracking.setProjectionVector(osim.Vec3(0, 0, 1));
        # contactTracking.setDivideByDuration(True)
        contactTracking.setDivideByMass(True)
        problem.addGoal(contactTracking);

    
    
    #Configure the solver
    # % ====================
    solver = study.initCasADiSolver();
    solver.set_optim_finite_difference_scheme('forward')
    solver.set_parameters_require_initsystem(False);
    duration = gait_end - gait_start
    num_mesh = round(duration/stepsize);
    solver.set_num_mesh_intervals(num_mesh);
    solver.set_verbosity(2);
    solver.set_optim_solver('ipopt');
    solver.set_optim_convergence_tolerance(1e-2);
    solver.set_optim_constraint_tolerance(1e-4);
    # solver.set_optim_max_iterations(6000);
    solver.set_scale_variables_using_bounds(True);
    # solver.set_minimize_implicit_auxiliary_derivatives(True);
    # solver.set_implicit_auxiliary_derivatives_weight(implicitWeight);

    if not trackGRF: 
        solver.set_optim_max_iterations(3000)
    else: 
        solver.set_optim_max_iterations(6000)    

    wantguess = True
    # set an initial guess up
    ### for now lets see if it can come up with anything on its own... and how long it takes. 
    if wantguess:
        if halfcycle:
            lastguess = osim.MocoTrajectory('torque_statetrack_grftrack_solution_redoarms_halfcycle_py.sto')
        else:
            lastguess = osim.MocoTrajectory('torque_statetrack_grftrack_solution_redoarms_py.sto')
        solver.setGuess(lastguess)


    # solve and visualize
    try:
        solution = study.solve()
        if halfcycle:
            solution.write('torque_statetrack_grftrack_solution_redoarms_halfcycle_py.sto')
        else:
            solution.write('torque_statetrack_grftrack_solution_redoarms_py.sto')
        print('ran the base')
    except:
        print(os.getcwd())
        print('could not solve the problem')
        whatfailed[subjectname + '_' + conditionname + '_' + trialname] = os.getcwd()
        solution = solution.unseal()
        if halfcycle:
            solution.write('torque_statetrack_grftrack_solution_unseal_redoarms_halfcycle_py.sto')
        else: 
            solution.write('torque_statetrack_grftrack_solution_unseal_redoarms_py.sto')
        return whatfailed
    
    osim.STOFileAdapter.write(solution.exportToControlsTable(), 'torque_statetrack_grftrack_controls_py.sto')
    osim.STOFileAdapter.write(solution.exportToStatesTable(), 'torque_statetrack_grftrack_states_py.sto')
    
    # do some post analysis
    analyzeStrings_forcePaths = osim.StdVectorString();
    analyzeStrings_forcePaths.append('.*externalloads.*');
    # analyzeStrings_forcePaths.append('/contactHeel_r');
    # analyzeStrings_forcePaths.append('/contactLateralMidfoot_r');
    # analyzeStrings_forcePaths.append('/contactMedialToe_r');
    # analyzeStrings_forcePaths.append('/contactMedialMidfoot_r');
    # analyzeStrings_forcePaths.append('/contactHeel_l');
    # analyzeStrings_forcePaths.append('/contactLateralMidfoot_l');
    # analyzeStrings_forcePaths.append('/contactMedialToe_l');
    # analyzeStrings_forcePaths.append('/contactMedialMidfoot_l');
    analyzeStrings_forcePaths.append('.*contact.*');
    table_jointMoments = study.calcGeneralizedForces(solution, analyzeStrings_forcePaths);
    osim.STOFileAdapter.write(table_jointMoments, 'torque_statetrack_grftrack_redo_moments_py.sto');
    ouf.IDplotter(osim.TimeSeriesTable('torque_statetrack_grftrack_redo_moments_py.sto'), 'torque_statetrack_grftrack', False, [subjectname, conditionname, trialname])

    # # grab anything else that might be useful
    # analyzeStrings_all = osim.StdVectorString()
    # analyzeStrings_all.append('.*')
    # table_all = study.analyze(solution, analyzeStrings_all)
    # osim.STOFileAdapter.write(table_all, 'torque_statetrack_grftrack_redo_all_py.sto')

    pdb.set_trace()
    study.visualize(solution)

    # post analysis etc. 
    # solution1 = osim.MocoTrajectory('muscle_statetrack_grfprescribe_solution_100con.sto')
    # solution2 = osim.MocoTrajectory('muscle_statetrack_grfprescribe_solution_100con_py.sto')

    # likely just use the matlab infrastructure to do the post analysis.
    # otherwise have to update everything to python - not worth the time likely for post analysis. s
    return whatfailed

###################################################################################################
# sets up moco track problem with a torque driven model, and prescribed GRF.
def torqueStateTrackGRFPrescribe(repodir, subjectname, conditionname, trialname, whatfailed, trackGRF, halfcycle):
    # establish a few weights for the problem

    return whatfailed

###################################################################################################
# sets up moco track problem with a torque driven model, and prescribed GRF.
def muscleInverse(repodir, subjectname, conditionname, trialname, whatfailed, trackGRF, halfcycle):
    # establish a few weights for the problem

    return whatfailed

def replaceMusclePaths(modelProcessor, kinematicsFile, newfit, subjdir):
    import time
    start_time = time.time()
    
    print('here we are, trying to get some new paths... ')
    if newfit:
        # create and fit the paths in the osim model with the polynomial path fitting tool
        fitter = osim.PolynomialPathFitter()
        # Set the model.
        # The model should contain path-based force objects (e.g., Muscles) that use
        # geometry-based paths (e.g., GeometryPath) to model path lengths and moment
        # arms. The fitter will create a set of FunctionBasedPaths that use
        # MultivariatePolynomialFunctions to model the path lengths and moment arms
        # # of the original model.
        model = modelProcessor.process()
        model.printToXML('pre_polyfit.osim')
        fitter.setModel(modelProcessor)
        model.initSystem()
        # set the coordinate values table
        # The fitter will randomly sample around the coordinate values provided in the
        # table to generate model configurations for which to compute path lengths and
        # moment arms. This table has many more rows than are needed for the fitter to
        # generate a good fit, so we will remove some of the rows to speed up the
        # fitting process.
        # values = osim.TimeSeriesTable('coordinates_updated.mot')
        values = osim.TimeSeriesTable(kinematicsFile)
        times = values.getIndependentColumn()
        for i in range(len(times)):
            # if i % 2 != 0:
            #     values.removeRow(times[i])
            # test way to cut anything after certain index
            if i > 150: 
                values.removeRow(times[i])
        fitter.setCoordinateValues(osim.TableProcessor(values))

        # Configure optional settings
        # ---------------------------
        # Use these settings to modify the default settings to tailor the fitting
        # process to your model and motion. See the documentation for
        # PolynomialPathFitter for more options.

        # Set a (relative) directory to where the fitting results will be saved.
        # Files printed to this directory include the set of FunctionBasedPaths
        # created by the fitter, the path lengths and moment arms computed for each
        # model configuration, and the path lengths and moment arms fitted to the
        # polynomial functions. File names will be prepended with the name of the
        # model.
        results_dir = subjdir + '/welknatural/trial01/' + 'pathresults'
        fitter.setOutputDirectory(results_dir)
        
        # Set the maximum order of the polynomials used to fit the path lengths
        # and moment arms. Higher order polynomials might lead to a better fit,
        # but could increase the computational time required to evaluate the
        # path length functions.
        fitter.setMaximumPolynomialOrder(5)

        # Set the number of random samples taken at each frame around the nominal 
        # coordinate values.
        fitter.setNumSamplesPerFrame(30)

        # By default, coordinate values are sample around the nominal coordinate
        # values using bounds of [-10, 10] degrees. You can set custom bounds for
        # individual coordinates using the appendCoordinateSamplingBounds() method.
        # fitter.appendCoordinateSamplingBounds('/jointset/hip_r/hip_flexion_r', osim.Vec2(-50, 50))
        # fitter.appendCoordinateSamplingBounds('/jointset/hip_l/hip_flexion_l', osim.Vec2(-50, 50))
        # fitter.appendCoordinateSamplingBounds('/jointset/walker_knee_r/knee_angle_r', osim.Vec2(-50, 50))
        # fitter.appendCoordinateSamplingBounds('/jointset/walker_knee_l/knee_angle_l', osim.Vec2(-50, 50))
        # fitter.appendCoordinateSamplingBounds('/jointset/ankle_r/ankle_angle_r', osim.Vec2(-50, 50))
        # fitter.appendCoordinateSamplingBounds('/jointset/ankle_l/ankle_angle_l', osim.Vec2(-50, 50))

        # Set the global coordinate sampling bounds. This will be used for any
        # coordinates that do not have custom bounds set. We use reasonably
        # large bounds here to sample a wide range of the model's coordinate space
        # around the reference trajectory.
        fitter.setGlobalCoordinateSamplingBounds(osim.Vec2(-30, 30))

        # Use stepwise regression to fit the path lengths and moment arms. This
        # setting evaluates the fit after adding polynomial terms one at a time
        # to determine the minimum number of coefficients needed to achieve the
        # path length and moment arm tolerances. Stepwise regression includes 
        # polynomial terms up to the maximum order set by 
        # setMaximumPolynomialOrder().
        fitter.setUseStepwiseRegression(True)

        # Set the path length and moment arm tolerances. When the RMS errors 
        # between the original path lengths and moment arms and the fitted 
        # polynomial paths are below these tolerances, the fitting process will
        # stop for a given path. Tighter tolerances may result in a better fit,
        # but at the expense of higher polynomial orders (or more polynomial 
        # terms, if using stepwise regression).
        fitter.setPathLengthTolerance(1e-3)
        fitter.setMomentArmTolerance(1e-3)

        # Run the fitter
        # --------------
        # Information about each step fitting process will be printed to the
        # console including the path length and moment arm RMS error for
        # each force object and averaged across all force objects.        
        fitter.run()

        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"Execution time: {elapsed_time:.4f} seconds")

        # # Plot the results
        # # ----------------
        # # Use the plotting helper functions to visualize the results of the
        # # fitting process and determine if the fits are good enough for your needs,
        # # or if the model or fitting settings need to be modified.

        # # Plot the sampled coordinate values used to generate the path lengths
        # # and moment arms.
        # plot_coordinate_samples(results_dir, model.getName())

        # # Plot the path lengths and moment arms computed from the original model
        # # paths (blue) and the fitted polynomial paths (orange).
        # #
        # # For most muscles the fit is very good, but there are noticeable fitting
        # # errors in a few muscles (e.g., /forceset/gaslat_r and /forceset/glmax1_r).
        # # Errors like these usually arise from the fitting process struggling with
        # # discontinuities due to wrapping geometry issues in the original model.
        # # Depending on size of the errors, you may want to adjust the wrapping
        # # geometry in the original model and re-run the fitter.
        # plot_path_lengths(results_dir, model.getName())
        # plot_moment_arms(results_dir, model.getName())

        # Evaluate the fitted functions on a 'new' trajectory
        # ---------------------------------------------------
        # You can use PolynomialPathFitter to evaluate a set of previously fitted
        # FunctionBasedPaths on a new trajectory. This can be useful if you want to
        # evaluate the path lengths and moment arms of a model on a trajectory that
        # was not used to fit the functions. This example uses the same trajectory
        # used to fit the functions, but you can replace it with any trajectory
        # that has a set of coordinate values that are consistent with the model.
        # Note that we do not need to use the 'fitter' object to call this function
        # (i.e., it is a static function).
        functionBasedPathsFile = os.path.join(
            results_dir, f'{model.getName()}_FunctionBasedPathSet.xml')
        osim.PolynomialPathFitter.evaluateFunctionBasedPaths(
            model, osim.TableProcessor(values), functionBasedPathsFile)

    ####################################################    
    # Replacing the original paths with the fitted paths
    # --------------------------------------------------
    # You can use a ModelProcessor to replace the original paths in the model
    # with the fitted paths. This can be useful if you want to use the fitted
    # paths in a simulation or analysis tool but keep the original paths in the
    # model file.
    if not newfit:
        print('just replacing the paths with last fit...')
        results_dir = os.path.join(subjdir, 'welknatural\\trial01\\pathresults')
        model = modelProcessor.process()
        functionBasedPathsFile = os.path.join(
            results_dir, f'{model.getName()}_FunctionBasedPathSet.xml')
    try: 
        modelProcessor.append(osim.ModOpReplacePathsWithFunctionBasedPaths(
            functionBasedPathsFile))
        model = modelProcessor.process()
        model.initSystem()
        model.printToXML('model_polyfit_redo.osim')
    except: 
        print('could not replace the paths, likely the file does not exist')
        pdb.set_trace()
        return 

    return model
