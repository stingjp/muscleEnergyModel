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


# start with the blanket analyze subject function to call other simulations. 
def analyzeSubject(subject, condition, trial, whatfailed):
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
    whatfailed = muscleStateTrackGRFPrescribe_thirdpass(repodir, subjectname, condname, trialname, whatfailed)
    return whatfailed


# analyze mimic to just do some post processing
def analyzeSubject_post(subject, condition, trial):
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
    # muscleStateTrackGRFPrescribe_thirdpass(repodir, subjectname, condname, trialname)
    # solution1 = osim.MocoTrajectory('muscle_statetrack_grfprescribe_solution_100con_py.sto')
    solution2 = osim.MocoTrajectory('muscle_statetrack_grfprescribe_solution_redoarms_py.sto')

    # ID plotter
    # ouf.IDplotter(osim.TimeSeriesTable('muscletrack_moments_py.sto'), 'muscletrack', True)
    ouf.IDplotter(osim.TimeSeriesTable('muscletrack_redo_moments_py.sto'), 'muscletrack_redo', True)


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
def muscleStateTrackGRFPrescribe_thirdpass(repodir, subjectname, conditionname, trialname, whatfailed):
    # create the tracking problem
    track = osim.MocoTrack()
    track.setName("muscle_statetrack_grfprescribe")
    # construct a ModelProcessor and add it to the tool.
    modelProcessor = osim.ModelProcessor("simple_model_all_the_probes.osim")
    modelProcessor.append(osim.ModOpAddExternalLoads("grf_walk.xml"))
    weldem = osim.StdVectorString()
    weldem.append('mtp_r')
    weldem.append('mtp_l')
    weldem.append('subtalar_r')
    weldem.append('subtalar_l')
    weldem.append('radius_hand_r')
    weldem.append('radius_hand_l')
    modelProcessor.append(osim.ModOpReplaceJointsWithWelds(weldem))
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
    tableProcessor.append(osim.TabOpLowPassFilter(6))
    tableProcessor.append(osim.TabOpUseAbsoluteStateNames())
    track.setStatesReference(tableProcessor)
    # prescribeTable = osim.TableProcessor('muscleprescribe_states.sto')
    tempkintable = osim.TimeSeriesTable('results_IK_redoarms.mot')
    track.set_states_global_tracking_weight(10)
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
    model.printToXML('post_simple_model_all_the_probes_muscletrack_redo.osim')
    model.initSystem()
    forceSet = model.getForceSet()
    for i in range(forceSet.getSize()):
        forcePath = forceSet.get(i).getAbsolutePathString()
        if 'pelvis' in forcePath:
            print('need to dial in the pelvis actuators...')
            effort.setWeightForControl(forcePath, 1e-2)
            # if 'pelvis_ty' in forcePath:
            #     effort.setWeightForControl(forcePath, 1e8)
            # if 'hip_rotation' in forcePath:
            #     effort.setWeightForControl(forcePath, 1e4)
        elif 'reserve' in forcePath and 'subtalar' in forcePath:
            effort.setWeightForControl(forcePath, 10)
        elif 'reserve' in forcePath and 'hip_rotation' in forcePath:
            effort.setWeightForControl(forcePath, 10)
        # if 'hip_rotation' in forcePath:
        #    effort.setWeightForControl(forcePath, 10)

    # set up the moment tracking goal
    # test a moment tracking goal from the id moments
    # Add a joint moment tracking goal to the problem.
    jointMomentTracking = osim.MocoGeneralizedForceTrackingGoal('joint_moment_tracking', 50) # type: ignore
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

    wantguess = True
    # set an initial guess up
    if wantguess:
        ## the try does the most recent solution, the except does the 100con solution from previous study.
        
        # twosteptraj = osim.MocoTrajectory('muscle_stateprescribe_grfprescribe_solution.sto')
        try:
            twosteptraj = osim.MocoTrajectory('muscle_statetrack_grfprescribe_solution_redoarms_py.sto')
        except:
            twosteptraj = osim.MocoTrajectory('muscle_statetrack_grfprescribe_solution_100con.sto')
        
        # twosteptraj = osim.MocoTrajectory('thirdpass_IG.sto')    
        steps = twosteptraj.getNumTimes()
        # solver changes. 
        solver = osim.MocoCasADiSolver.safeDownCast(study.updSolver())
        solver.resetProblem(problem)
        solver.set_optim_convergence_tolerance(1e-2)
        solver.set_optim_constraint_tolerance(1e-4)
        # solver.set_optim_max_iterations(3)
        # solver.set_optim_finite_difference_scheme('forward')
        # solver.set_optim_finite_difference_scheme('central')

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


    # solve and visualize
    try:
        solution = study.solve()
        solution.write('muscle_statetrack_grfprescribe_solution_redoarms_py.sto')
        print('ran the base')
    except:
        print(os.getcwd())
        print('could not solve the problem')
        whatfailed[subjectname + '_' + conditionname + '_' + trialname] = os.getcwd()
        solution = solution.unseal()
        solution.write('muscle_statetrack_grfprescribe_solution_unseal_redoarms_py.sto')
        return whatfailed
    
    osim.STOFileAdapter.write(solution.exportToControlsTable(), 'muscletrack_redo_controls_py.sto')
    osim.STOFileAdapter.write(solution.exportToStatesTable(), 'muscletrack_redo_states_py.sto')

    # try:
    #     solution = osim.MocoTrajectory('muscle_statetrack_grfprescribe_solution_redoarms_py.sto')
    # except:
    #     print('could not read the solution')
    #     whatfailed[subjectname + '_' + conditionname + '_' + trialname] = os.getcwd()
    #     return whatfailed
    
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
    ouf.IDplotter(osim.TimeSeriesTable('muscletrack_redo_moments_py.sto'), 'muscletrack_redo', False)

    # add muscle forces to the post analysis
    analyzeStrings_muscleForces = osim.StdVectorString();
    analyzeStrings_muscleForces.append('.*active_fiber_force');
    analyzeStrings_muscleForces.append('.*passive_fiber_force');
    analyzeStrings_muscleForces.append('.*fiber_force_along_tendon');
    table_muscleForces = study.analyze(solution, analyzeStrings_muscleForces);
    osim.STOFileAdapter.write(table_muscleForces, 'muscletrack_redo_muscleforces_py.sto');

    # pdb.set_trace()

    # analyzeStrings_probe = osim.StdVectorString()
    # analyzeStrings_probe.append('/probeset/.*')
    # table_probe = study.analyze(solution, analyzeStrings_probe)
    # osim.STOFileAdapter.write(table_probe, 'muscletrack_redo_probes_py.sto')

    # pdb.set_trace()



    analyzeStrings_all = osim.StdVectorString()
    analyzeStrings_all.append('.*')
    table_all = study.analyze(solution, analyzeStrings_all)
    osim.STOFileAdapter.write(table_all, 'muscletrack_redo_all_py.sto')

    # pdb.set_trace()

    # study.visualize(solution)

    # post analysis etc. 
    # solution1 = osim.MocoTrajectory('muscle_statetrack_grfprescribe_solution_100con.sto')
    # solution2 = osim.MocoTrajectory('muscle_statetrack_grfprescribe_solution_100con_py.sto')

    # likely just use the matlab infrastructure to do the post analysis.
    # otherwise have to update everything to python - not worth the time likely for post analysis. s
    return whatfailed
