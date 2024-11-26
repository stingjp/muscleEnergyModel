import os
os.add_dll_directory('C:/Users/jonstingel/opensim-core-4.5.1-2024-08-23-cf3ef35/bin')
import opensim as osim
import OsimUtilityfunctions as ouf
os.chdir('C:/Users/jonstingel/code/musclemodel/results/welk008/welknatural/trial01/')

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

# get the subject name and gait timings
    # % load 'G:\Shared drives\Exotendon\muscleModel\muscleEnergyModel\subjectgaitcycles.mat';
    # load 'C:\Users\jonstingel\code\muscleModel\muscleEnergyModel\subjectgaitcycles.mat';
    # workdir = pwd;
    # [~,trialname,~] = fileparts(pwd);
    # cd ../
    # [~,conditionname,~] = fileparts(pwd);
    # cd ../
    # [~,subjectname,~] = fileparts(pwd);
    # cd(workdir);

    # gait_start = subjectgaitcycles.(genvarname(subjectname)).(genvarname(conditionname)).(genvarname(trialname)).initial;
    # gait_end = subjectgaitcycles.(genvarname(subjectname)).(genvarname(conditionname)).(genvarname(trialname)).final;

gait_start = 67.197
gait_end = 67.965

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
twosteptraj = osim.MocoTrajectory('muscle_statetrack_grfprescribe_solution.sto')    
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
ouf.IDplotter(osim.TimeSeriesTable('muscletrack_moments_py.sto'), 'muscletrack')


study.visualize(solution)

# post analysis etc. 
# solution1 = osim.MocoTrajectory('muscle_statetrack_grfprescribe_solution_100con.sto')
# solution2 = osim.MocoTrajectory('muscle_statetrack_grfprescribe_solution_100con_py.sto')


# likely just use the matlab infrastructure to do the post analysis.
# otherwise have to update everything to python - not worth the time likely for post analysis. s
