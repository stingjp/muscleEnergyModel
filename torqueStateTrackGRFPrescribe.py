import os
os.add_dll_directory('C:/Users/jonstingel/opensim-core-4.5.1-2024-08-23-cf3ef35/bin')
import opensim as osim

os.chdir('C:/Users/jonstingel/code/musclemodel/results/welk008/welknatural/trial01/')

# create the tracking problem
track = osim.MocoTrack()
track.setName("torque_markertrack_grfprescribe")

# construct a ModelProcessor and add it to the tool.
modelProcessor = osim.ModelProcessor("simple_model_all_the_probes.osim")
weldem = osim.StdVectorString()
weldem.append('subtalar_r')
weldem.append('mtp_r')
weldem.append('subtalar_l')
weldem.append('mtp_l')
weldem.append('radius_hand_r')
weldem.append('radius_hand_l')
modelProcessor.append(osim.ModOpReplaceJointsWithWelds(weldem))
# add grf
modelProcessor.append(osim.ModOpAddExternalLoads("grf_walk.xml"))
# remove all muscles for torque driven analysis
modelProcessor.append(osim.ModOpRemoveMuscles())
# add CoordinateActuators to the model DOF.
# ignores pelvis coordinates with already have.
modelProcessor.append(osim.ModOpAddReserves(250))
track.setModel(modelProcessor)
# set the states to be tracked
tableProcessor = osim.TableProcessor('results_IK_redoarms.mot')
tableProcessor.append(osim.TabOpLowPassFilter(15))
tableProcessor.append(osim.TabOpUseAbsoluteStateNames())
track.setStatesReference(tableProcessor)
track.set_states_global_tracking_weight(40)
# avoid exceptions if markers in file are no longer in the model (arms removed)
track.set_allow_unused_references(True)
# since there is only coordinate position data in the states references,
# this fills in the missing coordinate speed data using
# the derivative of splined position data
track.set_track_reference_position_derivatives(True)
# set specific weights for the individual weight set
coordinateweights = osim.MocoWeightSet()
coordinateweights.cloneAndAppend(osim.MocoWeight("pelvis_tx", 0.01))
coordinateweights.cloneAndAppend(osim.MocoWeight("pelvis_ty", 0.01))
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

# get the gait timings
gait_start = 67.197
gait_end = 67.965

# gait_start = 2.071
# gait_end = 2.754

# % get the subject name and gait timings
# load 'C:\Users\jonstingel\code\musclemodel\muscleEnergyModel\subjectgaitcycles.mat';
# workdir = pwd;
# [~,trialname,~] = fileparts(pwd);
# cd ../
# [~,conditionname,~] = fileparts(pwd);
# cd ../
# [~,subjectname,~] = fileparts(pwd);
# cd(workdir);
# gait_start = subjectgaitcycles.(genvarname(subjectname)).(genvarname(conditionname)).(genvarname(trialname)).initial;
# gait_end = subjectgaitcycles.(genvarname(subjectname)).(genvarname(conditionname)).(genvarname(trialname)).final;

# set the times and mesh interval, mesh points are computed internally.
track.set_initial_time(gait_start)
track.set_final_time(gait_end)
track.set_mesh_interval(0.01)
# set the study and problem
study = track.initialize()
problem = study.updProblem()
# get reference to the MocoControlGoal that is added to every MocoTrack problem
effort = osim.MocoControlGoal.safeDownCast(problem.updGoal("control_effort"))
effort.setWeight(0.01)
initactivationgoal = osim.MocoInitialActivationGoal("init_activation")
initactivationgoal.setWeight(10)
problem.addGoal(initactivationgoal)

# put large weight on the pelvis CoordinateActuators, which act as the
model = modelProcessor.process()
model.initSystem()
forceSet = model.getForceSet()
for i in range(forceSet.getSize()):
    forcePath = forceSet.get(i).getAbsolutePathString()
    if 'pelvis' in forcePath:
        effort.setWeightForControl(forcePath, 150)
        # if 'pelvis_ty' in forcePath:
        #     effort.setWeightForControl(forcePath, 1e8)
        # if 'hip_rotation' in forcePath:
        #     effort.setWeightForControl(forcePath, 1e4)

# constrain with some periodicity
periodicityGoal = osim.MocoPeriodicityGoal("periodicity")
periodicityGoal.setMode('endpoint_constraint')
for i in range(model.getNumStateVariables()):
    currentStateName = model.getStateVariableNames().getitem(i)
    if 'pelvis_tx/value' not in currentStateName:
        periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair(currentStateName))
forceSet = model.getForceSet()
for i in range(forceSet.getSize()):
    forcePath = forceSet.get(i).getAbsolutePathString()
    periodicityGoal.addControlPair(osim.MocoPeriodicityGoalPair(forcePath))
problem.addGoal(periodicityGoal)    

# test a moment tracking goal from the id moments
# Add a joint moment tracking goal to the problem.
jointMomentTracking = osim.MocoGeneralizedForceTrackingGoal('joint_moment_tracking', 10) # type: ignore
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
jointMomentTracking.setWeightForGeneralizedForcePattern('.*knee.*', 100)
jointMomentTracking.setWeightForGeneralizedForcePattern('.*beta.*', 0)
jointMomentTracking.setWeightForGeneralizedForcePattern('.*ankle.*', 100)
jointMomentTracking.setWeightForGeneralizedForcePattern('.*hip.*', 0)
jointMomentTracking.setWeightForGeneralizedForcePattern('.*lumbar.*', 0)
jointMomentTracking.setWeightForGeneralizedForcePattern('.*arm.*', 0)
jointMomentTracking.setWeightForGeneralizedForcePattern('.*elbow.*', 0)
jointMomentTracking.setWeightForGeneralizedForcePattern('.*pro.*', 0)
jointMomentTracking.setWeightForGeneralizedForcePattern('.*wrist.*', 0)
problem.addGoal(jointMomentTracking)

# solver changes
solver = osim.MocoCasADiSolver.safeDownCast(study.updSolver())
# solver.resetProblem(problem);
solver.set_optim_convergence_tolerance(1e-2) # 1e-2
solver.set_optim_constraint_tolerance(1e-2) # 1e-2
# solver.set_minimize_implicit_auxiliary_derivatives(true);
# solver.set_implicit_auxiliary_derivatives_weight(1e-8);
solver.set_optim_finite_difference_scheme('forward');
solver.set_parameters_require_initsystem(False);

# solve
solution = study.solve()
solution.write('torque_statetrack_grfprescribe_solution_py.sto')
import pdb; pdb.set_trace()
study.visualize(solution)
# end or visualize