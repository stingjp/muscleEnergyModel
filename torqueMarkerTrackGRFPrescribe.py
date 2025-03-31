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
# set the mocotrack markers reference
track.setMarkersReferenceFromTRC("motion_capture.trc")
track.set_allow_unused_references(True)
track.set_markers_global_tracking_weight(100)
# set individual markers weights
markerWeights = osim.MocoWeightSet()
markerWeights.cloneAndAppend(osim.MocoWeight("R.ASIS", 20))
markerWeights.cloneAndAppend(osim.MocoWeight("L.ASIS", 20))
markerWeights.cloneAndAppend(osim.MocoWeight("R.PSIS", 20))
markerWeights.cloneAndAppend(osim.MocoWeight("L.PSIS", 20))
markerWeights.cloneAndAppend(osim.MocoWeight("R.Knee", 15))
markerWeights.cloneAndAppend(osim.MocoWeight("R.Ankle", 100))
markerWeights.cloneAndAppend(osim.MocoWeight("R.Heel", 100))
markerWeights.cloneAndAppend(osim.MocoWeight("R.MT5", 100))
markerWeights.cloneAndAppend(osim.MocoWeight("R.Toe", 100))
markerWeights.cloneAndAppend(osim.MocoWeight("L.Knee", 15))
markerWeights.cloneAndAppend(osim.MocoWeight("L.Ankle", 100))
markerWeights.cloneAndAppend(osim.MocoWeight("L.Heel", 100))
markerWeights.cloneAndAppend(osim.MocoWeight("L.MT5", 100))
markerWeights.cloneAndAppend(osim.MocoWeight("L.Toe", 100))
track.set_markers_weight_set(markerWeights)
# get the subject name and gait timings
gait_start = 67.197
gait_end = 67.965

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
# add control effort to the problem
effort = osim.MocoControlGoal.safeDownCast(problem.updGoal("control_effort"))
effort.setWeight(1)
# add initial activation goal
initactivationgoal = osim.MocoInitialActivationGoal("init_activation")
initactivationgoal.setWeight(10)
problem.addGoal(initactivationgoal)
# put large weight on the pelvis CoordinateActuators, which act as the
# residual, or 'hand-of-god' forces which we would like to keep small
model = modelProcessor.process()
model.initSystem()
forceSet = model.getForceSet()
for i in range(forceSet.getSize()):
    forcePath = forceSet.get(i).getAbsolutePathString()
    if 'pelvis' in forcePath:
        effort.setWeightForControl(forcePath, 5000)

# solver changes
solver = osim.MocoCasADiSolver.safeDownCast(study.updSolver())
# solver.resetProblem(problem)
solver.set_optim_convergence_tolerance(1e-2) # 1e-2
solver.set_optim_constraint_tolerance(1e-2) # 1e-2
# solver.set_minimize_implicit_auxiliary_derivatives(true)
# solver.set_implicit_auxiliary_derivatives_weight(1e-8)
solver.set_optim_finite_difference_scheme('forward')
solver.set_parameters_require_initsystem(False)
# solve - the bool indicates to visualize the solution
solution = study.solve() # to visualize
solution.write('torque_markertrack_grfprescribe_solution_py.sto')

import pdb; pdb.set_trace()
study.visualize(solution)
# end or visualize