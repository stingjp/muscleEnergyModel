import opensim as osim
import numpy as np

model = osim.Model('G:/Shared drives/Exotendon/muscleModel/results/welk005/welknatural/trial01/post_simple_model_all_the_probes_muscletrack_python.osim')
solution = osim.TimeSeriesTable('G:/Shared drives/Exotendon/muscleModel/results/welk005/welknatural/trial01/muscle_statetrack_grfprescribe_solution.sto')
statestable = osim.TimeSeriesTable('G:/Shared drives/Exotendon/muscleModel/results/welk005/welknatural/trial01/muscletrack_states.sto')
controlstable = osim.TimeSeriesTable('G:/Shared drives/Exotendon/muscleModel/results/welk005/welknatural/trial01/muscletrack_controls.sto')
jr = osim.analyzeSpatialVec(model, statestable, controlstable, ['.*walker_knee.*reaction_on_parent.*'])
jr = jr.flatten(['_mx', '_my', '_mz', '_fx', '_fy', '_fz'])
traj = np.empty(jr.getNumRows())
idk = -np.inf
for itime in range(jr.getNumRows()):
    for irxn in range(int(jr.getNumColumns() / 6)):
        fx = jr.getDependentColumnAtIndex(6 * irxn + 3)[itime]
        fy = jr.getDependentColumnAtIndex(6 * irxn + 4)[itime]
        fz = jr.getDependentColumnAtIndex(6 * irxn + 5)[itime]
        norm = np.sqrt(fx**2 + fy**2 + fz**2)
        traj[itime] = norm
        idk = np.max([norm, idk])



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

print(idk)
print(weight)
print(avg)
print(idk/weight)
print(avg/weight)


# print(max / weight, avg / weight)


welk009 natural 1 norm by body weight
idk - 11.459
avg - 3.54
welk009 exo 1 norm buy body weight
idk - 11.66
avg - 4.31



welk005 exo 1
idk - 9.534
avg - 3.29
welk005 nat 1
idk - 10.79
avg - 3.25
