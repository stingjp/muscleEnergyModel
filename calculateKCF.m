import org.opensim.modeling.*
opensimSimulation.analyzeSpatialVec()
model = Model('post_simple_model_all_the_probes_muscletrack_python.osim');
statestable = TimeSeriesTable('muscletrack_states.sto');
controlstable = TimeSeriesTable('muscletrack_controls.sto');
% don't know how to figure out the output paths

test = StdVectorString();
test.add(java.lang.String('/jointset/walker_knee_r|reaction_on_parent_1'))
test.add(java.lang.String('/jointset/walker_knee_r|reaction_on_parent_2'))
test.add(java.lang.String('/jointset/walker_knee_r|reaction_on_parent_3'))
test.add(java.lang.String('/jointset/walker_knee_r|reaction_on_parent_4'))
test.add(java.lang.String('/jointset/walker_knee_r|reaction_on_parent_5'))
test.add(java.lang.String('/jointset/walker_knee_r|reaction_on_parent_6'))
test.add(java.lang.String('/jointset/walker_knee_l|reaction_on_parent_1'))
test.add(java.lang.String('/jointset/walker_knee_l|reaction_on_parent_2'))
test.add(java.lang.String('/jointset/walker_knee_l|reaction_on_parent_3'))
test.add(java.lang.String('/jointset/walker_knee_l|reaction_on_parent_4'))
test.add(java.lang.String('/jointset/walker_knee_l|reaction_on_parent_5'))
test.add(java.lang.String('/jointset/walker_knee_l|reaction_on_parent_6'))

jr = analyzeSpatialVec(model, statestable, controlstable, test)

def calc_knee_reaction_force(self, root_dir, solution):
        modelProcessor = self.create_model_processor(root_dir)
        model = modelProcessor.process()
        jr = osim.analyzeSpatialVec(model, solution,
                                    ['.*walker_knee.*reaction_on_parent.*'])
                                    
        StdVectorString(java.lang.String('/jointset/walker_knee_r/reaction_on_parent'))
                                    
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