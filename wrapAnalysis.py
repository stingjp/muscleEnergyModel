import os
os.add_dll_directory('C:/Users/jonstingel/opensim-core-4.5.1-2024-08-23-cf3ef35/bin')
import opensim as osim
import OsimUtilityfunctions as ouf
# load the solution
os.chdir('C:/Users/jonstingel/code/muscleModel/results/welk008/welknatural/trial01/')
solution = osim.MocoTrajectory('muscle_statetrack_grfprescribe_solution_redoarms_py.sto')

# call the ID function
ouf.IDplotter(solution, 'muscletrack_redo')



