# How connecting the legs with a spring improves human running economy
Experimental data that was collected of runners with and without an exotendon while running at 2.67 m/s. The 3D musculoskeletal simulations utilize an optimal control problem implementation with OpenSim Moco.

## Publications
https://pubmed.ncbi.nlm.nih.gov/37745177/

## SimTK project
https://simtk.org/projects/exotendon_sims

There are 2 data share repositories here 1) with raw/processed experimental data, and 2) the results directory and files used when executing the code in this repository.

The experimental data includes 7 subjects with the following types of data: 
- Motion capture
- ground reaction forces
- Electromyography
- metabolic energy expenditure

All the needed experimental data files are contained within the code results directory as well. 
Within each subject folder the subdirectories are:
- /welkexo - Exotendon running condition
- /welknatural - Natural running condition

Within each of these are subdirectories for each trial, or gait cycle simulated. Each of those then contains the results and files for the simulation of that gait cycle, /RRAfiles (which contains the RRA setup and results), /expdata (which contains all necessary experimental data and model files for that subject and gait cycle. 


## Github
https://github.com/stingjp/muscleEnergyModel/tree/paperSubmission
Note: branch: paperSubmission is specifically for the project on simulating the exotendon. 


## Running the code
The majority of the codebase was written in matlab 2021b. The python files in the repository were written in python 3.7.3.

### Setting up the file directories (optional - unnecessary when downloading the results from SimTK):
1. To set up the results directories call ./createResultsDirectories.py
2. Call ./analyzeSubject_setup.m - be sure to uncomment renameExperimentalData(), metabolicsModelSetup(). 
3. At this point I manually ran RRA in Opensim (https://opensim.stanford.edu/) in order to adjust model masses. Each subject's manually adjusted model is ./simple_model_all_the_probes_adjusted.osim. The adjustments/averages can be found in ./RRAMassMods.xlsx. 


### Running simulations
#### Within a subject-condition-trial: 

analyzeSubject.m allows you to run multiple simulations/analyses on a single subject. You can customize what simulations and analyses you would like to run. Options and desctiptions follow, starting with the ones used for paper results:
- torqueStateTrackGRFPrescribe.m -> formulates a torque-driven state tracking (IK results) problem with prescribed GRFs. 
- muscleStateTrackGRFPrescribe.m -> formulates a muscle-driven state tracking problem with prescribed GRFs. 
- muscleStateTrackGRFPrescribe_firstpass.m, muscleStateTrackGRFPrescribe_secondpass.m, muscleStateTrackGRFPrescribe_thirdpass.m -> these are the same as muscleStateTrackGRFPrescribe.m ^ but have approximate weights that may aid in generating progressively better results. The secondpass script is currently set up to provide the final results from the simulations. 

#### To batch process subjects, conditions, or trials, call runsubjects1.m (edit for which subjects, conditions, trials to run). This wrapper script can call analyzeSubject.m or analyzeSubject_setup.m - be sure to specify and set up that file for what analyses and simulations you would like to complete. 


#### Post-simulation analysis can be run after the simulation solves within the same analyzeSubject script. (below)
- analyzeMetabolicCost.m -> takes the simulated solution and performs all the metabolics computations. Computes average metaboalic cost of the gait cycle, stance and swing costs, as well as individual muscle costs. 
- computeIDFromResults.m -> takes the solution, and performs an analysis of the dynamics. It will return instances where reserves and residuals exceed the recommended thresholds by (Hicks et al. 2015). 
- computeKinematicDifferences.m -> takes the solution and plots differences between the dynamically consistent solution and the input kinematics. 
- computeMarkerRMSE.m -> takes the solution and computes the output marker trajectory RMSE from the experimental markers. 
- Throughout the analysis there is a data structure called 'Issues' that keeps track of any flags that the analysis files returns, and what they are. This is helpful when batching many simulations and subjects.


#### Other post-simulation analyses (often run on multiple subject results):
- plottingJoinMoments.m -> searches through subject results and plots combined and individual net joint moments for both natural and exotendon conditions. 
- plottingJointCoordinates.m -> searches through subject results and plots combined and individual joint coordinate trajectories for both conditions.
- plottingGRF.m -> plots the GRF and moments, COP, and external forces including the force in the exotendon. 
- plottingCOM.m -> plots the COM trajectory for each subject in both conditions. 
- activityPlots.m -> plots the muscle excitations and activations from the simulations. 
- combiningEMGandActivity.m -> plots weighted averaged muscle activations from simulations in comparison with EMG signals that were collected. 
- findMuscleSavers.m -> computes per-muscle metabolic savings and determines which are the greatest savesrs. 
- muscleFiberForcesPlotting.m -> plots active and passive muscle fiber forces throughout the gait cycle in both conditions. 
- muscleLengthsPlotting.m -> plots the muscle fiber length and velocity in both conditions across the gait cycle. 
- muscleMetabolicContinuousPlotting.m -> Options to plot the full or breakdown the per-muscle metabolic rate throughout the gait cycle for both conditions. 
- muscleStatePlotting.m -> in-depth breakdown of the per-muscle metabolic rate. 
- plotting_stancevsswing.ipynb -> Produces paper figure of breakdown between metabolic rate in stance and swing. (requires running collectMetabolicResults.py)
- jointCoordinateVariance.m -> Computes the variance in each of the joint coordinate values in the model based on the kinematics results. 
- muscleStatePlotting.m -> does a low-level breakdown of plots for each subject and each muscle (time costly). 
- collectMetabolicsResults.py -> searches the directory and copies all the metabolic results to a single location. 
- plotting_stancevsswing.ipynb -> creates the paper figure splitting metabolics between stance and swing phases. 
- combiningMuscleGroups_metabolics.xlsx -> contains all the average muscle metabolic breakdowns across subjects. Combines muscles into functional groups and contains the paper plot of the muscle group costs. 
metaboliccostcomparisons_rev2.xlsx -> contians the final whole body metabolic cost comparisons. Also contains tabs for statistics performed on the whole body, stance, swing, etc. 

    
#### There are other scripts and code in the repository, but were not used in the final analysis of this project. Therefore, their functionality is not tested. 
