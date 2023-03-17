# How connecting the legs with a spring improves human running economy
Experimental data that was collected of runners with and without an exotendon while running at 2.67 m/s. 

##Publications
(pending)

## SimTK project
https://simtk.org/projects/exotendon_sims

## Github
https://github.com/stingjp/muscleEnergyModel/tree/paperSubmission
Note: branch: paperSubmission is specifically for the project on simulating the exotendon. 

The experimental data includes 7 subjects with the following types of data: 
- Motion capture
- ground reaction forces
- Electromyography
- metabolic energy expenditure

The file SubjectInfo_MarkerEditing.xlsx includes several tabs with information on subjects, as well as info that was used in processing. Of note: 
    1) markerediting - includes the frames that were processed from the raw files for each subject and trial. 
    2) SpringStiffness/length - incluedes the leg length and exotendon parameters for each subject. These parameters are computed from a baseline characterization of the material and the method of fabrication. In green are the values used in the OpenSim model. 
    3) physicalparams - includes all demographic info for the participants. 
    - "simulations" and "grfFileExporting" were used more as checklists in processing of data. 

Within each subject folder there are 4 subdirectories:
/MarkerData
    This includes all the processed and edited marker trajectory files for the running and calibration trials that were collected on the 3rd day of the protocol. (.trc)
/metabolics
    This includes the output raw data from the metabolic energy measurments. These files also contain computations of trial averages and comparisons between conditions. There is a seperate file for each of the days in the experimental protocol. (.xlsx)
/OpenSim
    Within this directory there are an additional set of subdirectories:
        /EmgControls - contains subdirectories of processed EMG signals for each of the different trials on the 3rd day of the experimental protocol. these files were generated from the raw EMG signals. (.sto, .mat)
        /ID - contains subdirectories of processed ground reaction force data for each of the trials on the 3rd day of the experimental protocol. (.mot)
        /IK - contains subdirectories of inverse kinematics results and supporting files for each of the motion capture trials on the 3rd day of the experimental protocol, including the static trial. (.mot, .sto, .xml)
        /Scaling - contains a setup file for scaling a generic model in OpenSim, as well as two scaled models for each respective subject - the second scaled model is identical to the first, but with an exotendon added. (.osim, .xml)
/RawAnalog
    This directory contains all of unprocessed analog output files from the 3rd day of the experimental protocol. (.anc) 


