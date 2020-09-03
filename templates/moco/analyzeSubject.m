%{ 
jon stingel - 08/06/2020
Ensure that all the functions here are the main file
directory that is added to matlab path. That way 
they can be called from each subject condition 
results directory. The scripts utilize the file structure
to extract data. 

TODO:
making a script to set up the file structure...
%}

% edit experimental data for simulations
renameExperimentalData();

% run simulations of the subject, and get metabolic cost of motion
close all;
metabolicsModelSetup();
close all;
% torqueMarkerTrackGRFPrescribe();
close all;
torqueStateTrackGRFPrescribe();
close all;
muscleStatePrescribeGRFPrescribe();
close all;
muscleStatePrescribeGRFPrescribeWithEMG();
close all;

disp('success');
