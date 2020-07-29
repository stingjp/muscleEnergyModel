% file to simulate and analyze th metabolics of the motion condition and trial here

metabolicsModelSetup('subject_walk_armless_noprobe.osim');
torqueStateTrackGRFPrescribe();
close all
muscleStatePrescribeGRFPrescribe();
close all
muscleStatePrescribeGRFPrescribeWithEMG();
close all
% analyze and save the metabolics
analyzeMetabolicCost();
close all
disp('success')
