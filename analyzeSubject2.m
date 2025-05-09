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
% print out issues to look at for each subject:

tic

import org.opensim.modeling.*
workingdir = pwd;
[~,trialname,~] = fileparts(pwd);
cd ../
[~,condname,~] = fileparts(pwd);
cd ../
[~,subjectname,~] = fileparts(pwd);
experimentname = subjectname(1:4);
cd(workingdir)


if exist('Issues','var')
    Issues = [Issues; [java.lang.String(subjectname), java.lang.String(strcat(condname,trialname))]];
else
%     global Issues
    Issues = [java.lang.String('coordinate actuator'), java.lang.String('ratio to net')];
end

% edit experimental data for simulations
% renameExperimentalData();

% run simulations of the subject, and get metabolic cost of motion
% close all;
% metabolicsModelSetup();
% close all;
% torqueStateTrackGRFPrescribe();
% close all;
% torqueStatePrescribeGRFPrescribe();
% close all;
torqueMarkerTrackGRFPrescribe();
close all;
muscleStatePrescribeGRFPrescribe_2();
close all;
% torqueStateTrackGRFTrack();
% close all;
% torqueMarkerTrackGRFTrack();
% close all;
% Issues = muscleStatePrescribeGRFPrescribe(Issues);
% close all;
% Issues = muscleStatePrescribeGRFPrescribeWithEMG(Issues);
% close all;
% Issues = muscleStateTrackGRFPrescribe_firstPass(Issues);
% close all;
% Issues = muscleStateTrackGRFPrescribe_secondpass(Issues);
% close all;
% Issues = muscleStateTrackGRFPrescribe_thirdpass(Issues);
% close all;
% Issues = muscleStateTrackGRFPrescribe(Issues); % set up as secondpass with different mesh
% close all;
% Issues = muscleStateTrackGRFTrack(Issues)s
close all;

%%% only uncomment this if the above simulations are commented out
% this will load the existing solutions and perform the post analyses
solution1 = MocoTrajectory('muscle_stateprescribe_grfprescribe_solution.sto');
analyzeMetabolicCost(solution1, 'muscleprescribe');
% Issues = computeIDFromResult(Issues, solution1, tag);

% solution2 = MocoTrajectory('muscle_stateprescribe_grfprescribe_withemg_solution.sto');
% solution1 = MocoTrajectory('muscle_statetrack_grfprescribe_solution_100con_rra.sto');

% 
% Issues = [Issues; [java.lang.String('muscledrivensim'), java.lang.String('inverseproblem')]];

% analyzeMetabolicCost(solution1);
% trackorprescribe = 'prescribe';
% computeKinematicDifferences(solution1, trackorprescribe);
% Issues = [Issues; [java.lang.String('muscledrivensimwithEMG'); java.lang.String('inverseproblem')]];
% Issues = computeIDFromResult(Issues, solution2);
% analyzeMetabolicCostWithEMG(solution2);

toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% finisher
disp('finished the subject-condition-trial');
% Issues
save('issuesfile.mat','Issues');
disp('end this subject-condition-trial')
