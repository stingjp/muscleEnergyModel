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
if ~exist('errnames','var')
    errnames = [];
    trannames = [];
    marknames = [];
end
if ~exist('errvalus','var')
    errvalus = [];
    tranvalus = [];
    markvalus = [];
end
if ~exist('errstruct','var')
    errstruct = struct;
    transtruct = struct;
    markstruct = struct;
end
if ~exist('residualnames','var')
    residualnames = [];
    maxresidualvalus = [];
    avgresidualvalus = [];
    maxresidualpercvalus = [];
    avgresidualpercvalus = [];
end
if ~exist('residualmomnames','var')
    residualmomnames = [];
    maxresidualmomvalus = [];
    avgresidualmomvalus = [];
    maxresidualmompercvalus = [];
    avgresidualmompercvalus = [];
end
if ~exist('reservenames','var')
    reservenames = [];
    maxreservevalus = [];
    avgreservevalus = [];
    maxreservepercvalus = [];
    avgreservepercvalus = [];
end



% edit experimental data for simulations
% renameExperimentalData();

% run simulations of the subject, and get metabolic cost of motion
% close all;
% metabolicsModelSetup();
% close all;
% runRRA_1('./RRAfiles/RRA_Setup_1.xml');  % this is not tested, run manually
% close all;
% runRRA_2('./RRAfiles/RRA_Setup_2.xml');  % this is not tested, run manually
% close all;
% torqueMarkerTrackGRFPrescribe();
% close all;
% torqueStateTrackGRFPrescribe();
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
% Issues = muscleStateTrackGRFPrescribe(Issues);
% close all;

%%% only uncomment this if the above simulations are commented out
% this will load the existing solutions and perform the post analyses
% solution1 = MocoTrajectory('muscle_stateprescribe_grfprescribe_solution.sto');
% solution2 = MocoTrajectory('muscle_stateprescribe_grfprescribe_withemg_solution.sto');
% solution1 = MocoTrajectory('muscle_statetrack_grfprescribe_solution.sto');

disp('only run the rest of the file when the simulations have completed. Code following will analyze the solutions.')
keyboard

tag = 'muscleprescribe';
tag = 'muscletrack';

if strcmp(tag, 'muscletrack')
    solution2 = MocoTrajectory('muscle_statetrack_grfprescribe_solution.sto');
    conornot = true;
else
    solution2 = MocoTrajectory('muscle_stateprescribe_grfprescribe_solution.sto');
    conornot = false;
end

% try
%     solution2 = MocoTrajectory('muscle_statetrack_grfprescribe_solution_100con.sto');
%     conornot = true;
% catch
%     disp('NO 100 CON SOLUTION');
%     solution2 = MocoTrajectory('muscle_statetrack_grfprescribe_solution.sto');
%     conornot = false;
% end

% 
% if strcmp(subjectname, 'welk002') %|| strcmp(subjectname, 'welk003')
%     solution2 = MocoTrajectory('muscle_statetrack_grfprescribe_solution.sto');
%     conornot = false;
% else
%     solution2 = MocoTrajectory('muscle_statetrack_grfprescribe_solution_100con.sto');
%     conornot = true;
% end
    
disp('switching to the 100 con solution')
% analyze the tracking simulations 
Issues = [Issues; [java.lang.String('muscledrivensim'), java.lang.String('inverseproblem')]];




% main tracking solution analysis
%
% analyzeMetabolicCost(solution1, 'muscletrack');
% here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analyzeMetabolicCostSecond(solution2, 'muscletrack');
[Issues, maxreservepercvalus_new, avgreservepercvalus_new, maxreservevalus_new, avgreservevalus_new, reservenames_new, ...
    maxresidualpercvalus_new, avgresidualpercvalus_new, maxresidualvalus_new, avgresidualvalus_new, residualnames_new, ...
    maxresidualmompercvalus_new, avgresidualmompercvalus_new, maxresidualmomvalus_new, avgresidualmomvalus_new, residualmomnames_new...
    ] = computeIDFromResult(Issues, solution2, 'muscletrack');

% add to combined arrays 
% reserves first
reservenames = [reservenames, reservenames_new];
maxreservevalus = [maxreservevalus, maxreservevalus_new];
avgreservevalus = [avgreservevalus, avgreservevalus_new];
maxreservepercvalus = [maxreservepercvalus, maxreservepercvalus_new];
avgreservepercvalus = [avgreservepercvalus, avgreservepercvalus_new];
% now residuals
residualnames = [residualnames, residualnames_new];
maxresidualvalus = [maxresidualvalus, maxresidualvalus_new];
avgresidualvalus = [avgresidualvalus, avgresidualvalus_new];
maxresidualpercvalus = [maxresidualpercvalus, maxresidualpercvalus_new];
avgresidualpercvalus = [avgresidualpercvalus, avgresidualpercvalus_new];
% now residual moments
residualmomnames = [residualmomnames, residualmomnames_new];
maxresidualmomvalus = [maxresidualmomvalus, maxresidualmomvalus_new];
avgresidualmomvalus = [avgresidualmomvalus, avgresidualmomvalus_new];
maxresidualmompercvalus = [maxresidualmompercvalus, maxresidualmompercvalus_new];
avgresidualmompercvalus = [avgresidualmompercvalus, avgresidualmompercvalus_new];


% computing the kinematic differences
trackorprescribe = 'track';
% computeKinematicDifferences(solution2, trackorprescribe);

% joint coordinate errors
[errnames_new, errvalus_new, errstruct_new, trannames_new, tranvalus_new, transtruct_new] = computeKinematicRMSE( ... 
    solution2, trackorprescribe);
errnames = [errnames, errnames_new];
errvalus = [errvalus, errvalus_new];
errfields = fields(errstruct_new);
for i=1:length(fields(errstruct_new))
    if ~isfield(errstruct, errfields{i})
        errstruct.(errfields{i}) = [errstruct_new.(errfields{i})];
    else
        errstruct.(errfields{i}) = [errstruct.(errfields{i}), errstruct_new.(errfields{i})];
    end
end
% pelvis translation errors
trannames = [trannames, trannames_new];
tranvalus = [tranvalus, tranvalus_new];
tranfields = fields(transtruct_new);
for i=1:length(fields(transtruct_new))
    if ~isfield(transtruct, tranfields{i})
        transtruct.(tranfields{i}) = [transtruct_new.(tranfields{i})];
    else
        transtruct.(tranfields{i}) = [transtruct.(tranfields{i}), transtruct_new.(tranfields{i})];
    end
end
% get the marker errors
[marknames_new, markvalus_new, markstruct_new] = computeMarkerRMSE( ... 
    solution2, trackorprescribe, conornot);
if strcmp(subjectname, 'welk002') || strcmp(subjectname, 'welk003')
    % don't want the two sternum markers
    % dont want the exotendon markers. 
    marknames_new = marknames_new(1:37);
    markvalus_new = markvalus_new(1:37);

    temp1 = marknames_new(1:3);
    temp1 = [temp1; 'Sternum'];
    temp10 = markvalus_new(1:3);
    temp10 = [temp10; mean(markvalus_new(4:5))];

    temp2 = marknames_new(6:end);
    temp2(1) = string("r.ASIS");
    temp2(2) = string('L.ASIS');
    temp2(3) = string('r.PSIS');
    temp2(4) = string('L.PSIS');
    marknames_new = [temp1;temp2];

    temp20 = markvalus_new(6:end);
    markvalus_new = [temp10; temp20];
end


marknames = [marknames, marknames_new];
markvalus = [markvalus, markvalus_new];
markfields = fields(markstruct_new);
for i=1:length(fields(markstruct_new))
    if ~isfield(markstruct, markfields{i})
        markstruct.(markfields{i}) = [markstruct_new.(markfields{i})];
    else
        markstruct.(markfields{i}) = [markstruct.(markfields{i}), markstruct_new.(markfields{i})];
    end
end


%}

% try
%     solution2 = MocoTrajectory('muscle_statetrack_grfprescribe_solution_100con.sto');
%     % 100 con tracking solution analysis
%     analyzeMetabolicCostSecond(solution2,'muscletrack');
%     % Issues = computeIDFromResult(Issues, solution2, 'muscletrack');
% catch
%     disp('must be subj 2 or 3')
% end
% solution3 = MocoTrajectory('muscle_stateprescribe_grfprescribe_solution.sto');
% prescribe inverse solution analysis
% analyzeMetabolicCost(solution3, 'muscleprescribe');
% Issues = computeIDFromResult(Issues, solution3, 'muscleprescribe');


%{
% analyzeMetabolicCost(solution1, 'muscletrack');
trackorprescribe = 'track';
% joint coordinate errors
[errnames_new, errvalus_new, errstruct_new, trannames_new, tranvalus_new, transtruct_new] = computeKinematicRMSE( ... 
    solution1, trackorprescribe);
errnames = [errnames, errnames_new];
errvalus = [errvalus, errvalus_new];
errfields = fields(errstruct_new);
for i=1:length(fields(errstruct_new))
    if ~isfield(errstruct, errfields{i})
        errstruct.(errfields{i}) = [errstruct_new.(errfields{i})];
    else
        errstruct.(errfields{i}) = [errstruct.(errfields{i}), errstruct_new.(errfields{i})];
    end
end
% pelvis translation errors
trannames = [trannames, trannames_new];
tranvalus = [tranvalus, tranvalus_new];
tranfields = fields(transtruct_new);
for i=1:length(fields(transtruct_new))
    if ~isfield(transtruct, tranfields{i})
        transtruct.(tranfields{i}) = [transtruct_new.(tranfields{i})];
    else
        transtruct.(tranfields{i}) = [transtruct.(tranfields{i}), transtruct_new.(tranfields{i})];
    end
end
% get the marker errors
[marknames_new, markvalus_new, markstruct_new] = computeMarkerRMSE( ... 
    solution1, trackorprescribe);
marknames = [marknames, marknames_new];
markvalus = [markvalus, markvalus_new];
markfields = fields(markstruct_new);
for i=1:length(fields(markstruct_new))
    if ~isfield(markstruct, markfields{i})
        markstruct.(markfields{i}) = [markstruct_new.(markfields{i})];
    else
        markstruct.(markfields{i}) = [markstruct.(markfields{i}), markstruct_new.(markfields{i})];
    end
end

% compute some average values for the different RMSE values
%}


% Issues = [Issues; [java.lang.String('muscledrivensimwithEMG'); java.lang.String('inverseproblem')]];
% Issues = computeIDFromResult(Issues, solution2);
% analyzeMetabolicCostWithEMG(solution2);

% Issues = [Issues; [java.lang.String('muscledrivensim'), java.lang.String('inverseproblem')]];
% analyzeMetabolicCost(solution2, 'muscleprescribe');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% finisher
disp('finished the subject-condition-trial');
% Issues
save('issuesfile.mat','Issues');
disp('end this subject-condition-trial')