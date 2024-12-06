clear all
clc

% repodir = 'G:\Shared drives\Exotendon\muscleModel\muscleEnergyModel';
repodir = 'C:\Users\jonstingel\code\musclemodel\muscleEnergyModel';
resultsdir = strcat(repodir, '/../results');
cd(resultsdir)


% conditions
% walsconditions = ['walsslack','walslow','walsmed','walshigh','walsmax']
% jackconditions = ['jackpower1','jackpower2','jackpower3','jackpower4','jackpower5','jackpower6',
%                   'jacktau1','jacktau2','jacktau3','jacktau4','jacktau5']
% dembconditions = ['dembnoloadfree','dembnoloadslow','dembloadedfree','dembloadedmatched']
% sildconditions = ['sildbw0','sildbw5','sildbw10','sild10w0','sild10w5','sild10w10',
%                   'sild20w0','sild20w5','sild20w10','sild30w0','sild30w5','sild30w10',
%                   'sildbwrun0','sild10wrun0','sild20wrun0','sild30wrun0']

%%%%%
% dembconditions = {'dembnoloadfree', 'dembloadedfree'}; %
% dembsubjects = {'demb010','demb011','demb012','demb014', 'demb005','demb007','demb009'}; %
welkconditions = {'welknatural','welkexo'}; % ,'welknaturalslow','welknaturalnatural', ...
                  % 'welknaturalexo','welkexonatural','welkexoexo','welkexofast'};
% welksubjects = {'welk002','welk003','welk005','welk008','welk009','welk010','welk013'};
welksubjects = {'welk002','welk003','welk005','welk008','welk009','welk010','welk013'};
% load 'G:\Shared drives\Exotendon\muscleModel\muscleEnergyModel\subjectGaitCycles.mat';
load 'C:\Users\jonstingel\code\musclemodel\muscleEnergyModel\subjectGaitCycles.mat';
% an issues holder for this script
% global Issues
Issues = [[java.lang.String('running multiple subjects'); java.lang.String('here we go')]];
% joint coord errors
errnames = [];
errvalus = [];
errstruct = struct;
% pelvis position errors
trannames = [];
tranvalus = [];
transtruct = struct;
% marker errors
marknames = [];
markvalus = [];
markstruct = struct;
% reserve errors
reservenames = [];
maxreservevalus = [];
avgreservevalus = [];
maxreservepercvalus = [];
avgreservepercvalus = [];
% residual errors
residualnames = [];
maxresidualvalus = [];
avgresidualvalus = [];
maxresidualpercvalus = [];
avgresidualpercvalus = [];
% residual moment errors
residualmomnames = [];
maxresidualmomvalus = [];
avgresidualmomvalus = [];
maxresidualmompercvalus = [];
avgresidualmompercvalus = [];


for subj=1:length(welksubjects)
    subject = char(welksubjects(subj));
    subjdir = strcat(resultsdir, strcat('/',subject));
    
    for cond=1:length(welkconditions)
       condition = char(welkconditions(cond));
       conddir = strcat(subjdir, strcat('/',condition));
       trials = fieldnames(subjectgaitcycles.(genvarname(subject)).(genvarname(condition)));
        for trial=1:length(trials)
            test = char(trials(trial));
            trialdir = strcat(conddir, strcat('/',test));
            cd(trialdir)
            % run the analysis
            try
                analyzeSubject()
                % analyzeSubject_setup()
                disp('ran');
            catch
                disp('issue');
                disp(trialdir);
            end
        end
    end
end

% possibly do all the error stuff here
% compute some average values for the different RMSE values
% TODO add in code for averaging different things
% start with marker errors
[max_markerr, idxmax_markerr] = max(markvalus(:));
marknameslong = marknames(:);
max_markerrname = marknameslong(idxmax_markerr)
disp(max_markerr)
mean_markerr = mean(markvalus,2);
avg_markerr = mean(mean_markerr)
% get rid of upper body and redo
marknamesabbrev = marknames(5:end,:);
markvalusabbrev = markvalus(5:end,:);
[max_markerrabbrev, idxmax_markerrabbrev] = max(markvalusabbrev(:));
marknamesabbrevlong = marknamesabbrev(:);
max_markerrnameabbrev = marknamesabbrevlong(idxmax_markerrabbrev)
disp(max_markerrabbrev)
mean_markerrabbrev = mean(markvalusabbrev,2);
avg_markerrabbrev = mean(mean_markerrabbrev)

% now for the joint kinematic errors
[max_coorderr, idxmax_coorderr] = max(errvalus(:));
coordnameslong = errnames(:);
max_coordname = coordnameslong(idxmax_coorderr)
disp(max_coorderr)
mean_coorderr = mean(errvalus,2);
avg_coorderr = mean(mean_coorderr)
% remove lumbar and try again
errnamesabbrev = errnames(1:17,:);
errvalusabbrev = errvalus(1:17,:);
[max_coorderrabbrev, idxmax_coorderrabbrev] = max(errvalusabbrev(:));
coordnamesabbrevlong = errnamesabbrev(:);
max_coordnameabbrev = coordnamesabbrevlong(idxmax_coorderrabbrev)
disp(max_coorderrabbrev)
mean_coorderrabbrev = mean(errvalusabbrev,2);
avg_coorderrabbrev = mean(mean_coorderrabbrev)

% next is the translational coordinates
[max_poserr, idxmax_poserr] = max(tranvalus(:));
posnameslong = trannames(:);
max_posname = posnameslong(idxmax_poserr)
disp(max_poserr)
mean_poserr = mean(tranvalus,2);
avg_poserr = mean(mean_poserr)

% next is reserve errors - raw values max
% have to remove the lumbar actuators
maxreservevalusabbrev = maxreservevalus(1:12,:);
reservenamesabbrev = reservenames(1:12,:);
[max_maxreserveerr, idxmax_maxreserveerr] = max(maxreservevalusabbrev(:));
reservenameslong = reservenamesabbrev(:);
max_maxreservename = reservenameslong(idxmax_maxreserveerr)
disp(max_maxreserveerr)
% reserves errors - percent of net joint moment max
maxreservepercvalusabbrev = maxreservepercvalus(1:12,:);
[max_maxreserveperc, idxmax_maxreserveperc] = max(maxreservepercvalusabbrev(:));
max_maxreservepercname = reservenameslong(idxmax_maxreserveperc)
disp(max_maxreserveperc)
% reserve errors - RMS raw valus
avgreservevalusabbrev = avgreservevalus(1:12,:);
mean_avgreserveerr = mean(avgreservevalusabbrev,2);
avg_avgreserveerr = mean(mean_avgreserveerr)
% reserve errors - RMS perc valus
avgreservepercvalusabbrev = avgreservepercvalus(1:12,:);
mean_avgreservepercerr = mean(avgreservepercvalusabbrev,2);
avg_avgreservepercerr = mean(mean_avgreservepercerr)

% residual errors raw max
[max_maxresidualerr, idxmax_maxresidualerr] = max(maxresidualvalus(:));
residualnameslong = residualnames(:);
max_maxresidualname = residualnameslong(idxmax_maxresidualerr)
disp(max_maxresidualerr)
% residual errors percent max
[max_maxresidualperc, idxmax_maxresidualperc] = max(maxresidualpercvalus(:));
max_maxresidualpercname = residualnameslong(idxmax_maxresidualperc)
disp(max_maxresidualperc)
% residual errors - RMS raw valus
mean_avgresidualerr = mean(avgresidualvalus,2);
avg_avgresidualerr = mean(mean_avgresidualerr)
% residual errors - RMS perc valus
mean_avgresidualpercerr = mean(avgresidualpercvalus,2);
avg_avgresidualpercerr = mean(mean_avgresidualpercerr)

% residual moment errors raw max
[max_maxresidualmomerr, idxmax_maxresidualmomerr] = max(maxresidualmomvalus(:));
residualmomnameslong = residualmomnames(:);
max_maxresidualmomname = residualmomnameslong(idxmax_maxresidualmomerr)
disp(max_maxresidualmomerr)
% residual moment errors percent max
[max_maxresidualmomperc, idxmax_maxresidualmomperc] = max(maxresidualmompercvalus(:));
max_maxresidualmompercname = residualmomnameslong(idxmax_maxresidualmomperc)
disp(max_maxresidualmomperc)
% residual moment errors - RMS raw valus
mean_avgresidualmomerr = mean(avgresidualmomvalus,2);
avg_avgresidualmomerr = mean(mean_avgresidualmomerr)
% residual moment errors - RMS perc valus
mean_avgresidualmompercerr = mean(avgresidualmompercvalus,2);
avg_avgresidualmompercerr = mean(mean_avgresidualmompercerr)


disp('finished all the sims!')
cd(resultsdir)
save('bigissuesfile.mat','Issues');
save('kinematicErrors.mat','errstruct');
disp('end')