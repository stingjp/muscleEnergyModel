% written by Jon Stingel
% 20210323
import org.opensim.modeling.*
repodir = 'G:\Shared drives\Exotendon\muscleModel\muscleEnergyModel';
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

%%%%% - remember to only put in the exo conditions that you are looking to see the reductions from
% dembconditions = {'dembnoloadfree', 'dembloadedfree'}; %
% dembsubjects = {'demb010','demb011','demb012','demb014', 'demb005','demb007','demb009'}; %
welkexoconditions = {'welkexo'}; % ,'welkexoexo'}; % ,'welknaturalslow','welknaturalnatural', ...
                  % 'welknaturalexo','welkexonatural','welkexoexo','welkexofast'};
welknaturalconditions = {'welknatural'};% ,'welknaturalnatural'};
welksubjects = {'welk002','welk003','welk005','welk008','welk009','welk010','welk013'};
tag = 'muscletrack';

thingstoplot = {'PassiveFiberForce','ActiveFiberForce'};

load 'G:\Shared drives\Exotendon\muscleModel\muscleEnergyModel\subjectgaitcycles.mat';




% loop through the subjects
for subj=1:length(welksubjects)
    subject = char(welksubjects(subj));
    subjdir = strcat(resultsdir, strcat('/',subject));
    
    % loop through each of the things we want to plot
    for thing=1:length(thingstoplot)
        tempthing = char(thingstoplot(thing))
        welknaturalstruct = struct();
        welkexostruct = struct();

        % loop through conditions - exo first
        for cond=1:length(welkexoconditions)
           condition = char(welkexoconditions(cond));
           conddir = strcat(subjdir, strcat('/',condition));
           trials = fieldnames(subjectgaitcycles.(genvarname(subject)).(genvarname(condition)));
            % loop the trials
            for trial=1:length(trials)
                % what do we actually want to do here
                test = char(trials(trial));
                trialdir = strcat(conddir, strcat('/',test));
                cd(trialdir)
                disp(trialdir)

                % now figure out how to get and plot the signal i want
                % have all the muscle analysis files already
                % do I want to do average or individual?
                tempfile = strcat(trialdir, strcat('/analyzemuscles',tag,'_MuscleAnalysis_', strcat(tempthing, '.sto')));
                tempTimeSeriesTable = TimeSeriesTable(tempfile);
                temptime = tempTimeSeriesTable.getIndependentColumn();
                times = zeros(temptime.size(),1);
                for i=0:temptime.size()-1
                    times(i+1) = temptime.get(i);
                end
                timespercent = (times - times(1)) / (times(end) - times(1)) *100;
                timespercent101 = [0:1:100]';
                welkexostruct.time = timespercent101;

                % now for each of the things
                numCols = tempTimeSeriesTable.getNumColumns(); % including time
                labels = tempTimeSeriesTable.getColumnLabels();            

                for i=0:labels.size()-1
                    muscle = char(labels.get(i));
                    tempcol = tempTimeSeriesTable.getDependentColumn(java.lang.String(muscle)).getAsMat();
                    tempcolinterp = interp1(timespercent, tempcol, timespercent101);

                    if ~isfield(welkexostruct, muscle)
                        welkexostruct.(genvarname(muscle)) = [];
                    end
                    welkexostruct.(genvarname(muscle)) = [welkexostruct.(genvarname(muscle)), tempcolinterp];
                end
            end            
        end
        % done with the exo conditions

        % loop through conditions - now for the natural
        for cond=1:length(welknaturalconditions)
           condition = char(welknaturalconditions(cond));
           conddir = strcat(subjdir, strcat('/',condition));
           trials = fieldnames(subjectgaitcycles.(genvarname(subject)).(genvarname(condition)));
            % loop the trials
            for trial=1:length(trials)
                % what do we actually want to do here
                test = char(trials(trial));
                trialdir = strcat(conddir, strcat('/',test));
                cd(trialdir)
                disp(trialdir)

                % now figure out how to get and plot the signal i want
                % have all the muscle analysis files already
                % do I want to do average or individual?
                tempfile = strcat(trialdir, strcat('/analyzemuscles',tag,'_MuscleAnalysis_', strcat(tempthing, '.sto')));
                tempTimeSeriesTable = TimeSeriesTable(tempfile);
                temptime = tempTimeSeriesTable.getIndependentColumn();
                times = zeros(temptime.size(),1);
                for i=0:temptime.size()-1
                    times(i+1) = temptime.get(i);
                end
                timespercent = (times - times(1)) / (times(end) - times(1)) *100;
                timespercent101 = [0:1:100]';
                welknaturalstruct.time = timespercent101;

                % now for each of the things
                numCols = tempTimeSeriesTable.getNumColumns(); % including time
                labels = tempTimeSeriesTable.getColumnLabels();            

                for i=0:labels.size()-1
                    muscle = char(labels.get(i));
                    tempcol = tempTimeSeriesTable.getDependentColumn(java.lang.String(muscle)).getAsMat();
                    tempcolinterp = interp1(timespercent, tempcol, timespercent101);

                    if ~isfield(welknaturalstruct, muscle)
                        welknaturalstruct.(genvarname(muscle)) = [];
                    end
                    welknaturalstruct.(genvarname(muscle)) = [welknaturalstruct.(genvarname(muscle)), tempcolinterp];
                end
            end
        end
        
        tempfig = figure('Position',[1,1,1920,1080]);
        % do more stuff
        % averaging and whatnot
        for i=0:(labels.size()/2)-1
            subplot(5,8,i+1);
            templabel = char(labels.get(i));
            muscleplot = welknaturalstruct.(genvarname(char(templabel)));
            plot(welknaturalstruct.time, muscleplot, ':')
            hold on;
            plot(welknaturalstruct.time, mean(muscleplot,2), 'k-', 'LineWidth', 2)
            title(templabel)
            xlabel('% gait cycle')
            ylabel(tempthing)
            grid on;
        end
        print(tempfig, ...
            strcat(strcat('G:\Shared drives\Exotendon\muscleModel\analysis\', strcat(subject,'\')), strcat(strcat(tempthing, '_natural'), '.png')),...
            '-dpng', '-r500')
        disp('print 1')
        
        
        tempfig2 = figure('Position',[1,1,1920,1080]);
        % do more stuff
        % averaging and whatnot
        for i=0:(labels.size()/2)-1
            subplot(5,8,i+1);
            templabel = char(labels.get(i));
            muscleplot = welkexostruct.(genvarname(char(templabel)));
            plot(welkexostruct.time, muscleplot, ':')
            hold on;
            plot(welkexostruct.time, mean(muscleplot,2), 'k-', 'LineWidth', 2)
            title(templabel)
            xlabel('% gait cycle')
            ylabel(tempthing)
            grid on;
        end
        print(tempfig2, ...
            strcat(strcat('G:\Shared drives\Exotendon\muscleModel\analysis\', strcat(subject,'\')), strcat(strcat(tempthing, '_exo'), '.png')),...
            '-dpng', '-r500')
        disp('print 2')
        
    end
    
end



