% written by Jon Stingel
% 20210329
% gather and plot all the muscle activities from simulations. 
import org.opensim.modeling.*
repodir = 'C:\Users\JP\code\repos\Stanford\delplab\projects\muscleModel\muscleEnergyModel';
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
welkexoconditions = {'welkexo','welkexoexo'}; % ,'welknaturalslow','welknaturalnatural', ...
                  % 'welknaturalexo','welkexonatural','welkexoexo','welkexofast'};
welknaturalconditions = {'welknatural','welknaturalnatural'};
welksubjects = {'welk001'}; % welk002

thingstoplot = {'excitation','activation'};

load 'C:\Users\JP\code\repos\Stanford\delplab\projects\muscleModel\muscleEnergyModel\subjectgaitcycles.mat';




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
                % disp(trialdir)

                % now figure out how to get and plot the signal i want
                % have all the muscle analysis files already
                % do I want to do average or individual?
                
                if tempthing == 'activation'
                    % do something for the activations
                    % disp('getting activations...')

                    tempfile = strcat(trialdir, '/muscleprescribe_states.sto');
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
                    
                    % loop through everything and get the activations.
                    for i=0:labels.size()-1
                        templab = char(labels.get(i));
                        % check if it is an activation for the right leg
                        if string(templab(length(templab)-9:end)) == 'activation' && templab(length(templab)-11) == 'r'
                            % we want the activations - doing right leg for
                            % simplicity
                            tempcol = tempTimeSeriesTable.getDependentColumn(java.lang.String(templab)).getAsMat();
                            tempcolinterp = interp1(timespercent, tempcol, timespercent101);
                            
                            if ~isfield(welkexostruct, templab(11:length(templab)-11))
                                % fix the naming
                                welkexostruct.(genvarname(templab(11:length(templab)-11))) = [];
                            end
                            welkexostruct.(genvarname(templab(11:length(templab)-11))) = [welkexostruct.(genvarname(templab(11:length(templab)-11))), tempcolinterp]; 
                        end
                    end
                end

                % now for the case of excitations
                if tempthing == 'excitation'
                    % do something for the excitations
                    % disp('getting excitations...')
                    
                    tempfile = strcat(trialdir, '/muscleprescribe_controls.sto');
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
                    
                    % loop through everything and get the activations.
                    for i=0:labels.size()-1
                        templab = char(labels.get(i));
                        % check if it is an activation for the right leg                        
                        if string(templab(length(templab)-1:end)) ~= '_l'
                            % we want the controls - doing right leg for
                            % simplicity
                            tempcol = tempTimeSeriesTable.getDependentColumn(java.lang.String(templab)).getAsMat();
                            tempcolinterp = interp1(timespercent, tempcol, timespercent101);
                            
                            if ~isfield(welkexostruct, templab(11:end))
                                % fix the naming
                                welkexostruct.(genvarname(templab(11:end))) = [];
                            end
                            welkexostruct.(genvarname(templab(11:end))) = [welkexostruct.(genvarname(templab(11:end))), tempcolinterp]; 
                        end
                    end                    
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
                % disp(trialdir)

                % now figure out how to get and plot the signal i want
                % have all the muscle analysis files already
                % do I want to do average or individual?
                
                if tempthing == 'activation'
                    % do something for the activations
                    % disp('getting activations...')

                    tempfile = strcat(trialdir, '/muscleprescribe_states.sto');
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
                    
                    % loop through everything and get the activations.
                    for i=0:labels.size()-1
                        templab = char(labels.get(i));
                        % check if it is an activation for the right leg
                        if string(templab(length(templab)-9:end)) == 'activation' && templab(length(templab)-11) == 'r'
                            % we want the activations - doing right leg for
                            % simplicity
                            tempcol = tempTimeSeriesTable.getDependentColumn(java.lang.String(templab)).getAsMat();
                            tempcolinterp = interp1(timespercent, tempcol, timespercent101);
                            
                            if ~isfield(welknaturalstruct, templab(11:length(templab)-11))
                                % fix the naming
                                welknaturalstruct.(genvarname(templab(11:length(templab)-11))) = [];
                            end
                            welknaturalstruct.(genvarname(templab(11:length(templab)-11))) = [welknaturalstruct.(genvarname(templab(11:length(templab)-11))), tempcolinterp]; 
                        end
                    end
                end

                % now for the case of excitations
                if tempthing == 'excitation'
                    % do something for the excitations
                    % disp('getting excitations...')
                    
                    tempfile = strcat(trialdir, '/muscleprescribe_controls.sto');
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
                    
                    % loop through everything and get the activations.
                    for i=0:labels.size()-1
                        templab = char(labels.get(i));
                        % check if it is an activation for the right leg                        
                        if string(templab(length(templab)-1:end)) ~= '_l'
                            % we want the controls - doing right leg for
                            % simplicity
                            tempcol = tempTimeSeriesTable.getDependentColumn(java.lang.String(templab)).getAsMat();
                            tempcolinterp = interp1(timespercent, tempcol, timespercent101);
                            
                            if ~isfield(welknaturalstruct, templab(11:end))
                                % fix the naming
                                welknaturalstruct.(genvarname(templab(11:end))) = [];
                            end
                            welknaturalstruct.(genvarname(templab(11:end))) = [welknaturalstruct.(genvarname(templab(11:end))), tempcolinterp]; 
                        end
                    end                    
                end
            end
        end
        
        
        tempfig = figure('Position',[1,1,1920,1080]);
        % do more stuff
        % averaging and whatnot
        labels = fields(welknaturalstruct);
        for i=2:length(labels)
            subplot(7,8,i-1);
            templabel = char(labels(i));
            muscleplot = welknaturalstruct.(genvarname(templabel));
            plot(welknaturalstruct.time, muscleplot, ':')
            hold on;
            plot(welknaturalstruct.time, mean(muscleplot,2), 'k-', 'LineWidth', 2)
            title(templabel)
            xlabel('% gait cycle')
            ylabel(tempthing)
            grid on;
        end
        print(tempfig, ...
            strcat(strcat('C:\Users\JP\code\repos\Stanford\delplab\projects\muscleModel\analysis\', strcat(subject,'\')), strcat(strcat(tempthing, '_natural'), '.png')),...
            '-dpng', '-r500')
        disp('print 1')
        
        
        tempfig2 = figure('Position',[1,1,1920,1080]);
        % do more stuff
        % averaging and whatnot
        labels = fields(welkexostruct);
        for i=2:length(labels)
            subplot(7,8,i-1);
            templabel = char(labels(i));
            muscleplot = welkexostruct.(genvarname(templabel));
            plot(welkexostruct.time, muscleplot, ':')
            hold on;
            plot(welkexostruct.time, mean(muscleplot,2), 'k-', 'LineWidth', 2)
            title(templabel)
            xlabel('% gait cycle')
            ylabel(tempthing)
            grid on;
        end
        print(tempfig2, ...
            strcat(strcat('C:\Users\JP\code\repos\Stanford\delplab\projects\muscleModel\analysis\', strcat(subject,'\')), strcat(strcat(tempthing, '_exo'), '.png')),...
            '-dpng', '-r500')
        disp('print 2')
        
    end
    
end
