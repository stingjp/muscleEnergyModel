% written by Jon Stingel
% 20211021
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
welknaturalconditions = {'welknatural'};% ,'welknaturalnatural'};exit
% welksubjects = {'welk001','welk002','welk003','welk004'};
% welksubjects = {'welk005','welk007','welk008','welk009','welk010','welk013'};
welksubjects = {'welk002','welk003','welk005','welk008','welk009','welk010','welk013'};
% welksubjects = {'welk005'};
thingstoplot = {'externalforces'}; % 'probes', 'shortening', 'mechanical', 'activation'

load 'G:\Shared drives\Exotendon\muscleModel\muscleEnergyModel\subjectgaitcycles.mat';
load 'G:\Shared drives\Exotendon\muscleModel\muscleEnergyModel\subjectmass.mat';
tag = 'muscleprescribe'


% loop through each of the things we want to plot
for thing=1:length(thingstoplot)
    tempthing = char(thingstoplot(thing))

    % create stucture for combined subject figures
    welknaturalstruct_combine = struct();
    welkexostruct_combine = struct();

    allsubjmeans_nat = struct();
    allsubjmeans_exo = struct();

    meansall_exo = struct();
    meansall_nat = struct();
    stdsall_exo = struct();
    stdsall_nat = struct();
    
    allsubj_nat = struct();
    allsubj_exo = struct();

    % loop through the subjects
    for subj=1:length(welksubjects)
        subject = char(welksubjects(subj));
        subjdir = strcat(resultsdir, strcat('/',subject));
        
        % create the struct for individual figures
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
                tempfile = strcat(trialdir, '/analyzemuscles',tag,'_ForceReporter_forces', '.sto');
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
                    coord = char(labels.get(i));
                    
                    % need to screen only the things that we want
                    if contains(char(coord), 'GRF_') % 'value' denotes a coordinate value
                        % tempsplit = split(coord,'/');
                        % coordshort = string(tempsplit(4));
                        
                        tempcol = tempTimeSeriesTable.getDependentColumn(java.lang.String(coord)).getAsMat();
                        tempcolinterp = interp1(timespercent, tempcol, timespercent101);
                        if ~isfield(welkexostruct, coord)
                            welkexostruct.(genvarname(coord)) = [];
                        end
                        welkexostruct.(genvarname(coord)) = [welkexostruct.(genvarname(coord)), tempcolinterp];
                    elseif contains(char(coord), 'HOBL')
                        tempcol = tempTimeSeriesTable.getDependentColumn(java.lang.String(coord)).getAsMat();
                        tempcolinterp = interp1(timespercent, tempcol, timespercent101);
                        if ~isfield(welkexostruct, coord)
                             welkexostruct.(genvarname(coord)) = [];
                        end
                        welkexostruct.(genvarname(coord)) = [welkexostruct.(genvarname(coord)), tempcolinterp];
                    
                    
                    % elseif contains(char(muscle), 'metabolics_combined')
                    %     % we want these measures
                    %     tempcol = tempTimeSeriesTable.getDependentColumn(java.lang.String(muscle)).getAsMat();
                    %     tempcolinterp = interp1(timespercent, tempcol, timespercent101);
                    %     if ~isfield(welkexostruct, muscle)
                    %         welkexostruct.(genvarname(muscle)) = [];
                    %     end
                    %     welkexostruct.(genvarname(muscle)) = [welkexostruct.(genvarname(muscle)), tempcolinterp];
                    % elseif contains(char(muscle), 'activation_maintenance_rate_')
                    %     % we want these measures
                    %     tempcol = tempTimeSeriesTable.getDependentColumn(java.lang.String(muscle)).getAsMat();
                    %     tempcolinterp = interp1(timespercent, tempcol, timespercent101);
                    %     if ~isfield(welkexostruct, muscle)
                    %         welkexostruct.(genvarname(muscle)) = [];
                    %     end
                    %     welkexostruct.(genvarname(muscle)) = [welkexostruct.(genvarname(muscle)), tempcolinterp];
                    % elseif contains(char(muscle), 'shortening_rate')
                    %     % we want these measures
                    %     tempcol = tempTimeSeriesTable.getDependentColumn(java.lang.String(muscle)).getAsMat();
                    %     tempcolinterp = interp1(timespercent, tempcol, timespercent101);
                    %     if ~isfield(welkexostruct, muscle)
                    %         welkexostruct.(genvarname(muscle)) = [];
                    %     end
                    %     welkexostruct.(genvarname(muscle)) = [welkexostruct.(genvarname(muscle)), tempcolinterp];
                    % elseif contains(char(muscle), 'mechanical_work_rate')
                    %     % we want these measures
                    %     tempcol = tempTimeSeriesTable.getDependentColumn(java.lang.String(muscle)).getAsMat();
                    %     tempcolinterp = interp1(timespercent, tempcol, timespercent101);
                    %     if ~isfield(welkexostruct, muscle)
                    %         welkexostruct.(genvarname(muscle)) = [];
                    %     end
                    %     welkexostruct.(genvarname(muscle)) = [welkexostruct.(genvarname(muscle)), tempcolinterp];
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
                disp(trialdir)

                % now figure out how to get and plot the signal i want
                % have all the muscle analysis files already
                % do I want to do average or individual?
                tempfile = strcat(trialdir, '/analyzemuscles',tag,'_ForceReporter_forces', '.sto');
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
                    coord = char(labels.get(i));
                    
                    % need to screen only the things that we want
                    if contains(char(coord), 'GRF_')
                        % tempsplit = split(coord,'/');
                        % coordshort = string(tempsplit(4));
                        
                        % we also want the whole body measure
                        tempcol = tempTimeSeriesTable.getDependentColumn(java.lang.String(coord)).getAsMat();
                        tempcolinterp = interp1(timespercent, tempcol, timespercent101);
                        if ~isfield(welknaturalstruct, coord)
                            welknaturalstruct.(genvarname(coord)) = [];
                        end
                        welknaturalstruct.(genvarname(coord)) = [welknaturalstruct.(genvarname(coord)), tempcolinterp];
                    elseif contains(char(coord), 'HOBL')
                        tempcol = tempTimeSeriesTable.getDependentColumn(java.lang.String(coord)).getAsMat();
                        tempcolinterp = interp1(timespercent, tempcol, timespercent101);
                        if ~isfield(welknaturalstruct, coord)
                             welknaturalstruct.(genvarname(coord)) = [];
                        end
                        welknaturalstruct.(genvarname(coord)) = [welknaturalstruct.(genvarname(coord)), tempcolinterp];

                    % elseif contains(char(muscle), 'metabolics_combined')
                    %     % we want these measures
                    %     tempcol = tempTimeSeriesTable.getDependentColumn(java.lang.String(muscle)).getAsMat();
                    %     tempcolinterp = interp1(timespercent, tempcol, timespercent101);
                    %     if ~isfield(welknaturalstruct, muscle)
                    %         welknaturalstruct.(genvarname(muscle)) = [];
                    %     end
                    %     welknaturalstruct.(genvarname(muscle)) = [welknaturalstruct.(genvarname(muscle)), tempcolinterp];
                    % elseif contains(char(muscle), 'activation_maintenance_rate')
                    %     % we want these measures
                    %     tempcol = tempTimeSeriesTable.getDependentColumn(java.lang.String(muscle)).getAsMat();
                    %     tempcolinterp = interp1(timespercent, tempcol, timespercent101);
                    %     if ~isfield(welknaturalstruct, muscle)
                    %         welknaturalstruct.(genvarname(muscle)) = [];
                    %     end
                    %     welknaturalstruct.(genvarname(muscle)) = [welknaturalstruct.(genvarname(muscle)), tempcolinterp];
                    % elseif contains(char(muscle), 'shortening_rate')
                    %     % we want these measures
                    %     tempcol = tempTimeSeriesTable.getDependentColumn(java.lang.String(muscle)).getAsMat();
                    %     tempcolinterp = interp1(timespercent, tempcol, timespercent101);
                    %     if ~isfield(welknaturalstruct, muscle)
                    %         welknaturalstruct.(genvarname(muscle)) = [];
                    %     end
                    %     welknaturalstruct.(genvarname(muscle)) = [welknaturalstruct.(genvarname(muscle)), tempcolinterp];
                    % elseif contains(char(muscle), 'mechanical_work_rate')
                    %     % we want these measures
                    %     tempcol = tempTimeSeriesTable.getDependentColumn(java.lang.String(muscle)).getAsMat();
                    %     tempcolinterp = interp1(timespercent, tempcol, timespercent101);
                    %     if ~isfield(welknaturalstruct, muscle)
                    %         welknaturalstruct.(genvarname(muscle)) = [];
                    %     end
                        % welknaturalstruct.(genvarname(muscle)) = [welknaturalstruct.(genvarname(muscle)), tempcolinterp];
                    end
                end
            end
        end

        % this is only per subject one
        stdNat = struct();
        stdExo = struct();
        meanNat = struct();
        meanExo = struct();
    
        meanboth = struct();
        stdboth = struct();

        % okay now to plot etc. 
        % keyboard

        newlabels = fields(welkexostruct);
        % need to redo the labels
        tempfig = figure('Position',[1,1,1920,1080]);
        % do more stuff
        % averaging and whatnot
        for i=2:length(newlabels)
            subplot(4,5,i-1);
            templabel = newlabels(i);
            templabel = char(templabel);
                        
            % plot each of the gait cycles
            if contains(templabel, 'HOBL')
                % muscleplot_nat = welknaturalstruct.(genvarname(char(templabel)));
                muscleplot_exo = welkexostruct.(genvarname(char(templabel)));
                % plot(welknaturalstruct.time, muscleplot_nat, 'r:')
                plot(welkexostruct.time, muscleplot_exo, 'b:')
                % plot the subject average
                % plot(welknaturalstruct.time, mean(muscleplot_nat,2), 'r-', 'LineWidth', 1)
                hold on
                plot(welkexostruct.time, mean(muscleplot_exo,2), 'b-', 'LineWidth', 1)
            else
                muscleplot_nat = welknaturalstruct.(genvarname(char(templabel)));
                muscleplot_exo = welkexostruct.(genvarname(char(templabel)));
                plot(welknaturalstruct.time, muscleplot_nat, 'r:')
                hold on;
                plot(welkexostruct.time, muscleplot_exo, 'b:')
                % plot the subject average
                plot(welknaturalstruct.time, mean(muscleplot_nat,2), 'r-', 'LineWidth', 1)
                plot(welkexostruct.time, mean(muscleplot_exo,2), 'b-', 'LineWidth', 1)

                % add to the structs
                if ~isfield(stdNat, templabel)
                    stdNat.(genvarname(char(templabel))) = [];
                    stdExo.(genvarname(char(templabel))) = [];
                    meanNat.(genvarname(char(templabel))) = [];
                    meanExo.(genvarname(char(templabel))) = [];
                    meanboth.(genvarname(char(templabel))) = [];
                    stdboth.(genvarname(char(templabel))) = [];
                end
                stdNat.(genvarname(char(templabel))) = [stdNat.(genvarname(char(templabel))), 2.*std(muscleplot_nat,0,2)];
                stdExo.(genvarname(char(templabel))) = [stdExo.(genvarname(char(templabel))), 2.*std(muscleplot_exo,0,2)];
                meanNat.(genvarname(char(templabel))) = [meanNat.(genvarname(char(templabel))), mean(muscleplot_nat,2)];
                meanExo.(genvarname(char(templabel))) = [meanExo.(genvarname(char(templabel))), mean(muscleplot_exo,2)];
                meanboth.(genvarname(char(templabel))) = [meanboth.(genvarname(char(templabel))), mean([muscleplot_nat,muscleplot_exo],2)];
                stdboth.(genvarname(char(templabel))) = [stdboth.(genvarname(char(templabel))), 2.*std([muscleplot_nat,muscleplot_exo],0,2)];
                % now the full subjects as well
                if ~isfield(allsubj_nat, templabel)
                    allsubj_nat.(genvarname(char(templabel))) = [];
                    allsubj_exo.(genvarname(char(templabel))) = [];
                end
                allsubj_nat.(genvarname(char(templabel))) = [allsubj_nat.(genvarname(char(templabel))), muscleplot_nat];
                allsubj_exo.(genvarname(char(templabel))) = [allsubj_exo.(genvarname(char(templabel))), muscleplot_exo];

            end
            
            ylabel('Force [N]');
            % title(templabel)
            xlabel('% gait cycle')
            % select the name out 
            templabel2 = strrep(templabel,'_',' ');
            % for activation maintenance rate
            % templabel2 = templabel()
            title(templabel2);
            grid on;
        end
        
        % print(tempfig, ...
        %     strcat(strcat('G:\Shared drives\Exotendon\muscleModel\analysis\', strcat(subject,'\')), strcat(strcat(tempthing, '_acrossconditions'), '.png')),...
        %     '-dpng', '-r500')
        disp('print 1')
        
        % keyboard

        % % save the individual subject if you want to
        % meannattable = struct2table(meanNat); writetable(meannattable, 'mean_externalForces_nat.csv')
        % meanexotable = struct2table(meanExo); writetable(meanexotable, 'mean_externalForces_exo.csv')
        % stdexotable = struct2table(stdExo); writetable(stdexotable, 'std_externalForces_exo.csv')
        % stdnattable = struct2table(stdNat); writetable(stdnattable, 'std_externalForces_nat.csv')
        % meanbothtable = struct2table(meanboth); writetable(meanbothtable, 'mean_externalForces_both.csv')
        % stdbothtable = struct2table(stdboth); writetable(stdbothtable, 'std_externalForces_both.csv')

        % add the subject average to the combined struct?
        welknaturalstruct_combine.(genvarname(subject)) = welknaturalstruct;
        welkexostruct_combine.(genvarname(subject)) = welkexostruct;
        
        
        
    end

    % loop through the subjects again?
    markr = {'r:','r--'};
    markb = {'b:','b--'};
    

    % now plot across subjects
    tempfig2 = figure('Position',[1,1,1920,1080]);
        % then loop through the muscles inside each subject
    for i=2:length(newlabels)
        subplot(4,5,i-1);
        templabel = newlabels(i);
        templabel = char(templabel);
        % loop through the subjects
        for subj=1:length(welksubjects)
            subject = char(welksubjects(subj));
            indsubjectmass = subjectmass.(genvarname(subject))
            if contains(templabel,'HOBL')
                % muscleplot_nat = welknaturalstruct_combine.(genvarname(subject)).(genvarname(char(templabel)));
                muscleplot_exo = welkexostruct_combine.(genvarname(subject)).(genvarname(char(templabel)));
                % have all of them, want the average plotted for each subject
                % plot(welknaturalstruct.time, mean(muscleplot_nat,2), char(markr(subj)))
                hold on;
                plot(welkexostruct.time, mean(muscleplot_exo,2), 'b--'); %char(markb(subj)))
                hold on

                % add to the full all means struct
                if ~isfield(allsubjmeans_exo, templabel)
                    allsubjmeans_exo.(genvarname(templabel)) = [];
                end
                allsubjmeans_exo.(genvarname(templabel)) = [allsubjmeans_exo.(genvarname(templabel)), mean(muscleplot_exo,2)];
            
            elseif contains(templabel, 'GRF_F')
                muscleplot_nat = welknaturalstruct_combine.(genvarname(subject)).(genvarname(char(templabel)));
                muscleplot_exo = welkexostruct_combine.(genvarname(subject)).(genvarname(char(templabel)));
                % have all of them, want the average plotted for each subject   
                plot(welknaturalstruct.time, mean(muscleplot_nat,2)./(indsubjectmass*9.81), 'r-'); %char(markr(subj)))
                hold on;
                plot(welkexostruct.time, mean(muscleplot_exo,2)./(indsubjectmass*9.81), 'b-'); %char(markb(subj)))

                % add to the full means structs
                 if ~isfield(allsubjmeans_exo, templabel)
                    allsubjmeans_exo.(genvarname(templabel)) = [];
                    allsubjmeans_nat.(genvarname(templabel)) = [];
                end
                allsubjmeans_exo.(genvarname(templabel)) = [allsubjmeans_exo.(genvarname(templabel)), mean(muscleplot_exo,2)./(indsubjectmass*9.81)];
                allsubjmeans_nat.(genvarname(templabel)) = [allsubjmeans_exo.(genvarname(templabel)), mean(muscleplot_nat,2)./(indsubjectmass*9.81)];
            else
                muscleplot_nat = welknaturalstruct_combine.(genvarname(subject)).(genvarname(char(templabel)));
                muscleplot_exo = welkexostruct_combine.(genvarname(subject)).(genvarname(char(templabel)));
                % have all of them, want the average plotted for each subject   
                plot(welknaturalstruct.time, mean(muscleplot_nat,2), 'r-'); %char(markr(subj)))
                hold on;
                plot(welkexostruct.time, mean(muscleplot_exo,2), 'b-'); %char(markb(subj)))

                % add to the full means structs
                 if ~isfield(allsubjmeans_exo, templabel)
                    allsubjmeans_exo.(genvarname(templabel)) = [];
                    allsubjmeans_nat.(genvarname(templabel)) = [];
                end
                allsubjmeans_exo.(genvarname(templabel)) = [allsubjmeans_exo.(genvarname(templabel)), mean(muscleplot_exo,2)];
                allsubjmeans_nat.(genvarname(templabel)) = [allsubjmeans_exo.(genvarname(templabel)), mean(muscleplot_nat,2)];

            end
        end
        
        templabel2 = strrep(templabel,'_',' ');
        title(templabel2)
        
        xlabel('% gait cycle')
        ylabel('Force [N]')
        grid on;
    end
    subplot(4,5,i);
    plot(welknaturalstruct.time, mean(muscleplot_nat,2), 'r-'); %char(markr(subj-1)))
    hold on;
    plot(welkexostruct.time, mean(muscleplot_exo,2), 'b-'); %char(markb(subj-1)))
    plot(welknaturalstruct.time, mean(muscleplot_nat,2), 'r-'); %char(markr(subj)))
    hold on;
    plot(welkexostruct.time, mean(muscleplot_exo,2), 'b-'); % char(markb(subj)))
    % legend('Subj 1 natural','Subj 1 exo','Subj 2 natural','Subj 2 exo')
    % 
    % print(tempfig2, ...
    %     strcat('G:\Shared drives\Exotendon\muscleModel\analysis\', tempthing, '_combined_newsubjs', '.png'),...
    %     '-dpng', '-r500')
    disp('print 2')

        

    % make the fig for a single line for each nat and exo for all subjects
    % and gait cycles
    % now plot across subjects
    tempfig3 = figure('Position',[1,1,1920,1080]);
        % then loop through the muscles inside each subject
    for i=2:length(newlabels)
        subplot(4,5,i-1);
        templabel = newlabels(i);
        templabel = char(templabel);
                    
        if contains(templabel,'HOBL')
                % muscleplot_nat = welknaturalstruct_combine.(genvarname(subject)).(genvarname(char(templabel)));
            muscleplot_exo = allsubjmeans_exo.(genvarname(templabel));
            % have all of them, want the average plotted for each subject
                % plot(welknaturalstruct.time, mean(muscleplot_nat,2), char(markr(subj)))
            hold on;
            plot(welkexostruct.time, mean(muscleplot_exo,2), 'b'); %char(markb(subj)))
            hold on
            % take an average and std across all subjects
            if ~isfield(meansall_exo, templabel)
                meansall_exo.(genvarname(templabel)) = [];
                stdsall_exo.(genvarname(templabel)) = [];
            end
            meansall_exo.(genvarname(templabel)) = [meansall_exo.(genvarname(templabel)), mean(muscleplot_exo, 2)];
            stdsall_exo.(genvarname(templabel)) = [stdsall_exo.(genvarname(templabel)), std(muscleplot_exo, 0, 2)];
        else
            muscleplot_nat = allsubjmeans_nat.(genvarname(templabel));
            muscleplot_exo = allsubjmeans_exo.(genvarname(templabel));
            % have all of them, want the average plotted for each subject   
            plot(welknaturalstruct.time, mean(muscleplot_nat,2), 'r'); %char(markr(subj)))
            hold on;
            plot(welkexostruct.time, mean(muscleplot_exo,2), 'b'); %char(markb(subj)))
            % take an average and std across all subjects
            if ~isfield(meansall_exo, templabel)
                meansall_exo.(genvarname(templabel)) = [];
                stdsall_exo.(genvarname(templabel)) = [];
                meansall_nat.(genvarname(templabel)) = [];
                stdsall_nat.(genvarname(templabel)) = [];
            end
            meansall_exo.(genvarname(templabel)) = [meansall_exo.(genvarname(templabel)), mean(muscleplot_exo, 2)];
            stdsall_exo.(genvarname(templabel)) = [stdsall_exo.(genvarname(templabel)), std(muscleplot_exo, 0, 2)];
            meansall_nat.(genvarname(templabel)) = [meansall_nat.(genvarname(templabel)), mean(muscleplot_nat, 2)];
            stdsall_nat.(genvarname(templabel)) = [stdsall_nat.(genvarname(templabel)), std(muscleplot_nat, 0, 2)];
        end
    

    templabel2 = strrep(templabel,'_',' ');
    title(templabel2)
    
    xlabel('% gait cycle')
    ylabel('Force [N]')
    grid on;
    end
    subplot(4,5,i);
    plot(welknaturalstruct.time, mean(muscleplot_nat,2), 'r-'); %char(markr(subj-1)))
    hold on;
    plot(welkexostruct.time, mean(muscleplot_exo,2), 'b-'); %char(markb(subj-1)))
    plot(welknaturalstruct.time, mean(muscleplot_nat,2), 'r-'); %char(markr(subj)))
    hold on;
    plot(welkexostruct.time, mean(muscleplot_exo,2), 'b-'); % char(markb(subj)))
    % legend('Subj 1 natural','Subj 1 exo','Subj 2 natural','Subj 2 exo')
    % 
    % print(tempfig3, ...
    %     strcat('G:\Shared drives\Exotendon\muscleModel\analysis\', tempthing, '_combined_newsubjs_singleavg', '.png'),...
    %     '-dpng', '-r500')
    disp('print 3')

    
    cd("G:\Shared drives\Exotendon\muscleModel\muscleEnergyModel\")

    % save the data from the means and stds
    meansall_nat_table = struct2table(meansall_nat);
    meansall_exo_table = struct2table(meansall_exo);
    stdsall_nat_table = struct2table(stdsall_nat);
    stdsall_exo_table = struct2table(stdsall_exo);

    writetable(meansall_nat_table, 'meansall_GRF_nat.csv');
    writetable(meansall_exo_table, 'meansall_GRF_exo.csv');
    writetable(stdsall_nat_table, 'stdsall_GRF_nat.csv');
    writetable(stdsall_exo_table, 'stdsall_GRF_exo.csv');

    % now plot across subjects
    tempfig4 = figure('Position',[1,1,1920,1080]);
        % then loop through the muscles inside each subject
    for i=2:3
        subplot(3,1,i-1);
        templabel = newlabels(i);
        templabel = char(templabel);
                    
        if contains(templabel,'HOBL')
                % muscleplot_nat = welknaturalstruct_combine.(genvarname(subject)).(genvarname(char(templabel)));
            muscleplot_exo = allsubjmeans_exo.(genvarname(templabel));
            % have all of them, want the average plotted for each subject
                % plot(welknaturalstruct.time, mean(muscleplot_nat,2), char(markr(subj)))
            hold on;
            plot(welkexostruct.time, mean(muscleplot_exo,2), 'c' ,'LineWidth',4); %char(markb(subj)))
            hold on

        else
            muscleplot_nat = allsubjmeans_nat.(genvarname(templabel));
            muscleplot_exo = allsubjmeans_exo.(genvarname(templabel));
            % have all of them, want the average plotted for each subject   
            plot(welknaturalstruct.time, mean(muscleplot_nat,2), 'r'); %char(markr(subj)))
            hold on;
            plot(welkexostruct.time, mean(muscleplot_exo,2), 'b'); %char(markb(subj)))
        end
    

        templabel2 = strrep(templabel,'_',' ');
        % title(templabel2)
        % 
        % xlabel('% gait cycle')
        ylabel('Force [N]')
        fontsize(26)
        grid off;
    end
    subplot(3,1,i);
    plot(welknaturalstruct.time, mean(muscleplot_nat,2), 'r-'); %char(markr(subj-1)))
    hold on;
    plot(welkexostruct.time, mean(muscleplot_exo,2), 'b-'); %char(markb(subj-1)))
    plot(welknaturalstruct.time, mean(muscleplot_nat,2), 'r-'); %char(markr(subj)))
    hold on;
    plot(welkexostruct.time, mean(muscleplot_exo,2), 'b-'); % char(markb(subj)))
    % legend('Subj 1 natural','Subj 1 exo','Subj 2 natural','Subj 2 exo')
    
    print(tempfig4, ...
        strcat('G:\Shared drives\Exotendon\muscleModel\analysis\', tempthing, 'only_exoforce', '.png'),...
        '-dpng', '-r500')
    disp('print 4')


end
