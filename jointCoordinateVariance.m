% written by Jon Stingel
% 20220901
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
thingstoplot = {'coordinates'}; % 'probes', 'shortening', 'mechanical', 'activation'
load 'G:\Shared drives\Exotendon\muscleModel\muscleEnergyModel\subjectgaitcycles.mat';

% loop through each of the things we want to plot
for thing=1:length(thingstoplot)
    tempthing = char(thingstoplot(thing))

    % create stucture for combined subject figures
    welknaturalstruct_combine = struct();
    welkexostruct_combine = struct();

    subjectexovar_combine = {};
    subjectnatvar_combine = {};

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
                % tempfile = strcat(trialdir, '/muscle_statetrack_grfprescribe_solution', '.sto');
                % tempfile = strcat(trialdir, '/torque_statetrack_grfprescribe_tracked_states','.sto');
                tempfile = strcat(trialdir, '/coordinates_updated.mot');
                
                tempTimeSeriesTable = TimeSeriesTable(tempfile);
                temptime = tempTimeSeriesTable.getIndependentColumn();

                gait_start = subjectgaitcycles.(genvarname(subject)).(genvarname(condition)).(genvarname(test)).initial;
                gait_end = subjectgaitcycles.(genvarname(subject)).(genvarname(condition)).(genvarname(test)).final;
                
                time = [];
                for t=0:temptime.size()-1
                    try
                        tempix = temptime.get(t).doubleValue();
                    catch
                        tempix = temptime.get(t);
                    end
                    if tempix >= gait_start && tempix <= gait_end
                        time = [time, tempix];
                    end
                end
                        
%                 times = zeros(temptime.size(),1);
%                 for i=0:temptime.size()-1
%                     times(i+1) = temptime.get(i);
%                 end
%                 timespercent = (times - times(1)) ./ (times(end) - times(1)) .*100;
                
                timespercent = (time - time(1))./(time(end)-time(1)).*100;
                timespercent101 = [0:1:100]';
                welkexostruct.time = timespercent101;

                % grab the start and finish index for the rows in the
                % table
                if temptime.get(tempTimeSeriesTable.getRowIndexAfterTime(gait_start)) == time(1)
                    row_idx_start = tempTimeSeriesTable.getRowIndexAfterTime(gait_start);
                else
                    row_idx_start = tempTimeSeriesTable.getRowIndexBeforeTime(gait_start);
                end
                if temptime.get(tempTimeSeriesTable.getRowIndexAfterTime(gait_end)) == time(end)
                    row_idx_end = tempTimeSeriesTable.getRowIndexAfterTime(gait_end);
                else
                    row_idx_end = tempTimeSeriesTable.getRowIndexBeforeTime(gait_end);
                end

                % now for each of the things
                numCols = tempTimeSeriesTable.getNumColumns(); % including time
                labels = tempTimeSeriesTable.getColumnLabels();            

                for i=0:labels.size()-1
                    coord = char(labels.get(i));
                    
                   
                    % we want it all
                    % tempsplit = split(coord,'/');
                    % coordshort = string(tempsplit(4));
                    tempcol = tempTimeSeriesTable.getDependentColumn(java.lang.String(coord)).getAsMat();
                    % now take the rows that we want
                    col = tempcol(row_idx_start:row_idx_end);

                    tempcolinterp = interp1(timespercent, col, timespercent101);
                    if ~isfield(welkexostruct, coord)
                        welkexostruct.(genvarname(coord)) = [];
                    end
                    welkexostruct.(genvarname(coord)) = [welkexostruct.(genvarname(coord)), tempcolinterp];

                    
                    end
                end
            end
        % end
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
%                 tempfile = strcat(trialdir, '/muscle_statetrack_grfprescribe_solution', '.sto');
                tempfile = strcat(trialdir, '/coordinates_updated.mot');
                tempTimeSeriesTable = TimeSeriesTable(tempfile);
                temptime = tempTimeSeriesTable.getIndependentColumn();
                
                gait_start = subjectgaitcycles.(genvarname(subject)).(genvarname(condition)).(genvarname(test)).initial;
                gait_end = subjectgaitcycles.(genvarname(subject)).(genvarname(condition)).(genvarname(test)).final;
                
                time = [];
                for t=0:temptime.size()-1
                   try
                        tempix = temptime.get(t).doubleValue();
                    catch
                        tempix = temptime.get(t);
                    end
                    if tempix >= gait_start && tempix <= gait_end
                        time = [time, tempix];
                    end
                end
                        
%                 times = zeros(temptime.size(),1);
%                 for i=0:temptime.size()-1
%                     times(i+1) = temptime.get(i);
%                 end
%                 timespercent = (times - times(1)) ./ (times(end) - times(1)) .*100;
                
                timespercent = (time - time(1))./(time(end)-time(1)).*100;
                timespercent101 = [0:1:100]';
                welknaturalstruct.time = timespercent101;

                % grab the start and finish index for the rows in the
                % table
                if temptime.get(tempTimeSeriesTable.getRowIndexAfterTime(gait_start)) == time(1)
                    row_idx_start = tempTimeSeriesTable.getRowIndexAfterTime(gait_start);
                else
                    row_idx_start = tempTimeSeriesTable.getRowIndexBeforeTime(gait_start);
                end
                if temptime.get(tempTimeSeriesTable.getRowIndexAfterTime(gait_end)) == time(end)
                    row_idx_end = tempTimeSeriesTable.getRowIndexAfterTime(gait_end);
                else
                    row_idx_end = tempTimeSeriesTable.getRowIndexBeforeTime(gait_end);
                end
                %{
                times = zeros(temptime.size(),1);
                for i=0:temptime.size()-1
                    times(i+1) = temptime.get(i);
                end
                timespercent = (times - times(1)) / (times(end) - times(1)) *100;
                timespercent101 = [0:1:100]';
                welknaturalstruct.time = timespercent101;
                %}
                % now for each of the things
                numCols = tempTimeSeriesTable.getNumColumns(); % including time
                labels = tempTimeSeriesTable.getColumnLabels();            

                
                for i=0:labels.size()-1
                    coord = char(labels.get(i));
                        
                        % we also want the whole body measure
                    tempcol = tempTimeSeriesTable.getDependentColumn(java.lang.String(coord)).getAsMat();
                    col = tempcol(row_idx_start:row_idx_end);
                    tempcolinterp = interp1(timespercent, col, timespercent101);
                    if ~isfield(welknaturalstruct, coord)                         
                        welknaturalstruct.(genvarname(coord)) = [];
                    end
                    welknaturalstruct.(genvarname(coord)) = [welknaturalstruct.(genvarname(coord)), tempcolinterp];
                   
                end
            end
        end


        % okay now to compute the variance metrics that I want for each subject and each condition
        
        % metrics: range, mean, std dev. 

        subjectexovar = {};
        subjectnatvar = {};

        % loop through all the coords
        newlabels = fields(welkexostruct);
        for idk=2:length(newlabels)
            coord = newlabels(idk);
            coord = char(coord);
            coorddata = welkexostruct.(genvarname(char(coord)));
            coorddata2 = welknaturalstruct.(genvarname(char(coord)));

            coordtimevar = [];
            coordtimevar2 = [];
            % now to loop through time
            for t=1:length(coorddata)
                % do the calcs and add them to the vector we want
                temprange = range(coorddata(t,:));
                tempavg = mean(coorddata(t,:));
                tempstd = std(coorddata(t,:));
                coordtimevar = [coordtimevar; [temprange,tempavg,tempstd]];

                % now natural 
                temprange2 = range(coorddata2(t,:));
                tempavg2 = mean(coorddata2(t,:));
                tempstd2 = std(coorddata2(t,:));
                coordtimevar2 = [coordtimevar2; [temprange2, tempavg2,tempstd2]];

            end

            



            % okay have one coord through time
            subjectexovar.(genvarname(char(coord))) = mean(coordtimevar(:,1));
            subjectnatvar.(genvarname(char(coord))) = mean(coordtimevar2(:,1));
        
            
        end

        % add to the combined
        subjectexovar_combine.(genvarname(subject)) = subjectexovar;
        subjectnatvar_combine.(genvarname(subject)) = subjectnatvar;
        

        
        % okay now to plot etc. 
        newlabels = fields(welkexostruct);
        % need to redo the labels
        tempfig = figure('Position',[1,1,1920,1080]);
        % do more stuff
        % averaging and whatnot
        for i=2:length(newlabels)
            subplot(6,7,i-1);
            templabel = newlabels(i);
            templabel = char(templabel);
            
            % plot each of the gait cycles
            subjectexovar_plot = subjectexovar.(genvarname(char(templabel)));
            subjectnatvar_plot = subjectnatvar.(genvarname(char(templabel)));

            bar([subjectexovar_plot,subjectnatvar_plot]);
            xticklabels({strcat('exo:',' ',string(subjectexovar_plot)),strcat('natural:',' ', string(subjectnatvar_plot))});
            xtickangle(30);
            ylabel('Range');
            % select the name out 
            templabel2 = strrep(templabel,'_',' ');
            % for activation maintenance rate
            % templabel2 = templabel()
            title(templabel2);

            
        end
        
        print(tempfig, ...
            strcat(strcat('G:\Shared drives\Exotendon\muscleModel\analysis\', strcat(subject,'\')), strcat(strcat('coordvariance', '_acrossconditions'), '.png')),...
            '-dpng', '-r500')
        disp('print 1')
        


        % add the subject average to the combined struct?
        welknaturalstruct_combine.(genvarname(subject)) = welknaturalstruct;
        welkexostruct_combine.(genvarname(subject)) = welkexostruct;
        
    end

    
%     % loop through the subjects again?
%     markr = {'r:','r--'};
%     markb = {'b:','b--'};
%     
%     % now plot across subjects
%     tempfig2 = figure('Position',[1,1,1920,1080]);
%         % then loop through the muscles inside each subject
%     for i=2:length(newlabels)
%         subplot(6,7,i-1);
%         templabel = newlabels(i);
%         templabel = char(templabel);
%         % loop through the subjects
%         for subj=1:length(welksubjects)
%             subject = char(welksubjects(subj));
%             
%             subjectexovar_plot = subjectexovar_combine.(genvarname(subject)).(genvarname(char(templabel)));
%             subjectnatvar_plot = subjectnatvar_combine.(genvarname(subject)).(genvarname(char(templabel)));
% 
% 
%             % have all of them, want the average plotted for each subject
%             bar(subjectnatvar_plot, 'r-');
%             hold on;
%             bar(subjectexovar_plot, 'b-');%char(markb(subj)))
%         end
%         
%         templabel2 = strrep(templabel,'_',' ');
%         title(templabel2)
%         
%         xlabel('% gait cycle')
%         ylabel('Coordinate Value [rad/m]')
%         grid on;
%     end
% %     subplot(5,7,i);
% %     plot(welknaturalstruct.time, mean(muscleplot_nat,2), char(markr(subj-1)))
% %     hold on;
% %     plot(welkexostruct.time, mean(muscleplot_exo,2), char(markb(subj-1)))
% %     plot(welknaturalstruct.time, mean(muscleplot_nat,2), char(markr(subj)))
% %     hold on;
% %     plot(welkexostruct.time, mean(muscleplot_exo,2), char(markb(subj)))
% %     legend('Subj 1 natural','Subj 1 exo','Subj 2 natural','Subj 2 exo')
% 
%     print(tempfig2, ...
%         strcat('G:\Shared drives\Exotendon\muscleModel\analysis\', tempthing, '_combined', '.png'),...
%         '-dpng', '-r500')
%     disp('print 2')
end
