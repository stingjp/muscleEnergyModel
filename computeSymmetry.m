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

    naturalstruct_combine = struct();
    exostruct_combine = struct();


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
                keyboard
                tempfile = strcat(trialdir, '/coordinates_updated.mot');
                
                tempTimeSeriesTable = TimeSeriesTable(tempfile);
                temptime = tempTimeSeriesTable.getIndependentColumn();

                gait_start = subjectgaitcycles.(genvarname(subject)).(genvarname(condition)).(genvarname(test)).initial;
                gait_end = subjectgaitcycles.(genvarname(subject)).(genvarname(condition)).(genvarname(test)).final;
                
                time = [];
                for t=0:temptime.size()-1
                    tempix = temptime.get(t);
                    try
                        if tempix >= gait_start && tempix <= gait_end
                            time = [time, tempix];
                        end
                    catch
                        if tempix.doubleValue() >= gait_start && tempix.doubleValue() <= gait_end
                             time = [time, tempix.doubleValue()];
                        end
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
                    
                    % get the range, as well as average value
                    temprange = range(tempcolinterp);
                    tempavg = mean(tempcolinterp);

                    if ~isfield(welkexostruct, coord)
                        welkexostruct.(genvarname(coord)) = [];
                    end
                    if ~isfield(exostruct_combine, coord)
                        exostruct_combine.(genvarname(coord)) = [];
                    end
                    welkexostruct.(genvarname(coord)) = [welkexostruct.(genvarname(coord)), [temprange, tempavg]];
                    exostruct_combine.(genvarname(coord)) = [exostruct_combine.(genvarname(coord)); [temprange, tempavg]];

                    
                    
                    
                    %{
                    % need to screen only the things that we want
                    if contains(char(coord), 'value') % 'value' denotes a coordinate value
                        tempsplit = split(coord,'/');
                        coordshort = string(tempsplit(4));
                        
                        tempcol = tempTimeSeriesTable.getDependentColumn(java.lang.String(coord)).getAsMat();
                        tempcolinterp = interp1(timespercent, tempcol, timespercent101);
                        if ~isfield(welkexostruct, coordshort)
                            welkexostruct.(genvarname(coordshort)) = [];
                        end
                        welkexostruct.(genvarname(coordshort)) = [welkexostruct.(genvarname(coordshort)), tempcolinterp];
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
                    %}
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
                    tempix = temptime.get(t);
                    try
                        if tempix >= gait_start && tempix <= gait_end
                            time = [time, tempix];
                        end
                    catch
                        if tempix.doubleValue() >= gait_start && tempix.doubleValue() <= gait_end
                            time = [time, tempix];
                        end
                    end
                end
                        
%                 times = zeros(temptime.size(),1);
%                 for i=0:temptime.size()-1
%                     times(i+1) = temptime.get(i);
%                 end
%                 timespercent = (times - times(1)) ./ (times(end) - times(1)) .*100;
                testtime = double(time);
                timespercent = (testtime - testtime(1))./(testtime(end)-testtime(1)).*100;
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

                    % get range and avg
                    temprange = range(tempcolinterp);
                    tempavg = mean(tempcolinterp);


                    if ~isfield(welknaturalstruct, coord)                         
                        welknaturalstruct.(genvarname(coord)) = [];
                    end
                    if ~isfield(naturalstruct_combine, coord)
                        naturalstruct_combine.(genvarname(coord)) = [];
                    end
                    welknaturalstruct.(genvarname(coord)) = [welknaturalstruct.(genvarname(coord)), [temprange, tempavg]];
                    naturalstruct_combine.(genvarname(coord)) = [naturalstruct_combine.(genvarname(coord)); [temprange, tempavg]];
                   
                    %{
                    % need to screen only the things that we want
                    if contains(char(coord), 'value')
                        tempsplit = split(coord,'/');
                        coordshort = string(tempsplit(4));
                        
                        % we also want the whole body measure
                        tempcol = tempTimeSeriesTable.getDependentColumn(java.lang.String(coord)).getAsMat();
                        tempcolinterp = interp1(timespercent, tempcol, timespercent101);
                        if ~isfield(welknaturalstruct, coordshort)
                            welknaturalstruct.(genvarname(coordshort)) = [];
                        end
                        
                        welknaturalstruct.(genvarname(coordshort)) = [welknaturalstruct.(genvarname(coordshort)), tempcolinterp];
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
                    %}
                end
            end
        end


        
        % okay now to plot etc. 
        
        % newlabels = fields(welkexostruct);
        % % need to redo the labels
        % tempfig = figure('Position',[1,1,1920,1080]);
        % % do more stuff
        % % averaging and whatnot
        % for i=2:length(newlabels)
        %     subplot(6,7,i-1);
        %     templabel = newlabels(i);
        %     templabel = char(templabel);
        %     % plot each of the gait cycles
        %     muscleplot_nat = welknaturalstruct.(genvarname(char(templabel)));
        %     muscleplot_exo = welkexostruct.(genvarname(char(templabel)));
        %     plot(welknaturalstruct.time, muscleplot_nat, 'r:')
        %     hold on;
        %     plot(welkexostruct.time, muscleplot_exo, 'b:')
        %     % plot the subject average
        %     plot(welknaturalstruct.time, mean(muscleplot_nat,2), 'r-', 'LineWidth', 1)
        %     plot(welkexostruct.time, mean(muscleplot_exo,2), 'b-', 'LineWidth', 1)
        %     ylabel('Coordinate Value [rad/m]');
        %     % title(templabel)
        %     xlabel('% gait cycle')
        %     % select the name out 
        %     templabel2 = strrep(templabel,'_',' ');
        %     % for activation maintenance rate
        %     % templabel2 = templabel()
        %     title(templabel2);
        %     grid on;
        % end
        
        % print(tempfig, ...
        %     strcat(strcat('G:\Shared drives\Exotendon\muscleModel\analysis\', strcat(subject,'\')), strcat(strcat(tempthing, '_acrossconditions'), '.png')),...
        %     '-dpng', '-r500')
        % disp('print 1')
        


        % add the subject average to the combined struct?
        welknaturalstruct_combine.(genvarname(subject)) = welknaturalstruct;
        welkexostruct_combine.(genvarname(subject)) = welkexostruct;
        
    end

    keyboard
    newlabels = fields(exostruct_combine);
    
    exo_stds = exostruct_combine;
    nat_stds = naturalstruct_combine;
    
    exo_means = exostruct_combine;
    nat_means = naturalstruct_combine;

    % loop through the subjects again?
%     markr = {'r:','r--'};
%     markb = {'b:','b--'};
    
    % now plot across subjects
%     tempfig2 = figure('Position',[1,1,1280,1920]);
        % then loop through the muscles inside each subject
    for i=1:length(newlabels)
%         subplot(9,3,i-1);
        templabel = newlabels(i);
        templabel = char(templabel)
%         temp1 = [];
%         temp2 = [];

        % output the mean range
        exo_means.(genvarname(templabel)) = mean(exostruct_combine.(genvarname(templabel)),1);
        nat_means.(genvarname(templabel)) = mean(naturalstruct_combine.(genvarname(templabel)),1);
        exo_stds.(genvarname(templabel)) = std(exostruct_combine.(genvarname(templabel)),1);
        nat_stds.(genvarname(templabel)) = std(naturalstruct_combine.(genvarname(templabel)),1);

        % output the mean mean
    end
    % first is range, second is avg

    



% 
% %         % loop through the subjects
% %         for subj=1:length(welksubjects)
% %             subject = char(welksubjects(subj));
% %             muscleplot_nat = welknaturalstruct_combine.(genvarname(subject)).(genvarname(char(templabel)));
% %             muscleplot_exo = welkexostruct_combine.(genvarname(subject)).(genvarname(char(templabel)));
% %             
% %             temp1 = [temp1, mean(muscleplot_nat,2)];
% %             temp2 = [temp2, mean(muscleplot_exo,2)];
% 
%             
%             % have all of them, want the average plotted for each subject
% %             plot(welknaturalstruct.time, mean(muscleplot_nat,2), 'r:','LineWidth',0.4);%char(markr(subj)))
% %             hold on;
% %             plot(welkexostruct.time, mean(muscleplot_exo,2), 'b:','LineWidth',0.4);%char(markb(subj)))
% %         end
%         
% %         plot(mean(temp1,2), 'r', 'LineWidth', 2)
% %         plot(mean(temp2,2), 'b', 'LineWidth', 2)
% %         legend(num2str(min(mean(temp1,2))), num2str(min(mean(temp2,2))));
%         
%         
%         % note that we can do max or min for flexions/extensions
%         % doing ROM now too
%         disp(strcat('nat: ', templabel))
%         maxnat = max(temp1) - min(temp1)
%         natmaxavg = mean(maxnat)
%         natsd = std(maxnat)
%         natse = natsd/sqrt(length(maxnat))
% 
%         disp(strcat('exo: ',templabel))
%         maxexo = max(temp2) - min(temp2)
%         exomaxavg = mean(maxexo)
%         exosd = std(maxexo)
%         exose = exosd/sqrt(length(maxexo))
% 
% 
% 
% 
% 
% 
% 
%         templabel2 = strrep(templabel,'_',' ');
%         title(templabel2)
%         
%         xlabel('% gait cycle')
%         ylabel('Coordinate Value')
%         grid on;
%     end
%     subplot(5,7,i);
%     plot(welknaturalstruct.time, mean(muscleplot_nat,2), char(markr(subj-1)))
%     hold on;
%     plot(welkexostruct.time, mean(muscleplot_exo,2), char(markb(subj-1)))
%     plot(welknaturalstruct.time, mean(muscleplot_nat,2), char(markr(subj)))
%     hold on;
%     plot(welkexostruct.time, mean(muscleplot_exo,2), char(markb(subj)))
%     legend('Subj 1 natural','Subj 1 exo','Subj 2 natural','Subj 2 exo')

    print(tempfig2, ...
        strcat('G:\Shared drives\Exotendon\muscleModel\analysis\', tempthing, '_combined', '.png'),...
        '-dpng', '-r500')
    disp('print 2')
end
