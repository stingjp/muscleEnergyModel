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
    totalsubjsstruct_combine = struct();

    naturalstruct_combine = struct();
    exostruct_combine = struct();
    totalstruct_combine = struct();


    % loop through the subjects
    for subj=1:length(welksubjects)
        subject = char(welksubjects(subj));
        subjdir = strcat(resultsdir, strcat('/',subject));
        
        % create the struct for individual figures
        welknaturalstruct = struct();
        welkexostruct = struct();
    

        % add the subject to the total struct
        totalsubjsstruct_combine.(genvarname(subject)) = [];


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
                
%                 tempfile = strcat(trialdir, '/coordinates_updated.mot');
                tempfile = strcat(trialdir, '/muscle_coordinates_short.sto');
                
                tempTimeSeriesTable = TimeSeriesTable(tempfile);
                temptime = tempTimeSeriesTable.getIndependentColumn();

                gait_start = subjectgaitcycles.(genvarname(subject)).(genvarname(condition)).(genvarname(test)).initial;
                gait_end = subjectgaitcycles.(genvarname(subject)).(genvarname(condition)).(genvarname(test)).final;
                
                
                % % need a way to match right and left curves. 
                % kneer = tempTimeSeriesTable.getDependentColumn('knee_angle_r').getAsMat();
                % kneel = tempTimeSeriesTable.getDependentColumn('knee_angle_l').getAsMat();
                % [rm, ixr] = max(kneer);
                % [lm, ixl] = max(kneel);
                % newkneer = [kneer(ixr:end); kneer(1:ixr-1)];
                % newkneel = [kneel(ixl:end); kneel(1:ixl-1)];


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

                
                % need a way to match right and left curves. 
                kneer = tempTimeSeriesTable.getDependentColumn('knee_angle_r').getAsMat();
                kneel = tempTimeSeriesTable.getDependentColumn('knee_angle_l').getAsMat();
                
                % now take the rows that we want
                % col = tempcol(row_idx_start:row_idx_end);
                tcol_r = kneer;
                ttempcolinterp_r = interp1(timespercent, tcol_r, timespercent101);
                tcol_l = kneel;
                ttempcolinterp_l = interp1(timespercent, tcol_l, timespercent101);

                
                [rm, ixr] = max(ttempcolinterp_r);
                [lm, ixl] = max(ttempcolinterp_l);
                newkneer = [ttempcolinterp_r(ixr:end); ttempcolinterp_r(1:ixr-1)];
                newkneel = [ttempcolinterp_l(ixl:end); ttempcolinterp_l(1:ixl-1)];


                % now for each of the things
                numCols = tempTimeSeriesTable.getNumColumns(); % including time
                labels = tempTimeSeriesTable.getColumnLabels();            

                for i=0:labels.size()-1
                    coord = char(labels.get(i));
                    
                    if endsWith(coord, '_r')
                        coordshort = coord(1:end-2);
                        coordleft = strcat(coordshort,'_l');

                        
                       
                        tempcol_r = tempTimeSeriesTable.getDependentColumn(java.lang.String(coord)).getAsMat();
                        tempcol_l = tempTimeSeriesTable.getDependentColumn(java.lang.String(coordleft)).getAsMat();
                        
                        % now take the rows that we want
    %                     col = tempcol(row_idx_start:row_idx_end);
                        col_r = tempcol_r;
                        tempcolinterp_r = interp1(timespercent, col_r, timespercent101);
                        col_l = tempcol_l;
                        tempcolinterp_l = interp1(timespercent, col_l, timespercent101);
                        

                        % match the right and lefts
                        % [rm, ixr] = max(tempcolinterp_r);
                        % [lm, ixl] = max(tempcolinterp_l);
                        newtempcolinterp_r = [tempcolinterp_r(ixr:end); tempcolinterp_r(1:ixr-1)];
                        newtempcolinterp_l = [tempcolinterp_l(ixl:end); tempcolinterp_l(1:ixl-1)];
    
                        % figure();
                        % hold on;
                        % plot(newtempcolinterp_r); plot(newtempcolinterp_l);

                        maerror = zeros(length(newtempcolinterp_l), 1);
                        % get the MAE
                        for h=1:length(newtempcolinterp_l)
                            % get the mae at each time point
%                             maerror(h) = abs(newtempcolinterp_l(h) - newtempcolinterp_r(h));
                            maerror(h) = (newtempcolinterp_l(h) - newtempcolinterp_r(h))^2;
                        end
                        
                        


                        % add to the struct
                        if ~isfield(welkexostruct, coordshort)
                            welkexostruct.(genvarname(coordshort)) = [];
                        end
                        if ~isfield(exostruct_combine, coordshort)
                            exostruct_combine.(genvarname(coordshort)) = [];
                        end
                        if ~isfield(totalstruct_combine, coordshort)
                            totalstruct_combine.(genvarname(coordshort)) = [];
                        end
                        welkexostruct.(genvarname(coordshort)) = [welkexostruct.(genvarname(coordshort)), maerror];
                        exostruct_combine.(genvarname(coordshort)) = [exostruct_combine.(genvarname(coordshort)), maerror];
                        totalstruct_combine.(genvarname(coordshort)) = [totalstruct_combine.(genvarname(coordshort)), maerror];

                        totalsubjsstruct_combine.(genvarname(subject)) = [totalsubjsstruct_combine.(genvarname(subject)), maerror];
                    end

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
%                 tempfile = strcat(trialdir, '/coordinates_updated.mot');
                tempfile = strcat(trialdir, '/muscle_coordinates_short.sto');
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

                
                % need a way to match right and left curves. 
                kneer = tempTimeSeriesTable.getDependentColumn('knee_angle_r').getAsMat();
                kneel = tempTimeSeriesTable.getDependentColumn('knee_angle_l').getAsMat();
                

                % now take the rows that we want
                % col = tempcol(row_idx_start:row_idx_end);
                tcol_r = kneer;
                ttempcolinterp_r = interp1(timespercent, tcol_r, timespercent101);
                tcol_l = kneel;
                ttempcolinterp_l = interp1(timespercent, tcol_l, timespercent101);

                [rm, ixr] = max(ttempcolinterp_r);
                [lm, ixl] = max(ttempcolinterp_l);
                newkneer = [ttempcolinterp_r(ixr:end); ttempcolinterp_r(1:ixr-1)];
                newkneel = [ttempcolinterp_l(ixl:end); ttempcolinterp_l(1:ixl-1)];


                % now for each of the things
                numCols = tempTimeSeriesTable.getNumColumns(); % including time
                labels = tempTimeSeriesTable.getColumnLabels();            

                
                for i=0:labels.size()-1
                    coord = char(labels.get(i));
                    if endsWith(coord, '_r')
                        coordshort = coord(1:end-2);
                        coordleft = strcat(coordshort,'_l');

                        
                       
                        tempcol_r = tempTimeSeriesTable.getDependentColumn(java.lang.String(coord)).getAsMat();
                        tempcol_l = tempTimeSeriesTable.getDependentColumn(java.lang.String(coordleft)).getAsMat();
                        
                        % now take the rows that we want
    %                     col = tempcol(row_idx_start:row_idx_end);
                        col_r = tempcol_r;
                        tempcolinterp_r = interp1(timespercent, col_r, timespercent101);
                        col_l = tempcol_l;
                        tempcolinterp_l = interp1(timespercent, col_l, timespercent101);
                        

                        % match the right and lefts
                        % [rm, ixr] = max(tempcolinterp_r);
                        % [lm, ixl] = max(tempcolinterp_l);
                        newtempcolinterp_r = [tempcolinterp_r(ixr:end); tempcolinterp_r(1:ixr-1)];
                        newtempcolinterp_l = [tempcolinterp_l(ixl:end); tempcolinterp_l(1:ixl-1)];

                        % figure();
                        % hold on;
                        % plot(newtempcolinterp_r); plot(newtempcolinterp_l);

                        maerror = zeros(length(newtempcolinterp_l), 1);
                        % get the MAE
                        for h=1:length(newtempcolinterp_l)
                            % get the mae at each time point
%                             maerror(h) = abs(newtempcolinterp_l(h) - newtempcolinterp_r(h));
                            maerror(h) = (newtempcolinterp_l(h) - newtempcolinterp_r(h))^2;
                        end

                        % add to the struct
                        if ~isfield(welknaturalstruct, coordshort)
                            welknaturalstruct.(genvarname(coordshort)) = [];
                        end
                        if ~isfield(naturalstruct_combine, coordshort)
                            naturalstruct_combine.(genvarname(coordshort)) = [];
                        end
                        if ~isfield(totalstruct_combine, coordshort)
                            totalstruct_combine.(genvarname(coordshort)) = [];
                        end
                        welknaturalstruct.(genvarname(coordshort)) = [welknaturalstruct.(genvarname(coordshort)), maerror];
                        naturalstruct_combine.(genvarname(coordshort)) = [naturalstruct_combine.(genvarname(coordshort)), maerror];
                        totalstruct_combine.(genvarname(coordshort)) = [totalstruct_combine.(genvarname(coordshort)), maerror];

                        totalsubjsstruct_combine.(genvarname(subject)) = [totalsubjsstruct_combine.(genvarname(subject)), maerror];
                    end
    
%                     end
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
    % switching to RMS - that's what I have throughout the paper.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get MAE across all time for each subject across all coords
    welkallcoordsexo_combine = struct();
    subjs = fields(welkexostruct_combine)
    for subj=1:length(subjs)
        subject = char(subjs(subj));
        % loop through coords and get the MAE
        coordinates = fields(welkexostruct_combine.(genvarname(subject)));
        welkallcoordsexo_combine.(genvarname(subject)) = mean([ ...
            mean(mean(welkexostruct_combine.(genvarname(subject)).hip_flexion,2)), ...
            mean(mean(welkexostruct_combine.(genvarname(subject)).hip_adduction,2)), ...
            mean(mean(welkexostruct_combine.(genvarname(subject)).hip_rotation,2)), ...
            mean(mean(welkexostruct_combine.(genvarname(subject)).knee_angle,2)), ...
            mean(mean(welkexostruct_combine.(genvarname(subject)).ankle_angle,2)), ...
            mean(mean(welkexostruct_combine.(genvarname(subject)).subtalar_angle,2)), ...
            ]);
    end
    % same for natural
    welkallcoordsnat_combine = struct();
    subjs = fields(welknaturalstruct_combine)
    for subj=1:length(subjs)
        subject = char(subjs(subj));
        % loop through coords and get the MAE
        coordinates = fields(welknaturalstruct_combine.(genvarname(subject)));
        welkallcoordsnat_combine.(genvarname(subject)) = mean([ ...
            mean(mean(welknaturalstruct_combine.(genvarname(subject)).hip_flexion,2)), ...
            mean(mean(welknaturalstruct_combine.(genvarname(subject)).hip_adduction,2)), ...
            mean(mean(welknaturalstruct_combine.(genvarname(subject)).hip_rotation,2)), ...
            mean(mean(welknaturalstruct_combine.(genvarname(subject)).knee_angle,2)), ...
            mean(mean(welknaturalstruct_combine.(genvarname(subject)).ankle_angle,2)), ...
            mean(mean(welknaturalstruct_combine.(genvarname(subject)).subtalar_angle,2)), ...
            ]);
    end

    % now find the average across all conditions
    totalmeanacrosscoords = struct();
    for subj=1:length(subjs)
        subject = char(subjs(subj));
        % find the conditions average for each subject
        totalmeanacrosscoords.(genvarname(subject)) = mean([ ...
            welkallcoordsexo_combine.(genvarname(subject)), ...
            welkallcoordsnat_combine.(genvarname(subject))
            ]);
    end

    % figure out which subject has the highest MAE across all coordinates and time and conditions
    % max is for welk009 with MAE of .0685 rad or 3.9 deg. 
    keyboard

    % figure out how to get the std for these as well 
    totalmeansacrosscoords = struct();
    totalstdsacrosscoords = struct();
    for subj=1:length(subjs)
        subject = char(subjs(subj));
        % find the MAE across all time and coordinates and conditions for
        % each subject
        temp = mean(totalsubjsstruct_combine.(genvarname(subject)),1);
        temp2 = sqrt(temp);

        meanRMS = mean(temp2);
        stdRMS = std(temp2);
        totalmeansacrosscoords.(genvarname(subject)) = meanRMS;
        totalstdsacrosscoords.(genvarname(subject)) = stdRMS;


%         totalmeansacrosscoords.(genvarname(subject)) = sqrt(mean(totalsubjsstruct_combine.(genvarname(subject))(:)));
%         totalstdsacrosscoords.(genvarname(subject)) = sqrt(std(totalsubjsstruct_combine.(genvarname(subject))(:)));
    end

    % RMS 
    % welk009 -mean- .0779 = 4.5 deg, 
    % welk009 -std-  .0554 = 3.1742


    % manually find the biggest and get the stds as well
    % .0685 for welk009 or 3.9 degrees
    % std is .0662 or 3.7 degrees
    keyboard


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get MAE across all time all subjects for each coord
    totalstdacrosssubjs = struct(); % figure this out

    exoacrosssubjs = struct();
    exoacrosssubjsstd = struct();
    exoacrosssubjsfields = fields(exostruct_combine);
    for each=1:length(exoacrosssubjsfields)
        exoacrosssubjs.(genvarname(char(exoacrosssubjsfields(each)))) = mean(mean(exostruct_combine.(genvarname(char(exoacrosssubjsfields(each)))),2));
        exoacrosssubjsstd.(genvarname(char(exoacrosssubjsfields(each)))) = std(exostruct_combine.(genvarname(char(exoacrosssubjsfields(each))))(:));
    end
    % same for natural 
    natacrosssubjs = struct();
    natacrosssubjsstd = struct();
    natacrosssubjsfields = fields(naturalstruct_combine);
    for each=1:length(natacrosssubjsfields)
        natacrosssubjs.(genvarname(char(natacrosssubjsfields(each)))) = mean(mean(naturalstruct_combine.(genvarname(char(natacrosssubjsfields(each)))),2));
        natacrosssubjsstd.(genvarname(char(natacrosssubjsfields(each)))) = std(naturalstruct_combine.(genvarname(char(natacrosssubjsfields(each))))(:));
    end
    % get average across all conditions
    totalmeanacrosssubjs = struct();
    for each=1:length(exoacrosssubjsfields)
        totalmeanacrosssubjs.(genvarname(char(exoacrosssubjsfields(each)))) = mean([exoacrosssubjs.(genvarname(char(exoacrosssubjsfields(each)))), natacrosssubjs.(genvarname(char(natacrosssubjsfields(each))))]);
    end

    % find the coordinate with max MAE across all subjects and conditions
    % it is the ankle angle at .0753 rad or 4.3 deg
    keyboard
    % compare the total struct
    totalmeanacrosssubj = struct();
    totalstdacrosssubj = struct();
    totalfields = fields(totalstruct_combine);
    for each=1:length(totalfields)
        temp = mean(totalstruct_combine.(genvarname(char(totalfields(each)))),1);
        temp2 = sqrt(temp);

        meanRMS = mean(temp2);
        stdRMS = std(temp2);
        totalmeanacrosssubj.(genvarname(char(totalfields(each)))) = meanRMS;
        totalstdacrosssubj.(genvarname(char(totalfields(each)))) = stdRMS;
        
        
%         totalmeanacrosssubj.(genvarname(char(totalfields(each)))) = mean(totalstruct_combine.(genvarname(char(totalfields(each))))(:))
%         totalstdacrosssubj.(genvarname(char(totalfields(each)))) = std(totalstruct_combine.(genvarname(char(totalfields(each))))(:))
    end

    % RMS
    % mean ankle - .0883 = 5.0592 deg
    % std ankle  - .0304 = 1.7418 deg
    




    % find the max one - still ankle with .0753 or 4.3 degrees
    % for the ankle angle the std is .0551 or 3.157 degrees
    keyboard



    % should consider the standard dev. 
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
