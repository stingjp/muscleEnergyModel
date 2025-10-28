% written by Jon Stingel
% 20211002
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
thingstoplot = {'NormalizedFiberLength','NormFiberVelocity'};

load 'G:\Shared drives\Exotendon\muscleModel\muscleEnergyModel\subjectgaitcycles.mat';

welkexo_peaklength = struct();
welknat_peaklength = struct();
welkexo_peakvel = struct();
welknat_peakvel = struct();




% loop through each of the things we want to plot
for thing=1:length(thingstoplot)
    tempthing = char(thingstoplot(thing))

    % create stucture for combined subject figures
    welknaturalstruct_combine = struct();
    welkexostruct_combine = struct();
    


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
                
                tempfile = strcat(trialdir, strcat('/analyzemuscles',tag,'_MuscleAnalysis_', tempthing, '.sto'));

                % if strcmp(subject, 'welk002') || strcmp(subject, 'welk003')
                %     tempfile = strcat(trialdir, strcat('/analyzemuscles_MuscleAnalysis_', strcat(tempthing, '.sto')));
                % else
                %     tempfile = strcat(trialdir, strcat('/analyzemuscles',tag,'_MuscleAnalysis_', tempthing, '.sto'));
                % end
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
                
                tempfile = strcat(trialdir, strcat('/analyzemuscles',tag,'_MuscleAnalysis_', tempthing, '.sto'));

                % if strcmp(subject, 'welk002') || strcmp(subject, 'welk003')
                %     tempfile = strcat(trialdir, strcat('/analyzemuscles_MuscleAnalysis_', strcat(tempthing, '.sto')));
                % else
                %     tempfile = strcat(trialdir, strcat('/analyzemuscles',tag,'_MuscleAnalysis_', tempthing, '.sto'));
                % end
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
        
        

        % tempfig = figure('Position',[1,1,1920,1080]);
        % % do more stuff
        % % averaging and whatnot
        % for i=0:(labels.size()/2)-1
        %     subplot(5,8,i+1);
        %     templabel = char(labels.get(i));
        %     muscleplot_nat = welknaturalstruct.(genvarname(char(templabel)));
        %     muscleplot_exo = welkexostruct.(genvarname(char(templabel)));
        %     plot(welknaturalstruct.time, muscleplot_nat, 'r:')
        %     hold on;
        %     plot(welkexostruct.time, muscleplot_exo, 'b:')
        %     plot(welknaturalstruct.time, mean(muscleplot_nat,2), 'r-', 'LineWidth', 1)
        %     plot(welkexostruct.time, mean(muscleplot_exo,2), 'b-', 'LineWidth', 1)
        %     title(templabel)
        %     xlabel('% gait cycle')
        %     ylabel(tempthing)
        %     grid on;
        % end
        % print(tempfig, ...
        %     strcat('G:\Shared drives\Exotendon\muscleModel\analysis\',subject,'\',tempthing, tag,'_acrossconditions', '.png'),...
        %     '-dpng', '-r500')
        % disp('print 1')
        

        % now need to loop through both natural and exo to find the 3 glutes
        labels_nat = fields(welknaturalstruct);
        glutemax = {'glmax1_r','glmax2_r','glmax3_r'};
        glutemed = {'glmed1_r','glmed2_r','glmed3_r'};
        glutemin = {'glmin1_r','glmin2_r','glmin3_r'};
        
        glutemax_data_nat = [];
        glutemed_data_nat = [];
        glutemin_data_nat = [];
        
        glutemax_data_exo = [];
        glutemed_data_exo= [];
        glutemin_data_exo = [];
        
        % loop the naturals first
        for i=1:length(labels_nat)
            templabel_nat = string(labels_nat(i));
            if any(contains(templabel_nat,glutemax))
                tempglute = welknaturalstruct.(genvarname(templabel_nat));
                glutemax_data_nat = [glutemax_data_nat, tempglute];
            end
            if any(contains(templabel_nat,glutemed))
                tempglute = welknaturalstruct.(genvarname(templabel_nat));
                glutemed_data_nat = [glutemed_data_nat, tempglute];
            end
            if any(contains(templabel_nat,glutemin))
                tempglute = welknaturalstruct.(genvarname(templabel_nat));
                glutemin_data_nat = [glutemin_data_nat, tempglute];
            end
        end
        glutemax_data_nat = mean(glutemax_data_nat, 2);
        glutemed_data_nat = mean(glutemed_data_nat, 2);
        glutemin_data_nat = mean(glutemin_data_nat, 2);
        
        labels_exo = fields(welkexostruct);
        % loop the exos now
        for i=1:length(labels_exo)
            templabel_exo = string(labels_exo(i));
            if any(contains(templabel_exo,glutemax))
                tempglute = welkexostruct.(genvarname(templabel_exo));
                glutemax_data_exo = [glutemax_data_exo, tempglute];
            end
            if any(contains(templabel_exo,glutemed))
                tempglute = welkexostruct.(genvarname(templabel_exo));
                glutemed_data_exo = [glutemed_data_exo, tempglute];
            end
            if any(contains(templabel_exo,glutemin))
                tempglute = welkexostruct.(genvarname(templabel_exo));
                glutemin_data_exo = [glutemin_data_exo, tempglute];
            end
        end
        glutemax_data_exo = mean(glutemax_data_exo, 2);
        glutemed_data_exo = mean(glutemed_data_exo, 2);
        glutemin_data_exo = mean(glutemin_data_exo, 2);

        
        % make sure the new averaged will get into figure
        welknaturalstruct.glmax_avg_r = glutemax_data_nat;
        welknaturalstruct.glmed_avg_r = glutemed_data_nat;
        welknaturalstruct.glmin_avg_r = glutemin_data_nat;
        
        welkexostruct.glmax_avg_r = glutemax_data_exo;
        welkexostruct.glmed_avg_r = glutemed_data_exo;
        welkexostruct.glmin_avg_r = glutemin_data_exo;
        
        % need to get new total labels
        testlabels_nat = fields(welknaturalstruct);
        testlabels_exo = fields(welkexostruct);              
        

        % add the subject average to the combined struct?
        welknaturalstruct_combine.(genvarname(subject)) = welknaturalstruct;
        welkexostruct_combine.(genvarname(subject)) = welkexostruct;
        
        if strcmp(char(tempthing), 'NormalizedFiberLength')
            welknat_peaklength.(genvarname(subject)) = [];
            welkexo_peaklength.(genvarname(subject)) = [];
        elseif strcmp(char(tempthing), 'NormFiberVelocity')
            welknat_peakvel.(genvarname(subject)) = [];
            welkexo_peakvel.(genvarname(subject)) = [];
        end

    end


    
    
    labels = fields(welkexostruct);
    % loop through the subjects again?
    % now plot across subjects
    tempfig2 = figure('Position',[1,1,1920,1080]);
        % then loop through the muscles inside each subject
    j = 2;
    for i=2:length(labels)
        templabel = labels(i);
        if ~contains(templabel,'_l')
            
            testnat = [];
            testexo = [];

            subplot(5,9,j-1);
            j = j+1;
            % loop through the subjects
            for subj=1:length(welksubjects)

                subject = char(welksubjects(subj));
                muscleplot_nat = welknaturalstruct_combine.(genvarname(subject)).(genvarname(char(templabel)));
                muscleplot_exo = welkexostruct_combine.(genvarname(subject)).(genvarname(char(templabel)));
                

                % add some stuff for post-process analysis
                if strcmp(char(tempthing), 'NormalizedFiberLength')
                    welkexo_peaklength.(genvarname(subject)) = [welkexo_peaklength.(genvarname(subject)), max(mean(muscleplot_exo,2))];
                    welknat_peaklength.(genvarname(subject)) = [welknat_peaklength.(genvarname(subject)), max(mean(muscleplot_nat,2))];
                elseif strcmp(char(tempthing), 'NormFiberVelocity')
                    welkexo_peakvel.(genvarname(subject)) = [welkexo_peakvel.(genvarname(subject)), min(mean(muscleplot_exo,2))];
                    welknat_peakvel.(genvarname(subject)) = [welknat_peakvel.(genvarname(subject)), min(mean(muscleplot_nat,2))];
                end
                
                
                testnat = [testnat, muscleplot_nat];
                testexo = [testexo, muscleplot_exo];
                % have all of them, want the average plotted for each subject
                plot(welknaturalstruct.time, mean(muscleplot_nat,2), 'r:','LineWidth',0.4)
                hold on;
                plot(welkexostruct.time, mean(muscleplot_exo,2), 'b:','LineWidth',0.4)
            end

            % plot the means here.
            plot(mean(testnat,2),'r','LineWidth',2)
            plot(mean(testexo,2),'b','LineWidth',2)

            title(templabel)
            xlabel('% gait cycle')
            ylabel(tempthing)
%             grid on;

            % try to get legend with peaks
%             legend(strcat('nat min: ',num2str(min(mean(testnat,2)))), ...
%                 strcat('exo min: ',num2str(min(mean(testexo,2)))), ...
%                 strcat('nat max: ',num2str(max(mean(testnat,2)))), ...
%                 strcat('exo max: ',num2str(max(mean(testexo,2)))));



        end
    end

    % _withlegend
    print(tempfig2, ...
        strcat('G:\Shared drives\Exotendon\muscleModel\analysis\', tempthing,tag, 'nolegend_combined', '.png'),...
        '-dpng', '-r500')
    disp('print 2')

    
end

