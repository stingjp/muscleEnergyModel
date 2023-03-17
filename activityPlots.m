% written by Jon Stingel
% 20210329
% gather and plot all the muscle activities from simulations. 
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
welkexoconditions = {'welkexo'};%,'welkexoexo'}; % ,'welknaturalslow','welknaturalnatural', ...
                  % 'welknaturalexo','welkexonatural','welkexoexo','welkexofast'};
welknaturalconditions = {'welknatural'};%,'welknaturalnatural'};
% welksubjects = {'welk002','welk003','welk005','welk007','welk008','welk009','welk010','welk013'};
welksubjects = {'welk002','welk003','welk005','welk008','welk009','welk010','welk013'}; % 'welk008
tag = 'muscletrack';
thingstoplot = {'excitation','activation'};

load 'G:\Shared drives\Exotendon\muscleModel\muscleEnergyModel\subjectgaitcycles.mat';



exomeans_excitation = struct();
naturalmeans_excitation = struct();
exomeans_activation = struct();
naturalmeans_activation = struct();

exopeaks_excitation = struct();
exopeaks_activation = struct();
naturalpeaks_excitation = struct();
naturalpeaks_activation = struct();

 
% loop through the subjects
for subj=1:length(welksubjects)
    subject = char(welksubjects(subj));
    subjdir = strcat(resultsdir, strcat('/',subject));
    
    exomeans_excitation.(genvarname(subject)) = [];
    naturalmeans_excitation.(genvarname(subject)) = [];
    exomeans_activation.(genvarname(subject)) = [];
    naturalmeans_activation.(genvarname(subject)) = [];
    
    exopeaks_excitation.(genvarname(subject)) = [];
    exopeaks_activation.(genvarname(subject)) = [];
    naturalpeaks_excitation.(genvarname(subject)) = [];
    naturalpeaks_activation.(genvarname(subject)) = [];

    
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
                    
                    tempfile = strcat(trialdir, strcat('/',tag,'_states_100con.sto'));
                    disp(tempfile)
                    % if strcmp(subject,'welk002') || strcmp(subject,'welk003')
                    %     tempfile = strcat(trialdir,'/muscleprescribe_states.sto');
                    % else
                    %     tempfile = strcat(trialdir, strcat('/',tag,'_states.sto'));
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
                    
                    tempfile = strcat(trialdir, strcat('/',tag,'_controls_100con.sto'));
                    disp(tempfile)
                    % if strcmp(subject,'welk002') || strcmp(subject,'welk003')
                    %     tempfile = strcat(trialdir,'/muscleprescribe_controls.sto');
                    % else
                    %     tempfile = strcat(trialdir, strcat('/',tag,'_controls.sto'));
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
                    
                    tempfile = strcat(trialdir, strcat('/',tag,'_states_100con.sto'));

                    % if strcmp(subject,'welk002') || strcmp(subject,'welk003')
                    %     tempfile = strcat(trialdir,'/muscleprescribe_states.sto');
                    % else
                    %     tempfile = strcat(trialdir, strcat('/',tag,'_states.sto'));
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
                    
                    tempfile = strcat(trialdir, strcat('/',tag,'_controls_100con.sto'));

                    % if strcmp(subject,'welk002') || strcmp(subject,'welk003')
                    %     tempfile = strcat(trialdir,'/muscleprescribe_controls.sto');
                    % else
                    %     tempfile = strcat(trialdir, strcat('/',tag,'_controls.sto'));
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
        
        % now need to loop through both natural and exo to find the 3 glutes
        disp('can do weighted avgs here')
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
            if any(strcmp(glutemax, templabel_nat))
                tempglute = welknaturalstruct.(genvarname(templabel_nat));
                glutemax_data_nat = [glutemax_data_nat, tempglute];
            end
            if any(strcmp(glutemed, templabel_nat))
                tempglute = welknaturalstruct.(genvarname(templabel_nat));
                glutemed_data_nat = [glutemed_data_nat, tempglute];
            end
            if any(strcmp(glutemin, templabel_nat))
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
            if any(strcmp(glutemax, templabel_exo))
                tempglute = welkexostruct.(genvarname(templabel_exo));
                glutemax_data_exo = [glutemax_data_exo, tempglute];
            end
            if any(strcmp(glutemed, templabel_exo))
                tempglute = welkexostruct.(genvarname(templabel_exo));
                glutemed_data_exo = [glutemed_data_exo, tempglute];
            end
            if any(strcmp(glutemin, templabel_exo))
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
        
        
        % now create a figure 
        % tempfig = figure('Position',[1,1,1920,1080]);
        % for i=2:length(testlabels_nat)
        %     subplot(7,9,i-1);
        %     templabel = char(testlabels_nat(i));
        %     muscleplot1 = welknaturalstruct.(genvarname(templabel));
        %     plot(welknaturalstruct.time, muscleplot1, ':')
        %     hold on;
        %     plot(welknaturalstruct.time, mean(muscleplot1,2), 'k-', 'LineWidth', 2)
        %     title(templabel)
        %     xlabel('% gait cycle')
        %     ylabel(tempthing)
        %     grid on;
        % end
        % print(tempfig, ...
        %     strcat(repodir,'\..\analysis\',subject,'\', tempthing, '_natural', '.png'),...
        %     '-dpng', '-r500')
        % disp('print 1')
        
        
        % tempfig2 = figure('Position',[1,1,1920,1080]);
        % for i=2:length(testlabels_exo)
        %     subplot(7,9,i-1);
        %     templabel = char(testlabels_exo(i));
        %     muscleplot2 = welkexostruct.(genvarname(templabel));
        %     plot(welkexostruct.time, muscleplot2, ':')
        %     hold on;
        %     plot(welkexostruct.time, mean(muscleplot2,2), 'k-', 'LineWidth', 2)
        %     title(templabel)
        %     xlabel('% gait cycle')
        %     ylabel(tempthing)
        %     grid on;
        % end
        % print(tempfig2, ...
        %     strcat(repodir,'\..\analysis\',subject,'\', tempthing, '_exo', '.png'),...
        %     '-dpng', '-r500')
        % disp('print 2')
        
        
        % combined exo and natural fig
        % combineexonaturalfig = figure('Position',[1,1,1920,1080]);
        % title('red=natural, blue=exo');
        labels = fields(welkexostruct);
        for i=2:length(labels)
            % subplot(7,9,i-1);
            templabel = char(labels(i));
            muscleplot2 = welkexostruct.(genvarname(templabel));
            muscleplot1 = welknaturalstruct.(genvarname(templabel));
            % % plot(welkexostruct.time, muscleplot1, ':')
            % hold on;
            % plot(welkexostruct.time, mean(muscleplot2,2), 'b-', 'LineWidth', 2)
            % plot(welknaturalstruct.time, mean(muscleplot1,2), 'r-', 'LineWidth', 2)
            % title(templabel)
            % xlabel('% gait cycle')
            % ylabel(tempthing)
            % grid on;
            % % legend('exo','natural')
            if tempthing == 'excitation'
                exomeans_excitation.(genvarname(subject)) = [exomeans_excitation.(genvarname(subject)), mean(muscleplot2, 2)];
                naturalmeans_excitation.(genvarname(subject)) = [naturalmeans_excitation.(genvarname(subject)), mean(muscleplot1, 2)];
                excitelabels = fields(welkexostruct);
                % exomeans_excitation. = [exomeans_excitation, mean(muscleplot2,2)];
                % naturalmeans_excitation = [naturalmeans_excitation, mean(muscleplot1, 2)];

                exopeaks_excitation.(genvarname(subject)) = [exopeaks_excitation.(genvarname(subject)), max(mean(muscleplot2, 2))];
                naturalpeaks_excitation.(genvarname(subject)) = [naturalpeaks_excitation.(genvarname(subject)), max(mean(muscleplot1, 2))];

            end
            if tempthing == 'activation'
                exomeans_activation.(genvarname(subject)) = [exomeans_activation.(genvarname(subject)), mean(muscleplot2, 2)];
                naturalmeans_activation.(genvarname(subject)) = [naturalmeans_activation.(genvarname(subject)), mean(muscleplot1, 2)];
                activelabels = fields(welkexostruct);
                % exomeans_activation = [exomeans_activation, mean(muscleplot2, 2)];
                % naturalmeans_activation = [naturalmeans_activation, mean(muscleplot1, 2)];

                exopeaks_activation.(genvarname(subject)) = [exopeaks_activation.(genvarname(subject)), max(mean(muscleplot2, 2))];
                naturalpeaks_activation.(genvarname(subject)) = [naturalpeaks_activation.(genvarname(subject)), max(mean(muscleplot1, 2))];

            end
        end
        % print(combineexonaturalfig, ...
        %     strcat(repodir,'\..\analysis\',subject,'\', tempthing, '_combined', '.png'),...
        %     '-dpng', '-r500')
        disp('print combined')
    end
end


% keyboard
% combined exo and natural fig - averages for subjects and mean
subjectcombineexonaturalfig1 = figure('Position',[1,1,1200,1920]);
title('red=natural, blue=exo');
% do more stuff
% averaging and whatnot
    

for i=1:length(excitelabels)-1    
    % make a subplot
    subplot(16,4,i);
    templabel = char(excitelabels(i+1));
    hold on;
    tempsubjavgs1 = [];
    tempsubjavgs2 = [];
    
    % for each muscle go through each subject
    for subj=1:length(welksubjects)
        subject = char(welksubjects(subj));
        plot(welkexostruct.time, exomeans_excitation.(genvarname(subject))(:,i), 'b:', 'LineWidth', 0.4)
        plot(welknaturalstruct.time, naturalmeans_excitation.(genvarname(subject))(:,i), 'r:', 'LineWidth', 0.4)
        
        % add each signal to the temp vector to get means before moving to
        % the next muscle
        tempsubjavgs1 = [tempsubjavgs1, naturalmeans_excitation.(genvarname(subject))(:,i)];
        tempsubjavgs2 = [tempsubjavgs2, exomeans_excitation.(genvarname(subject))(:,i)];
    end
    
    % now plot the means from the tempsubjavgs1/2
    plot(welkexostruct.time, mean(tempsubjavgs2, 2), 'b-', 'LineWidth', 2)
    plot(welknaturalstruct.time, mean(tempsubjavgs1, 2), 'r-', 'LineWidth', 2)
    
    
    title(templabel)
    xlabel('% gait cycle')
    ylabel(tempthing)
%     grid on;
    % legend('exo','natural')
%     legend(strcat('exo peak: ', num2str(max(mean(tempsubjavgs2, 2)))),strcat('nat peak: ',num2str(max(mean(tempsubjavgs1, 2)))))

end

print(subjectcombineexonaturalfig1, ...
    strcat(repodir, '\..\analysis\', 'excitation_all_subjects_combined_nolegend.png'),...
    '-dpng', '-r500')
disp('print combined')


% now for the activations
% combined exo and natural fig - averages for subjects and mean
subjectcombineexonaturalfig2 = figure(100); %'Number',100,'Position',[1,1,1920,1080]);
set(gcf,'WindowStyle','Docked','Position',[1,1,2480,3508])
% set(gcf,'Position',[1,1,700,4508])

title('red=natural, blue=exo');
% do more stuff
% averaging and whatnot
    
for i=1:length(activelabels)-1   
    % make a subplot
    subplot(11,4,i);
%     axis('square')
%     axis('tight')
    templabel = char(activelabels(i+1));
    hold on;
    tempsubjavgs1 = [];
    tempsubjavgs2 = [];
    
    % for each muscle go through each subject
    for subj=1:length(welksubjects)
        subject = char(welksubjects(subj));
        plot(welkexostruct.time, exomeans_activation.(genvarname(subject))(:,i), 'b:', 'LineWidth', 0.4)
        plot(welknaturalstruct.time, naturalmeans_activation.(genvarname(subject))(:,i), 'r:', 'LineWidth', 0.4)
        
        % add each signal to the temp vector to get means before moving to
        % the next muscle
        tempsubjavgs1 = [tempsubjavgs1, naturalmeans_activation.(genvarname(subject))(:,i)];
        tempsubjavgs2 = [tempsubjavgs2, exomeans_activation.(genvarname(subject))(:,i)];
    end
    
    disp(templabel)
    % spm tests for simulated activities
    addpath(genpath('G:\Shared drives\Exotendon\Matlab\common\')) ;
    spm = spm1d.stats.nonparam.ttest_paired(squeeze(tempsubjavgs1)', squeeze(tempsubjavgs2)');
    spmi   = spm.inference(0.05, 'two_tailed',true, 'interp',true);
    disp(spmi)
    tempfig = figure;
    spmi.plot();
    spmi.plot_threshold_label();
    spmi.plot_p_values();
    title(templabel)
    print(tempfig, strcat(repodir, '\..\analysis\activitySPM\',templabel,'.png'), ...
        '-dpng', '-r500')
    % now plot the means from the tempsubjavgs1/2
    figure(100);
    plot(welkexostruct.time, mean(tempsubjavgs2, 2), 'b-', 'LineWidth', 2)
    plot(welknaturalstruct.time, mean(tempsubjavgs1, 2), 'r-', 'LineWidth', 2)
    
    
    title(templabel)
    xlabel('% gait cycle')
    ylabel(tempthing)
%     grid on;
    % legend('exo','natural')
%     legend(strcat('exo peak: ', num2str(max(mean(tempsubjavgs2, 2)))),strcat('nat peak: ',num2str(max(mean(tempsubjavgs1, 2)))))
end

print(subjectcombineexonaturalfig2, ...
    strcat(repodir, '\..\analysis\', 'activation_all_subjects_combined_nolegend.png'),...
    '-dpng', '-r500')
disp('print combined')
