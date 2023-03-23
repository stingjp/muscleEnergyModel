% written by Jon Stingel

repodir = 'G:\Shared drives\Exotendon\muscleModel\muscleEnergyModel';
resultsdir = strcat(repodir, '/../results');
cd(resultsdir)


welkconditions = {'welkexo'}; % ,'welknaturalslow','welknaturalnatural', ...
                  % 'welknaturalexo','welkexonatural','welkexoexo','welkexofast'};
welkaltconditions = {'welknatural'};
welksubjects = {'welk002','welk003','welk005','welk008','welk009','welk010','welk013'};
load 'G:\Shared drives\Exotendon\muscleModel\muscleEnergyModel\subjectgaitcycles.mat';



exoMetabolicsAvg_new = [];
naturalMetabolicsAvg_new = [];

exoMetabolicsInd_allsubj = [];
naturalMetabolicsInd_allsubj = [];

holdingsubjects_natural = [];
holdingsubjects_exo = [];


% loop through subjects
for subj=1:length(welksubjects)
    subject = char(welksubjects(subj));
    subjdir = strcat(resultsdir, strcat('/',subject));
    
    exoMetabolicsInd = [];
    exoNames = [];
    
    % loop through conditions
    for cond=1:length(welkconditions)
        condition = char(welkconditions(cond));
        conddir = strcat(subjdir, strcat('/',condition));
        trials = fieldnames(subjectgaitcycles.(genvarname(subject)).(genvarname(condition)));
        % loop the trials
        for trial=1:length(trials)
            test = char(trials(trial));
            trialdir = strcat(conddir, strcat('/',test));
            cd(trialdir)
            disp(trialdir)
            temptable = readtable('muscleMetabolicsAll.csv');
            exoNames = temptable{:,'Var2'};
            exoMetabolicsInd = [exoMetabolicsInd, temptable{:,'Var1'}];            
        end
    end
    % done with the exo conditions
    
    % do the same thing for the natural cases
    naturalMetabolicsInd = [];
    naturalNames = [];

    % loop through conditions
    for cond=1:length(welkaltconditions)
        condition = char(welkaltconditions(cond));
        conddir = strcat(subjdir, strcat('/',condition));
        trials = fieldnames(subjectgaitcycles.(genvarname(subject)).(genvarname(condition)));
        % loop the trials
        for trial=1:length(trials)
            test = char(trials(trial));
            trialdir = strcat(conddir, strcat('/',test));
            cd(trialdir)
            disp(trialdir)
            temptable = readtable('muscleMetabolicsAll.csv');
            naturalNames = temptable{:,'Var2'};
            naturalMetabolicsInd = [naturalMetabolicsInd, temptable{:,'Var1'}];
        end
    end

    % add the individual subject gait cycles to the full matrix
    exoMetabolicsInd_allsubj = [exoMetabolicsInd_allsubj, exoMetabolicsInd(:,1:4)];
    naturalMetabolicsInd_allsubj = [naturalMetabolicsInd_allsubj, naturalMetabolicsInd(:,1:4)];

    
    % need to average for the two different conditions across the trials
    % for each 
    % then compute the differences between the exo and natural
    % then plot somehow
    
    % take averages 
    exoMetabolicsAvg = [];
    naturalMetabolicsAvg = [];
    for i=1:size(exoMetabolicsInd,1)
        exoMetabolicsAvg = [exoMetabolicsAvg, mean(exoMetabolicsInd(i,:))];
        naturalMetabolicsAvg = [naturalMetabolicsAvg, mean(naturalMetabolicsInd(i,:))];
    end
    % if you want to look at individual values and cycles. 
    % exoNames or natural names
    % exoMetabolicsInd and naturalMetabolicsInd will give all the different values for 4 cycles
    % exoMetabolicsAvg and naturalMetabolicsAvg will give the values averaged for 4 cycles

    
    
    
    % here is where we are going to make the figures for all the muscles
    % for each of the individual subjects
    % have 4 different things for each muscle, 40 muscles, 4 gait cycles, 2
    % conditions
    numMuscles = 40;
    f = 1;
    target = size(exoMetabolicsInd,1);
%     while f<=target
% 
%         % loop each muscle and make a figure
%         tempfig = figure(f); 
%         set(tempfig, 'Position', [1,1,1920,1080]);
%         % x axis is going to be the 4 gait cycles
%         % will do subplots for the different metabolic components
%         % also have multiple conditions (double bar)
%         
%         % y axis is going to be metabolic cost
%         tempvec_exo = exoMetabolicsInd(f:f+3,:);
%         tempvec_nat = naturalMetabolicsInd(f:f+3,:);
%         
%         subplot(2,2,1);
%         % combined metabolics
%         bar([tempvec_nat(1,:); tempvec_exo(1,:)]');
% %         yl = ylim;
%         xlabel('Gait cycle #')
%         ylabel('Total Metabolic Cost [W/kg]')
% %         ylim([0,max([tempvec_nat(1,:), tempvec_exo(1,:)])+.03]);
%         legend('Natural Running', 'Exotendon Running','location','northoutside')
%         
%         
%         % now for the activation cost
%         subplot(2,2,2);
%         bar([tempvec_nat(2,:); tempvec_exo(2,:)]')
%         xlabel('Gait cycle #');
%         ylabel('Activation Cost [W/kg]');
% %         ylim(yl)
% %         ylim([0,max([tempvec_nat(2,:), tempvec_exo(2,:)])+.03]);
%         legend('Natural Running', 'Exotendon Running','location','northoutside')
% 
%         % now for the shortening and lengthening cost
%         subplot(2,2,3);
%         bar([tempvec_nat(3,:); tempvec_exo(3,:)]')
%         xlabel('Gait cycle #');
%         ylabel('Shortening/Lengthening Cost [W/kg]');
% %         ylim(yl)
% %         ylim([0,max([tempvec_nat(3,:), tempvec_exo(3,:)])+.03]);
%         legend('Natural Running', 'Exotendon Running','location','northoutside')
% 
%         % now for the mechanical work rate
%         subplot(2,2,4);
%         bar([tempvec_nat(4,:); tempvec_exo(4,:)]')
%         xlabel('Gait cycle #');
%         ylabel('Mechanical Work rate [W/kg]');
% %         ylim(yl)
% %         ylim([0,max([tempvec_nat(4,:), tempvec_exo(4,:)])+.03]);
%         legend('Natural Running', 'Exotendon Running','location','northoutside')
%         
%         musc = char(exoNames(f));
%         sgtitle([subject, ' ', musc(21:end-6)]);
%         savefig(strcat('G:\Shared drives\Exotendon\muscleModel\analysis\',subject,'\',musc(21:end-6), '_metabolicBreakdown.fig'));
%         print(tempfig, ...
%             strcat('G:\Shared drives\Exotendon\muscleModel\analysis\',subject,'\',musc(21:end-6), '_metabolicBreakdown.png'),...
%             '-dpng', '-r500')
%         %%% oh crap need to do this for each of the muscles... need the
%         %%% title to be the muscle... 
%         
%         %%% TODO fix 4 indexing
%         f = f+4;
%         f
%         close
%     end
    
    
    % now to do something with both subjects
    % maybe need to add the things to another data structure for safe keep
    holdingsubjects_natural = [holdingsubjects_natural, naturalMetabolicsAvg'];
    holdingsubjects_exo = [holdingsubjects_exo, exoMetabolicsAvg'];
    
end


% now make a figure with all the subjects for each muscle

% probably want to average across gait cycles for each subject, 
% and then plot each of the subjects on x axis

keyboard
numMuscles = 40;
f = 1;
target = size(exoMetabolicsInd,1);
while f<=target

    % loop each muscle and make a figure
    tempfig = figure(f); 
    set(tempfig, 'Position', [1,1,1920,1080]);
    % x axis is going to be the 4 gait cycles
    % will do subplots for the different metabolic components
    % also have multiple conditions (double bar)

    % y axis is going to be metabolic cost
    tempvec_exo = holdingsubjects_exo(f:f+3,:);
    tempvec_nat = holdingsubjects_natural(f:f+3,:);

    subplot(2,2,1);
    % combined metabolics
    bar([tempvec_nat(1,:); tempvec_exo(1,:)]');
%     yl = ylim;
    xlabel('Subject #')
    ylabel('Total Metabolic Cost [W/kg]')
%         ylim([0,max([tempvec_nat(1,:), tempvec_exo(1,:)])+.03]);
    legend('Natural Running', 'Exotendon Running','location','northoutside')
    grid on;
    title(strcat('avg nat: ', string(mean(tempvec_nat(1,:),2)), ' avg exo: ', string(mean(tempvec_exo(1,:),2)), ...
        ' diff: ', string(mean(tempvec_nat(1,:),2)-mean(tempvec_exo(1,:),2))));

    % now for the activation cost
    subplot(2,2,2);
    bar([tempvec_nat(2,:); tempvec_exo(2,:)]')
    xlabel('Subject #');
    ylabel('Activation Cost [W/kg]');
%     ylim(yl);
%         ylim([0,max([tempvec_nat(2,:), tempvec_exo(2,:)])+.03]);
    legend('Natural Running', 'Exotendon Running','location','northoutside')
    grid on;
    title(strcat('avg nat: ', string(mean(tempvec_nat(2,:),2)), ' avg exo: ', string(mean(tempvec_exo(2,:),2)), ...
        ' diff: ', string(mean(tempvec_nat(2,:),2)-mean(tempvec_exo(2,:),2))));

    % now for the shortening and lengthening cost
    subplot(2,2,3);
    bar([tempvec_nat(3,:); tempvec_exo(3,:)]')
    xlabel('Subject #');
    ylabel('Shortening/Lengthening Cost [W/kg]');
%     ylim(yl);
%         ylim([0,max([tempvec_nat(3,:), tempvec_exo(3,:)])+.03]);
    legend('Natural Running', 'Exotendon Running','location','northoutside')
    grid on;
    title(strcat('avg nat: ', string(mean(tempvec_nat(3,:),2)), ' avg exo: ', string(mean(tempvec_exo(3,:),2)), ...
        ' diff: ', string(mean(tempvec_nat(3,:),2)-mean(tempvec_exo(3,:),2))));
    % now for the mechanical work rate
    subplot(2,2,4);
    bar([tempvec_nat(4,:); tempvec_exo(4,:)]')
    xlabel('Subject #');
    ylabel('Mechanical Work rate [W/kg]');
%     ylim(yl);
%         ylim([0,max([tempvec_nat(4,:), tempvec_exo(4,:)])+.03]);
    legend('Natural Running', 'Exotendon Running','location','northoutside')
    grid on;
    title(strcat('avg nat: ', string(mean(tempvec_nat(4,:),2)), ' avg exo: ', string(mean(tempvec_exo(4,:),2)), ...
        ' diff: ', string(mean(tempvec_nat(4,:),2)-mean(tempvec_exo(4,:),2))));
    
    musc = char(exoNames(f));
    sgtitle(['Subjects Average ', ' ', musc(21:end-6)]);
    savefig(strcat('G:\Shared drives\Exotendon\muscleModel\analysis\muscleMetabolics','\',musc(21:end-6), '_metabolicBreakdown_subjects.fig'));
    print(tempfig, ...
        strcat('G:\Shared drives\Exotendon\muscleModel\analysis\muscleMetabolics','\',musc(21:end-6), '_metabolicBreakdown_subjects.png'),...
        '-dpng', '-r500')
   
    %%% TODO fix 4 indexing
    f = f+4;
    f
    close
end



%% have all the values, now have to compute differences



% do we want to plot differences of any kind?

% % first take the difference between exo and natural for each subject
% muscleDifferences_raw = [];
% muscleDifferences_perc = [];
% 
% for s = 1:length(welksubjects)
%     alldifferences_raw = exoMetabolicsAvg_new(s,:) - naturalMetabolicsAvg_new(s,:);
%     alldifferences_perc = (exoMetabolicsAvg_new(s,:) - naturalMetabolicsAvg_new(s,:)) ./ (naturalMetabolicsAvg_new(s,:)) .*100;
%     
%     muscleDifferences_raw = [muscleDifferences_raw; alldifferences_raw];
%     muscleDifferences_perc = [muscleDifferences_perc; alldifferences_perc];
% end
% 
% 
% 
% % then average the differences
% muscleDifferencesAvg_raw = mean(muscleDifferences_raw, 1);
% muscleDifferencesAvg_perc = mean(muscleDifferences_perc, 1);
% 
% muscleDifferencesSTD_raw = std(muscleDifferences_raw, 0, 1);
% muscleDifferencesSTD_perc = std(muscleDifferences_perc, 0, 1);
% 
% 
% % take individual condition averages for plots
% naturalMetabolicsAvg_new_avg = mean(naturalMetabolicsAvg_new,1);
% exoMetabolicsAvg_new_avg = mean(exoMetabolicsAvg_new, 1);
% 
% %% not sure if I want to use the differences or just the raw values. 
% % might be good to do a figure for each subject that is raw values and then
% % from there we can do averaged raw values for the subject comparison. Am I
% % going a figure for each muscle??
% 
% % .... yes?
% 
% keyboard
% % here is where we make the figure for multiple subject comparisons
