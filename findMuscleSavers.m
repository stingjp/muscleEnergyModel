% written by Jon Stingel

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
welkconditions = {'welkexo'}; % ,'welknaturalslow','welknaturalnatural', ...
                  % 'welknaturalexo','welkexonatural','welkexoexo','welkexofast'};
welkaltconditions = {'welknatural'};
welksubjects = {'welk002','welk003','welk005','welk008','welk009','welk010','welk013'};
load 'G:\Shared drives\Exotendon\muscleModel\muscleEnergyModel\subjectgaitcycles.mat';



exoMetabolicsAvg_new = [];
naturalMetabolicsAvg_new = [];

exoMetabolicsInd_allsubj = [];
naturalMetabolicsInd_allsubj = [];

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

%     keyboard
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


    % lets loop through and get rid of the ones that we don't want right
    % now
    
    names_new = [];
    tempnew = [];
    tempnew2 = [];
    for i=1:length(exoNames)
        % this loop is grabbing all the values that are the full individual
        % muscle metabolic costs (not individual heat terms etc.)
        if contains(char(exoNames(i)), 'combined')
            % exoMetabolicsAvg_new = [exoMetabolicsAvg_new, exoMetabolicsAvg(i)];
            % naturalMetabolicsAvg_new = [naturalMetabolicsAvg_new, naturalMetabolicsAvg(i)];
            % names_new = [names_new, exoNames(i)];
            tempnew = [tempnew, exoMetabolicsAvg(i)];
            tempnew2 = [tempnew2, naturalMetabolicsAvg(i)];
            names_new = [names_new, exoNames(i)];    
        end
    end
    
    
    
    if size(exoMetabolicsAvg_new, 1) == 0
        exoMetabolicsAvg_new = tempnew;
        naturalMetabolicsAvg_new = tempnew2;
    else
        exoMetabolicsAvg_new = [exoMetabolicsAvg_new; tempnew];
        naturalMetabolicsAvg_new = [naturalMetabolicsAvg_new; tempnew2];
    end
end

%% have all the values, now have to compute differences
keyboard
% going to make a new structure that is the combined muscle paths/groups
disp('NOTE: These values are for single muscle, single leg.')




% first take the difference between exo and natural for each subject
muscleDifferences_raw = [];
muscleDifferences_perc = [];

for s = 1:length(welksubjects)
    alldifferences_raw = exoMetabolicsAvg_new(s,:) - naturalMetabolicsAvg_new(s,:);
    alldifferences_perc = (exoMetabolicsAvg_new(s,:) - naturalMetabolicsAvg_new(s,:)) ./ (naturalMetabolicsAvg_new(s,:)) .*100;
    
    muscleDifferences_raw = [muscleDifferences_raw; alldifferences_raw];
    muscleDifferences_perc = [muscleDifferences_perc; alldifferences_perc];
end



% then average the differences
muscleDifferencesAvg_raw = mean(muscleDifferences_raw, 1);
muscleDifferencesAvg_perc = mean(muscleDifferences_perc, 1);

muscleDifferencesSTD_raw = std(muscleDifferences_raw, 0, 1);
muscleDifferencesSTD_perc = std(muscleDifferences_perc, 0, 1);


% take individual condition averages for plots
naturalMetabolicsAvg_new_avg = mean(naturalMetabolicsAvg_new,1);
exoMetabolicsAvg_new_avg = mean(exoMetabolicsAvg_new, 1);


% now to figure out sorting of the values based on savings
[smallissaving, savingix] = sort(muscleDifferencesAvg_raw);
names_ordered = [];
for i=1:length(savingix)
    names_ordered = [names_ordered, names_new(savingix(i))];
end


names_ordered_2 = string(names_ordered);
names_ordered_3 = categorical(names_ordered_2);
names_ordered_4 = reordercats(names_ordered_3, names_ordered);
figure();
bar(names_ordered_4, smallissaving)

% do the same for the percent change
[smallissaving_perc, savingix_perc] = sort(muscleDifferencesAvg_perc);
names_ordered_perc = [];
for i=1:length(savingix)
    names_ordered_perc = [names_ordered_perc, names_new(savingix_perc(i))];
end
names_ordered_perc_2 = string(names_ordered_perc);
names_ordered_perc_3 = categorical(names_ordered_perc_2);
names_ordered_perc_4 = reordercats(names_ordered_perc_3, names_ordered_perc);
figure();
bar(names_ordered_perc_4, smallissaving_perc)




%% this figure doesn't really work, I think the averages make it skewed
% make a figure with the biggest savers - raw savings
savecompare1 = [];
savecompare2 = [];
savecomparenames = [];
% grab the biggest savers (differences)
for i=1:4
    % for the savers
    savecompare1 = [savecompare1, naturalMetabolicsAvg_new_avg(savingix(i))];
    savecompare2 = [savecompare2, exoMetabolicsAvg_new_avg(savingix(i))];
    savecomparenames = [savecomparenames, names_new(savingix(i))];
    
    % for the spenders
%     savecompare1 = [savecompare1, naturalMetabolicsAvg_new_avg(savingix(end-i+1))];
%     savecompare2 = [savecompare2, exoMetabolicsAvg_new_avg(savingix(end-i+1))];
%     savecomparenames = [savecomparenames, names_new(savingix(end-i+1))];
end

figure();
% for savers
bar(savecompare1, 'r');
hold on;
bar(savecompare2, 'b');

% for spenders
% bar(savecompare2, 'b');
% hold on;
% bar(savecompare1, 'r');



% >>>>>>> 709e47e1a1398634845ed6ca8a7074acca7255b8
xticklabels(savecomparenames)
ylabel('Metabolic Cost [W/kg]')
legend('Exotendon Running','Natural Running')

