% written by Jon Stingel

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
welkconditions = {'welkexo','welkexoexo'}; % ,'welknaturalslow','welknaturalnatural', ...
                  % 'welknaturalexo','welkexonatural','welkexoexo','welkexofast'};
welkaltconditions = {'welknatural','welknaturalnatural'};
welksubjects = {'welk001'};
load 'C:\Users\JP\code\repos\Stanford\delplab\projects\muscleModel\muscleEnergyModel\subjectgaitcycles.mat';


exoMetabolicsInd = [];
exoNames = [];

% loop through subjects
for subj=1:length(welksubjects)
    subject = char(welksubjects(subj));
    subjdir = strcat(resultsdir, strcat('/',subject));
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
    
    % lets loop through and get rid of the ones that we don't want right
    % now
    exoMetabolicsAvg_new = [];
    naturalMetabolicsAvg_new = [];
    names_new = [];
    for i=1:length(exoNames)
       if contains(char(exoNames(i)), 'combined')
           exoMetabolicsAvg_new = [exoMetabolicsAvg_new, exoMetabolicsAvg(i)];
           naturalMetabolicsAvg_new = [naturalMetabolicsAvg_new, naturalMetabolicsAvg(i)];
           names_new = [names_new, exoNames(i)];
       end
    end
    keyboard
    % look at differences for the values that we care about and some how
    % save them for printing
    exoMinusNatural = exoMetabolicsAvg_new - naturalMetabolicsAvg_new;
    exoMinusNatural_percDiff = (exoMetabolicsAvg_new - naturalMetabolicsAvg_new) ./ (naturalMetabolicsAvg_new) .*100;
    
    [smallissaving, savingix] = sort(exoMinusNatural);
    names_ordered = [];
    for i=1:length(savingix)
        names_ordered = [names_ordered, names_new(savingix(i))];
    end
    names_ordered_2 = string(names_ordered);
    names_ordered_3 = categorical(names_ordered_2);
    names_ordered_4 = reordercats(names_ordered_3, names_ordered);
    bar(names_ordered_4, smallissaving)
    
    % do the same for the percent change
    [smallissaving_perc, savingix_perc] = sort(exoMinusNatural_percDiff);
    names_ordered_perc = [];
    for i=1:length(savingix)
        names_ordered_perc = [names_ordered_perc, names_new(savingix(i))];
    end
    names_ordered_perc_2 = string(names_ordered_perc);
    names_ordered_perc_3 = categorical(names_ordered_perc_2);
    names_ordered_perc_4 = reordercats(names_ordered_perc_3, names_ordered_perc);
    bar(names_ordered_perc_4, smallissaving_perc)

    % make a figure with the biggest savers but the whole values
    savecompare1 = [];
    savecompare2 = [];
    savecomparenames = [];
    for i=1:4
        savecompare1 = [savecompare1, naturalMetabolicsAvg_new(savingix(i))];
        savecompare2 = [savecompare2, exoMetabolicsAvg_new(savingix(i))];
        savecomparenames = [savecomparenames, names_new(savingix(i))];
    end
    
    y = [savecompare1;savecompare2];
    x = [savecomparenames;savecomparenames];
    bar(y.','stacked')
    xticklabels({'Vas. Lat.','Gas. Med.','Add. Long.','Semimem.'})
    ylabel('Metabolic Cost [W/kg]')
    legend('Exotendon Running','Natural Running')

end

disp('finished all the sims!')
cd(resultsdir)
save('bigissuesfile.mat','Issues');
disp('end')