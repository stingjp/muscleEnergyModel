clear all
clc



repodir = 'C:\Users\JP\code\repos\Stanford\delplab\projects\muscleModel\muscleEnergyModel';
resultsdir = strcat(repodir, '/../results');
cd(resultsdir)


% conditions
% walsconditions = ['walsslack','walslow','walsmed','walshigh','walsmax']
% welkconditions = ['welknatural','welkexo','welknaturalslow','welknaturalnatural',
%                   'welknaturalexo','welkexonatural','welkexoexo','welkexofast']
% jackconditions = ['jackpower1','jackpower2','jackpower3','jackpower4','jackpower5','jackpower6',
%                   'jacktau1','jacktau2','jacktau3','jacktau4','jacktau5']
% dembconditions = ['dembnoloadfree','dembnoloadslow','dembloadedfree','dembloadedmatched']
% sildconditions = ['sildbw0','sildbw5','sildbw10','sild10w0','sild10w5','sild10w10',
%                   'sild20w0','sild20w5','sild20w10','sild30w0','sild30w5','sild30w10',
%                   'sildbwrun0','sild10wrun0','sild20wrun0','sild30wrun0']

dembconditions = {'dembloadedfree',};%'dembnoloadfree'};
dembsubjects = {'demb005','demb007','demb009','demb010'};

load 'C:\Users\JP\code\repos\Stanford\delplab\projects\muscleModel\muscleEnergyModel\subjectgaitcycles.mat';

for subj=1:length(dembsubjects)
    subject = char(dembsubjects(subj));
    subjdir = strcat(resultsdir, strcat('/',subject));
    
    for cond=1:length(dembconditions)
       condition = char(dembconditions(cond));
       conddir = strcat(subjdir, strcat('/',condition));
       trials = fieldnames(subjectgaitcycles.(genvarname(subject)).(genvarname(condition)));
        for trial=1:length(trials)
            test = char(trials(trial));
            trialdir = strcat(conddir, strcat('/',test));
            cd(trialdir)
            % run the analysis
            try
                analyzeSubject();
                disp('ran');
            catch
                disp('issue');
                disp(trialdir);
            end
        end
    end
end

