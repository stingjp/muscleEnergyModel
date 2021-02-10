function [temptable] = tabletrimming(filename)
    
    % get the subject name and gait timings
    load 'C:\Users\JP\code\repos\Stanford\delplab\projects\muscleModel\muscleEnergyModel\subjectgaitcycles.mat';
    workdir = pwd;
    [~,trialname,~] = fileparts(pwd);
    cd ../
    [~,conditionname,~] = fileparts(pwd);
    cd ../
    [~,subjectname,~] = fileparts(pwd);
    cd(workdir);

    gait_start = subjectgaitcycles.(genvarname(subjectname)).(genvarname(conditionname)).(genvarname(trialname)).initial;
    gait_end = subjectgaitcycles.(genvarname(subjectname)).(genvarname(conditionname)).(genvarname(trialname)).final;
    gait_start = gait_start-0.02;
    gait_end = gait_end+0.02;

    % use the gait timings to load and save only a specific part of the files
    import org.opensim.modeling.*
    temptable = TimeSeriesTable(filename);

    temptable.trimTo(gait_end);
    temptable.trimFrom(gait_start);
end

