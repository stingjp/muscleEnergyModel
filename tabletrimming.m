function [temptable,starttime,renumber] = tabletrimming(filename, renumber)
    
    % get the subject name and gait timings
    load 'G:\Shared drives\Exotendon\muscleModel\muscleEnergyModel\subjectgaitcycles.mat';
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
    starttime = gait_start;

    if renumber==true
        gait_end = gait_end-gait_start;
        gait_start = 0;
        subjectgaitcycles.(genvarname(subjectname)).(genvarname(conditionname)).(genvarname(trialname)).initial = gait_start;
        subjectgaitcycles.(genvarname(subjectname)).(genvarname(conditionname)).(genvarname(trialname)).final = gait_end;
        save('G:\Shared drives\Exotendon\muscleModel\muscleEnergyModel\subjectgaitcycles.mat',"subjectgaitcycles");
    end
    renumber = ~renumber;
end

