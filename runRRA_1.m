function runRRA_1(filename)
    import org.opensim.modeling.*
    
    rra = RRATool(filename);

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

    % set the gait timings
    rra.setInitialTime(gait_start-0.02);
    rra.setFinalTime(gait_end+0.02);

    rra.run();

    disp('RRA run successful... files in results.')

end