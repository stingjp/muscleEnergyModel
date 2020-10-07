function analyzeMetabolicCost(solution)
    import org.opensim.modeling.*

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for the unconstrained solution
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Conduct an analysis using MuscleAnalysis and ProbeReporter.
%     solution = MocoTrajectory('muscle_stateprescribe_grfprescribe_solution.sto');
    Time = solution.getTimeMat();
    numColPoints = solution.getNumTimes();
    

    keyboard
    computeIDFromResult(solution);

    % get the subject name and mass
    load 'C:\Users\JP\code\repos\Stanford\delplab\projects\muscleModel\muscleEnergyModel\subjectmass.mat';
    workdir = pwd;
    cd ../../
    tempdir = pwd;
    [~,subjectname,~] = fileparts(pwd);
    cd(workdir);
    model_mass = subjectmass.(genvarname(subjectname)); % kg
    
    % full moco method
    analyze = AnalyzeTool();
    analyze.setName("analyzemuscles");
    analyze.setModelFilename("post_simple_model_all_the_probes.osim");
    analyze.setStatesFileName("muscleprescribe_states.sto");
    analyze.updAnalysisSet().cloneAndAppend(MuscleAnalysis());
    analyze.updAnalysisSet().cloneAndAppend(ProbeReporter());
    analyze.updControllerSet().cloneAndAppend(PrescribedController("muscleprescribe_controls.sto"));
    analyze.setInitialTime(Time(1));
    analyze.setFinalTime(Time(end));
    analyze.print("testing_AnalyzeTool_setup.xml");
    % Run the analysis.
    analyze = AnalyzeTool("testing_AnalyzeTool_setup.xml");
    analyze.run();


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % tables
    table_activefiberforce = TimeSeriesTable("analyzemuscles_MuscleAnalysis_ActiveFiberForce.sto");
    table_fibervelocity = TimeSeriesTable("analyzemuscles_MuscleAnalysis_FiberVelocity.sto");
    table_metabolics = TimeSeriesTable('analyzemuscles_ProbeReporter_probes.sto');
    table_lMT = TimeSeriesTable('analyzemuscles_MuscleAnalysis_Length.sto');
    table_fiberlength = TimeSeriesTable('analyzemuscles_MuscleAnalysis_FiberLength.sto');
    
    
    % get time
    time_os = table_activefiberforce.getIndependentColumn();
    time = [];
    for i=0:time_os.size()-1
        time = [time; time_os.get(i)];
    end
    time_met_os = table_metabolics.getIndependentColumn();
    time_met = [];
    for i=0:time_met_os.size()-1
        time_met = [time_met; time_met_os.get(i)];
    end

    % get all the muscles
    muscles = [];
    labels = table_activefiberforce.getColumnLabels();
    for i=0:labels.size()-1
        muscles = [muscles, labels.get(i)];
    end
    numMuscles = length(muscles);
    
    % get active fiber force and fiber velocity
    activefiberforce = [];
    fibervelocity = [];
    fiberlength = [];
    for i=1:length(muscles)
        temp_activefiberforce = table_activefiberforce.getDependentColumn(muscles(i)).getAsMat();
        activefiberforce = [activefiberforce, temp_activefiberforce];
        temp_fibervelocity = table_fibervelocity.getDependentColumn(muscles(i)).getAsMat();
        fibervelocity = [fibervelocity, temp_fibervelocity];
        temp_fiberlength = table_fiberlength.getDependentColumn(muscles(i)).getAsMat();
        fiberlength = [fiberlength, temp_fiberlength];
    end
    

    % get excitations
    controlData = solution.getControlsTrajectoryMat();
    controlNames_os = solution.getControlNames();
    controlNames = [];
    for i = 0:controlNames_os.size()-1
        controlNames= [controlNames, controlNames_os.get(i)];
    end

    % get activations
    stateData = solution.getStatesTrajectoryMat();
    stateNames_os = solution.getStateNames();
    stateNames = [];
    for i=0:stateNames_os.size()-1
        stateNames = [stateNames, stateNames_os.get(i)];
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % workspace - working on the typical metabolic outputs that we have
    
    % get metabolics probe information
    metabolics_all_os = table_metabolics.getDependentColumn('all_metabolics_TOTAL');
    metabolics_act_os = table_metabolics.getDependentColumn('all_activation_maintenance_rate_TOTAL');
    metabolics_short_os = table_metabolics.getDependentColumn('all_shortening_rate_TOTAL');
    metabolics_basal_os = table_metabolics.getDependentColumn('all_basal_rate_TOTAL');
    metabolics_mech_os = table_metabolics.getDependentColumn('all_mechanical_work_rate_TOTAL');
    
    metabolics_all = metabolics_all_os.getAsMat;
    metabolics_act = metabolics_act_os.getAsMat;
    metabolics_short = metabolics_short_os.getAsMat;
    metabolics_basal = metabolics_basal_os.getAsMat;
    metabolics_mech = metabolics_mech_os.getAsMat;

    metabolics_all_avg = 2*((trapz(time, metabolics_all)) / (time(end)-time(1))) / model_mass;
    metabolics_act_avg = 2*((trapz(time, metabolics_act)) / (time(end)-time(1))) / model_mass;
    metabolics_short_avg = 2*((trapz(time, metabolics_short)) / (time(end)-time(1))) / model_mass;
    metabolics_basal_avg = 2*((trapz(time, metabolics_basal)) / (time(end)-time(1))) / model_mass;
    metabolics_mech_avg = 2*((trapz(time, metabolics_mech)) / (time(end)-time(1))) / model_mass;

    temp_reg = metabolics_all_avg

    met_rows = {'trial'};
    met_table = table(metabolics_all_avg, metabolics_act_avg, metabolics_short_avg,...
                metabolics_basal_avg, metabolics_mech_avg, 'RowNames', met_rows);
            
    writetable(met_table, 'metabolicsTable.csv','WriteRowNames',true);

end
