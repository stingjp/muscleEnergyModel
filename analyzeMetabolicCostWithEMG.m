function analyzeMetabolicCostWithEMG(solution)
    import org.opensim.modeling.*
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for the EMG constrained solution
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Conduct an analysis using MuscleAnalysis and ProbeReporter.
    % solution = MocoTrajectory('muscle_stateprescribe_grfprescribe_withemg_solution.sto');
    Time = solution.getTimeMat();
    numColPoints = solution.getNumTimes();
    
    keyboard
    %%%% TODO include the new stance and swing things for this script


    % get the subject name and mass
    load 'C:\Users\JP\code\repos\Stanford\delplab\projects\muscleModel\muscleEnergyModel\subjectmass.mat';
    
    workdir = pwd;
    [~,trialname,~] = fileparts(pwd);
    cd ../
    [~,condname,~] = fileparts(pwd);
    cd ../
    [~,subjectname,~] = fileparts(pwd);
    experimentname = subjectname(1:4);
    cd(workdir);
    model_mass = subjectmass.(genvarname(subjectname)); % kg

    % full moco method
    analyze = AnalyzeTool();
    analyze.setName("analyzemusclesEMG");
    analyze.setModelFilename("post_simple_model_all_the_probes.osim");
    analyze.setStatesFileName("muscleprescribewithemg_states.sto");
    analyze.updAnalysisSet().cloneAndAppend(MuscleAnalysis());
    analyze.updAnalysisSet().cloneAndAppend(ProbeReporter());
    analyze.updControllerSet().cloneAndAppend(PrescribedController("muscleprescribewithemg_controls.sto"));
    analyze.updAnalysisSet().cloneAndAppend(ForceReporter());
    analyze.updAnalysisSet().cloneAndAppend(BodyKinematics());
    analyze.setInitialTime(Time(1));
    analyze.setFinalTime(Time(end));
    analyze.print("testing_AnalyzeTool_setup.xml");
    % Run the analysis.
    analyze = AnalyzeTool("testing_AnalyzeTool_setup.xml");
    analyze.run();


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % tables
    table_activefiberforce = TimeSeriesTable("analyzemusclesEMG_MuscleAnalysis_ActiveFiberForce.sto");
    table_fibervelocity = TimeSeriesTable("analyzemusclesEMG_MuscleAnalysis_FiberVelocity.sto");
    table_metabolics = TimeSeriesTable('analyzemusclesEMG_ProbeReporter_probes.sto');
    table_lMT = TimeSeriesTable('analyzemusclesEMG_MuscleAnalysis_Length.sto');
    table_fiberlength = TimeSeriesTable('analyzemusclesEMG_MuscleAnalysis_FiberLength.sto');
    
    
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
    % TODO: figure out good save points/methods for:
    % whole body, each muscle, through time and averaged
    
    % get metabolics probe information - TODO: figure out what is useful from here 
    metabolics_all_os = table_metabolics.getDependentColumn('all_metabolics_TOTAL');
    metabolics_act_os = table_metabolics.getDependentColumn('all_activation_maintenance_rate_TOTAL');
    metabolics_short_os = table_metabolics.getDependentColumn('all_shortening_rate_TOTAL');
    metabolics_basal_os = table_metabolics.getDependentColumn('all_basal_rate_TOTAL');
    metabolics_mech_os = table_metabolics.getDependentColumn('all_mechanical_work_rate_TOTAL');
    % individual muscles
    metabolics_gas_os = table_metabolics.getDependentColumn('gastroc_metabolics_TOTAL');
    metabolics_sol_os = table_metabolics.getDependentColumn('soleus_metabolics_TOTAL');
    metabolics_bifemlh_os = table_metabolics.getDependentColumn('bifemlh_metabolics_TOTAL');
    metabolics_recfem_os = table_metabolics.getDependentColumn('recfem_metabolics_TOTAL');

    metabolics_all = metabolics_all_os.getAsMat;
    metabolics_act = metabolics_act_os.getAsMat;
    metabolics_short = metabolics_short_os.getAsMat;
    metabolics_basal = metabolics_basal_os.getAsMat;
    metabolics_mech = metabolics_mech_os.getAsMat;

    % fix the basal stuff
    metabolics_basal_old = metabolics_basal;
    basal_coef = 2;
    basal_exp = 1;
    for i=1:length(metabolics_basal)
        metabolics_all(i) = metabolics_all(i) - metabolics_basal(i);
        metabolics_basal(i) = basal_coef*(model_mass^basal_exp);
        % metabolics_all(i) = metabolics_all(i) + metabolics_basal(i)
    end

    % individual muscles
    metabolics_gas = metabolics_gas_os.getAsMat;
    metabolics_sol = metabolics_sol_os.getAsMat;
    metabolics_bifemlh = metabolics_bifemlh_os.getAsMat;
    metabolics_recfem = metabolics_recfem_os.getAsMat;

    metabolics_basal_avg = ((trapz(time, metabolics_basal)) / (time(end)-time(1))) / model_mass;
    metabolics_act_avg = 2*((trapz(time, metabolics_act)) / (time(end)-time(1))) / model_mass;
    metabolics_short_avg = 2*((trapz(time, metabolics_short)) / (time(end)-time(1))) / model_mass;
    metabolics_mech_avg = 2*((trapz(time, metabolics_mech)) / (time(end)-time(1))) / model_mass;
    metabolics_all_avg = 2*((trapz(time, metabolics_all)) / (time(end)-time(1))) / model_mass;
    metabolics_all_avg = metabolics_all_avg + metabolics_basal_avg;

    % individual muscles
    metabolics_gas_avg = ((trapz(time, metabolics_gas)) / (time(end)-time(1))); %/ model_mass;
    metabolics_sol_avg = ((trapz(time, metabolics_sol)) / (time(end)-time(1))); %/ model_mass;
    metabolics_bifemlh_avg = ((trapz(time, metabolics_bifemlh)) / (time(end)-time(1))); %/ model_mass;
    metabolics_recfem_avg = ((trapz(time, metabolics_recfem)) / (time(end)-time(1))); %/ model_mass;

    % temp_reg
    temp_emg = metabolics_all_avg
    % difference = temp_reg - temp_emg

    met_rows = {'trial'};
    met_table = table(metabolics_all_avg, metabolics_act_avg, metabolics_short_avg,...
                metabolics_basal_avg, metabolics_mech_avg,...
                metabolics_gas_avg, metabolics_sol_avg, metabolics_bifemlh_avg, metabolics_recfem_avg,... 
                model_mass, subjectname,condname,...
                experimentname,trialname,'RowNames', met_rows);
                        
    writetable(met_table, 'metabolicsTable_withemg.csv','WriteRowNames',true);

end