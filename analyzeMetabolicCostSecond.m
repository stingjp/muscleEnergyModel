function analyzeMetabolicCost_second(solution, tag)
    import org.opensim.modeling.*
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for the unconstrained solution
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Conduct an analysis using MuscleAnalysis and ProbeReporter.
    % solution = MocoTrajectory('muscle_stateprescribe_grfprescribe_solution.sto');
    Time = solution.getTimeMat();
    numColPoints = solution.getNumTimes();
    
    % get the subject name and mass
    load 'G:\Shared drives\Exotendon\muscleModel\muscleEnergyModel\subjectmass.mat';
    workdir = pwd;
    [~,trialname,~] = fileparts(pwd);
    cd ../
    [~,condname,~] = fileparts(pwd);
    cd ../
    [~,subjectname,~] = fileparts(pwd);
    experimentname = subjectname(1:4);
    cd(workdir);
    model_mass = subjectmass.(genvarname(subjectname)); % kg
    
    % Redo the states and controls files needed for inverse analysis
    STOFileAdapter.write(solution.exportToControlsTable(), strcat(tag, "_controls_100con.sto"));
    STOFileAdapter.write(solution.exportToStatesTable(), strcat(tag, "_states_100con.sto"));


    % full moco method
    analyze = AnalyzeTool();
    analyze.setName(strcat("analyzemuscles_100con",tag));
    analyze.setModelFilename("post_simple_model_all_the_probes_muscletrack.osim");
    

    % if strcmp(subjectname,'welk002') || strcmp(subjectname,'welk003')
    %     analyze.setStatesFileName("muscleprescribe_states.sto");
    %     analyze.updControllerSet().cloneAndAppend(PrescribedController("muscleprescribe_controls.sto"));
    % else
    %     analyze.setStatesFileName(strcat(tag, "_states_100con.sto"));
    %     analyze.updControllerSet().cloneAndAppend(PrescribedController(strcat(tag,"_controls_100con.sto")));
    % end
    analyze.setStatesFileName(strcat(tag, "_states_100con.sto"));
    analyze.updControllerSet().cloneAndAppend(PrescribedController(strcat(tag,"_controls_100con.sto")));
   


    analyze.updAnalysisSet().cloneAndAppend(MuscleAnalysis());
    analyze.updAnalysisSet().cloneAndAppend(ProbeReporter());
    analyze.updAnalysisSet().cloneAndAppend(ForceReporter());
    analyze.updAnalysisSet().cloneAndAppend(BodyKinematics());
    analyze.setInitialTime(Time(1));
    analyze.setFinalTime(Time(end));
    analyze.print(strcat(tag,"_AnalyzeTool_setup_100con.xml"));
    % Run the analysis.
    analyze = AnalyzeTool(strcat(tag,"_AnalyzeTool_setup_100con.xml"));
    analyze.run();

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % tables
    table_activefiberforce = TimeSeriesTable(strcat("analyzemuscles_100con",tag,"_MuscleAnalysis_ActiveFiberForce.sto"));
    table_fibervelocity = TimeSeriesTable(strcat("analyzemuscles_100con",tag,"_MuscleAnalysis_FiberVelocity.sto"));
    table_metabolics = TimeSeriesTable(strcat("analyzemuscles_100con",tag,'_ProbeReporter_probes.sto'));
    table_lMT = TimeSeriesTable(strcat("analyzemuscles_100con",tag,'_MuscleAnalysis_Length.sto'));
    table_fiberlength = TimeSeriesTable(strcat("analyzemuscles_100con",tag,'_MuscleAnalysis_FiberLength.sto'));
    
    % get time
    time_os = table_activefiberforce.getIndependentColumn();
    time = zeros(time_os.size(),1);
    for i=0:time_os.size()-1
        time(i+1) = time_os.get(i);
    end
    time_met_os = table_metabolics.getIndependentColumn();
    time_met = zeros(time_met_os.size(), 1);
    for i=0:time_met_os.size()-1
        time_met(i+1) = time_met_os.get(i);
    end

    % get all the muscles
    muscles = [];
    labels = table_activefiberforce.getColumnLabels();
    for i=0:labels.size()-1
        muscles = [muscles; {labels.get(i)}];
    end
    numMuscles = length(muscles);
    
    % get active fiber force and fiber velocity
    activefiberforce = [];
    fibervelocity = [];
    fiberlength = [];
    for i=1:length(muscles)
        temp_activefiberforce = table_activefiberforce.getDependentColumn(muscles{i}).getAsMat();
        activefiberforce = [activefiberforce, temp_activefiberforce];
        temp_fibervelocity = table_fibervelocity.getDependentColumn(muscles{i}).getAsMat();
        fibervelocity = [fibervelocity, temp_fibervelocity];
        temp_fiberlength = table_fiberlength.getDependentColumn(muscles{i}).getAsMat();
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
    basal_coef = 1.2;
    basal_exp = 1;
    for i=1:length(metabolics_basal)
        metabolics_all(i) = metabolics_all(i) - metabolics_basal(i);
        metabolics_basal(i) = basal_coef*(model_mass^basal_exp);
        % metabolics_all(i) = metabolics_all(i) + metabolics_basal(i)
    end

    %%% workspace %%%
    % TODO figure out how to get all the muscles averages
    
    table_musc_metabolics = table_metabolics;
    table_musc_metabolics.removeColumn('all_metabolics_TOTAL');
    table_musc_metabolics.removeColumn('all_activation_maintenance_rate_TOTAL');
    table_musc_metabolics.removeColumn('all_shortening_rate_TOTAL');
    table_musc_metabolics.removeColumn('all_basal_rate_TOTAL');
    table_musc_metabolics.removeColumn('all_mechanical_work_rate_TOTAL');
    table_musc_metabolics.removeColumn('soleus_metabolics_TOTAL');
    table_musc_metabolics.removeColumn('gastroc_metabolics_TOTAL');
    table_musc_metabolics.removeColumn('bifemlh_metabolics_TOTAL');
    table_musc_metabolics.removeColumn('recfem_metabolics_TOTAL');
    
    
    
    % now it is each probe type for each muscle - hella probes
    nummuscmet = table_musc_metabolics.getNumColumns();
    muscmetlabels = table_musc_metabolics.getColumnLabels();
    muscMetabolicsMat = [];
    muscMetabolicsLabels = {};
    
    muscMetTime = table_musc_metabolics.getIndependentColumn();
    for i=0:nummuscmet-1
        templabel = muscmetlabels.get(i);
        tempcolumn = table_musc_metabolics.getDependentColumn(templabel);
        muscMetabolicsMat = [muscMetabolicsMat, tempcolumn.getAsMat()];
        muscMetabolicsLabels{i+1} = char(templabel);
    end
    
    muscMetabolicsMat;
    avgMuscMetMat = [];
    
    % loop through to average each over the gait cycle
    for i = 1:nummuscmet
        tempinteg = ((trapz(time, muscMetabolicsMat(:,i))) / (time(end)-time(1))) / model_mass;
        avgMuscMetMat = [avgMuscMetMat, tempinteg];
    end
    
    % write them all to a file that I can pull later to get differences
    % get everything set up for the table printout
    % met_rows = {'trial'};

%     avgMuscMetMat2 = num2cell(avgMuscMetMat);
%     musc_table = cell2table(avgMuscMetMat2);
%     musc_table.Properties.VariableNames = muscMetabolicsLabels;
    musc_table = table((avgMuscMetMat)', (muscMetabolicsLabels)');        
    writetable(musc_table, 'muscleMetabolicsALL_100con.csv');% ,'WriteRowNames',true);
    
    
    % look through the GRF file?
    % get grf for residual comparisons
    table_grf = TimeSeriesTable(strcat('analyzemuscles_100con',tag,'_ForceReporter_forces.sto'));
    grf_r_Fx = table_grf.getDependentColumn('calcn_r_Right_GRF_Fx').getAsMat();
    grf_r_Fy = table_grf.getDependentColumn('calcn_r_Right_GRF_Fy').getAsMat();
    grf_r_Fz = table_grf.getDependentColumn('calcn_r_Right_GRF_Fz').getAsMat();
    grf_l_Fx = table_grf.getDependentColumn('calcn_l_Left_GRF_Fx').getAsMat();
    grf_l_Fy = table_grf.getDependentColumn('calcn_l_Left_GRF_Fy').getAsMat();
    grf_l_Fz = table_grf.getDependentColumn('calcn_l_Left_GRF_Fz').getAsMat();
    grf_r_Tx = table_grf.getDependentColumn('calcn_r_Right_GRF_Tx').getAsMat();
    grf_r_Ty = table_grf.getDependentColumn('calcn_r_Right_GRF_Ty').getAsMat();
    grf_r_Tz = table_grf.getDependentColumn('calcn_r_Right_GRF_Tz').getAsMat();
    grf_l_Tx = table_grf.getDependentColumn('calcn_l_Left_GRF_Tx').getAsMat();
    grf_l_Ty = table_grf.getDependentColumn('calcn_l_Left_GRF_Ty').getAsMat();
    grf_l_Tz = table_grf.getDependentColumn('calcn_l_Left_GRF_Tz').getAsMat();

    
    % grab a list of indices that are nonzero in the Fy direction, indicating stance. 
    % get the corresponding others for swing of that leg. 
    tempstanceix = find(grf_r_Fy);
    % loop through these and grab the metabolic cost values
    met_stance = [];
    met_swing = metabolics_all;
    popped = 0;


    % loop through the number of steps we have in stance phase
    for i = 1:length(tempstanceix)
        % get the index corresponding to the next stance phase timestep
        tempix = tempstanceix(i);
        % get the metabolic rate at that time step - again in stance
        tempmet = metabolics_all(tempix);
        % add the metabolic rate to the vector
        met_stance = [met_stance, tempmet];
        % pop out a value from the swing vector that is at the same index as the stance one we just specified
        met_swing(tempix-popped) = [];
        % add one to popped
        popped = popped + 1;
    end

    % going to need time vecs for the integration: actual time does not matter, only amount
    % get the end time
    time_max = time(end);
    % get initial time
    time_min = time(1);
    % time difference
    time_diff = time_max - time_min;
    % get the time step size
    time_step = time_diff / length(time);
    % create a time vector that is the same length as the stance and swing vectors,
    % multiply by step size to get actual time spacing normalized
    stance_time = time_step*linspace(1, length(met_stance),length(met_stance));
    swing_time = time_step*linspace(1, length(met_swing), length(met_swing));


    % now actually get average values for each stance and swing for one leg
    metabolics_stance_avg = ((trapz(stance_time, met_stance)) / (stance_time(end)-stance_time(1))) / model_mass;
    metabolics_swing_avg = ((trapz(swing_time, met_swing)) / (swing_time(end)-swing_time(1))) / model_mass;
    % TODO do this analysis for stance and swing for all the types of metabolics, not just full
    


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


    temp_reg = metabolics_all_avg

    % get everything set up for the table printout
    met_rows = {'trial'};
    met_table = table(metabolics_all_avg, metabolics_act_avg, metabolics_short_avg,...
                metabolics_basal_avg, metabolics_mech_avg,...
                metabolics_gas_avg, metabolics_sol_avg, metabolics_bifemlh_avg, metabolics_recfem_avg,... 
                metabolics_swing_avg, metabolics_stance_avg,...
                model_mass, {subjectname},{condname},...
                {experimentname},{trialname},'RowNames', met_rows);
    met_table.Properties.VariableNames{'Var13'} = 'subjectname';
    met_table.Properties.VariableNames{'Var14'} = 'condname';
    met_table.Properties.VariableNames{'Var15'} = 'experimentname';
    met_table.Properties.VariableNames{'Var16'} = 'trialname';
    
    writetable(met_table, 'metabolicsTable_100con.csv','WriteRowNames',true);

end
