function torqueStateTrackGRFPrescribeStrict()
    import org.opensim.modeling.*;
    
    % create and name an instance of the MocoTrack tool
    track = MocoTrack();
    track.setName("torque_statetrack_grfprescribe_strict");
    
    % construct a ModelProcessor and add it to the tool.
    modelProcessor = ModelProcessor("simple_model_all_the_probes_adjusted.osim");
    
%     modelProcessor = ModelProcessor("simple_model_all_the_probes.osim");
    % now to do stuff with the model
    % modelProcessor = ModelProcessor(model);
    % need to adjust some of the joints - weld them
    weldem = StdVectorString();
    % weldem.add('subtalar_r');
    weldem.add('mtp_r');
    % weldem.add('subtalar_l');
    weldem.add('mtp_l');
    % weldem.add('radius_hand_r');
    % weldem.add('radius_hand_l');
    
    modelProcessor.append(ModOpReplaceJointsWithWelds(weldem));
    % model = modelProcessor.process();
    % add ground reaction external loads in lieu of ground-contact model. 
    modelProcessor.append(ModOpAddExternalLoads("grf_walk.xml"));
    % remove all muscles for torque driven analysis
    modelProcessor.append(ModOpRemoveMuscles());
    % add CoordinateActuators to the model DOF. 
    % ignores pelvis coordinates with already have. 
    modelProcessor.append(ModOpAddReserves(250));
    torquemodel = modelProcessor.process();
    torquemodel.print('torquemodel.osim');
    track.setModel(modelProcessor);
    
    % construct a TableProcessor of the coordinate data and pass it to the tracking tool. 
    % 1
    % track.setStatesReference(TableProcessor('torque_markertrack_grfprescribe_solution.sto'));
    % 2 
    % tableProcessor = TableProcessor('coordinates_updated.mot');
    % 3
    tableProcessor = TableProcessor(tabletrimming('coordinates_updated.mot')); %***
    tableProcessor.append(TabOpLowPassFilter(6));
    % 4 
    % tempTable = TimeSeriesTable('./ResultsRRA_2/subject01_walk1_RRA_Kinematics_q.sto');
    % tableProcessor = TableProcessor(tempTable);
    
    tableProcessor.append(TabOpUseAbsoluteStateNames());
    
    track.setStatesReference(tableProcessor);
    track.set_states_global_tracking_weight(1e2);
    % avoid exceptions if markers in file are no longer in the model (arms removed)
    track.set_allow_unused_references(true);
    % since there is only coordinate position data in the states references, 
    % this fills in the missing coordinate speed data using 
    % the derivative of splined position data
    track.set_track_reference_position_derivatives(true);
    
    % set specific weights for the individual weight set
    coordinateweights = MocoWeightSet();
    % coordinateweights.cloneAndAppend(MocoWeight("pelvis_tx", 0));
    % coordinateweights.cloneAndAppend(MocoWeight("pelvis_ty", 0));
    % coordinateweights.cloneAndAppend(MocoWeight("pelvis_tz", 0));
    % coordinateweights.cloneAndAppend(MocoWeight("pelvis_list", 0));
    % coordinateweights.cloneAndAppend(MocoWeight("pelvis_rotation", 0));
    % coordinateweights.cloneAndAppend(MocoWeight("pelvis_tilt", 0));
    coordinateweights.cloneAndAppend(MocoWeight("hip_rotation_r", 1e-14));
    coordinateweights.cloneAndAppend(MocoWeight("hip_rotation_l", 1e-14));
    coordinateweights.cloneAndAppend(MocoWeight("subtalar_angle_r", 1e-14));
    coordinateweights.cloneAndAppend(MocoWeight("subtalar_angle_l", 1e-14));
    % coordinateweights.cloneAndAppend(MocoWeight("ankle_angle_r", 0));
    % coordinateweights.cloneAndAppend(MocoWeight("ankle_angle_l", 0));
    
    track.set_states_weight_set(coordinateweights);
    
    

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

    % set the times and mesh interval, mesh points are computed internally. 
    track.set_initial_time(gait_start);
    track.set_final_time(gait_end);
    track.set_mesh_interval(0.01); %.01% may regret later

    % instead of calling solve, call initialize to get pre-configured
    % MocoStudy object, that can be further customized
    study = track.initialize();

    % get reference to the MocoControlGoal that is added to every MocoTrack problem
    problem = study.updProblem();
    effort = MocoControlGoal.safeDownCast(problem.updGoal('control_effort'));
    effort.setWeight(1e-16);
    % initactivationgoal = MocoInitialActivationGoal('init_activation');
    % initactivationgoal.setWeight(10);
    % problem.addGoal(initactivationgoal);



    % put large weight on the pelvis CoordinateActuators, which act as the 
    % residual, or 'hand-of-god' forces which we would like to keep small
    
    model = modelProcessor.process();
    model.initSystem();
    forceSet = model.getForceSet();
    for i=0:forceSet.getSize()-1
        forcePath = forceSet.get(i).getAbsolutePathString();
        if contains(string(forcePath), 'pelvis')
            % effort.setWeightForControl(forcePath, 1e-4); % here
            % if contains(string(forcePath), 'pelvis_ty')
            %     effort.setWeightForControl(forcePath, 1e8);
            % end
        end
        % if contains(string(forcePath), 'hip_rotation')
        %    effort.setWeightForControl(forcePath, 1e4);
        % end
    end
    
    % solver changes
%     solver = MocoCasADiSolver.safeDownCast(study.updSolver());
%     solver.resetProblem(problem);
%     solver.set_optim_convergence_tolerance(1e-5); % 1e-2
%     solver.set_optim_constraint_tolerance(1e-5); % 1e-2
    
    
    % solve and visualize
    solution = study.solve();
    % study.visualize(solution);
    % generate a report and save
    solution.write('torque_statetrack_grfprescribe_strict_solution.sto');
    % study.visualize(MocoTrajectory("torque_statetrack_grfprescribe_solution.sto"));
        
    report = osimMocoTrajectoryReport(model, 'torque_statetrack_grfprescribe_strict_solution.sto');
    reportFilePath = report.generate();
    pdfFilePath = reportFilePath(1:end-2);
    pdfFilePath = strcat(pdfFilePath, 'pdf');
    ps2pdf('psfile',reportFilePath,'pdffile',pdfFilePath, ...
        'gscommand','C:\Program Files\gs\gs9.54.0\bin\gswin64.exe', ...
        'gsfontpath','C:\Program Files\gs\gs9.54.0\Resource\Font', ...
        'gslibpath','C:\Program Files\gs\gs9.54.0\lib');
    % open(pdfFilePath);
    % save('torque_statetrack_grfprescribe.mat');
    disp('end state torque track')
end