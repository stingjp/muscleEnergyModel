function torquePrescribeTrackGRFPrescribe()
    import org.opensim.modeling.*;
    
    % create and name an instance of the MocoTrack tool
    inverse = MocoInverse();
    inverse.setName("torque_stateprescribe_grfprescribe");
    
    % construct a ModelProcessor and add it to the tool.
    % modelProcessor = ModelProcessor("simple_model_all_the_probes_adjusted.osim");
    % modelProcessor = ModelProcessor("subject_redoarms.osim");
    modelProcessor = ModelProcessor("simple_model_all_the_probes.osim");
    % need to adjust some of the joints - weld them
    weldem = StdVectorString();
    weldem.add('subtalar_r');
    weldem.add('mtp_r');
    weldem.add('subtalar_l');
    weldem.add('mtp_l');
    weldem.add('radius_hand_r');
    weldem.add('radius_hand_l');
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
    inverse.setModel(modelProcessor);
    
    % tableProcessor = TableProcessor('results_IK_redoarms.mot');
    tempkintable = TableProcessor('torque_statetrack_grfprescribe_solution.sto').process();
    templabels_os = tempkintable.getColumnLabels();
    for i=0:templabels_os.size()-1
        temp = templabels_os.get(i);
        if ~startsWith(temp, '/jointset')
            tempkintable.removeColumn(temp);
        end
    end
    tableprocessorkin = TableProcessor(tempkintable);
    tableprocessorkin.append(TabOpUseAbsoluteStateNames());
    inverse.setKinematics(tableprocessorkin);
    
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
    inverse.set_initial_time(gait_start);
    inverse.set_final_time(gait_end);
    inverse.set_mesh_interval(0.03); %.01% may regret later

    % instead of calling solve, call initialize to get pre-configured
    % MocoStudy object, that can be further customized
    study = inverse.initialize();

    % get reference to the MocoControlGoal that is added to every MocoTrack problem
    problem = study.updProblem();
    effort = MocoControlGoal.safeDownCast(problem.updGoal('excitation_effort'));
%     effort.setWeight(0.1);
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
            effort.setWeightForControl(forcePath, 10000); % here
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
    solution.write('torque_stateprescribe_grfprescribe_solution_nokinematics.sto');
    solution.insertStatesTrajectory(tempkintable);
    solution.write('torque_stateprescribe_grfprescribe_solution.sto');
    % study.visualize(MocoTrajectory("torque_statetrack_grfprescribe_solution.sto"));
        
    report = osimMocoTrajectoryReport(model, 'torque_stateprescribe_grfprescribe_solution.sto');
    reportFilePath = report.generate();
    pdfFilePath = reportFilePath(1:end-2);
    pdfFilePath = strcat(pdfFilePath, 'pdf');
    ps2pdf('psfile',reportFilePath,'pdffile',pdfFilePath, ...
        'gscommand','C:\Program Files\gs\gs9.54.0\bin\gswin64.exe', ...
        'gsfontpath','C:\Program Files\gs\gs9.54.0\Resource\Font', ...
        'gslibpath','C:\Program Files\gs\gs9.54.0\lib');
    % open(pdfFilePath);
    % save('torque_statetrack_grfprescribe.mat');
    disp('end state torque inverse')
end