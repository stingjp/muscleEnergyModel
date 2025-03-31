function torqueMarkerTrackGRFTrack()
    
    import org.opensim.modeling.*;

    % create and name an instance of the MocoTrack tool
    track = MocoTrack();
    track.setName("torque_markertrack_grftrack");

    % construct a ModelProcessor and add it to the tool.
    modelProcessor = ModelProcessor("simple_model_all_the_probes_contact.osim");
    weldem = StdVectorString();
    weldem.add('subtalar_r');
    % weldem.add('mtp_r');
    weldem.add('subtalar_l');
    % weldem.add('mtp_l');
    weldem.add('radius_hand_r');
    weldem.add('radius_hand_l');
    modelProcessor.append(ModOpReplaceJointsWithWelds(weldem));
    % add ground reaction external loads in lieu of ground-contact model. 
    % modelProcessor.append(ModOpAddExternalLoads("grf_walk.xml"));
    % remove all muscles for torque driven analysis
    modelProcessor.append(ModOpRemoveMuscles());
    % add CoordinateActuators to the model DOF. 
    % ignores pelvis coordinates with already have. 
    modelProcessor.append(ModOpAddReserves(250));
    track.setModel(modelProcessor);

    % set the MocoTrack markers reference directly from trc file. 
    % data is filtered at 6Hz and if in millimeters, converted to m. 
    track.setMarkersReferenceFromTRC("motion_capture.trc");
    % avoid exceptions if markers in file are no longer in the model (arms removed)
    track.set_allow_unused_references(true);
    % increase global marker tracking weight, which is associated with internal 
    % MocoMarkerTrackingGoal term
    track.set_markers_global_tracking_weight(10);
    % increase tracking weights for individual markers in the data set
    markerWeights = MocoWeightSet();
    markerWeights.cloneAndAppend(MocoWeight("R.ASIS", 20));
    markerWeights.cloneAndAppend(MocoWeight("L.ASIS", 20));
    markerWeights.cloneAndAppend(MocoWeight("R.PSIS", 20));
    markerWeights.cloneAndAppend(MocoWeight("L.PSIS", 20));
    markerWeights.cloneAndAppend(MocoWeight("R.Knee", 10));
    markerWeights.cloneAndAppend(MocoWeight("R.Ankle", 10));
    markerWeights.cloneAndAppend(MocoWeight("R.Heel", 10));
    markerWeights.cloneAndAppend(MocoWeight("R.MT5", 5));
    markerWeights.cloneAndAppend(MocoWeight("R.Toe", 2));
    % markerWeights.cloneAndAppend(MocoWeight("L.Knee", 10));
    % markerWeights.cloneAndAppend(MocoWeight("L.Ankle", 10));
    % markerWeights.cloneAndAppend(MocoWeight("L.Heel", 10));
    % markerWeights.cloneAndAppend(MocoWeight("L.MT5", 5));
    % markerWeights.cloneAndAppend(MocoWeight("L.Toe", 2));
    track.set_markers_weight_set(markerWeights);


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
    track.set_mesh_interval(0.01);

    % instead of calling solve, call initialize to get pre-configured
    % MocoStudy object, that can be further customized
    study = track.initialize();

    % get reference to the MocoControlGoal that is added to every MocoTrack problem
    problem = study.updProblem();
    effort = MocoControlGoal.safeDownCast(problem.updGoal('control_effort'));
%     effort.setWeight(0.1);
    initactivationgoal = MocoInitialActivationGoal('init_activation');
    initactivationgoal.setWeight(10);
    problem.addGoal(initactivationgoal);

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
    
%     keyboard
    contactTracking = MocoContactTrackingGoal('contact', 1e0)    
    contactTracking.setExternalLoadsFile('grf_walk.xml');
    
    forceNamesRightFoot = StdVectorString();
    forceNamesRightFoot.add('/contactHeel_r');
    % # forceNamesRightFoot.add('/forceset/contactLateralRearfoot_r');
    forceNamesRightFoot.add('/contactLateralMidfoot_r');
    % # forceNamesRightFoot.add('/contactLateralToe_r');
    forceNamesRightFoot.add('/contactMedialToe_r');
    forceNamesRightFoot.add('/contactMedialMidfoot_r');
    % # contactTracking.addContactGroup(forceNamesRightFoot, 'Right_GRF');
    contactTrackingSplitRight = MocoContactTrackingGoalGroup(forceNamesRightFoot, 'Right_GRF');
    contactTrackingSplitRight.append_alternative_frame_paths('/bodyset/toes_r')
    contactTracking.addContactGroup(contactTrackingSplitRight);

    forceNamesLeftFoot = StdVectorString();
    forceNamesLeftFoot.add('/contactHeel_l');
    % # forceNamesLeftFoot.add('/forceset/contactLateralRearfoot_l');
    forceNamesLeftFoot.add('/contactLateralMidfoot_l');
    % # forceNamesLeftFoot.add('/contactLateralToe_l');
    forceNamesLeftFoot.add('/contactMedialToe_l');
    forceNamesLeftFoot.add('/contactMedialMidfoot_l');
    % # contactTracking.addContactGroup(forceNamesLeftFoot, 'Left_GRF');
    contactTrackingSplitLeft = MocoContactTrackingGoalGroup(forceNamesLeftFoot, 'Left_GRF');
    contactTrackingSplitLeft.append_alternative_frame_paths('/bodyset/toes_l')
    contactTracking.addContactGroup(contactTrackingSplitLeft);
    
    contactTracking.setProjection('plane');
    contactTracking.setProjectionVector(Vec3(0, 0, 1));

    contactTracking.setDivideByDisplacement(false)
    contactTracking.setDivideByMass(true)
    problem.addGoal(contactTracking);


    % solve - the bool indicates to visualize the solution
    % solution = track.solve(true); % to visualize
    solution = study.solve();
    solution.write('torque_markertrack_grftrack_solution.sto');

    % generate a pdf report containing plots of the variables in the solution. 
    % for details see osimMocoTrajectoryReport.m in moco resource/code/matlab/utilities 
    model = modelProcessor.process();
    report = osimMocoTrajectoryReport(model, 'torque_markertrack_grftrack_solution.sto');
    reportFilePath = report.generate();
    pdfFilePath = reportFilePath(1:end-2);
    pdfFilePath = strcat(pdfFilePath, 'pdf');
    ps2pdf('psfile',reportFilePath,'pdffile',pdfFilePath, ...
        'gscommand','C:\Program Files\gs\gs9.54.0\bin\gswin64.exe', ...
        'gsfontpath','C:\Program Files\gs\gs9.54.0\Resource\Font', ...
        'gslibpath','C:\Program Files\gs\gs9.54.0\lib');
    % open(pdfFilePath);
    % save('torque_markertrack_grfprescribe.mat');
    disp('end marker track');
end