function torqueMarkerTrackGRFPrescribe()
    
    import org.opensim.modeling.*;

    % create and name an instance of the MocoTrack tool
    track = MocoTrack();
    track.setName("torque_markertrack_grfprescribe");

    % construct a ModelProcessor and add it to the tool.
    modelProcessor = ModelProcessor("simple_model_all_the_probes.osim");
    % add ground reaction external loads in lieu of ground-contact model. 
    modelProcessor.append(ModOpAddExternalLoads("grf_walk.xml"));
    % remove all muscles for torque driven analysis
    modelProcessor.append(ModOpRemoveMuscles());
    % add CoordinateActuators to the model DOF. 
    % ignores pelvis coordinates with already have. 
    modelProcessor.append(ModOpAddReserves(250));
    track.setModel(modelProcessor);

    % set the MocoTrack markers reference directly from trc file. 
    % data is filtered at 6Hz and if in millimeters, converted to m. 
    track.setMarkersReferenceFromTRC("marker_trajectories.trc");
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
    markerWeights.cloneAndAppend(MocoWeight("L.Knee", 10));
    markerWeights.cloneAndAppend(MocoWeight("L.Ankle", 10));
    markerWeights.cloneAndAppend(MocoWeight("L.Heel", 10));
    markerWeights.cloneAndAppend(MocoWeight("L.MT5", 5));
    markerWeights.cloneAndAppend(MocoWeight("L.Toe", 2));
    track.set_markers_weight_set(markerWeights);

    % set the times and mesh interval, mesh points are computed internally. 
    track.set_initial_time(0.81);
    track.set_final_time(1.65);
    track.set_mesh_interval(0.05);

    % solve - the bool indicates to visualize the solution
    % solution = track.solve(true); % to visualize
    solution = track.solve();
    solution.write('torque_markertrack_grfprescribe_solution.sto');

    % generate a pdf report containing plots of the variables in the solution. 
    % for details see osimMocoTrajectoryReport.m in moco resource/code/matlab/utilities 
    model = modelProcessor.process();
    report = osimMocoTrajectoryReport(model, 'torque_markertrack_grfprescribe_solution.sto');
    reportFilePath = report.generate();
    pdfFilePath = reportFilePath(1:end-2);
    pdfFilePath = strcat(pdfFilePath, 'pdf');
    ps2pdf('psfile',reportFilePath,'pdffile',pdfFilePath, ...
        'gscommand','C:\Program Files\gs\gs9.52\bin\gswin64.exe', ...
        'gsfontpath','C:\Program Files\gs\gs9.52\Resource\Font', ...
        'gslibpath','C:\Program Files\gs\gs9.52\lib');
    % open(pdfFilePath);
    % save('torque_markertrack_grfprescribe.mat');
end