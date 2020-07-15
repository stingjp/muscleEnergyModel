function somthing()

    % run a marker tracking to get joint angles etc. 
    % torqueMarkerTrackGRFPrescribe();

    % run a inverse problem with the motion created above, and solve for muscles
    % torqueStateTrackGRFPrescribe();

    % run an inverse problem and with prescribed states and grf for muscle activations
%     muscleStatePrescribeGRFPrescribe();
    
    % run an inverse problem with prescribed states and grf, and
    % tracking experimental emg
    muscleStatePrescribeGRFPrescribeWithEMG();
    
end



function torqueMarkerTrackGRFPrescribe()
    
    import org.opensim.modeling.*;

    % create and name an instance of the MocoTrack tool
    track = MocoTrack();
    track.setName("torque_markertrack_grfprescribe");

    % construct a ModelProcessor and add it to the tool.
    modelProcessor = ModelProcessor("subject_walk_armless.osim");
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



function torqueStateTrackGRFPrescribe()

    import org.opensim.modeling.*;

    % create and name an instance of the MocoTrack tool
    track = MocoTrack();
    track.setName("torque_statetrack_grfprescribe");

    % construct a ModelProcessor and add it to the tool.
    modelProcessor = ModelProcessor("subject_walk_armless.osim");
    % add ground reaction external loads in lieu of ground-contact model. 
    modelProcessor.append(ModOpAddExternalLoads("grf_walk.xml"));
    % remove all muscles for torque driven analysis
    modelProcessor.append(ModOpRemoveMuscles());
    % add CoordinateActuators to the model DOF. 
    % ignores pelvis coordinates with already have. 
    modelProcessor.append(ModOpAddReserves(250));
    track.setModel(modelProcessor);

    % construct a TableProcessor of the coordinate data and pass it to the tracking tool. 
    track.setStatesReference(TableProcessor('coordinates.sto'));
    track.set_states_global_tracking_weight(10);
    % avoid exceptions if markers in file are no longer in the model (arms removed)
    track.set_allow_unused_references(true);
    % since there is only coordinate position data in the states references, 
    % this fills in the missing coordinate speed data using 
    % the derivative of splined position data
    track.set_track_reference_position_derivatives(true);

    % set the times and mesh interval, mesh points are computed internally. 
    track.set_initial_time(0.81);
    track.set_final_time(1.65);
    track.set_mesh_interval(0.05);

    % instead of calling solve, call initialize to get pre-configured
    % MocoStudy object, that can be further customized
    study = track.initialize();

    % get reference to the MocoControlGoal that is added to every MocoTrack problem
    problem = study.updProblem();
    effort = MocoControlGoal.safeDownCast(problem.updGoal('control_effort'));

    % put large weight on the pelvis CoordinateActuators, which act as the 
    % residual, or 'hand-of-god' forces which we would like to keep small
    model = modelProcessor.process()
    model.initSystem();
    forceSet = model.getForceSet();
    for i=0:forceSet.getSize()-1
        forcePath = forceSet.get(i).getAbsolutePathString();
        if contains(string(forcePath), 'pelvis')
            effort.setWeightForControl(forcePath, 10);
        end
    end

    % solve and visualize
    solution = study.solve();
    % study.visualize(solution);

    % generate a report and save
    solution.write('torque_statetrack_grfprescribe_solution.sto')
    report = osimMocoTrajectoryReport(model, 'torque_statetrack_grfprescribe_solution.sto');
    reportFilePath = report.generate();
    pdfFilePath = reportFilePath(1:end-2);
    pdfFilePath = strcat(pdfFilePath, 'pdf');
    ps2pdf('psfile',reportFilePath,'pdffile',pdfFilePath, ...
        'gscommand','C:\Program Files\gs\gs9.52\bin\gswin64.exe', ...
        'gsfontpath','C:\Program Files\gs\gs9.52\Resource\Font', ...
        'gslibpath','C:\Program Files\gs\gs9.52\lib');
    % open(pdfFilePath);
    % save('torque_statetrack_grfprescribe.mat');
end



function muscleStatePrescribeGRFPrescribe()
    import org.opensim.modeling.*;

    % construct MocoInverse tool
    inverse = MocoInverse();

    % construct ModelProcessor and sit it on the tool. 
    % replace default muscles with degrootefregly 2016 muscles, and adjust params
    modelProcessor = ModelProcessor('subject_walk_armless.osim');
    modelProcessor.append(ModOpAddExternalLoads('grf_walk.xml'));
    modelProcessor.append(ModOpIgnoreTendonCompliance());
    modelProcessor.append(ModOpReplaceMusclesWithDeGrooteFregly2016());
    % only valid for degroote
    modelProcessor.append(ModOpIgnorePassiveFiberForcesDGF());
    % only valid for degroote
    modelProcessor.append(ModOpScaleActiveFiberForceCurveWidthDGF(1.5));
    modelProcessor.append(ModOpAddReserves(1.0));
    inverse.setModel(modelProcessor);

    % construct TableProcessor of the coordinate data and pass it to the inverse tool
    % if no operators, it returns the base table
    inverse.setKinematics(TableProcessor('torque_statetrack_grfprescribe_solution.sto'));
    
    % set time and intervals
    inverse.set_initial_time(0.81);
    inverse.set_final_time(1.65);
    inverse.set_mesh_interval(0.02);
    % By default, Moco gives an error if the kinematics contains extra columns.
    % Here, we tell Moco to allow (and ignore) those extra columns.
    inverse.set_kinematics_allow_extra_columns(true);

    % Solve the problem and write the solution to a Storage file.
    % solution = inverse.solve(true); % to visualize
    solution = inverse.solve();
    solution.getMocoSolution().write('muscle_stateprescribe_grfprescribe_solution.sto');

    % Generate a report with plots for the solution trajectory.
    model = modelProcessor.process();
    report = osimMocoTrajectoryReport(model, ...
            'muscle_stateprescribe_grfprescribe_solution.sto', 'bilateral', true);
    % The report is saved to the working directory.
    reportFilePath = report.generate();
    pdfFilePath = reportFilePath(1:end-2);
    pdfFilePath = strcat(pdfFilePath, 'pdf');
    ps2pdf('psfile',reportFilePath,'pdffile',pdfFilePath, ...
        'gscommand','C:\Program Files\gs\gs9.52\bin\gswin64.exe', ...
        'gsfontpath','C:\Program Files\gs\gs9.52\Resource\Font', ...
        'gslibpath','C:\Program Files\gs\gs9.52\lib');
    open(pdfFilePath);
end



function muscleStatePrescribeGRFPrescribeWithEMG()

    import org.opensim.modeling.*;

    % setting up model and the kinematics
    inverse = MocoInverse();
    modelProcessor = ModelProcessor('subject_walk_armless.osim');
    modelProcessor.append(ModOpAddExternalLoads('grf_walk.xml'));
    modelProcessor.append(ModOpIgnoreTendonCompliance());
    modelProcessor.append(ModOpReplaceMusclesWithDeGrooteFregly2016());
    modelProcessor.append(ModOpIgnorePassiveFiberForcesDGF());
    modelProcessor.append(ModOpScaleActiveFiberForceCurveWidthDGF(1.5));
    modelProcessor.append(ModOpAddReserves(1.0));
    inverse.setModel(modelProcessor);

    inverse.setKinematics(TableProcessor('torque_statetrack_grfprescribe_solution.sto'));
    inverse.set_initial_time(0.81);
    inverse.set_final_time(1.65);
    inverse.set_mesh_interval(0.02);
    inverse.set_kinematics_allow_extra_columns(true);

    study = inverse.initialize();
    problem = study.updProblem();

    % add EMG tracking
    emgTracking = MocoControlTrackingGoal('emg_tracking');
    emgTracking.setWeight(50.0);
    % each column in electromyography.sto is normalized so the maximum value is 1 for each
    controlsRef = TimeSeriesTable('electromyography.sto');

    % Scalde down the tracked muscle activity based on the peak levels from 
    % 'Gait analysis: normal and pathological function' by Perry and Burnfield, 2010 
    soleus = controlsRef.updDependentColumn('soleus');
    gasmed = controlsRef.updDependentColumn('gastrocnemius');
    tibant = controlsRef.updDependentColumn('tibialis_anterior');
    medham = controlsRef.updDependentColumn('medial_hamstrings');
    bicfem = controlsRef.updDependentColumn('biceps_femoris');
    vaslat = controlsRef.updDependentColumn('vastus_lateralis');
    vasmed = controlsRef.updDependentColumn('vastus_medius');
    recfem = controlsRef.updDependentColumn('rectus_femoris')
    % glumax = controlsRef.updDependentColumn('gluteus_maximus');
    % glumed = controlsRef.updDependentColumn('gluteus_medius');

    for t = 0:controlsRef.getNumRows() - 1
        soleus.set(t, 0.77*soleus.get(t));
        gasmed.set(t, 0.87*gasmed.get(t));
        tibant.set(t, 0.37*tibant.get(t));
        medham.set(t, 0.40*medham.get(t));
        bicfem.set(t, 0.25*bicfem.get(t));
        vaslat.set(t, 0.30*vaslat.get(t));
        vasmed.set(t, 0.39*vasmed.get(t));
        recfem.set(t, 0.20*recfem.get(t));
        % glumax.set(t, 0.25*glumax.get(t));
        % glumed.set(t, )
    end

    emgTracking.setReference(TableProcessor(controlsRef));
    % associate actuators in the model with columns in electromyography.sto
    emgTracking.setReferenceLabel('/forceset/soleus_r','soleus');
    emgTracking.setReferenceLabel('/forceset/gasmed_r','gastrocnemius');
    emgTracking.setReferenceLabel('/forceset/gaslat_r','gastrocnemius');
    emgTracking.setReferenceLabel('/forceset/tibant_r','tibialis_anterior');
    emgTracking.setReferenceLabel('/forceset/semimem_r','medial_hamstrings');
    emgTracking.setReferenceLabel('/forceset/bflh_r','biceps_femoris');
    emgTracking.setReferenceLabel('/forceset/vaslat_r','vastus_lateralis');
    emgTracking.setReferenceLabel('/forceset/vasmed_r','vastus_medius');
    emgTracking.setReferenceLabel('/forceset/recfem_r','rectus_femoris');
    % emgTracking.setReferenceLabel('/forceset/glmax1_r','gluteus_maximus');
    % emgTracking.setReferenceLabel('/forceset/glmax2_r','gluteus_maximus');
    % emgTracking.setReferenceLabel('/forceset/glmax2_r','gluteus_maximus');
    problem.addGoal(emgTracking)

    % can we add in the rest of the actuators from the file to reference, or are they ignored

    solution = study.solve();
    solution.write('muscle_stateprescribe_grfprescribe_withemg_solution.sto')

    % keyboard



    % write the reference data in a way that's easy to compare to the solution.
    % controlsRef.removeColumn('medial_hamstrings');
    % controlsRef.removeColumn('biceps_femoris');
    % controlsRef.removeColumn('vastus_lateralis');
    % controlsRef.removeColumn('vastus_medius');
    % controlsRef.removeColumn('rectus_femoris');
    controlsRef.removeColumn('gluteus_maximus');
    controlsRef.removeColumn('gluteus_medius');
    columnLabels = StdVectorString();
    columnLabels.add('/forceset/soleus_r');
    columnLabels.add('/forceset/gasmed_r');
    columnLabels.add('/forceset/tibant_r');
    columnLabels.add('/forceset/semimem_r');
    columnLabels.add('/forceset/bflh_r');
    columnLabels.add('/forceset/vaslat_r');
    columnLabels.add('/forceset/vaslmed_r');
    columnLabels.add('/forceset/recfem_r');
    controlsRef.setColumnLabels(columnLabels);
    controlsRef.appendColumn('/forceset/gaslat_r', gasmed);
    STOFileAdapter.write(controlsRef, 'controls_reference.sto');

    % generate a report comparing MocoInverse solutions without and with EMG tracking
    model = modelProcessor.process();
    report = osimMocoTrajectoryReport(model, ...
             'muscle_stateprescribe_grfprescribe_solution.sto', ...
             'outputFilepath','muscle_stateprescribe_grfprescribe_withemg_solution_report.pdf', ...
             'bilateral', true, ...
             'refFiles', {'muscle_stateprescribe_grfprescribe_withemg_solution.sto', ...
                          'controls_reference.sto'});
    % the report is saved to the working directory
    reportFilePath = report.generate();
    pdfFilePath = reportFilePath(1:end-2);
    pdfFilePath = strcat(pdfFilePath, 'pdf');
    ps2pdf('psfile',reportFilePath,'pdffile',pdfFilePath, ...
        'gscommand','C:\Program Files\gs\gs9.52\bin\gswin64.exe', ...
        'gsfontpath','C:\Program Files\gs\gs9.52\Resource\Font', ...
        'gslibpath','C:\Program Files\gs\gs9.52\lib');
    % open(pdfFilePath);

    keyboard

    
    
    
    
    % need to then figure out how to calculate and save metabolics stuff

    % function that sets up all the parameters for the metabolics probe
    % calcWholeBodyMetabolicRate()
    % function that sets up all the params for each muscle to calculate probe
    % calcIndividualMetabolicRate()

    % this is a matlab function that actually calculates the umberger cost
    % calcUmbergetProbe()
    

    %{
    model = org.opensim.modeling.Model(model_path);
    
    mat.Time = Time;                % Time = res.time;
    mat.DatStore = DatStore;        % DatStore = lots
    mat.OptInfo = OptInfo;          % OptInfo = output (from gpops)
    mat.MuscleNames = MuscleNames;  % MuscleNames = DatStore.MuscleNames;

    
    wholeBodyEnergyCost = calcWholeBodyMetabolicRate(model, mat);
    muscleEnergyCost = calcIndividualMetabolicRate(model, mat);
    %}
    
end




