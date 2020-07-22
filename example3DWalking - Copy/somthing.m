function somthing()

    % modify a model for the metabolics study
    metabolicsModelSetup('subject_walk_armless_noprobe.osim');

    % run a marker tracking to get joint angles etc. 
    torqueMarkerTrackGRFPrescribe();

    % run a inverse problem with the motion created above, and solve for muscles
    torqueStateTrackGRFPrescribe();

    % run an inverse problem and with prescribed states and grf for muscle activations
    muscleStatePrescribeGRFPrescribe();
    
    % run an inverse problem with prescribed states and grf, and
    % tracking experimental emg
    muscleStatePrescribeGRFPrescribeWithEMG();

    % run an analysis to calculate the metabolic cost
    analyzeMetabolicCost();
    
end

function metabolicsModelSetup(modelFilename)
    import org.opensim.modeling.*
        
    muscleKeys = {'addbrev_r', ...
                    'addlong_r', ...
                    'addmagDist_r', ...
                    'addmagIsch_r', ...
                    'addmagMid_r', ...
                    'addmagProx_r', ...
                    'bflh_r', ...
                    'bfsh_r', ...
                    'edl_r', ...
                    'ehl_r', ...
                    'fdl_r', ...
                    'fhl_r', ...
                    'gaslat_r', ...
                    'gasmed_r', ...
                    'glmax1_r', ...
                    'glmax2_r', ...
                    'glmax3_r', ...
                    'glmed1_r', ...
                    'glmed2_r', ...
                    'glmed3_r', ...
                    'glmin1_r', ...
                    'glmin2_r', ...
                    'glmin3_r', ...
                    'grac_r', ...
                    'iliacus_r', ...
                    'perbrev_r', ...
                    'perlong_r', ...
                    'piri_r', ...
                    'psoas_r', ...
                    'recfem_r', ...
                    'sart_r', ...
                    'semimem_r', ...
                    'semiten_r', ...
                    'soleus_r', ...
                    'tfl_r', ...
                    'tibant_r', ...
                    'tibpost_r', ...
                    'vasint_r', ...
                    'vaslat_r', ...
                    'vasmed_r', ...
                    'addbrev_l', ...
                    'addlong_l', ...
                    'addmagDist_l', ...
                    'addmagIsch_l', ...
                    'addmagMid_l', ...
                    'addmagProx_l', ...
                    'bflh_l', ...
                    'bfsh_l', ...
                    'edl_l', ...
                    'ehl_l', ...
                    'fdl_l', ...
                    'fhl_l', ...
                    'gaslat_l', ...
                    'gasmed_l', ...
                    'glmax1_l', ...
                    'glmax2_l', ...
                    'glmax3_l', ...
                    'glmed1_l', ...
                    'glmed2_l', ...
                    'glmed3_l', ...
                    'glmin1_l', ...
                    'glmin2_l', ...
                    'glmin3_l', ...
                    'grac_l', ...
                    'iliacus_l', ...
                    'perbrev_l', ...
                    'perlong_l', ...
                    'piri_l', ...
                    'psoas_l', ...
                    'recfem_l', ...
                    'sart_l', ...
                    'semimem_l', ...
                    'semiten_l', ...
                    'soleus_l', ...
                    'tfl_l', ...
                    'tibant_l', ...
                    'tibpost_l', ...
                    'vasint_l', ...
                    'vaslat_l', ...
                    'vasmed_l'};
    ratioVals = [0.5, ...
        0.5, ...
        0.55200000000000005, ...
        0.55200000000000005, ...
        0.55200000000000005, ...
        0.55200000000000005, ...
        0.54249999999999998, ...
        0.52900000000000003, ...
        0.75, ...
        0.75, ...
        0.59999999999999998, ...
        0.59999999999999998, ...
        0.50700000000000001, ...
        0.56599999999999995, ...
        0.55000000000000004, ...
        0.55000000000000004, ...
        0.55000000000000004, ...
        0.55000000000000004, ...
        0.55000000000000004, ...
        0.55000000000000004, ...
        0.55000000000000004, ...
        0.55000000000000004, ...
        0.55000000000000004, ...
        0.5, ...
        0.5, ...
        0.59999999999999998, ...
        0.59999999999999998, ...
        0.5, ...
        0.5, ...
        0.38650000000000001, ...
        0.5, ...
        0.49249999999999999, ...
        0.42499999999999999, ...
        0.80300000000000005, ...
        0.5, ...
        0.69999999999999996, ...
        0.59999999999999998, ...
        0.54300000000000004, ...
        0.45500000000000002, ...
        0.503, ...
        0.5, ...
        0.5, ...
        0.55200000000000005, ...
        0.55200000000000005, ...
        0.55200000000000005, ...
        0.55200000000000005, ...
        0.54249999999999998, ...
        0.52900000000000003, ...
        0.75, ...
        0.75, ...
        0.59999999999999998, ...
        0.59999999999999998, ...
        0.50700000000000001, ...
        0.56599999999999995, ...
        0.55000000000000004, ...
        0.55000000000000004, ...
        0.55000000000000004, ...
        0.55000000000000004, ...
        0.55000000000000004, ...
        0.55000000000000004, ...
        0.55000000000000004, ...
        0.55000000000000004, ...
        0.55000000000000004, ...
        0.5, ...
        0.5, ...
        0.59999999999999998, ...
        0.59999999999999998, ...
        0.5, ...
        0.5, ...
        0.38650000000000001, ...
        0.5, ...
        0.49249999999999999, ...
        0.42499999999999999, ...
        0.80300000000000005, ...
        0.5, ...
        0.69999999999999996, ...
        0.59999999999999998, ...
        0.54300000000000004, ...
        0.45500000000000002, ...
        0.503];
    muscleRatio = containers.Map(muscleKeys, ratioVals);
    
    
    modelProcessor = ModelProcessor(modelFilename);
    % modelProcessor.append(ModOpAddExternalLoads('grf_walk.xml'));
    % modelProcessor.append(ModOpIgnoreTendonCompliance());
    % modelProcessor.append(ModOpReplaceMusclesWithDeGrooteFregly2016());
    % modelProcessor.append(ModOpIgnorePassiveFiberForcesDGF());
    % modelProcessor.append(ModOpScaleActiveFiberForceCurveWidthDGF(1.5));
    % modelProcessor.append(ModOpAddReserves(1.0));

    model = modelProcessor.process();
    
    muscleSet = model.getMuscles();
    numMuscles = muscleSet.getSize();
    
    probe_all = Umberger2010MuscleMetabolicsProbe();
    probe_all.setName('all_metabolics');
    
    probe_act = Umberger2010MuscleMetabolicsProbe();
    probe_act.setName('all_activation_maintenance_rate');
    probe_act.set_activation_maintenance_rate_on(true);
    probe_act.set_shortening_rate_on(false);
    probe_act.set_basal_rate_on(false);
    probe_act.set_mechanical_work_rate_on(false);
    
    probe_short = Umberger2010MuscleMetabolicsProbe();
    probe_short.setName("all_shortening_rate");
    probe_short.set_activation_maintenance_rate_on(false);
    probe_short.set_shortening_rate_on(true);
    probe_short.set_basal_rate_on(false);
    probe_short.set_mechanical_work_rate_on(false);
    
    probe_basal = Umberger2010MuscleMetabolicsProbe();
    probe_basal.setName("all_basal_rate");
    probe_basal.set_activation_maintenance_rate_on(false);
    probe_basal.set_shortening_rate_on(false);
    probe_basal.set_basal_rate_on(true);
    probe_basal.set_mechanical_work_rate_on(false);
    
    probe_mech = Umberger2010MuscleMetabolicsProbe();
    probe_mech.setName("all_mechanical_work_rate");
    probe_mech.set_activation_maintenance_rate_on(false);
    probe_mech.set_shortening_rate_on(false);
    probe_mech.set_basal_rate_on(false);
    probe_mech.set_mechanical_work_rate_on(true);    
    
    for m = 0:numMuscles-1
        musc = muscleSet.get(m);
        muscName = musc.getName();
        muscName = char(muscName);
        probe_all.addMuscle(muscName, muscleRatio(muscName));
        probe_act.addMuscle(muscName, muscleRatio(muscName));
        probe_short.addMuscle(muscName, muscleRatio(muscName));
        probe_basal.addMuscle(muscName, muscleRatio(muscName));
        probe_mech.addMuscle(muscName, muscleRatio(muscName));
    end
    
    model.addProbe(probe_all);
    model.addProbe(probe_act);
    model.addProbe(probe_short);
    model.addProbe(probe_basal);
    model.addProbe(probe_mech);
    
    %{
        % save this for later when everything is actually working
        
        for m = 0:numMuscles-1
            musc = muscleSet.get(m);
            muscName = musc.getName();
            muscName = char(muscName);
            
            % add a probe for each muscle
            probe = Umberger2010MuscleMetabolicsProbe();
            probe.setName(strcat("metabolics_", muscName));
            probe.addMuscle(muscName, muscleRatio(muscName));
            model.addProbe(probe);
            
            probe = Umberger2010MuscleMetabolicsProbe();
            probe.setName(strcat('activation_maintenance_rate_',muscName));
            probe.set_activation_maintenance_rate_on(true);
            probe.set_shortening_rate_on(false);
            probe.set_basal_rate_on(false);
            probe.set_mechanical_work_rate_on(false);
            probe.addMuscle(muscName, muscleRatio(muscName));
            model.addProbe(probe);
            
            probe = Umberger2010MuscleMetabolicsProbe();
            probe.setName(strcat("shortening_rate_",muscName));
            probe.set_activation_maintenance_rate_on(false);
            probe.set_shortening_rate_on(true);
            probe.set_basal_rate_on(false);
            probe.set_mechanical_work_rate_on(false);
            probe.addMuscle(muscName, muscleRatio(muscName));
            model.addProbe(probe);
            
            probe = Umberger2010MuscleMetabolicsProbe();
            probe.setName(strcat("basal_rate_",muscName));
            probe.set_activation_maintenance_rate_on(false);
            probe.set_shortening_rate_on(false);
            probe.set_basal_rate_on(true);
            probe.set_mechanical_work_rate_on(false);
            probe.addMuscle(muscName, muscleRatio(muscName));
            model.addProbe(probe);
            
            probe = Umberger2010MuscleMetabolicsProbe();
            probe.setName(strcat("mechanical_work_rate_",muscName));
            probe.set_activation_maintenance_rate_on(false);
            probe.set_shortening_rate_on(false);
            probe.set_basal_rate_on(false);
            probe.set_mechanical_work_rate_on(true);
            probe.addMuscle(muscName, muscleRatio(muscName));
            model.addProbe(probe);
            
        end
    %}
    
    model.finalizeConnections();
    model.print('simple_model_all_the_probes.osim');  

end


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



function torqueStateTrackGRFPrescribe()

    import org.opensim.modeling.*;

    % create and name an instance of the MocoTrack tool
    track = MocoTrack();
    track.setName("torque_statetrack_grfprescribe");

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
    modelProcessor = ModelProcessor('simple_model_all_the_probes.osim');
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
    %  open(pdfFilePath);
end



function muscleStatePrescribeGRFPrescribeWithEMG()

    import org.opensim.modeling.*;
    
    % setting up model and the kinematics
    inverse = MocoInverse();
    modelProcessor = ModelProcessor('simple_model_all_the_probes.osim');
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
    recfem = controlsRef.updDependentColumn('rectus_femoris');
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
    problem.addGoal(emgTracking);

    % can we add in the rest of the actuators from the file to reference, or are they ignored

    solution = study.solve();
    % study.visualize(solution);
    keyboard
    solution.write('muscle_stateprescribe_grfprescribe_withemg_solution.sto')
    STOFileAdapter.write(solution.exportToControlsTable(), 'testing_controls.sto')
    STOFileAdapter.write(solution.exportToStatesTable(), 'testing_states.sto')

    keyboard



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
             'bilateral', false, ...
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
end



function analyzeMetabolicCost()
    
    keyboard
    import org.opensim.modeling.*
    % Conduct an analysis using MuscleAnalysis and ProbeReporter.
    solution = MocoTrajectory('muscle_stateprescribe_grfprescribe_withemg_solution.sto');
    Time = solution.getTimeMat();
    numColPoints = solution.getNumTimes();
    
    % full moco method
    analyze = AnalyzeTool();
    analyze.setName("analyze");
    analyze.setModelFilename("simple_model_all_the_probes.osim");
    analyze.setStatesFileName("testing_states.sto");
    analyze.updAnalysisSet().cloneAndAppend(MuscleAnalysis());
    analyze.updAnalysisSet().cloneAndAppend(ProbeReporter());
    analyze.updControllerSet().cloneAndAppend(PrescribedController("testing_controls.sto"));
    analyze.setInitialTime(Time(1));
    analyze.setFinalTime(Time(end));
    analyze.print("testing_AnalyzeTool_setup.xml");
    % Run the analysis.
    analyze = AnalyzeTool("testing_AnalyzeTool_setup.xml");
    analyze.run();

    % tables
    table_force = TimeSeriesTable("analyze_MuscleAnalysis_ActiveFiberForce.sto");
    table_velocity = TimeSeriesTable("analyze_MuscleAnalysis_FiberVelocity.sto");
    table_metabolics = TimeSeriesTable('analyze_ProbeReporter_probes.sto');

    % time
    time_os = table_force.getIndependentColumn();
    time = [];
    for i=0:time_os.size()-1
        time = [time; time_os.get(i)];
    end
    time_met_os = table_metabolics.getIndependentColumn();
    time_met = [];
    for i=0:time_met_os.size()-1
        time_met = [time_met; time_met_os.get(i)];
    end

    % work through all the muscles
    muscles = [];
    labels = table_force.getColumnLabels();
    for i=0:labels.size()-1
        muscles = [muscles, labels.get(i)];
    end
    numMuscles = length(muscles);
    
    % force and velocity
    force = [];
    velocity = [];
    for i=1:length(muscles)
        temp_force = table_force.getDependentColumn(muscles(i)).getAsMat();
        force = [force, temp_force];
        temp_velocity = table_velocity.getDependentColumn(muscles(i)).getAsMat();
        velocity = [velocity, temp_velocity];
    end
    
    keyboard
    figure(1);
    hold on
    for i = 1:numMuscles
        plot(time, force(:,i).*-velocity(:,i))
    end
    
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Create an AnalyzeTool setup file.
    
    
    % controls
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
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % not sure about this stuff
    model = Model('subject_walk_armless.osim');
    musclesApoorva = model.getMuscles();
    
    probeSet = model.getProbeSet();
    probe = probeSet.get('metabolic_power');
    probeUmberger = Umberger2010MuscleMetabolicsProbe.safeDownCast(probe);
    
    rho = 1059.7;
    muscleEnergyRate = NaN(numColPoints, numMuscles);
    keyboard
    
    
    
    
    
    
    %{
    % big loop through all the muscles
    for m = 1:numMuscles()
        musc = musclesApoorva.get(muscles(m));
        Fmax = musc.getMaxIsometricForce();
        Lceopt = musc.getOptimalFiberLength();
        maxFiberVel = musc.getMaxContractionVelocity();
        
        rST = probeUmberger.getRatioSlowTwitchFibers(muscles(m));        
    end
    %}
    
    
end


