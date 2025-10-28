function torqueMarkerTrackGRFPrescribe()
    
    import org.opensim.modeling.*;

    % create and name an instance of the MocoTrack tool
    track = MocoTrack();
    track.setName("torque_markertrack_grfprescribe");

    % construct a ModelProcessor and add it to the tool.
    modelProcessor = ModelProcessor("simple_model_all_the_probes.osim");
    weldem = StdVectorString();
    weldem.add('subtalar_r');
    weldem.add('mtp_r');
    weldem.add('subtalar_l');
    weldem.add('mtp_l');
    weldem.add('radius_hand_r');
    weldem.add('radius_hand_l');
    modelProcessor.append(ModOpReplaceJointsWithWelds(weldem));
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
    track.setMarkersReferenceFromTRC("motion_capture.trc");
    % avoid exceptions if markers in file are no longer in the model (arms removed)
    track.set_allow_unused_references(true);
    % increase global marker tracking weight, which is associated with internal 
    % MocoMarkerTrackingGoal term
    track.set_markers_global_tracking_weight(100);
    % increase tracking weights for individual markers in the data set
    markerWeights = MocoWeightSet();
    markerWeights.cloneAndAppend(MocoWeight("R.ASIS", 20));
    markerWeights.cloneAndAppend(MocoWeight("L.ASIS", 20));
    markerWeights.cloneAndAppend(MocoWeight("R.PSIS", 20));
    markerWeights.cloneAndAppend(MocoWeight("L.PSIS", 20));
    markerWeights.cloneAndAppend(MocoWeight("R.Knee", 15));
    markerWeights.cloneAndAppend(MocoWeight("R.Ankle", 100));
    markerWeights.cloneAndAppend(MocoWeight("R.Heel", 100));
    markerWeights.cloneAndAppend(MocoWeight("R.MT5", 100));
    markerWeights.cloneAndAppend(MocoWeight("R.Toe", 100));
    markerWeights.cloneAndAppend(MocoWeight("L.Knee", 15));
    markerWeights.cloneAndAppend(MocoWeight("L.Ankle", 100));
    markerWeights.cloneAndAppend(MocoWeight("L.Heel", 100));
    markerWeights.cloneAndAppend(MocoWeight("L.MT5", 100));
    markerWeights.cloneAndAppend(MocoWeight("L.Toe", 100));
    track.set_markers_weight_set(markerWeights);


    % get the subject name and gait timings
    load 'C:\Users\jonstingel\code\musclemodel\muscleEnergyModel\subjectgaitcycles.mat';
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
    effort.setWeight(10);
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
            effort.setWeightForControl(forcePath, 1000); % here
            % if contains(string(forcePath), 'pelvis_ty')
            %     effort.setWeightForControl(forcePath, 1e8);
            % end
        end
        % if contains(string(forcePath), 'hip_rotation')
        %    effort.setWeightForControl(forcePath, 1e4);
        % end
    end
    

    % solver changes
    solver = MocoCasADiSolver.safeDownCast(study.updSolver());
    % solver.resetProblem(problem);
    solver.set_optim_convergence_tolerance(1e-2); % 1e-2
    solver.set_optim_constraint_tolerance(1e-2); % 1e-2
    % solver.set_minimize_implicit_auxiliary_derivatives(true);
    % solver.set_implicit_auxiliary_derivatives_weight(1e-8); 
    solver.set_optim_finite_difference_scheme('forward');
    solver.set_parameters_require_initsystem(false);


    % % create an initial guess
    % twosteptraj = MocoTrajectory('torque_markertrack_grfprescribe_solution.sto');
    % steps = twosteptraj.getNumTimes();
    % guess = solver.createGuess('bounds'); % bounds or random  
    % guess.write('boundsguess.sto');
    % % solver.setGuess(guess);

    % randomguess = MocoTrajectory('boundsguess.sto');
    % randomguess.resampleWithNumTimes(steps);
    % % go through and overwrite the states first
    % randomstatenames = randomguess.getStateNames();
    % % this will cover joint values, speeds, muscle activations, and norm
    % % tendon force
    % for s = 0:randomstatenames.size()-1
    %     statename = randomstatenames.get(s);
    %     % temprandom = randomguess.getStateMat(statename);
    %     temp2step = twosteptraj.getStateMat(statename);
    %     randomguess.setState(statename,temp2step);       
    % end
    % % go through all the controls - excitations
    % randomcontrolnames = randomguess.getControlNames();
    % % this covers all excitations and reserves
    % for c = 0:randomcontrolnames.size()-1
    %     controlname = randomcontrolnames.get(c);
    %     % temprandom = randomguess.getControlMat(controlname);
    %     temp2step = twosteptraj.getControlMat(controlname);
    %     randomguess.setControl(controlname, temp2step);
    % end
    % % go through others??
    % % randomparamnames = randomguess.getParameterNames();
    % % this is empty in the normal condition 
    % % multipliers
    % randommultnames = randomguess.getMultiplierNames();
    % for m = 0:randommultnames.size()-1
    %     multname = randommultnames.get(m);
    %     % temprandom = randomguess.getMultiplierMat(multname)
    %     try
    %         temp2step = twosteptraj.getMultiplierMat(mutlname);
    %         randomguess.setMultiplier(multname, temp2step);
    %     catch
    %         disp('did not have the multiplier in the 2 step problem solution');
    %     end
    % end
    % % now for the implicit derivatives
    % randomderivnames = randomguess.getDerivativeNames();
    % for d = 0:randomderivnames.size()-1
    %     derivname = randomderivnames.get(d);
    %     % temprandom = randomguess.getDerivativeMat(derivname);
    %     temp2step = twosteptraj.getDerivativeMat(derivname);
    %     randomguess.setDerivative(derivname, temp2step);
    % end    
    % % now set the guess for the solver
    % solver.setGuess(randomguess);

    % solve - the bool indicates to visualize the solution
    % solution = track.solve(true); % to visualize
    solution = study.solve();
    solution.write('torque_markertrack_grfprescribe_solution.sto');

    % generate a pdf report containing plots of the variables in the solution. 
    % for details see osimMocoTrajectoryReport.m in moco resource/code/matlab/utilities 
    % model = modelProcessor.process();
    % report = osimMocoTrajectoryReport(model, 'torque_markertrack_grfprescribe_solution.sto');
    % reportFilePath = report.generate();
    % pdfFilePath = reportFilePath(1:end-2);
    % pdfFilePath = strcat(pdfFilePath, 'pdf');
    % ps2pdf('psfile',reportFilePath,'pdffile',pdfFilePath, ...
    %     'gscommand','C:\Program Files\gs\gs9.54.0\bin\gswin64.exe', ...
    %     'gsfontpath','C:\Program Files\gs\gs9.54.0\Resource\Font', ...
    %     'gslibpath','C:\Program Files\gs\gs9.54.0\lib');
    % open(pdfFilePath);
    % save('torque_markertrack_grfprescribe.mat');
    disp('end marker track');
end