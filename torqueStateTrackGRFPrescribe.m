function torqueStateTrackGRFPrescribe()
    import org.opensim.modeling.*;
    
    % create and name an instance of the MocoTrack tool
    track = MocoTrack();
    track.setName("torque_statetrack_grfprescribe");
    
    % construct a ModelProcessor and add it to the tool.
    % modelProcessor = ModelProcessor("simple_model_all_the_probes_adjusted.osim");
    % modelProcessor = ModelProcessor("subject_redoarms.osim");
    modelProcessor = ModelProcessor("simple_model_all_the_probes.osim");


    
    % modelProcessor = ModelProcessor("simple_model_all_the_probes.osim");
    % now to do stuff with the model
    % modelProcessor = ModelProcessor(model);
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
    track.setModel(modelProcessor);
    
    % construct a TableProcessor of the coordinate data and pass it to the tracking tool. 
    % 1
    % track.setStatesReference(TableProcessor('torque_markertrack_grfprescribe_solution.sto'));
    % 2 
%     tableProcessor = TableProcessor('coordinates_updated.mot');
    tableProcessor = TableProcessor('results_IK_redoarms.mot');
    % 3
%     tableProcessor = TableProcessor(tabletrimming('coordinates_updated.mot')); %***
    tableProcessor.append(TabOpLowPassFilter(6));
    % 4 
    % tempTable = TimeSeriesTable('./ResultsRRA_2/subject01_walk1_RRA_Kinematics_q.sto');
    % tableProcessor = TableProcessor(tempTable);
    
    tableProcessor.append(TabOpUseAbsoluteStateNames());
    
    track.setStatesReference(tableProcessor);
    track.set_states_global_tracking_weight(10);
    % avoid exceptions if markers in file are no longer in the model (arms removed)
    track.set_allow_unused_references(true);
    % since there is only coordinate position data in the states references, 
    % this fills in the missing coordinate speed data using 
    % the derivative of splined position data
    track.set_track_reference_position_derivatives(true);
    
    % % set specific weights for the individual weight set
    coordinateweights = MocoWeightSet();
    coordinateweights.cloneAndAppend(MocoWeight("pelvis_tx", 0.01));
    coordinateweights.cloneAndAppend(MocoWeight("pelvis_ty", 0.01));
    coordinateweights.cloneAndAppend(MocoWeight("pelvis_tz", 0.01));
    coordinateweights.cloneAndAppend(MocoWeight("pelvis_list", 0.01));
    coordinateweights.cloneAndAppend(MocoWeight("pelvis_rotation", 0.01));
    coordinateweights.cloneAndAppend(MocoWeight("pelvis_tilt", 0.01));
    % coordinateweights.cloneAndAppend(MocoWeight("hip_rotation_r", 0));
    % coordinateweights.cloneAndAppend(MocoWeight("hip_rotation_l", 0));
    % coordinateweights.cloneAndAppend(MocoWeight("hip_flexion_r", 4));
    % coordinateweights.cloneAndAppend(MocoWeight("hip_flexion_l", 4));
    coordinateweights.cloneAndAppend(MocoWeight("knee_angle_r", 5));
    coordinateweights.cloneAndAppend(MocoWeight("knee_angle_l", 5));
    coordinateweights.cloneAndAppend(MocoWeight("ankle_angle_r", 6));
    coordinateweights.cloneAndAppend(MocoWeight("ankle_angle_l", 6));
    
    track.set_states_weight_set(coordinateweights);
    
    

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
    track.set_mesh_interval(0.01); %.01% may regret later

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
            effort.setWeightForControl(forcePath, 500); % here
            % if contains(string(forcePath), 'pelvis_ty')
            %     effort.setWeightForControl(forcePath, 1e8);
            % end
        end
        % if contains(string(forcePath), 'hip_rotation')
        %    effort.setWeightForControl(forcePath, 1e4);
        % end
    end
    
    % Constrain the states and controls to be periodic.
    periodicityGoal = MocoPeriodicityGoal("periodicity");
    periodicityGoal.setMode('endpoint_constraint');
    for i = 0:model.getNumStateVariables()-1
        currentStateName = string(model.getStateVariableNames().getitem(i));
        if (~contains(currentStateName,'pelvis_tx/value'))
            periodicityGoal.addStatePair(MocoPeriodicityGoalPair(currentStateName));
        end
    end
    forceSet = model.getForceSet();
    for i = 0:forceSet.getSize()-1
        forcePath = forceSet.get(i).getAbsolutePathString();
        periodicityGoal.addControlPair(MocoPeriodicityGoalPair(forcePath));
    end
    problem.addGoal(periodicityGoal);
    


    % contactTracking = osim.MocoContactTrackingGoal('contact', GRFTrackingWeight)    
    % contactTracking.setExternalLoadsFile('grf_walk_nat_1.xml');
    
    % forceNamesRightFoot = osim.StdVectorString();
    % forceNamesRightFoot.append('/contactHeel_r');
    % % # forceNamesRightFoot.append('/forceset/contactLateralRearfoot_r');
    % forceNamesRightFoot.append('/contactLateralMidfoot_r');
    % % # forceNamesRightFoot.append('/contactLateralToe_r');
    % forceNamesRightFoot.append('/contactMedialToe_r');
    % forceNamesRightFoot.append('/contactMedialMidfoot_r');
    % % # contactTracking.addContactGroup(forceNamesRightFoot, 'Right_GRF');
    % contactTrackingSplitRight = osim.MocoContactTrackingGoalGroup(forceNamesRightFoot, 'Right_GRF');
    % contactTrackingSplitRight.append_alternative_frame_paths('/bodyset/toes_r')
    % contactTracking.addContactGroup(contactTrackingSplitRight);

    % forceNamesLeftFoot = osim.StdVectorString();
    % forceNamesLeftFoot.append('/contactHeel_l');
    % % # forceNamesLeftFoot.append('/forceset/contactLateralRearfoot_l');
    % forceNamesLeftFoot.append('/contactLateralMidfoot_l');
    % % # forceNamesLeftFoot.append('/contactLateralToe_l');
    % forceNamesLeftFoot.append('/contactMedialToe_l');
    % forceNamesLeftFoot.append('/contactMedialMidfoot_l');
    % % # contactTracking.addContactGroup(forceNamesLeftFoot, 'Left_GRF');
    % contactTrackingSplitLeft = osim.MocoContactTrackingGoalGroup(forceNamesLeftFoot, 'Left_GRF');
    % contactTrackingSplitLeft.append_alternative_frame_paths('/bodyset/toes_l')
    % contactTracking.addContactGroup(contactTrackingSplitLeft);
    
    % contactTracking.setProjection('plane');
    % contactTracking.setProjectionVector(osim.Vec3(0, 0, 1));

    % contactTracking.setDivideByDisplacement(False)
    % contactTracking.setDivideByMass(True)
    % problem.addGoal(contactTracking);


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
    % twosteptraj = MocoTrajectory('torque_statetrack_grfprescribe_solution.sto');
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
    
    % solve and visualize
    solution = study.solve();
    % study.visualize(solution);
    % generate a report and save
    solution.write('torque_statetrack_grfprescribe_solution.sto');
    % study.visualize(MocoTrajectory("torque_statetrack_grfprescribe_solution.sto"));
        
    % report = osimMocoTrajectoryReport(model, 'torque_statetrack_grfprescribe_solution.sto');
    % reportFilePath = report.generate();
    % pdfFilePath = reportFilePath(1:end-2);
    % pdfFilePath = strcat(pdfFilePath, 'pdf');
    % ps2pdf('psfile',reportFilePath,'pdffile',pdfFilePath, ...
    %     'gscommand','C:\Program Files\gs\gs9.54.0\bin\gswin64.exe', ...
    %     'gsfontpath','C:\Program Files\gs\gs9.54.0\Resource\Font', ...
    %     'gslibpath','C:\Program Files\gs\gs9.54.0\lib');
    % open(pdfFilePath);
    % save('torque_statetrack_grfprescribe.mat');
    disp('end state torque track')
end