function [Issues] = muscleStateTrackGRFPrescribe(Issues)
    import org.opensim.modeling.*;
    
    % create and name an instance of the MocoTrack tool
    track = MocoTrack();
    track.setName("muscle_statetrack_grfprescribe");
    
    % construct ModelProcessor and sit it on the tool. 
    % replace default muscles with degrootefregly 2016 muscles, and adjust params
    modelProcessor = ModelProcessor('simple_model_all_the_probes_adjusted.osim');
    modelProcessor.append(ModOpAddExternalLoads('grf_walk.xml'));
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

    % set up the base model
    % modelProcessor.append(ModOpIgnoreTendonCompliance());
    modelProcessor.append(ModOpReplaceMusclesWithDeGrooteFregly2016());
    % only valid for degroote
    % modelProcessor.append(ModOpIgnorePassiveFiberForcesDGF());
    % only valid for degroote
    modelProcessor.append(ModOpScaleActiveFiberForceCurveWidthDGF(1.5));
    modelProcessor.append(ModOpAddReserves(1.0));
    



    % now do tweaks to get tendon compliance
    basemodel = modelProcessor.process();
    basemodel.print('basemodel_simple_model_all_the_probes.osim');
    
    % turn on the probes for the study - I think RRA turns some off?
    basemodel = probeActivate(basemodel);

    % updates
    basemodel.initSystem();
    basemuscles = basemodel.updMuscles();
    numBaseMuscles = basemuscles.getSize();
    % for m = 0:numBaseMuscles-1
    %     set tendon compliance on for certain muscles
    %     if lopt > lst want stiff (ignore)

    %     get the muscle
    %     basemusc = basemuscles.get(m);
    %     get lopt
    %     baselopt = basemusc.getOptimalFiberLength();
    %     get lst
    %     baselst = basemusc.getTendonSlackLength();
    %     set compliance if lopt > lst
    %     if baselopt < baselst
    %         basemusc.set_ignore_tendon_compliance(false)
    %     end
    % end

    %% do more model processor stuff
    modelProcessorDC = ModelProcessor(basemodel);
    modelProcessorDC.append(ModOpFiberDampingDGF(0.01));
    % modelProcessorDC.append(ModOpAddReserves(1, 2.5, true));
    modelProcessorDC.append(ModOpTendonComplianceDynamicsModeDGF('implicit'));

    % need to add the bhargava metabolics probe for cost function
    % bhargmet = Bhargava2004SmoothedMuscleMetabolics();
    % bhargmet.setName("simmetabolics");
    % bhargmet.set_use_smoothing(true);
    % modelmet = modelProcessorDC.process();
    % modelmetMuscles = modelmet.getMuscles();
    % nummodelmetMuscles = modelmetMuscles.getSize();
    % for m = 0:nummodelmetMuscles-1
    %     musc = modelmetMuscles.get(m);
    %     muscName = musc.getName();
    %     muscName = char(muscName);
    %     bhargmet.addMuscle(muscName, musc);
    % end
    
    % modelmet.addComponent(bhargmet);
    % modelmet.finalizeConnections();
    % modelProcessorDC2 = ModelProcessor(modelmet);
    % inverse.setModel(modelProcessorDC2);
    
    track.setModel(modelProcessorDC)



    % construct a TableProcessor of the coordinate data and pass it to the tracking tool. 
    % 1
    % track.setStatesReference(TableProcessor('torque_markertrack_grfprescribe_solution.sto'));
    % 2 
%     tableProcessor = TableProcessor('coordinates_updated.mot');
    % 3
    tableProcessor = TableProcessor(tabletrimming('coordinates_updated.mot')); %***
    tableProcessor.append(TabOpLowPassFilter(6));
    % 4 
%     tempTable =
%     TimeSeriesTable('./ResultsRRA_1/subject01_walk1_RRA_Kinematics_q.sto');
%     % was using these!!!
%     tableProcessor = TableProcessor(tempTable);
%     tableProcessor.append(TabOpLowPassFilter(6));

    
    tableProcessor.append(TabOpUseAbsoluteStateNames());
    
    track.setStatesReference(tableProcessor);
    
%     track.set_kinematics_allow_extra_columns(true);
    track.set_states_global_tracking_weight(10); % 50 % need to weigh benefit of higher global vs specific coordinate
    % avoid exceptions if markers in file are no longer in the model (arms removed)
    track.set_allow_unused_references(true);
    % since there is only coordinate position data in the states references, 
    % this fills in the missing coordinate speed data using 
    % the derivative of splined position data
    track.set_track_reference_position_derivatives(true);
    
    % set specific weights for the individual weight set
    coordinateweights = MocoWeightSet();
%     coordinateweights.cloneAndAppend(MocoWeight("pelvis_tx", 1000000));
%     coordinateweights.cloneAndAppend(MocoWeight("pelvis_ty", 1000000));
%     coordinateweights.cloneAndAppend(MocoWeight("pelvis_tz", 1000000));
    coordinateweights.cloneAndAppend(MocoWeight("pelvis_list", 1000000));
    coordinateweights.cloneAndAppend(MocoWeight("pelvis_rotation", 1000000));
    coordinateweights.cloneAndAppend(MocoWeight("pelvis_tilt", 1000000));
%     coordinateweights.cloneAndAppend(MocoWeight("hip_rotation_r", 1000));
%     coordinateweights.cloneAndAppend(MocoWeight("hip_rotation_l", 1000));
%     coordinateweights.cloneAndAppend(MocoWeight("hip_adduction_r", 100000));
%     coordinateweights.cloneAndAppend(MocoWeight("hip_adduction_l", 100000));
    coordinateweights.cloneAndAppend(MocoWeight("ankle_angle_r", 1000000));
    coordinateweights.cloneAndAppend(MocoWeight("ankle_angle_l", 1000000));
    track.set_states_weight_set(coordinateweights);

    % get the subject name and gait timings
    load 'C:\Users\JP\code\repos\Stanford\delplab\projects\muscleModel\muscleEnergyModel\subjectgaitcycles.mat';
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
    track.set_mesh_interval(0.05); % 0.03 for all current subjects %.05 % .01% 
    
    % initialize and set goals
    study = track.initialize();    

    % get reference to the MocoControlGoal that is added to every MocoTrack problem
    problem = study.updProblem();

    % set a constraint so that the model doesnt overlap feet
    distance = MocoFrameDistanceConstraint();
    distance.setName('minimum_distance');
    distance.addFramePair(java.lang.String('/bodyset/calcn_l'), java.lang.String('/bodyset/calcn_r'), 0.15, Inf); % 0.20
    distance.addFramePair(java.lang.String('/bodyset/toes_l'), java.lang.String('/bodyset/toes_r'), 0.15, Inf); %0.20
    distance.addFramePair(java.lang.String('/bodyset/calcn_l'), java.lang.String('/bodyset/toes_r'), 0.15, Inf); %0.20
    distance.addFramePair(java.lang.String('/bodyset/toes_l'), java.lang.String('/bodyset/calcn_r'), 0.15, Inf); %0.20
    problem.addPathConstraint(distance);
    
    % effort goal
    effort = MocoControlGoal.safeDownCast(problem.updGoal('control_effort'));
    effort.setWeight(.5); %0.1 for the new

    initactivationgoal = MocoInitialActivationGoal('init_activation');
    initactivationgoal.setWeight(10);
    problem.addGoal(initactivationgoal);


    % put large weight on the pelvis CoordinateActuators, which act as the 
    % residual, or 'hand-of-god' forces which we would like to keep small
    
    model = modelProcessorDC.process();
    model.print('post_simple_model_all_the_probes_muscletrack.osim');
    model.initSystem();
    forceSet = model.getForceSet();
    for i=0:forceSet.getSize()-1
        forcePath = forceSet.get(i).getAbsolutePathString();
        if contains(string(forcePath), 'pelvis')
            effort.setWeightForControl(forcePath, 10); % here 1000
            % if contains(string(forcePath), 'pelvis_ty')
            %     effort.setWeightForControl(forcePath, 1e8);
            % end
        end
%         if contains(string(forcePath), 'hip_rotation')
%            effort.setWeightForControl(forcePath, 10);
%         end
    end
    

    % track.set_guess_file('muscle_stateprescribe_grfprescribe_solution.sto');
    % twosteptraj = MocoTrajectory('muscle_stateprescribe_grfprescribe_solution.sto');
    % solver.setGuess(twosteptraj);
    % solver.setGuessFile('muscle_stateprescribe_grfprescribe_solution.sto');
    % solver.set_optim_convergence_tolerance(1e-4);


    % solver = study.initCasADiSolver();
    % guess = solver.createGuess('bounds'); % bounds or random
    % solver.setGuess(guess);

    twosteptraj = MocoTrajectory('muscle_stateprescribe_grfprescribe_solution.sto');
    steps = twosteptraj.getNumTimes();

    solver = MocoCasADiSolver.safeDownCast(study.updSolver());
    solver.resetProblem(problem)
    
%     solver.set_optim_convergence_tolerance(1e-4); % 1e-2
%     solver.set_optim_constraint_tolerance(1e-4); % 1e-2
    % solver.set_num_mesh_intervals(steps);

    guess = solver.createGuess('bounds'); % bounds or random  
    guess.write('boundsguess.sto');
    % solver.setGuess(guess);

    % twosteptraj = MocoTrajectory('muscle_stateprescribe_grfprescribe_solution.sto');
    randomguess = MocoTrajectory('boundsguess.sto');
    randomguess.resampleWithNumTimes(steps);
    
    % go through and overwrite the states first
    randomstatenames = randomguess.getStateNames();
    % this will cover joint values, speeds, muscle activations, and norm
    % tendon force
    for s = 0:randomstatenames.size()-1
        statename = randomstatenames.get(s);
        % temprandom = randomguess.getStateMat(statename);
        temp2step = twosteptraj.getStateMat(statename);
        randomguess.setState(statename,temp2step);       
    end
    
    % go through all the controls - excitations
    randomcontrolnames = randomguess.getControlNames();
    % this covers all excitations and reserves
    for c = 0:randomcontrolnames.size()-1
        controlname = randomcontrolnames.get(c);
        % temprandom = randomguess.getControlMat(controlname);
        temp2step = twosteptraj.getControlMat(controlname);
        randomguess.setControl(controlname, temp2step);
    end
    
    % go through others??
    % randomparamnames = randomguess.getParameterNames();
    % this is empty in the normal condition
        
    % multipliers
    randommultnames = randomguess.getMultiplierNames();
    for m = 0:randommultnames.size()-1
        multname = randommultnames.get(m);
        % temprandom = randomguess.getMultiplierMat(multname)
        try
            temp2step = twosteptraj.getMultiplierMat(mutlname);
            randomguess.setMultiplier(multname, temp2step);
        catch
            disp('did not have the multiplier in the 2 step problem solution');
        end
    end
    
    % now for the implicit derivatives
    randomderivnames = randomguess.getDerivativeNames();
    for d = 0:randomderivnames.size()-1
        derivname = randomderivnames.get(d);
        % temprandom = randomguess.getDerivativeMat(derivname);
        temp2step = twosteptraj.getDerivativeMat(derivname);
        randomguess.setDerivative(derivname, temp2step);
    end    
    
    
    % now set the guess for the solver
    solver.setGuess(randomguess);

    
    % solve and visualize
    
    solution = study.solve();
%     solution = MocoTrajectory('muscle_statetrack_grfprescribe_solution.sto');
    % study.visualize(solution);
    % generate a report and save
    solution.write('muscle_statetrack_grfprescribe_solution.sto');
    % study.visualize(MocoTrajectory("torque_statetrack_grfprescribe_solution.sto"));
    
    STOFileAdapter.write(solution.exportToControlsTable(), 'muscleprescribe_controls.sto');
    STOFileAdapter.write(solution.exportToStatesTable(), 'muscleprescribe_states.sto');

        
    report = osimMocoTrajectoryReport(model, ...
                                    'muscle_statetrack_grfprescribe_solution.sto', ...
                                    'bilateral', true);
    reportFilePath = report.generate();
    pdfFilePath = reportFilePath(1:end-2);
    pdfFilePath = strcat(pdfFilePath, 'pdf');
    ps2pdf('psfile',reportFilePath,'pdffile',pdfFilePath, ...
        'gscommand','C:\Program Files\gs\gs9.54.0\bin\gswin64.exe', ...
        'gsfontpath','C:\Program Files\gs\gs9.54.0\Resource\Font', ...
        'gslibpath','C:\Program Files\gs\gs9.54.0\lib');
    % open(pdfFilePath);
    % save('torque_statetrack_grfprescribe.mat');
    disp('end state muscle track')


    % post analysis and validation
    
    Issues = [Issues; [java.lang.String('muscledrivensim'); java.lang.String('trackingproblem')]];
    analyzeMetabolicCost(solution);
    Issues = computeIDFromResult(Issues, solution);
    analyzeMetabolicCost(solution);
    computeKinematicDifferences(solution);

end