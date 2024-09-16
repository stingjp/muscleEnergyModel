function [Issues] = muscleStatePrescribeGRFPrescribeStrict(Issues)
    import org.opensim.modeling.*;

    tag = 'muscleprescribestrict';
    % construct MocoInverse tool
    inverse = MocoInverse();
    
    % need to get the geometry path for the pathwrapset
    
    % construct ModelProcessor and sit it on the tool. 
    % replace default muscles with degrootefregly 2016 muscles, and adjust params
    modelProcessor = ModelProcessor('simple_model_all_the_probes_adjusted.osim');

%     modelProcessor = ModelProcessor("simple_model_all_the_probes.osim");
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
%     modelProcessor.append(ModOpAddReserves(1.0));
    modelProcessor.append(ModOpAddReserves(1.0));
    



    % now do tweaks to get tendon compliance
    basemodel = modelProcessor.process();

    % turn on the probes for the study
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

    inverse.setModel(modelProcessorDC);


    tempkintable = TableProcessor('torque_statetrack_grfprescribe_strict_solution.sto').process();
    % tempkintable = TableProcessor('torque_statetrack_grfprescribe_solution.sto').process();
    % tempkintable = TableProcessor('./ResultsRRA_1/subject01_walk1_RRA_states.sto').process();
    % tempkintable = TableProcessor('./coordinates_updated.mot').process();
    templabels_os = tempkintable.getColumnLabels();
    % templabels = []
    for i=0:templabels_os.size()-1
        % templabels = [templabels, templabels_os.get(i)];
        temp = templabels_os.get(i);
        if ~startsWith(temp, '/jointset') %~temp.startsWith('/jointset')
            tempkintable.removeColumn(temp);
        end
    end

    inverse.setKinematics(TableProcessor(tempkintable));
    % inverse.setKinematics(TableProcessor('torque_statetrack_grfprescribe_solution.sto'));

    % construct TableProcessor of the coordinate data and pass it to the inverse tool
    % if no operators, it returns the base table
    % inverse.setKinematics(TableProcessor('torque_statetrack_grfprescribe_solution.sto'));    


    % get the subject name and gait timings
    load 'G:\Shared drives\Exotendon\muscleModel\muscleEnergyModel\subjectGaitCycles.mat';
    workdir = pwd;
    [~,trialname,~] = fileparts(pwd);
    cd ../
    [~,conditionname,~] = fileparts(pwd);
    cd ../
    [~,subjectname,~] = fileparts(pwd);
    cd(workdir);

    gait_start = subjectgaitcycles.(genvarname(subjectname)).(genvarname(conditionname)).(genvarname(trialname)).initial;
    gait_end = subjectgaitcycles.(genvarname(subjectname)).(genvarname(conditionname)).(genvarname(trialname)).final;

    % set time and intervals
    inverse.set_initial_time(gait_start);
    inverse.set_final_time(gait_end);
    inverse.set_mesh_interval(0.02); %.05 .02, .01% may need to adjust this
    % By default, Moco gives an error if the kinematics contains extra columns.
    % Here, we tell Moco to allow (and ignore) those extra columns.
    inverse.set_kinematics_allow_extra_columns(true);

    % set inverse goals
    inverse.set_minimize_sum_squared_activations(true);
    inverse.set_reserves_weight(1e-10);% 3e-2 30

    study = inverse.initialize();
    problem = study.updProblem();

    % TODO test
    % excitation_effort goal
    excitegoal = problem.updGoal('excitation_effort');
    excitegoal.setWeight(5e-4); % 9e-1 5e-4 2.5e-4
    % 'activation_effort' goal
    % activegoal = problem.updGoal('activation_effort');
    % activegoal.setWeight(5e-4);

    % add metabolic cost goal
   
    % metGoal = MocoOutputGoal("met", 0.01);
    % metGoal.setOutputPath("/simmetabolics|total_metabolic_rate");
    % metGoal.setDivideByDisplacement(true);
    % metGoal.setDivideByMass(true);
    % problem.addGoal(metGoal);


    % add goals to the problem and scale them to get close to ~1
    % effortgoal = MocoControlGoal('effort');
    % effortgoal.setWeight(5e-3);
    % problem.addGoal(effortgoal);

    initactivationgoal = MocoInitialActivationGoal('init_activation');
    initactivationgoal.setWeight(1); % 1
    problem.addGoal(initactivationgoal);
    
    % for post problem processing
    model = modelProcessorDC.process();
    model.print('post_simple_model_all_the_probes_muscleprescribe.osim');

    %%% moving on to solve
    % set up the solver and solve the problem
    solver = MocoCasADiSolver.safeDownCast(study.updSolver());
    solver.resetProblem(problem);
    solver.set_optim_convergence_tolerance(.001); % 1e-2
    solver.set_optim_constraint_tolerance(1e-4); % 1e-2
    
    solution = study.solve();
    solution.insertStatesTrajectory(tempkintable);
    
    % solution = MocoTrajectory('muscle_stateprescribe_grfprescribe_solution.sto');
    
    % solution.write('muscleguess.sto');
    % study.visualize(solution);


    % post processing
    % solution.write('muscle_statetrack_grfprescribe_solution.sto');
    solution.write('muscle_stateprescribe_grfprescribe_strict_solution.sto')
    STOFileAdapter.write(solution.exportToControlsTable(), 'muscleprescribe_strict_controls.sto');
    STOFileAdapter.write(solution.exportToStatesTable(), 'muscleprescribe_strict_states.sto');

    % STOFileAdapter.write(solution.exportToControlsTable(), 'muscleprescribe_controls_old.sto');
    % STOFileAdapter.write(solution.exportToStatesTable(), 'muscleprescribe_states_old.sto');


    % Solve the problem and write the solution to a Storage file.
    % solution = inverse.solve(true); % to visualize
    % solution = inverse.solve();
    % solution.getMocoSolution().write('muscle_stateprescribe_grfprescribe_solution.sto');
    % STOFileAdapter.write(solution.getMocoSolution().exportToControlsTable(), 'muscleprescribe_controls.sto')
    % STOFileAdapter.write(solution.getMocoSolution().exportToStatesTable(), 'muscleprescribe_states.sto')
    
    
    % Generate a report with plots for the solution trajectory.
    % model = modelProcessor.process();
    report = osimMocoTrajectoryReport(model, ...
            'muscle_stateprescribe_grfprescribe_strict_solution.sto', 'bilateral', true);
    % The report is saved to the working directory.
    reportFilePath = report.generate();
    pdfFilePath = reportFilePath(1:end-2);
    pdfFilePath = strcat(pdfFilePath, 'pdf');
    ps2pdf('psfile',reportFilePath,'pdffile',pdfFilePath, ...
        'gscommand','C:\Program Files\gs\gs9.54.0\bin\gswin64.exe', ...
        'gsfontpath','C:\Program Files\gs\gs9.54.0\Resource\Font', ...
        'gslibpath','C:\Program Files\gs\gs9.54.0\lib');
    %  open(pdfFilePath);
    
    
    
    % post analysis and validation
    Issues = [Issues; [java.lang.String('muscledrivensim'); java.lang.String('inverseproblem')]];
    analyzeMetabolicCost(solution, 'muscleprescribe');
    % Issues = computeIDFromResult(Issues, solution, tag);
    % analyzeMetabolicCost(solution);
    % trackorprescribe = 'prescribe';
    % computeKinematicDifferences(solution, trackorprescribe);
    
end