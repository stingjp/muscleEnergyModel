function muscleStatePrescribeGRFPrescribe()
    import org.opensim.modeling.*;

    % construct MocoInverse tool
    inverse = MocoInverse();

    % construct ModelProcessor and sit it on the tool. 
    % replace default muscles with degrootefregly 2016 muscles, and adjust params
    modelProcessor = ModelProcessor('simple_model_all_the_probes.osim');
    modelProcessor.append(ModOpAddExternalLoads('grf_walk.xml'));

    % set up the base model
    modelProcessor.append(ModOpIgnoreTendonCompliance());
    modelProcessor.append(ModOpReplaceMusclesWithDeGrooteFregly2016());
    % only valid for degroote
    modelProcessor.append(ModOpIgnorePassiveFiberForcesDGF());
    % only valid for degroote
    modelProcessor.append(ModOpScaleActiveFiberForceCurveWidthDGF(1.5));
    modelProcessor.append(ModOpAddReserves(1.0));
    
    % now do tweaks to get tendon compliance
    basemodel = modelProcessor.process();
    basemodel.initSystem();
    basemuscles = basemodel.updMuscles();
    numBaseMuscles = basemuscles.getSize();
    for m = 0:numBaseMuscles-1
        % set tendon compliance on for certain muscles
        % if lopt > lst want stiff (ignore)

        % get the muscle
        basemusc = basemuscles.get(m);
        % get lopt
        baselopt = basemusc.getOptimalFiberLength();
        % get lst
        baselst = basemusc.getTendonSlackLength();
        % set compliance if lopt > lst
        if baselopt < baselst
            basemusc.set_ignore_tendon_compliance(false)
        end
    end

    %% do more model processor stuff
    modelProcessorDC = ModelProcessor(basemodel);
    modelProcessorDC.append(ModOpFiberDampingDGF(0.01));
    % modelProcessorDC.append(ModOpAddReserves(1, 2.5, true));
    modelProcessorDC.append(ModOpTendonComplianceDynamicsModeDGF('implicit'));
    inverse.setModel(modelProcessorDC);


    tempkintable = TableProcessor('torque_statetrack_grfprescribe_solution.sto').process();
    templabels_os = tempkintable.getColumnLabels();
    % templabels = []
    for i=0:templabels_os.size()-1
        % templabels = [templabels, templabels_os.get(i)];
        temp = templabels_os.get(i);
        if ~temp.startsWith('/jointset')
            tempkintable.removeColumn(temp);
        end
    end

    inverse.setKinematics(TableProcessor(tempkintable));
    % inverse.setKinematics(TableProcessor('torque_statetrack_grfprescribe_solution.sto'));

    % construct TableProcessor of the coordinate data and pass it to the inverse tool
    % if no operators, it returns the base table
    % inverse.setKinematics(TableProcessor('torque_statetrack_grfprescribe_solution.sto'));    


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

    % set time and intervals
    inverse.set_initial_time(gait_start);
    inverse.set_final_time(gait_end);
    inverse.set_mesh_interval(0.02);
    % By default, Moco gives an error if the kinematics contains extra columns.
    % Here, we tell Moco to allow (and ignore) those extra columns.
    inverse.set_kinematics_allow_extra_columns(true);

    % set inverse goals
    inverse.set_minimize_sum_squared_activations(true);
    inverse.set_reserves_weight(20); % 10

    study = inverse.initialize();
    problem = study.updProblem();

    % add goals to the problem and scale them to get close to ~1
    % effortgoal = MocoControlGoal('effort');
    % effortgoal.setWeight(100);
    % problem.addGoal(effortgoal);

    % initactivationgoal = MocoInitialActivationGoal('init_activation');
    % initactivationgoal.setWeight(1);
    % problem.addGoal(initactivationgoal);
    

    % TODO test
    % excitation_effort goal
    excitegoal = problem.updGoal('excitation_effort');
    excitegoal.setWeight(1e-4); % 1e-4
    % 'activation_effort' goal
    % activegoal = problem.updGoal('activation_effort');
    % activegoal.setWeight(1e-4);
    


    %%% moving on to solve
    % set up the solver and solve the problem
    solver = MocoCasADiSolver.safeDownCast(study.updSolver());
    solver.resetProblem(problem);

    solution = study.solve();
    solution.insertStatesTrajectory(tempkintable);
    % study.visualize(solution);


    % post problem processing
    model = modelProcessorDC.process();
    model.print('post_simple_model_all_the_probes.osim');

    solution.write('muscle_stateprescribe_grfprescribe_solution.sto');
    STOFileAdapter.write(solution.exportToControlsTable(), 'muscleprescribe_controls.sto');
    STOFileAdapter.write(solution.exportToStatesTable(), 'muscleprescribe_states.sto');


    % Solve the problem and write the solution to a Storage file.
    % solution = inverse.solve(true); % to visualize
    % solution = inverse.solve();
    % solution.getMocoSolution().write('muscle_stateprescribe_grfprescribe_solution.sto');
    % STOFileAdapter.write(solution.getMocoSolution().exportToControlsTable(), 'muscleprescribe_controls.sto')
    % STOFileAdapter.write(solution.getMocoSolution().exportToStatesTable(), 'muscleprescribe_states.sto')
    
    
    % Generate a report with plots for the solution trajectory.
    % model = modelProcessor.process();
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


    keyboard
    computeIDFromResult(solution);

    % analyze the metabolic cost
    analyzeMetabolicCost(solution);
    % ID analysis and evaluate the reserves
%     keyboard
%     computeIDFromResult(model, solution);
end