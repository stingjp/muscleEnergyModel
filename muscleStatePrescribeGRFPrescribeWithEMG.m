function [Issues] = muscleStatePrescribeGRFPrescribeWithEMG(Issues)

    import org.opensim.modeling.*;
    
    % setting up model and the kinematics
    inverse = MocoInverse();
    modelProcessor = ModelProcessor('simple_model_all_the_probes_adjusted.osim');
    modelProcessor.append(ModOpAddExternalLoads('grf_walk.xml'));

    % now to do stuff with the model
    % modelProcessor = ModelProcessor(model);
    % need to adjust some of the joints - weld them
    weldem = StdVectorString();
    weldem.add('subtalar_r');
    weldem.add('mtp_r');
    weldem.add('subtalar_l');
    weldem.add('mtp_l');
    % weldem.add('radius_hand_r');
    % weldem.add('radius_hand_l');
    
    modelProcessor.append(ModOpReplaceJointsWithWelds(weldem));
    % model = modelProcessor.process();
    % set up the base model
    modelProcessor.append(ModOpIgnoreTendonCompliance());
    modelProcessor.append(ModOpReplaceMusclesWithDeGrooteFregly2016());
    % only for DGF
    modelProcessor.append(ModOpIgnorePassiveFiberForcesDGF());
    modelProcessor.append(ModOpScaleActiveFiberForceCurveWidthDGF(1.5));
    modelProcessor.append(ModOpAddReserves(1.0));


    % now do tweaks to get tendon compliance
    basemodel = modelProcessor.process();
    
    % turn on the probes for the study
    basemodel = probeActivate(basemodel);

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

    inverse.set_initial_time(gait_start);
    inverse.set_final_time(gait_end);
    inverse.set_mesh_interval(0.02);
    inverse.set_kinematics_allow_extra_columns(true);

    % goal act^2
    inverse.set_minimize_sum_squared_activations(true);
    inverse.set_reserves_weight(20); % 30 % 10

    study = inverse.initialize();
    problem = study.updProblem();

    % add goals to the problem and scale them to get close to ~1
    % effortgoal = MocoControlGoal('effort');
    % effortgoal.setWeight(0.0001);
    % problem.addGoal(effortgoal);

    % initactivationgoal = MocoInitialActivationGoal('init_activation');
    % initactivationgoal.setWeight(0.0001);
    % problem.addGoal(initactivationgoal);

    % testing
    % excitation_effort goal
    excitegoal = problem.updGoal('excitation_effort');
    excitegoal.setWeight(1e-4); % 5e-4  % 1e-4
    % 'activation_effort' goal
    % activegoal = problem.updGoal('activation_effort');
    % activegoal.setWeight(1e-4);
    

    % add emg tracking stuff
    model = modelProcessorDC.process();
    muscleset = model.getMuscles();
    firstmuscle = muscleset.get(0).getName();
        
    % add EMG tracking
    emgTracking = MocoControlTrackingGoal('emg_tracking');
    emgTracking.setWeight(10.0);
    % each column in electromyography.sto is normalized so the maximum value is 1 for each
    controlsRef = TimeSeriesTable('electromyography.sto');

    % Scalde down the tracked muscle activity based on the peak levels from 
    % 'Gait analysis: normal and pathological function' by Perry and Burnfield, 2010 
    soleus = controlsRef.updDependentColumn('SOL');
    gasmed = controlsRef.updDependentColumn('GAS');
    tibant = controlsRef.updDependentColumn('TA');
    medham = controlsRef.updDependentColumn('MH');
    bicfem = controlsRef.updDependentColumn('BF');
    vaslat = controlsRef.updDependentColumn('VL');
    vasmed = controlsRef.updDependentColumn('VM');
    recfem = controlsRef.updDependentColumn('RF');
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
    
    if strcmp(firstmuscle, 'addbrev_r')
        % associate actuators in the model with columns in electromyography.sto
        emgTracking.setReferenceLabel('/forceset/soleus_r','SOL');
        emgTracking.setReferenceLabel('/forceset/gasmed_r','GAS');
        emgTracking.setReferenceLabel('/forceset/gaslat_r','GAS');
        emgTracking.setReferenceLabel('/forceset/tibant_r','TA');
        emgTracking.setReferenceLabel('/forceset/semimem_r','MH');
        emgTracking.setReferenceLabel('/forceset/bflh_r','BF');
        emgTracking.setReferenceLabel('/forceset/vaslat_r','VL');
        emgTracking.setReferenceLabel('/forceset/vasmed_r','VM');
        emgTracking.setReferenceLabel('/forceset/recfem_r','RF');
        % emgTracking.setReferenceLabel('/forceset/glmax1_r','gluteus_maximus');
        % emgTracking.setReferenceLabel('/forceset/glmax2_r','gluteus_maximus');
        % emgTracking.setReferenceLabel('/forceset/glmax2_r','gluteus_maximus');
    else
        % associate actuators in the model with columns in electromyography.sto
        emgTracking.setReferenceLabel('/forceset/soleus_r','SOL');
        emgTracking.setReferenceLabel('/forceset/med_gas_r','GAS');
        emgTracking.setReferenceLabel('/forceset/lat_gas_r','GAS');
        emgTracking.setReferenceLabel('/forceset/tib_ant_r','TA');
        emgTracking.setReferenceLabel('/forceset/semimem_r','MH');
        emgTracking.setReferenceLabel('/forceset/bifemlh_r','BF');
        emgTracking.setReferenceLabel('/forceset/vas_lat_r','VL');
        emgTracking.setReferenceLabel('/forceset/vas_med_r','VM');
        emgTracking.setReferenceLabel('/forceset/rect_fem_r','RF');
        % emgTracking.setReferenceLabel('/forceset/glmax1_r','gluteus_maximus');
        % emgTracking.setReferenceLabel('/forceset/glmax2_r','gluteus_maximus');
        % emgTracking.setReferenceLabel('/forceset/glmax2_r','gluteus_maximus');
    end
        
    problem.addGoal(emgTracking);

    % can we add in the rest of the actuators from the file to reference, or are they ignored

    solver = MocoCasADiSolver.safeDownCast(study.updSolver());
    solver.resetProblem(problem);

    solution = study.solve();
    solution.insertStatesTrajectory(tempkintable);
    % study.visualize(solution);
    solution.write('muscle_stateprescribe_grfprescribe_withemg_solution.sto')
    STOFileAdapter.write(solution.exportToControlsTable(), 'muscleprescribewithemg_controls.sto')
    STOFileAdapter.write(solution.exportToStatesTable(), 'muscleprescribewithemg_states.sto')

    
    % write the reference data in a way that's easy to compare to the solution.
    % controlsRef.removeColumn('medial_hamstrings');
    % controlsRef.removeColumn('biceps_femoris');
    % controlsRef.removeColumn('vastus_lateralis');
    % controlsRef.removeColumn('vastus_medius');
    % controlsRef.removeColumn('rectus_femoris');
    controlsRef.removeColumn('GMAX');
    controlsRef.removeColumn('GMED');
    columnLabels = StdVectorString();
    
    if strcmp(firstmuscle, 'addbrev_r')
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
    else
        columnLabels.add('/forceset/soleus_r');
        columnLabels.add('/forceset/med_gas_r');
        columnLabels.add('/forceset/tib_ant_r');
        columnLabels.add('/forceset/semimem_r');
        columnLabels.add('/forceset/bifemlh_r');
        columnLabels.add('/forceset/vas_lat_r');
        columnLabels.add('/forceset/vas_lmed_r');
        columnLabels.add('/forceset/rect_fem_r');
        controlsRef.setColumnLabels(columnLabels);
        controlsRef.appendColumn('/forceset/lat_gas_r', gasmed);
    end
    
    STOFileAdapter.write(controlsRef, 'controls_reference.sto');

    % generate a report comparing MocoInverse solutions without and with EMG tracking
    % model = modelProcessor.process();
    model.print('post_simple_model_all_the_probes.osim');
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

    % post analysis and validation
    Issues = [Issues; [java.lang.String('muscledrivensimwithEMG'); java.lang.String('inverseproblem')]];
    analyzeMetabolicCostWithEMG(solution);
    Issues = computeIDFromResult(Issues, solution);
    analyzeMetabolicCostWithEMG(solution);
    
end