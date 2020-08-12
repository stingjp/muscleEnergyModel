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

    % Solve the problem and write the solution to a Storage file.
    % solution = inverse.solve(true); % to visualize
    solution = inverse.solve();
    solution.getMocoSolution().write('muscle_stateprescribe_grfprescribe_solution.sto');
    STOFileAdapter.write(solution.getMocoSolution().exportToControlsTable(), 'muscleprescribe_controls.sto')
    STOFileAdapter.write(solution.getMocoSolution().exportToStatesTable(), 'muscleprescribe_states.sto')
    
    
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