function computeIDFromResultMuscle(solution)
    import org.opensim.modeling.*

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% compute ID moments 
    % solution = MocoTrajectory(solution);
    Time = solution.getTimeMat();
    numColPoints = solution.getNumTimes();

    workingdir = pwd;

    idtool = InverseDynamicsTool();
    modelid = Model('post_simple_model_all_the_probes.osim');
    idtool.setModel(modelid);

    % % get the coordinates and all the states
    % idtable = solution.exportToStatesTable();
    % idlabels = idtable.getColumnLabels();
    % idlabelssize = idlabels.size();
    % idlabelslist = [];
    % for l=0:idlabelssize-1
    %     % get the actual label names in a list
    %     actuallabel = idlabels.get(l);
    %     idlabelslist = [idlabelslist, actuallabel]
    % end

    % try storage instead
    idstorage = solution.exportToStatesStorage();

    % set the tool up
    idtool.setCoordinateValues(idstorage);
    idresult = 'muscle_joint_moment_breakdown.sto';
    idtool.setResultsDir(workingdir);
    idtool.setOutputGenForceFileName(idresult);

    % run the tool
    idtool.run();

    % load the net joint moments
    idresultfile = strcat(workingdir, '\', idresult);
    netjointmoments = TimeSeriesTable(idresultfile);

    % get the states
    statestraj = solution.exportToStatesTrajectory(modelid);


    %% get some details from the model
    % muscles
    muscleset = modelid.getMuscles();
    nummuscles = muscleset.getSize();
    tendonforces = zeros(length(Time), nummuscles);
    musclepaths = [];
    for m=0:nummuscles-1
        tempmusc = muscleset.get(m);
        tempmuscpath = tempmusc.getAbsolutePathString();
        musclepaths = [musclepaths, tempmuscpath];
        for t=1:length(Time)
            tempstate = statestraj.get(t-1);
            modelid.realizeDynamics(tempstate);
            tendonforces(t,m+1) = tempmusc.getTendonForce(tempstate);
        end
    end

    % coordinates
    coordinateset = modelid.getCoordinateSet();
    numcoordinates = coordinateset.getSize();
    coordinatepaths = [];
    for c=0:numcoordinates-1
        tempcoord = coordinateset.get(c);
        tempcoordpath = tempcoord.getAbsolutePathString();
        coordinatepaths = [coordinatepaths, tempcoordpath];
    end
    
    
    % forces in the model
    forceset = modelid.getForceSet();
    numforces = forceset.getSize();
    muscleforcepaths = [];
    coordactforcepaths = [];
    coordinatecatch = {'reserve','lumbar'};
    for f=0:numforces-1
        tempforce = forceset.get(f);
        tempforcepath = tempforce.getAbsolutePathString();
        if any(contains(char(tempforcepath),coordinatecatch))
            coordactforcepaths = [coordactforcepaths, tempforcepath];
        else
            muscleforcepaths = [muscleforcepaths, tempforcepath];
        end
    end

    
    
    keyboard
    % coordinate actuator moments
    numcoordact = length(coordactforcepaths);
    coordactmoments = zeros(length(Time), numcoordact);
    for c=0:numcoordact-1
        tempcoordact = CoordinateActuator().safeDownCast(modelid.getComponent(coordactforcepaths(c+1)));
        for t=0:length(Time)-1
            tempstate = statestraj.get(t);
            modelid.realizeDynamics(tempstate);
            coordactmoments(t+1,c+1) = tempcoordact.getActuation(tempstate);
        end
    end




    keyboard
    % figure out how to visualize and plot all the stuff that I want to see, 
    % and throw a warning if the reserves are certain % of net moments
    % TODO
    %% now going to actually compute and plot some stuff?

    % get the tendon forces



    figure(1);

end
