function [Issues] = computeIDFromResultMuscle(Issues, solution)
    import org.opensim.modeling.*
%     Issues = [[java.lang.String('coordinate actuator'); java.lang.String('ratio to net')]];
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

    
    % get the control actuator moments
    controlstable = solution.exportToControlsTable();
    numcoordact = length(coordactforcepaths);
    tempmoments = zeros(length(Time), numcoordact);
    coordactmoments = struct();
    coordactmoments.time = Time;
    for c=0:numcoordact-1
        tempcoordact = CoordinateActuator().safeDownCast(modelid.getComponent(coordactforcepaths(c+1)));
        tempcoordact_forcemultiply = tempcoordact.getOptimalForce();
        tempcontrols = controlstable.getDependentColumn(coordactforcepaths(c+1));
        for t=0:length(Time)-1
            tempcontrol = tempcontrols.get(t);
            tempmoments(t+1,c+1) = tempcoordact_forcemultiply*tempcontrol;
        end
        tempcoordactname = coordactforcepaths(c+1);
        tempcoordactname = tempcoordactname.split('/');
        tempcoordactname = char(tempcoordactname(end));
        coordactmoments.(genvarname(tempcoordactname)) = tempmoments(:,c+1);
    end
    coordactmomentstable = osimTableFromStruct(coordactmoments);

    
    % get moment arms, and plot stuff
    % aim for 5% of net joint moment or less
    forcecheck = {'tx','ty','tz'};
    templabels = coordactmomentstable.getColumnLabels();
    for i=1:length(coordinatepaths);
        tempcoordinate = coordinatepaths(i);
        coord = Coordinate.safeDownCast(modelid.getComponent(tempcoordinate));
        % skip the patella coordinate, loop others
        if ~any(contains(char(tempcoordinate),'beta'))
            % compare the reserve for each coordinate to the net joint moment
            tempmomentname = tempcoordinate.split('/');
            tempmomentname = char(tempmomentname(end));
            % figure out how to get the reserve ones here
            matched = false;
            c = 0;
            while ~matched
                tempcoordact = templabels.get(c);
                if any(contains(char(tempcoordact), tempmomentname)) || any(contains(tempmomentname, char(tempcoordact))) 
                    matched = true;
                else
                    c = c+1;
                end
            end
            % need condition to get the full net moment or force name
            if any(contains(char(tempcoordinate),forcecheck))
                % this is a translational force
                tempmomentname = strcat(tempmomentname,'_force');
            else
                tempmomentname = strcat(tempmomentname,'_moment');
            end
            tempnet = netjointmoments.getDependentColumn(tempmomentname);
            tempnet = tempnet.getAsMat();
            tempind = coordactmomentstable.getDependentColumn(tempcoordact);
            tempind = tempind.getAsMat();
            % get the reserve moment
            % do a validation check on the reserve vs net
            tempnetmean = mean(tempnet);
            tempindmean = mean(tempind);
            ratio = tempindmean/tempnetmean;
            if ratio > .05
                Issues = [Issues; [tempcoordact, java.lang.String(string(ratio))]];
            end
            % loop through time to find the largest instantaneous ratio
            instmaxratio = 0;
            for n=1:length(Time)
                instnet = tempnet(n);
                instind = tempind(n);
                instratio = instind/instnet;
                if instratio > instmaxratio
                    instmaxratio = instratio;
                end
            end
            if instmaxratio > .05
                Issues = [Issues; [tempcoordact, java.lang.String(strcat('instant ratio:',string(instmaxratio),'reserve:',string(instind)))]];
            end
            % get the peak of each and compare:
            % maybe talk to Scott or Jen about this? Nick?
        end
        
        
        % now need moment arms I think and then we should be good
        momentarms = zeros(length(Time), nummuscles);
        for m=1:length(muscleforcepaths)
            if ~any(contains(char(muscleforcepaths(m)),'exotendon'))
                muscle = Muscle.safeDownCast(modelid.getComponent(muscleforcepaths(m)));
                
                for t=1:length(Time)
                    state = statestraj.get(t-1);
                    momentarms(t,m) = muscle.computeMomentArm(state, coord);
                end
            end
        end

        %% TODO
        % here is where I should do any plotting of specific things 
        % that I want to plot for moments


    end
    


%     figure(1);

end
