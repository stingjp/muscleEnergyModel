function [Issues, maxreservepercvalus, avgreservepercvalus, maxreservevalus, avgreservevalus, reservenames, ...
    maxresidualpercvalus, avgresidualpercvalus, maxresidualvalus, avgresidualvalus, residualnames, ...
    maxresidualmompercvalus, avgresidualmompercvalus, maxresidualmomvalus, avgresidualmomvalus, residualmomnames ...
    ] = computeIDFromResult(Issues, solution, tag)

    import org.opensim.modeling.*
%     Issues = [[java.lang.String('coordinate actuator'); java.lang.String('ratio to net')]];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    keyboard
    %% compute ID moments 
    % solution = MocoTrajectory(solution);
    Time = solution.getTimeMat();
    numColPoints = solution.getNumTimes();

    workingdir = pwd;

    idtool = InverseDynamicsTool(java.lang.String('idguitesting.xml'));
    if strcmp(tag, 'muscletrack')
        modelid1 = Model('post_simple_model_all_the_probes_muscletrack.osim');
        modelid2 = Model('simple_model_all_the_probes_adjusted.osim');
    elseif strcmp(tag, 'muscleprescribe')
        modelid1 = Model('post_simple_model_all_the_probes_muscleprescribe.osim');
        modelid2 = Model('simple_model_all_the_probes.osim');
    elseif strcmp(tag, 'muscletrack_redo')
        modelid1 = Model('post_simple_model_all_the_probes_muscletrack_redo.osim');
        modelid2 = Model('simple_model_all_the_probes.osim');
    elseif strcmp(tag, 'muscletrack_paths_redo')
        modelid1 = Model('post_simple_model_all_the_probes_muscletrack_paths_redo.osim');
        modelid2 = Model('simple_model_all_the_probes.osim');
    end
    % modelid1 = Model('post_simple_model_all_the_probes_muscletrack.osim');
    modelid1.initSystem();
    % modelid2 = Model('simple_model_all_the_probes_adjusted.osim');
    idtool.setModel(modelid2);

    if strcmp(tag, 'grftrack')
        idtool.setExternalLoadsFileName('grftrack_grf_walk.xml')
    else
        idtool.setExternalLoadsFileName('grf_walk.xml')
    end


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
    sto = Storage();
    solutionstatestable = solution.exportToStatesTable();
    labels = solutionstatestable.getColumnLabels();
    numlabels = labels.size();
    properlabels = org.opensim.modeling.ArrayStr();
    properlabels.set(0,"time");
    for i=0:numlabels-1
        templabel = labels.get(i);
        if contains(char(templabel),'jointset') && contains(char(templabel),'value')
            % we want this one
            templabel2 = char(templabel)
            tempsplit = split(templabel2,'/')
            templabel3 = tempsplit(4)
            
            properlabels.set(i+1, templabel3)
        else
            % get rid of this column
            solutionstatestable.removeColumn(templabel);
        % properlabels.set(i+1,labels.get(i));
        end
    end
    sto.setColumnLabels(properlabels);
    
    statetime = solutionstatestable.getIndependentColumn();
    try
        starttime = statetime.get(0).doubleValue();
        endtime = statetime.get(statetime.size()-1).doubleValue();
    catch
        starttime = statetime.get(0);
        endtime = statetime.get(statetime.size()-1);
    end
    timelength = statetime.size();
    
    for i=0:timelength-1
        temprow = solutionstatestable.getRowAtIndex(i).getAsMat();
        temprow2 = org.opensim.modeling.Vector().createFromMat(temprow);
        try
            sto.append(statetime.get(i).doubleValue(), temprow2);
        catch
            sto.append(statetime.get(i), temprow2);
        end
    end
    keyboard
    % save some files. 
    if strcmp(tag, 'muscletrack')
        sto.print('muscletrack_coordinates_short.sto');
        idtool.setCoordinatesFileName('muscletrack_coordinates_short.sto');
    elseif strcmp(tag, 'muscleprescribe')
        sto.print('muscleprescribe_coordinates_short.sto');
        idtool.setCoordinatesFileName('muscleprescribe_coordinates_short.sto');
    elseif strcmp(tag, 'muscletrack_redo')
        sto.print('muscletrack_redo_coordinates_short.sto');
        idtool.setCoordinatesFileName('muscletrack_redo_coordinates_short.sto');
    elseif strcmp(tag, 'muscletrack_paths_redo')
        sto.print('muscletrack_paths_redo_coordinates_short.sto');
        idtool.setCoordinatesFileName('muscletrack_paths_redo_coordinates_short.sto');
    end
    
    
    % % idstorage = solution.exportToStatesStorage();
    % % set the tool up
    % % idtool.setCoordinateValues(idstorage);
    % sto.print('muscle_coordinates_short.sto');
    keyboard
    % % idtool.setCoordinateValues(sto);
    % idtool.setCoordinatesFileName('muscle_coordinates_short.sto');
    idresult = strcat(tag, '_muscle_joint_moment_breakdown.sto');
    idtool.setResultsDir(workingdir);
    idtool.setOutputGenForceFileName(idresult);
    idtool.setEndTime(endtime);
    idtool.setStartTime(starttime);
    idtool.set_results_directory(java.lang.String(workingdir));
    
    % run the tool
    idtool.print(strcat(tag,'idtesting.xml'));
    idtool.run();

    % load the net joint moments
    idresultfile = strcat(workingdir, '\', idresult);
    
    
    netjointmoments = TimeSeriesTable(idresultfile);
    
    % get the states
    statestraj = solution.exportToStatesTrajectory(modelid1);

    %% get some details from the model
    % muscles
    muscleset = modelid1.getMuscles();
    nummuscles = muscleset.getSize();
    tendonforces = zeros(length(Time), nummuscles);
    musclepaths = [];
    for m=0:nummuscles-1
        tempmusc = muscleset.get(m);
        tempmuscpath = tempmusc.getAbsolutePathString();
        musclepaths = [musclepaths, tempmuscpath];
        for t=1:length(Time)
            tempstate = statestraj.get(t-1);
            modelid1.realizeDynamics(tempstate);
            tendonforces(t,m+1) = tempmusc.getTendonForce(tempstate);
        end
    end

    % coordinates
    coordinateset = modelid1.getCoordinateSet();
    numcoordinates = coordinateset.getSize();
    coordinatepaths = [];
    for c=0:numcoordinates-1
        tempcoord = coordinateset.get(c);
        tempcoordpath = tempcoord.getAbsolutePathString();
        coordinatepaths = [coordinatepaths, tempcoordpath];
    end
    
    
    % forces in the model
    forceset = modelid1.getForceSet();
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
        tempcoordact = CoordinateActuator().safeDownCast(modelid1.getComponent(coordactforcepaths(c+1)));
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

    keyboard
    % get moment arms, and plot stuff
    % aim for 5% of net joint moment or less for reserves
    % lots for residuals
    % get grf for residual comparisons
    table_grf = TimeSeriesTable(strcat('analyzemuscles',tag,'_ForceReporter_forces.sto'));
    grf_r_Fx = table_grf.getDependentColumn('calcn_r_Right_GRF_Fx').getAsMat();
    grf_r_Fy = table_grf.getDependentColumn('calcn_r_Right_GRF_Fy').getAsMat();
    grf_r_Fz = table_grf.getDependentColumn('calcn_r_Right_GRF_Fz').getAsMat();
    grf_l_Fx = table_grf.getDependentColumn('calcn_l_Left_GRF_Fx').getAsMat();
    grf_l_Fy = table_grf.getDependentColumn('calcn_l_Left_GRF_Fy').getAsMat();
    grf_l_Fz = table_grf.getDependentColumn('calcn_l_Left_GRF_Fz').getAsMat();
    grf_r_Tx = table_grf.getDependentColumn('calcn_r_Right_GRF_Tx').getAsMat();
    grf_r_Ty = table_grf.getDependentColumn('calcn_r_Right_GRF_Ty').getAsMat();
    grf_r_Tz = table_grf.getDependentColumn('calcn_r_Right_GRF_Tz').getAsMat();
    grf_l_Tx = table_grf.getDependentColumn('calcn_l_Left_GRF_Tx').getAsMat();
    grf_l_Ty = table_grf.getDependentColumn('calcn_l_Left_GRF_Ty').getAsMat();
    grf_l_Tz = table_grf.getDependentColumn('calcn_l_Left_GRF_Tz').getAsMat();

    table_com = TimeSeriesTable(strcat('analyzemuscles',tag,'_BodyKinematics_pos_global.sto'));
    com_x = table_com.getDependentColumn('center_of_mass_X').getAsMat();
    com_y = table_com.getDependentColumn('center_of_mass_Y').getAsMat();
    com_z = table_com.getDependentColumn('center_of_mass_Z').getAsMat();
    
    % check the coordactlabels
    % check tempcoordinate name
    
    % set up data structs for reserve or residual error
    % reserve errors
    maxreservevalus = [];
    maxreservepercvalus = [];
    avgreservevalus = [];
    avgreservepercvalus = [];
    reservenames = [];
    %residual forces
    maxresidualvalus = [];
    maxresidualpercvalus = [];
    avgresidualvalus = [];
    avgresidualpercvalus = [];
    residualnames = [];
    % residual moments
    maxresidualmomvalus = [];
    maxresidualmompercvalus = [];
    avgresidualmomvalus = [];
    avgresidualmompercvalus = [];
    residualmomnames = [];

    % TODO start here and fix. 
    forcecheck = {'tx','ty','tz'};
    templabels = coordactmomentstable.getColumnLabels();
    for i=1:length(coordinatepaths);
        tempcoordinate = coordinatepaths(i);
        coord = Coordinate.safeDownCast(modelid1.getComponent(tempcoordinate));
        
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

                if contains(char(tempmomentname),'pelvis')
                    % do the residual force stuff
                    tempnetexternal = sqrt(grf_r_Fx.^2 + grf_r_Fy.^2 + grf_r_Fz.^2);
                    tempnetexternal_peak = max(abs(tempnetexternal));
                    tempnetexternal_rms = rms(tempnetexternal);

                    tempind = coordactmomentstable.getDependentColumn(tempcoordact).getAsMat();
                    tempind_peak = max(abs(tempind));
                    tempind_rms = rms(tempind);

                    % ratio_peak = tempind_peak/tempnetexternal_peak;
                    % ratio_rms = tempind_rms/tempnetexternal_rms;
                    ratio_peak = max(abs(tempind)./tempnetexternal_peak);
                    ratio_rms = rms(abs(tempind)./tempnetexternalcom_peak);
                    
                    % residual forces collection
                    maxresidualvalus = [maxresidualvalus; max(abs(tempind))];
                    avgresidualvalus = [avgresidualvalus; rms(abs(tempind))];
                    maxresidualpercvalus = [maxresidualpercvalus; ratio_peak];
                    avgresidualpercvalus = [avgresidualpercvalus; ratio_rms];
                    residualnames = [residualnames; string(tempmomentname)];

                    if ratio_peak > 0.05
                        Issues = [Issues; [java.lang.String(tempcoordact), java.lang.String(strcat('peak ratio:',string(ratio_peak)))]];
                    end
                    if ratio_rms > 0.05
                       Issues = [Issues; [java.lang.String(tempcoordact), java.lang.String(strcat('rms ratio:',string(ratio_rms)))]];
                    end 
                
                else
                    % do the reserve force stuff
                    disp('need to figure this out');
                end
            else
                tempmomentname = strcat(tempmomentname,'_moment');
            
                if contains(char(tempmomentname),'pelvis')
                    % do the residual moment stuff
                    %residual moments that are less than 1% of COM height times the magni- tude of the measured net external force
                    tempnetexternal = sqrt(grf_r_Fx.^2 + grf_r_Fy.^2 + grf_r_Fz.^2);
                    tempnetexternalcom = com_y .* tempnetexternal;
                    tempnetexternalcom_peak = max(abs(tempnetexternalcom));
                    tempnetexternalcom_rms = rms(tempnetexternalcom);

                    tempind = coordactmomentstable.getDependentColumn(tempcoordact).getAsMat();
                    tempind_peak = max(abs(tempind));
                    tempind_rms = rms(tempind);
 
                    % ratio_peak = tempind_peak/tempnetexternalcom_peak;
                    % ratio_rms = tempind_rms/tempnetexternalcom_rms;
                    ratio_peak = max(abs(tempind)./tempnetexternalcom_peak);
                    ratio_rms = rms(abs(tempind)./tempnetexternalcom_peak);

                    % residual moment collection
                    maxresidualmomvalus = [maxresidualmomvalus; max(abs(tempind))];
                    maxresidualmompercvalus = [maxresidualmompercvalus; ratio_peak];
                    avgresidualmomvalus = [avgresidualmomvalus; rms(abs(tempind))];
                    avgresidualmompercvalus = [avgresidualmompercvalus; ratio_rms];
                    residualmomnames = [residualmomnames; string(tempmomentname)];

                    if ratio_peak > 0.01
                        Issues = [Issues; [java.lang.String(tempcoordact), java.lang.String(strcat('peak ratio:',string(ratio_peak)))]];
                    end
                    if ratio_rms > 0.01
                       Issues = [Issues; [java.lang.String(tempcoordact), java.lang.String(strcat('rms ratio:',string(ratio_rms)))]];
                    end 

                else
                    % do the reserve moment stuff
                    tempnet = netjointmoments.getDependentColumn(tempmomentname).getAsMat();
                    tempnet_peak = max(abs(tempnet));
                    tempnet_rms = rms(tempnet);
                    
                    tempind = coordactmomentstable.getDependentColumn(tempcoordact).getAsMat();
                    tempind_peak = max(abs(tempind));
                    tempind_rms = rms(tempind);

                    % ratio_peak = tempind_peak/tempnet_peak;
                    % ratio_rms = tempind_rms/tempnet_rms;
                    ratio_peak = max(abs(tempind)./tempnet_peak);
                    ratio_rms = rms(abs(tempind)./tempnet_peak);

                    % reserve value collection
                    maxreservevalus = [maxreservevalus; max(abs(tempind))];
                    maxreservepercvalus = [maxreservepercvalus; ratio_peak];
                    avgreservevalus = [avgreservevalus; rms(abs(tempind))];
                    avgreservepercvalus = [avgreservepercvalus; ratio_rms];
                    reservenames = [reservenames; string(tempmomentname)];

                    if ratio_peak > 0.05
                        Issues = [Issues; [java.lang.String(tempcoordact), java.lang.String(strcat('peak ratio:',string(ratio_peak)))]];
                    end
                    if ratio_rms > 0.05
                       Issues = [Issues; [java.lang.String(tempcoordact), java.lang.String(strcat('rms ratio:',string(ratio_rms)))]];
                    end
                end
            end
            

            % % loop through time to find the largest instantaneous ratio
            % instmaxratio = 0;
            % for n=1:length(Time)
            %     instnet = tempnet(n);
            %     instind = tempind(n);
            %     instratio = instind/instnet;
            %     if instratio > instmaxratio
            %         instmaxratio = instratio;
            %     end
            % end
            % if instmaxratio > .05
            %     Issues = [Issues; [java.lang.String(tempcoordact), java.lang.String(strcat('instant ratio:',string(instmaxratio),'reserve:',string(instind)))]];
            % end
            % get the peak of each and compare:
            % maybe talk to Scott or Jen about this? Nick?
        end
    end
        
        
    % now need moment arms I think and then we should be good
    momentarms = zeros(length(Time), nummuscles);
    for m=1:length(muscleforcepaths)
        if ~any(contains(char(muscleforcepaths(m)),'exotendon')) && ~any(contains(char(muscleforcepaths(m)),'HOBL'))
            muscle = Muscle.safeDownCast(modelid1.getComponent(muscleforcepaths(m)));
            
            for t=1:length(Time)
                state = statestraj.get(t-1);
                momentarms(t,m) = muscle.computeMomentArm(state, coord);
            end
        end
    end

    %% TODO
    % here is where I should do any plotting of specific things 
    % that I want to plot for moments


    


%     figure(1);

end
