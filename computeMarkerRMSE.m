function [marknames, markvalus, markstruct] = computeMarkerRMSE(solution, trackorprescribe, conornot)
    import org.opensim.modeling.*
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    outpaths = {'/markerset.*location'};
    pathsstd = StdVectorString();
    for p=1:length(outpaths)
        pathsstd.add(java.lang.String(outpaths(p)));
    end
    
    testmodel = Model('post_simple_model_all_the_probes_muscletrack.osim');
    if conornot
        teststates = TimeSeriesTable('muscletrack_states_100con.sto');
        testcontrols = TimeSeriesTable('muscletrack_controls_100con.sto');
    else
        teststates = TimeSeriesTable('muscletrack_states.sto');
        testcontrols = TimeSeriesTable('muscletrack_controls.sto');
    end

    simMarkers = opensimSimulation.analyzeVec3(testmodel, teststates, testcontrols, pathsstd);
    filepath = strcat(pwd, '\motion_capture.trc');
    expMarkers = TimeSeriesTableVec3(java.lang.String(filepath));

    % get the times while we can
    simtime = double(simMarkers.getIndependentColumn().toArray());
    exptime = double(expMarkers.getIndependentColumn().toArray());
    % need to load the start and end times of the gait cycle
    starttime = simtime(1);
    endtime = simtime(end);

    % get exp time start index and end index
    for i=1:length(exptime)
        temptime = exptime(i);
        if temptime < starttime
            startid = i;
        elseif temptime < endtime
            endid = i+1;
        end
    end
    
    % set up the data structures
    marknames = [];
    markvalus = [];
    markstruct = struct;

    expnames = expMarkers.getColumnLabels();
    simnames = simMarkers.getColumnLabels();
    % now need to get names, then loop through 
    for m=0:simnames.size()-1
        simm1 = simnames.get(m);
        temp = char(simm1);
        simm = temp(12:end-9);
        
        if expnames.contains(simm)
            % we have a marker in both
            % get both columns, and shorten to the right time as the sim
            % also resample to make sure the same length
            expcol = expMarkers.getDependentColumn(simm).getAsMat();
            expcolscale = expcol.*0.001;
            simcol = simMarkers.getDependentColumn(simm1).getAsMat();
            
            % need to resample the vector
            test1 = interp1(exptime(startid:endid),expcolscale(startid:endid,1),simtime);
            test2 = interp1(exptime(startid:endid),expcolscale(startid:endid,2),simtime);
            test3 = interp1(exptime(startid:endid),expcolscale(startid:endid,3),simtime);
            
            % compute marker errors
            dist = sqrt(((test1-simcol(:,1)).^2)+((test2-simcol(:,2)).^2)+((test3-simcol(:,3)).^2));
            % so the dist is the error in m
            
            % get the RMS error
            sqerr = dist.^2;
            msqerr = mean(sqerr);
            rmse = sqrt(msqerr); 
            
            % add to vectors
            % this set of variables is for translational errors. 
            marknames = [marknames; string(simm)];
            markvalus = [markvalus; rmse];
            markstruct.(genvarname(string(simm))) = rmse;
        end
    end
end
