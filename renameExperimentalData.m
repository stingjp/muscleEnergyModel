function renameExperimentalData()
    
    % set up all the new files for the study - from exp/ik results/scaling
    destination = pwd;
    subject_file = 'subject.osim';
    mocap_file = 'motion_capture.trc';
    grf_file = 'ground_reaction.mot';
    ik_file = 'coordinates.mot';
    emg_file = 'electromyography.sto';

    cd expdata;
    expfiles = dir();
    L = size(expfiles);
    L = L(1);

    for i=1:L
        tempfile = expfiles(i).name;
        if contains(tempfile, 'osim')
            % copy the osim file
            full_dest = strcat(destination, strcat('\',subject_file));
            full_source = strcat(pwd, strcat('\',tempfile));
            copyfile(full_source, full_dest);        
        elseif contains(tempfile, 'trc')
            % copy the mocap file
            full_dest = strcat(destination, strcat('\', mocap_file));
            full_source = strcat(pwd, strcat('\', tempfile));
            copyfile(full_source, full_dest);
        elseif contains(tempfile, 'electromyography')
            if ~contains(tempfile, 'anc')
                % copy the EMG file
                full_dest = strcat(destination, strcat('\', emg_file));
                full_source = strcat(pwd, strcat('\', tempfile));
                copyfile(full_source, full_dest);
            end
        elseif contains(tempfile, 'ground_reaction') || contains(tempfile, 'grf_walk')
            if ~contains(tempfile, 'svg')
                % copy the grf file
                full_dest = strcat(destination, strcat('\', grf_file));
                full_source = strcat(pwd, strcat('\', tempfile));
                copyfile(full_source, full_dest);
            end
        elseif contains(tempfile, 'ik_solution')
            % copy ik file
            full_dest = strcat(destination, strcat('\', ik_file));
            full_source = strcat(pwd, strcat('\', tempfile));
            copyfile(full_source, full_dest);
        end
    end
    cd ..
    
    cd ik;
    ikfiles = dir();
    L = size(ikfiles);
    L = L(1);
    for i=1:L
        tempfile = ikfiles(i).name();
        if contains(tempfile, 'ik_solution')
            % copy ik file
            full_dest = strcat(destination, strcat('\', ik_file));
            full_source = strcat(pwd, strcat('\', tempfile));
            copyfile(full_source, full_dest);
        end
    end
    cd ..
end
