function renameExperimentalData()
    
    % set up all the new files for the study - from exp/ik results/scaling
    destination = pwd;
    subject_file = 'subject.osim';
    old_subject_file = 'subject_old.osim';
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
            if contains(tempfile, 'old.osim')
                % copy over the old one
                full_dest = strcat(destination, strcat('\',old_subject_file));
                full_source = strcat(pwd, strcat('\',tempfile));
                copyfile(full_source, full_dest);
            else
                % copy the osim file
                full_dest = strcat(destination, strcat('\',subject_file));
                full_source = strcat(pwd, strcat('\',tempfile));
                copyfile(full_source, full_dest);
            end
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
        elseif contains(tempfile, 'ground_reaction') || contains(tempfile, 'grf_walk') || contains(tempfile, 'GRFs')
            if ~contains(tempfile, 'svg')
                % TODO: edit the grf file so that it has the new name as the source in line 1
                % we actually might be okay here...

                % copy the grf file
                full_dest = strcat(destination, strcat('\', grf_file));
                full_source = strcat(pwd, strcat('\', tempfile));
                copyfile(full_source, full_dest);
            end
        elseif contains(tempfile, 'ik_solution') || contains(tempfile, 'results_ik') || contains(tempfile, 'results_IK')
            % copy ik file
            full_dest = strcat(destination, strcat('\', ik_file));
            full_source = strcat(pwd, strcat('\', tempfile));
            copyfile(full_source, full_dest);
        end
    end
    cd ..
    try
        cd ik;
        ikfiles = dir();
        L = size(ikfiles);
        L = L(1);
        for i=1:L
            tempfile = ikfiles(i).name();
            if contains(tempfile, 'ik_solution') || contains(tempfile, 'results_IK') || contains(tempfile, 'results_ik')
                % copy ik file
                full_dest = strcat(destination, strcat('\', ik_file));
                full_source = strcat(pwd, strcat('\', tempfile));
                copyfile(full_source, full_dest);
            end
        end
        cd ..
    catch
        % there is no IK folder
    end 
end
