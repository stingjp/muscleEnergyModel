config = ReadYaml('config.yaml');
data_path = config.data_path;
repo_path = config.repo_path;
addpath(repo_path)

subjects = {'subject024', 'subject077', 'subject088', 'subject112', ...
            'subject126', 'subject127', 'subject128'};
filenames = {'static','stp000','stp025','stp050','stp075','stp100'};

import org.opensim.modeling.*
for s = 1:length(subjects)
    fprintf('%s: extracting files...\n', subjects{s})
    for f = 1:length(filenames)
       
        % Switch to "preprocess" directory and copy over c3d files
        preprocess_path = fullfile(data_path,'preprocess',subjects{s});
        cd(preprocess_path)
        c3d_path = fullfile(data_path, 'raw', 'c3d', subjects{s}, ...
            [filenames{f} '.c3d']); 
        copyfile(c3d_path, fullfile(preprocess_path, [filenames{f} '.c3d']))
        
        % Add virtual markers for hip joint center calculation
        % TODO
        
%         % Get C3D file object
%         c3d = osimC3D(c3d_path);
%         
%         % Rotate data
%         c3d.rotateData('x', 90)
%         c3d.rotateData('y', 90);
%         
%         % Get marker data and write TRC file
%         markers = c3d.getTable_markers();
%         trc_adapter = TRCFileAdapter();
%         trc_adapter.write(markers,[filenames{f} '.trc']);
%         
%         % Get force data and write to MOT file
%         keyboard
%         forces = c3d.getTable_forces();
%         postfix = StdVectorString();
%         postfix.add('_x');
%         postfix.add('_y');
%         postfix.add('_z');
%         forces_flattened = forces.flatten(postfix);
%         sto_adapter = STOFileAdapter();
%         sto_adapter.write(forces_flattened,[filenames{f} '.mot']);
        
        data = btk_loadc3d([filenames{f} '.c3d']);
        
        % Rotate marker and GRF data to match default OpenSim ground frame
        rotation = [{'x' 90 'y' 90}];
        
        markers_rotated = rotateCoordinateSys(data.marker_data.Markers,rotation);
        data.marker_data.Markers = markers_rotated;
        
        grf_names = {'P','F','M'};
        for j = 1:length(grf_names)
            for k = 1:2
                grf_rotated = rotateCoordinateSys(data.fp_data.GRF_data(k).(grf_names{j}),rotation);
                if strcmp('P',grf_names{j})
                    data.fp_data.GRF_data(k).(grf_names{j}) = grf_rotated/1000;
                else
                    data.fp_data.GRF_data(k).(grf_names{j}) = grf_rotated;
                end
            end
        end
        
        printMOT(data)
        printTRC(data.marker_data.Markers,data.marker_data.Info.frequency,data.marker_data.Filename)
                
       fprintf('--> %s', filenames{f})
       fprintf('\n')
    end
    fprintf('\n\n')
    fprintf('%s: estimating hip joint centers...\n', subjects{s})
    
    % Estimate hip joint centers
    static_cal_file = fullfile(preprocess_path, 'static.trc');
    right_leg_cal_file = fullfile(preprocess_path, 'stp000.trc');
    left_leg_cal_file = fullfile(preprocess_path, 'stp000.trc');
    to_add_HJC_files = {'static.trc'};
    static_ref_dir = preprocess_path;
    to_add_HJC_dir = preprocess_path;
    write_dir = preprocess_path;
    use_same_pelvic_marker_set = true;
    
    estimateHJC(static_cal_file , right_leg_cal_file, ...
        left_leg_cal_file, to_add_HJC_files, ...
        static_ref_dir, to_add_HJC_dir, write_dir, ...
        use_same_pelvic_marker_set)
    
    fprintf('\n\n')
end