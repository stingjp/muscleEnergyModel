function printMOT(structData)
% printMOT() Printes a structure of forces to .mot format

%   Export force data Into OpenSim ready format 
%   Author: J.J. Dunne, Thor Besier, C.J. Donnelly, S. Hamner. 

%% Pre- allocate some arrays and set some index values
nFrames      = length(structData.fp_data.GRF_data(1).F);
nForces      = length(structData.fp_data.GRF_data);
forceArray   = [];
torqueArray  = [];
[pathstr, name, ext] = fileparts(structData.marker_data.Filename);
analogRate   = structData.fp_data.Info(1).frequency;


%% Dump out all the force  Data
    for i = 1:nForces
          length( structData.fp_data.GRF_data(i).F);
          forceArray = [forceArray structData.fp_data.GRF_data(i).F];
          forceArray = [forceArray structData.fp_data.GRF_data(i).P];
          
    end

    for i = 1:nForces
          torqueArray = [torqueArray structData.fp_data.GRF_data(i).M];
    end
% Create Arrays for Frame and Time 
    frameArray  =   [0:nFrames-1]';                % Create frame Number array
    timeArray   =   (frameArray)/analogRate;       % Create time array

% Append to GRF Data
   forceData=[timeArray forceArray torqueArray];
    
%% create the headers for the file   
   Cord = {'x' 'y' 'z'};
   bodyForceHeader =[];
   bodyTorqueHeader=[];
   
    for i = 1:nForces
        ForceHeader = [];
        PointHeader = [];
        TorqueHeader= [];
        
        for u = 1:3
            forceHeader = [num2str(i) '_ground_force_v' char(Cord(u))];
            ForceHeader = [ForceHeader {forceHeader}];
            pointHeader = [num2str(i) '_ground_force_p' char(Cord(u))];
            PointHeader = [PointHeader {pointHeader}];
        end
            bodyForceHeader = [bodyForceHeader ForceHeader PointHeader];
        
        for u = 1:3
            torqueHeader= [num2str(i) '_ground_torque_' char(Cord(u))];
            TorqueHeader= [TorqueHeader {torqueHeader}];
        end
            bodyTorqueHeader = [bodyTorqueHeader TorqueHeader];
        
    end
       headers = [{'time'} bodyForceHeader bodyTorqueHeader];
       nHeaders = length(headers);
%% Print Data 

% Find the size of the Matrix    
    [m n]            =   size(forceData);    
% Create a Print Path    
    motFileName      = fullfile(pathstr, [name '.mot']);
    fid              = fopen(motFileName,'w');
    
    display('   Printing marker external loads file');
     
    % Print trial header
    fprintf(fid,'name %s\n',name);
    fprintf(fid,'datacolumns %d\n',n);
    fprintf(fid,'datarows %d\n',m);
    fprintf(fid,'range %f %f\n',timeArray(1),timeArray(end));
    fprintf(fid,'endheader\n');
    
        % Print array headers
    fmt1 = ['%s\t'];
    for i = 1 : nHeaders
        fprintf(fid,fmt1,char(headers(i)));
    end
    fprintf(fid,'\n');
    
    % Print data
    fmt3 = [repmat(' %2.6g\t',1,9*nForces+1) '\n']; 
    for i   = 1:nFrames
            fprintf(fid,fmt3,forceData(i,:));    % print a row's of Mkr Data
    end
    fclose(fid);
end

