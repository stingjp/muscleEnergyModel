function printTRC(mkrData,sampleRate,trialName)
% printTRC() Printes a structure of Mkr Data to TRC format
%   mkrData, 
%   sampleRate,
%   trialName,
%   dataPath

%   Export mkr data Into OpenSim ready format 
%   Author: J.J. Dunne, Thor Besier, C.J. Donnelly, S. Hamner.  

%% Pre- allocate some arrays and set some index values
mkNameArray  = [];
fields       = fieldnames(mkrData);
nMarkers     = length(fields);
markerArray  = [];
[pathstr, name, ext] = fileparts(trialName);

% Dump out all the mkr names and Data 
for i = 1 : nMarkers
        % name array
        mkNameArray = [mkNameArray {char(fields(i))}];
        % data array
        eval(['markerArray(:,(3*i)-2:(3*i)) = mkrData.' char(fields(i)) ';' ] )
end

[nFrames n] =   size(markerArray);

% Create Arrays for Frame and Time 
frameArray  =   [0:nFrames-1]';                % Create frame Number array
timeArray   =   (frameArray)/sampleRate;       % Create time array
    
% create the final array of data for printing 
trcData     =   [frameArray timeArray markerArray];


%% Print Data 

% Output the data to a csv file
    trcFileName     = fullfile(pathstr, [name '.trc']);
    fid             = fopen(trcFileName,'w');
    
    display('   Printing marker trajectory file');
    
    % Print header information
    fprintf(fid,'PathFileType\t4\t(X/Y/Z)\t%s\n',name);
    fprintf(fid,'dataRate\tCameraRate\tNumFrames\tNumMarkers\tUnits\tOrigdataRate\tOrigdataStartFrame\tOrigNumFrames\n');
    fprintf(fid,'%d\t%d\t%d\t%d\tmm\t%d\t%d\t%d\n',sampleRate,sampleRate,nFrames,nMarkers,sampleRate,0,n);
    fprintf(fid,'Frame#\tTime');
    
    % create a printing array for the headers
    fmt1=['\t%s\t\t'];
    % Print the headers to file
    for i = 1:nMarkers
        fprintf(fid,fmt1,char(mkNameArray(i)));% print a row of marker names
    end
        fprintf(fid,'\n\t\t');                 % \new line
    for i = 1:nMarkers
        fprintf(fid,'X%d\tY%d\tZ%d\t',i,i,i);  % print a row of XYZ's
    end
    fprintf(fid,'\n');                         % \new line
    
    % Print the data
    fmt3=[repmat(' %2.6g\t',1,3*nMarkers+2) '\n']; 
    for i   = 1:nFrames
            fprintf(fid,fmt3,trcData(i,:));    % print row of Mkr Data
    end
            
    fclose(fid);
  
end

