function data = btk_loadc3d(file, grw_threshold)

% function [marker_data, analog_data, fp_data, sub_info] = btk_loadc3d(file)
%
% Function to load the data from a c3d file into the structured array data.
% The file may be excluded if you wish to choose it from a windows dialog
% box
%
% INPUT -
%           file - the file path that you wish to load (leave blank to
%               choose from a dialog box)
%           threshold - the force threshold for calculating the GRW
%           (default 5 N)
%
% OUTPUT - all structured arrays containing the following data
%           marker_data - any calculated data from the reconstructed C3D
%               file including marker trajectories  and any calculated angles, 
%               moments, powers or GRF data
%           analog_data - analog data (often sampled at a higher rate) including
%               force plate data and EMG data that might be collected.
%           grw_data - structure with the position magnitude of ground 
%               reaction force vector and moments relative to the global
%               cooridinate system
%           fp_info - structure with the force outputs from the force 
%               plates including the ground reaction wrench (force vector 
%               calculated in the global axis frame) and relevant 
%               sampling and forceplate position information             
%           sub_info - extra data from the C3D file if it exists, inlcuding
%               height and weight etc.
%
%
% e.g. [marker_data] = btk_loadc3d('c:\data\Session_1\trial001.c3d') will only load the
% marker data from the C3D file into a structured array called 'marker_data'
%
% [marker_data, a_data] = btk_loadc3d() will open a windows dialog box allowing
% the user to find the file they want to open and load the reconstructed
% data 
%
% Author: Glen Lichtwark (The University of Queensland)
% Updated: 06/12/2011

%
% This function requires the installation of the BTK (Biomechanical
% Toolkit) Matlab wrapper which must be added to the Matlab path. Please
% ensure you install the correct wrapper for your operating platform (ie.
% Windows, Mac etc) - please see http://code.google.com/p/b-tk/ for more
% information and cite appropriately - 
% further info at http://b-tk.googlecode.com/svn/doc/Matlab/0.1/index.html

warning off

% if no file is given then prompt for file
if nargin < 1
    [filein, pathname] = uigetfile({'*.c3d','C3D file'}, 'C3D data file...');

    if isequal(filein,0) || isequal(pathname,0)
        disp('No motion file loaded');
       return
    end
    
    file = [pathname filein];

end

if nargin < 2
    grw_threshold = 5;
end

% if file is c3d then load using the BTK matlab wrapper
acq = btkReadAcquisition(file);
% save filename to structure
marker_data.Filename = file;

% get first and last frame data
marker_data.First_Frame = btkGetFirstFrame(acq);
marker_data.Last_Frame = btkGetLastFrame(acq);

% get marker data
[markers, markersInfo] = btkGetMarkers(acq);
marker_data.Markers = markers;
fnames = fieldnames(markers);

% Get the frames where the most marker data is available.
% first,create an empty array the same length as the total number of frames
countArray = zeros(length( marker_data.Markers.(fnames{1})(:,1)),1);
% Loop through each marker, adding a 1 to each frame when it is available
% at that frame. Thereby counting how markers are availble at that frame
% number and getting a type of histogram. 
for i = 1 : length(fnames)
    % find when the marker has a value 
    mkrFrames = find(marker_data.Markers.(fnames{i})(:,1) ~= 0);
    % for each frame that marker has a value, add a count. 
    countArray(mkrFrames) = countArray(mkrFrames) + 1 ;
end

% find the frame range where 75% of the markers are available and store the
% first and last value of those values
VldFrames(1) = min(find(countArray > 0.75 * max(countArray)));
VldFrames(2) = max(find(countArray > 0.75 * max(countArray)));

% convert data to millimeters if in meters - catch if units have been screwed with
if mean(markers.(fnames{1})(:,1))<2 && mean(markers.(fnames{1})(:,2))<2 && mean(markers.(fnames{1})(:,3))<2
   m_mm = 1;
else m_mm = 0;
end

if strcmp(markersInfo.units.ALLMARKERS,'m') || m_mm == 1  
    for i = 1:length(fnames)
        marker_data.Markers.(fnames{i}) = marker_data.Markers.(fnames{i}) * 1000;
    end
    markersInfo.units.ALLMARKERS = 'mm';
end

% add the marker data information to the structure
marker_data.Info = markersInfo;
marker_data.Info.NumFrames = btkGetPointFrameNumber(acq);
marker_data.Info.VldFrames = VldFrames;
% add a timeline field
marker_data.Time = (1/marker_data.Info.frequency:...
    1/marker_data.Info.frequency:...
    marker_data.Info.NumFrames/marker_data.Info.frequency)';

data.marker_data = marker_data;

% get analog data
[analogs, analogsInfo] = btkGetAnalogs(acq);
if ~isempty(fieldnames(analogs))
    analog_data.Channels = analogs;
    analog_data.Info = analogsInfo;
    analog.data.Info.NumFrames = btkGetAnalogFrameNumber(acq);
    analog_data.Time = (1/analog_data.Info.frequency:1/analog_data.Info.frequency:...
        analog.data.Info.NumFrames/analog_data.Info.frequency)';
    % write to data structure
    data.analog_data = analog_data;
end

% get the forceplate info (corners etc)
[FPs, FPsInfo] = btkGetForcePlatforms(acq);
if ~isempty(FPs)
    % get the ground reaction data
    grw_data = btkGetGroundReactionWrenches(acq,grw_threshold);
    %grw_data = btkGetForcePlatformWrenches(acq,grw_threshold);
    fp_data.GRF_data = grw_data;
    fp_data.Info = FPsInfo;
    fp_data.FP_data = FPs;
    % convert the forceplate corrdinates to a xyz vector. 
    for i = 1 : length(fp_data.FP_data)
        fp_data.FP_data(i).corners = fp_data.FP_data(i).corners';
    end
    fp_data.Time = (1/FPsInfo(1).frequency:1/FPsInfo(1).frequency:...
        length(grw_data(1).P)/FPsInfo(1).frequency)';
    % write to data structure
    data.fp_data = fp_data;
end

% get scalar data (e.g. ultrasound data is stored as scalar)
[scalars, scalarsInfo] = btkGetScalars(acq);
if ~isempty(scalars)
    scalar_data.Data = scalars;
    scalar_data.Info = scalarsInfo;
    scalar_data.Time = (1/scalar_data.Info.frequency:...
        1/scalar_data.Info.frequency:...
        marker_data.Info.NumFrames/marker_data.Info.frequency)';
        data.scalar_data = scalar_data;
end

% get some subject information
% md = btkGetMetaData(acq);
% sub_info.Filename = file;
% if isfield(md.children,'SUBJECTS')
%     sub_info.Name = md.children.SUBJECTS.children.NAMES.info.values{1};
%     sub_info.MarkerSet = md.children.SUBJECTS.children.MARKER_SETS.info.values{1};
% else sub_info.Name = 'UNKNOWN';
%     subinfo_MarkerSet = 'UNKNOWN';
% end
% if isfield(md.children,'PROCESSING')
%     proc_fields = fieldnames(md.children.PROCESSING.children);
%     for i = 1:length(proc_fields)
%         sub_info.Processing_Data.(proc_fields{i}) = md.children.PROCESSING.children.(proc_fields{i}).info.values;
%     end
%     if isfield(sub_info.Processing_Data,'Bodymass')
%         data.Mass = sub_info.Processing_Data.Bodymass;
%     end
%     if isfield(sub_info.Processing_Data,'Height')
%         data.Height = sub_info.Processing_Data.Height;
%     end
% end

% data.Name = sub_info.Name;
% data.sub_info = sub_info;




    

