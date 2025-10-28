function estimateHJC(static_cal_file , right_leg_cal_file, ...
              left_leg_cal_file, to_add_HJC_files, static_ref_dir, to_add_HJC_dir, write_dir, ...
              use_same_pelvic_marker_set)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FIRST STEP
%
%Use a static trial to average all marker locations.  These will be used
%to:
%
%   1.  Find a starting reference from which movements can be related
%   2.  Note marker locations with reference to each other, so that a
%       reference frame can be determined in HJC-written files in which
%       markers normally used for reference are missing.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%cd 'C:\Inversekin\Dysport'

% PelvMarks={'L.ASIS';'R.ASIS';'L.PSIS';'R.PSIS';'S2';'R.Iliac';'L.Iliac'};
% RTMarks={'R.TH1';'R.TH2';'R.TH3';'R.TH4';'R.Knee';'R.MKnee'};
% LTMarks={'L.TH1';'L.TH2';'L.TH3';'L.TH4';'L.Knee';'L.MKnee'};

PelvMarks={'LASIS';'RASIS';'LWST';'RWST';'RIC';'LIC';'LPSIS';'RPSIS'};
RTMarks={'RTC';'RILT';'RIMT';'RSLT';'RSMT';'RLK'};
LTMarks={'LLS1','LLS2','LLS3','LLS4'};

statfile=[]; statfiler=[]; statfilel=[];
RTMrkStat={}; LTMrkStat={};

display('Select Static Calibration Trial');
%[statfile,dir]=uigetfile('*.trc','Select Static Calibration Trial');
[path,name,ext] = fileparts(static_cal_file);
statfile = [name ext];
dir = path;

if statfile==0
    display('ERROR: No Static Trial Selected; Try Again')
    return;
end
[x2,tx2,sfx2,nsx2,nmrk2,mnames2,file2,inpath2]=load_trc(statfile,dir);
[rows,columns]=size(x2);
%Calculate the mean values of each marker
mx=mean(x2);

a=1;
for i=1:length(PelvMarks)
    for j=1:length(mnames2)
        mrk=strcmp(PelvMarks(i),mnames2(j));
        if mrk==1
            PMrkStat(a,:)=PelvMarks(i);
            a=a+1;
        end
    end
end

a=1;
for i=1:length(RTMarks)
    for j=1:length(mnames2)
        mrk=strcmp(RTMarks(i),mnames2(j));
        if mrk==1
            RTMrkStat(a,:)=RTMarks(i);
            a=a+1;
        end
    end
end

LTMrkStat={};
a=1;
for i=1:length(LTMarks)
    for j=1:length(mnames2)
        mrk=strcmp(LTMarks(i),mnames2(j));
        if mrk==1
            LTMrkStat(a,:)=LTMarks(i);
            a=a+1;
        end
    end
end

%The marker locations are averaged throughout the trial and arranged in
% a vector so that they can be used by soder.m
for j=1:length(PMrkStat)
    P_loc(j)=strmatch(PMrkStat(j),mnames2,'exact');
    PStatLoc(:,((j-1)*3+1):((j-1)*3+3))=x2(:,((P_loc(j)-1)*3+1):((P_loc(j)-1)*3+3));
end

RTStatLoc=[];
for j=1:length(RTMrkStat)
    RT_loc(j)=strmatch(RTMrkStat(j),mnames2,'exact');
    RTStatLoc(:,((j-1)*3+1):((j-1)*3+3))=x2(:,((RT_loc(j)-1)*3+1):((RT_loc(j)-1)*3+3));
end

LTStatLoc=[];
for j=1:length(LTMrkStat)
    LT_loc(j)=strmatch(LTMrkStat(j),mnames2,'exact');
    LTStatLoc(:,((j-1)*3+1):((j-1)*3+3))=x2(:,((LT_loc(j)-1)*3+1):((LT_loc(j)-1)*3+3));
end


if size(PStatLoc) > 1
    Pelv_Ref=mean(PStatLoc);
else Pelv_Ref=[];
end

if size(RTStatLoc) >1
    RT_Ref=mean(RTStatLoc);
else RT_Ref=[];
end

if size(LTStatLoc) >1
    LT_Ref=mean(LTStatLoc);
else LT_Ref=[];
end

%Write the average locations to a file for future reference

ref_data=[Pelv_Ref RT_Ref LT_Ref;Pelv_Ref RT_Ref LT_Ref];
[a,b]=size(ref_data);

time=[0;1];

ref_mrks=strvcat([PMrkStat; RTMrkStat; LTMrkStat;]);


done = writeTRCFile(time,ref_data,ref_mrks,static_ref_dir,'static_ref');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SECOND STEP
%
%Find the hip joint centers for each, or either, leg.
%
%It is found by an appended m.file, called genHipKinem2.m.  This file finds
%the hip joint center location using a least squares method to find where
%pelvis and thigh share a common point.  It finds the point relative to
%thigh and pelvic frames and averages the distance between the two.  It
%outputs the averaged HJC in the mid-asis pelvic frame, along with the
%average distance between the two HJC's calculated and the standard
%deviation of all locations.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%User selects the trial for right leg HJC calibration

display('Right Leg HJC Calibration Trial');
%[infile, inpath]=uigetfile('*.trc','Select Right Leg HJC Calibration Trial');
[path,name,ext] = fileparts(right_leg_cal_file);
infile = [name ext];
inpath = path;
display(['Loading file ' infile])

%CALCULATE THE RIGHT HJC relative to pelvis anatomical center, along with mean
%distance between estimated HJC's relative to pelvis and thigh and the
%standard deviation in those measurements

if infile~=0
    [RHJC,R_Ave,R_Std_Dev]=genHipKinem2(infile,inpath,1,Pelv_Ref,RT_Ref,PMrkStat,RTMrkStat,dir);
else
    RHJC=0;
    R_Ave=0;
    R_Std_Dev=0;
end


%User selects the trial for left leg HJC locating
display('Left Leg HJC Calibration Trial');
%[infile, inpath]=uigetfile('*.trc','Select Left Leg HJC Calibration Trial');
[path,name,ext] = fileparts(left_leg_cal_file);
infile = [name ext];
inpath = path;
display(['Loading file ' infile])

%CALCULATE THE LEFT HJC relative to pelvis anatomical center, along with mean
%distance between estimated HJC's relative to pelvis and thigh and the
%standard deviation in those measurements

if infile~=0
    [LHJC,L_Ave,L_Std_Dev]=genHipKinem2(infile,inpath,2,Pelv_Ref,LT_Ref,PMrkStat,LTMrkStat,dir);
else LHJC=0;
    L_Ave=0;
    L_Std_Dev=0;
    if RHJC==0;
        display('No HJC Calibration Trials Selected.  Please Try Again.')
        return;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%THIRD STEP
%
%Select files to write HJC's to, and calculate the HJC's
%
%User is prompted to select files to write HJC's.  Then, each file is
%analyzed, and a pelvic reference frame is determined based on available
%markers.  The HJC is then transformed into this reference frame, and from
%the marker data in the file, the HJC is calculated in the global frame and
%written into the file.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Select files to which the HJC locations should be added (multiple must be
%selected).  The same directory will be used to create a file for new files
%with HJC locations.




display('Select Files to Add Hip Joint Center Locations')
%[files,directory]=uigetfile('*.trc','Select Files to Add HJC Locations','multiselect','on');
files = to_add_HJC_files;
directory = to_add_HJC_dir;

files=strvcat(files(:,:));

if length(files)==1;
    display('No Files Selected to Write HJC.  Try Again.');
    return;
end

directory=strvcat(directory);
infile=files(1,:); inpath=directory;
%mkdir(directory,'Files_W_HJCs')


%Reformat file list and directory information

[a b]=size(files);
%write_directory=[directory char('Files_W_HJCs')];
write_directory = write_dir;

%Locate pelvic markers in each data file, find the pelvic center, locate
%the HJC's and marker data and HJC data into a new file
reset=0;
for j=1:a;
    infile=files(j,:);
    inpath=directory;
    x=0; mnames=0; tx=0; marks=0; time=0;
    [x,tx,sfx,nsx,nmrk,mnames,file,inpath]=load_trc(infile,inpath);
    
    %Figure out which markers in the static trial are available in the
    %files to which HJC's are being written
    c=1;
    for i=1:length(PMrkStat)
        for k=1:length(mnames)
            mrk=strcmp(PMrkStat(i),mnames(k));
            if mrk==1
                PMrkWritA(c,:)=PMrkStat(i);
                c=c+1;
            end
        end
    end
    
    
    %Locate where each marker locatiion is written into the file
    for i=1:length(PMrkWritA)
        PMrkWritANum(i,1)=strmatch(PMrkWritA(i),mnames,'exact');
    end
    
    
    %Allow the user to decide whether to use the same pelvic marker set
    %for all trials or to select each trial's marker set individually
    
    if reset==0
        %button=questdlg('Would you like to use the same pelvic marker set in every file to reference and write hip joint centers?','Pelvic Markers in HJC Write trials','Yes','No!','Yes');
        reset=1;
    end
    
    
    %If yes, allow the user to pick once and move on
    
    if use_same_pelvic_marker_set
        if reset==1
            %[PMrkWrit,d]=listdlg('ListString',mnames,'PromptString','Select Pelvic Markers to Find HJC for Recording','InitialValue',PMrkWritANum);
            PMrkWrit = sort(PMrkWritANum)';
            reset=2;
        end
        
    else  %If no, allow the user to pick markers for each trial being written

        % TODO
%         if reset==1;
%             [PMrkWrit,d]=listdlg('ListString',mnames,'PromptString','Select Pelvic Markers to Find HJC for Recording','InitialValue',PMrkWritANum);
%         end
    end
    
    %Locate Pelvic Markers in the static trial.
    
    for i=1:length(PMrkWrit)
        mrklist(i,1)=mnames(PMrkWrit(i));
        WrMrkStat(i)=strmatch(mrklist(i),ref_mrks,'exact');
        PelvRefWrite(1,((i-1)*3+1):((i-1)*3+3))=[ref_data(1,((WrMrkStat(i)-1)*3+1):((WrMrkStat(i)-1)*3+3))];
    end
    
    %Find the center of the markers and define the coordinate
    %system of the marker set as identical to the global frame.
    
    D_ref=[mean([reshape(PelvRefWrite,3,length(PMrkWrit))]')]';
    R_ref=[1 0 0;0 1 0;0 0 1];
    
    %For static data, locate the rotation and location of the
    %coordinate system that the known HJC is in.
    
    RasisNum=strmatch('RASIS',ref_mrks,'exact');
    LasisNum=strmatch('LASIS',ref_mrks,'exact');
    RICnum = strmatch('RIC',ref_mrks,'exact');
%     MidASISnum = strmatch('MidASIS',ref_mrks,'exact');
%     MidICnum = strmatch('MidIC',ref_mrks,'exact');
%     TopASISnum = strmatch('TopASIS',ref_mrks,'exact');
    SacNum=strmatch('S2',ref_mrks,'exact');
    
    rasis=[ref_data(1,((RasisNum-1)*3+1):((RasisNum-1)*3+3))]';
    lasis=[ref_data(1,((LasisNum-1)*3+1):((LasisNum-1)*3+3))]';
    ric = [ref_data(1,((RICnum-1)*3+1):((RICnum-1)*3+3))]';
%     MidASIS = [ref_data(1,((MidASISnum-1)*3+1):((MidASISnum-1)*3+3))]';
%     MidIC = [ref_data(1,((MidICnum-1)*3+1):((MidICnum-1)*3+3))]';
%     TopASIS = [ref_data(1,((TopASISnum-1)*3+1):((TopASISnum-1)*3+3))]';
    if isempty(SacNum)==1
        RPsisNum=strmatch('RPSIS',ref_mrks,'exact');
        rpsis=ref_data(1,((RPsisNum-1)*3+1):((RPsisNum-1)*3+3))';
        LPsisNum=strmatch('LPSIS',ref_mrks,'exact');
        lpsis=ref_data(1,((LPsisNum-1)*3+1):((LPsisNum-1)*3+3))';
        sacral=(rpsis+lpsis)/2;
    else
        sacral=[ref_data(1,((SacNum-1)*3+1):((SacNum-1)*3+3))]';
    end

%     y = TopASIS - MidASIS;
%     y = y/sqrt(sum(y.*y));
%     z = rasis - MidASIS;
%     z = z/sqrt(sum(z.*z));
%     X = cross(y,z);
%     X = X/sqrt(sum(X.*X));
%     R = [X y z];

    midasis=(lasis+rasis)/2;
    y=lasis-rasis;
    y=y/sqrt(sum(y.*y));
%     z = cross((ric-lasis),(rasis-lasis));
    z=cross((sacral-lasis),(rasis-lasis));
    z=z/sqrt(sum(z.*z));
    X=cross(y,z);
    R=[X y z];
    
    
    %Find the Transformation Matrix from the HJC system to the
    %marker-set system. Then, find the HJC in the m-s system.
    
    D=midasis-D_ref;
    T=[R D;0 0 0 1];
    
    rhjcms=T*RHJC;
    lhjcms=T*LHJC;
    
    %Locate the markers in the marker set throughout the trial being
    %written to, then use soder to find the transformations to each
    %time set of markers.
    
    [a,b]=size(x);
    for i=1:a;
        for k=1:length(PMrkWrit)
            marks(i,((k-1)*3+1):((k-1)*3+3))=[x(i,((PMrkWrit(k)-1)*3+1):((PMrkWrit(k)-1)*3+3))];
        end
        
        [T_pelv,Res]=soder([PelvRefWrite;marks(i,:)]);
        
        %From these T, find the HJC in the global frame.
        center(i,:)=[mean(reshape(marks(i,:),3,length(marks(1,:))/3)')];
        
        
        if RHJC~=0
            rr=[center(i,:)]'+[T_pelv(1:3,1:3)*rhjcms(1:3,1)];
            r_hjc(i,(1:3))=[rr(1:3)]';
        end
        if LHJC~=0
            ll=[center(i,:)]'+[T_pelv(1:3,1:3)*lhjcms(1:3,1)];
            l_hjc(i,(1:3))=[ll(1:3)]';
        end
        time(i,1)=i/sfx-1/sfx;
    end
    
    mrkdata=0; mrknames=0;
    if RHJC~=0
        if LHJC~=0
            mrkdata=[x r_hjc(1:length(x(:,1)),:) l_hjc(1:length(x(:,1)),:)];
            mrknames=char([mnames;cellstr('RHJC');cellstr('LHJC')]);
        else
            mrknames=char([mnames;cellstr('RHJC')]);
            mrkdata=[x r_hjc(1:length(x(:,1)),:)];
        end
    else
        if LHJC~=0
            mrkdata=[x l_hjc(1:length(x(:,1)),:)];
            mrknames=char([mnames;cellstr('LHJC')]);
        else
            disp('No HJC Calibration Trials Selected. Try Again.');
            return;
        end
    end
    nfile=strcat(files(j,:));
    filen=nfile(1,1:(length(nfile)-4));
     
    % rotate the marker data into an OpenSim model coordinate system
%     R = [1 0 0; 0 0 -1; 0 1 0];
%     for i=1:3:size(mrkdata,2)-2;
%         mrkdata(:,i:i+2)=mrkdata(:,i:i+2)*R;
%     end
    done = writeTRCFile(time,mrkdata,mrknames,write_directory,[filen]);
    if done==1;
        display(['File ' infile ' written with HJC locations']);
    end
    
end

%Write a text file that gives statistical information about HJC location
info=[R_Ave; R_Std_Dev; L_Ave; L_Std_Dev];
fid=fopen([write_directory,'/statistics.txt'],'w');
fprintf(fid,['Statistics of Calibration of HJC (locations relative to pelvis)      \n']);
fprintf(fid,['                                                                     \n']);
fprintf(fid,['Right Leg HJC                      Left Leg HJC                      \n']);
fprintf(fid,['X          Y           Z           X          Y          Z           \n']);
if RHJC~=0;
    fprintf(fid,'%-f',RHJC(1)); fprintf(fid,['   ']); fprintf(fid,'%-f',RHJC(2)); fprintf(fid,['   ']); fprintf(fid,'%-f',RHJC(3)); fprintf(fid,['   ']);
end
if LHJC~=0;
    fprintf(fid,'%-f',LHJC(1)); fprintf(fid,['   ']); fprintf(fid,'%-f',LHJC(2)); fprintf(fid,['   ']); fprintf(fid,'%-f',LHJC(3)); fprintf(fid,['   \n']);
end
fprintf(fid,['                                                                     \n']);
fprintf(fid,['Mean         Std. Dev.             Mean         Std. Dev.            \n']);
fprintf(fid,'%-f',R_Ave);fprintf(fid,['     ']);
fprintf(fid,'%-f',R_Std_Dev);fprintf(fid,['              ']);
fprintf(fid,'%-f',L_Ave);fprintf(fid,['     ']);
fprintf(fid,'%-f\n',L_Std_Dev);
fclose(fid);
end

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
function [pos,time,f,n,nmrk,mrk_names,file,inpath]=load_trc(infile,inpath)
%   [pos,time,f,n,nmrk,mrk_names]=load_trc(infile)
%   LOAD_TRC is used to open a data file from Motion Analysis Realtime
%   output (*.trc).
%
%   Inputs:
%       infile - trc file to be loaded
%                If infile is unspecified, the user is prompted to select the input file
%       inpath - directory of location where data file is located
%               when no path is specified, it defaults to current directory
%
%   Outputs:
%       pos     contains - the meaured marker positions in order of the markers
%               that is columns 1-3 are the x,y,z components of marker 1
%                       columns 4-6 are the x,y,z components of marker 2
%                          ....
%       time - column vector of time
%       f - sample frequency
%       n - number of data frames
%       nmrk - number of markers
%       mrk_names - marker names
%
%   Updated: Feb. 15, 2006 (JWF)
%
%   MATLAB Version 7.1

n = nargin;
if (n==0);
    [infile, inpath]=uigetfile('*.trc','Select input file');
    if infile==0;
        f='';
        n='';
        nmrk='';
        mrk_names='';
        data=[];
        return;
    end
    fid=fopen([inpath infile],'r');
    file = infile(1:length(infile)-4);
elseif (n==1);
    file = infile(1:length(infile)-4);
    fid=fopen(infile,'r');
else (n==1);
    file = infile(1:length(infile)-4);
    fid=fopen(fullfile(inpath,infile),'r');
end

if (fid==-1);
    disp('File not found');
    f='';
    n='';
    nmrk='';
    mrk_names='';
    data=[];
    return;
end

disp(['Loading file...' infile] );

%disregard header info
for h=1:2
    hdr=fgetl(fid);
end
file_info=fscanf(fid,'%f');
f=file_info(1);
nmrk=file_info(4);
hdr=fscanf(fid,'%s',4);
line=fgetl(fid);
line=fgetl(fid);
j=1;
jl=length(line);
for i=1:(nmrk+2)
    name=sscanf(line(j:jl),'%s',1);
    ii=findstr(line(j:jl),name);
    j=j+ii(1)+length(name);
    if i>2
        mrk_names(i-2,1)=cellstr(name);
    end
end
%mrk_names=fscanf(fid,'%s',nmrk+2);
for h=1:2
    hdr=fgetl(fid);
end

line=[];
data=[];
while(length(data)<((nmrk*3)+2))
    line=fgetl(fid);
    data=sscanf(line,'%f');
end
time(1,1)=data(2);
pos(1,:)=data(3:length(data));
i=1;
while feof(fid)==0
    i=i+1;
    line=fgetl(fid);
    data=sscanf(line,'%f');
    time(i,1)=data(2);
    for j=3:length(data)
        pos(i,j-2)=data(j);
        %%%%%%%%
    end
end

[n,nc]=size(pos);
if n==1
    time=1;
else
    time=time(1:n,1);
end
% Return the position data in m
pos=pos/1000;
end

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------

function[HJC,Ave,Std_Dev]=genHipKinem2(infile,inpath,leg,Pelv_Ref,T_Ref,PelvMarks,ThMarks,Ref_Dir);
% Description:
%       Finds the location of the hip joint center between a pelvis and a
%       thigh, given a calibration trial.  The HJC is located using a least
%       squares regression to find two points, on each of the pelvis and
%       thigh segments, that minimize the distance between them throughout
%       the trial.
%
% INPUTS:
%       infile: The .trc file to be used for locating the hip joint center.
%       inpath: The path of infile.
%       leg: For which leg the hip joint center is being located.
%          Right Leg: leg=1 ; Left Leg: leg=2
%       PelvMrks: The names of the pelvis markers available for calculating
%           the HJC in the rewritten files
% OUTPUTS:
%       HJC: The location of the hip joint center in the pelvis anatomical
%            reference frame.
%        Ave: The average distance between the two points throughout the
%             trial
%       Std_Dev: The standard deviation of the point distances for the
%             trial
%
% Reference: Piazza, S. J., A. Erdemir, et al. (2004). "Assessment of the
%       functional method of hip joint center location subject to reduced range
%       of hip motion." J Biomech 37(3): 349-56.
%
% AUTHOR: Joseph Farron, NMBL, University of Wisconsin-Madison
%
% DATE: November 2, 2005
% VERSION: 2.0
%
% MATLAB Version 7.1

%Load trial for analysis

[x,tx,sfx,nsx,nmrk,mnames,file,inpath]=load_trc(infile,inpath);

%set range of data to be used
[ll,ww]=size(x);
x=x((1:ll),:);

%Allow the user to select the pelvic markers to be used in
%calibration.  Then, load the markers.

for i=1:length(PelvMarks)
    PelvMrksCal(i,1)=strmatch(PelvMarks(i),mnames);
end
% 
% [Pmark_list,d]=listdlg('ListString',mnames,'PromptString','Select Pelvic Markers for Locating HJC','InitialValue',PelvMrksCal);
% 
% for i=1:length(Pmark_list)
%     PelvMrks(i,1)=mnames(Pmark_list(i));
% end

Pmark_list = PelvMrksCal;
PelvMrks = PelvMarks;

%Allow the user to select the thigh markers to be used in calibration.
% Then, load the markers.
 
for i=1:length(ThMarks)
    ThighMrksCal(i)=strmatch(ThMarks(i),mnames);
end
% 
% [THmark_list,d]=listdlg('ListString',mnames,'PromptString','Select Thigh Markers for Locating HJC','InitialValue',ThighMrksCal);
% 
% for i=1:length(THmark_list)
%     ThighMrks(i,1)=mnames(THmark_list(i));
% end

THmark_list = ThighMrksCal;
ThighMrks = ThMarks;

%Set the marker locations to refer to in calibration based on static
%trial averages

[Pelvis_Ref,Thigh_Ref]=ref_load(PelvMrks,ThighMrks,Ref_Dir);

%Find the initial vector to the center of the right thigh frame
pmat=reshape(Pelvis_Ref,3,length(PelvMrks));
d_p0=[mean([pmat]')]';
%Find the initial vector to the center of the right thigh frame
tmat=reshape(Thigh_Ref,3,length(ThighMrks));
d_t0=[mean([tmat]')]';

for i=1:length(x);
    
    
    %Identify which set of coordinates belongs to each requested pelvic
    %marker, then list the coordinates as "pcoords"
    
    for ii=1:length(PelvMrks)
        for jj=1:3
            pcoords(i,(ii-1)*3+jj)=x(i,(Pmark_list(ii)-1)*3+jj);
        end
    end
    
    %Identify which set of coordinates belongs to each requested right
    %thigh marker,then list the coordinates as "rtcoords"
    if length(ThighMrks)>2;
        for ii=1:length(ThighMrks)
            for jj=1:3
                tcoords(i,(ii-1)*3+jj)=x(i,(THmark_list(ii)-1)*3+jj);
            end
        end
    end
    
    
    %Find pelvis transfer function, rotation matrix, and coordinates
    [T_p,res_p]=soder([Pelvis_Ref;pcoords(i,:)]);
    R_p=T_p(1:3,1:3);
    d_p(:,i)=R_p*d_p0+T_p(1:3,4);
    %record R_p in horizontal vectors
    pp(i,:)=reshape(R_p',1,9);
    dp(:,i)=T_p(1:3,4);
    
    if length(ThighMrks)>2;
        
        %Find right thigh transfer function, rotation matrix, and vector
        [T_t,res_t]=soder([Thigh_Ref;tcoords(i,:)]);
        R_t=T_t(1:3,1:3);
        %record rotation matrices and location vectors
        d_t(:,i)=R_t*d_t0+T_t(1:3,4);
        r_t(i,:)=reshape(R_t',1,9);
        %Find directional vector and rotation matrix from pelvis to right
        %thigh and record them in rr and txr
        p_d_t=d_t(:,i)-d_p(:,i);
        [p_R_t]=(R_p)'*(R_t);
        tt(i,:)=reshape(p_R_t',1,9);
        txt(i,:)=R_p'*p_d_t;
    end
end
%Find matrix A and vector b so that A^-1*b=x, where x=[x y z u v w]', and
%(x,y,z) is the vector from the pelvis to the HJC, and (u,v,w) is the
%vector from the thigh to the HJC

if length(ThighMrks)>2;
    [A,b]=load_A(tt,txt);
    HJC_t=(A^(-1))*b;
end



% Locate Pelv Ref Marks in file
RASISnum=strmatch('RASIS',mnames);
LASISnum=strmatch('LASIS',mnames);
RICnum = strmatch('RIC',mnames);
% MidASISnum = strmatch('MidASIS',mnames);
% MidICnum = strmatch('MidIC',mnames);
% TopASISnum = strmatch('TopASIS',mnames);
SACnum=strmatch('S2',mnames);

% Find the pelvic and thigh rotation frames and vector locations to
% calculate the HJC in the global reference frame for the locations found
% by HJC_r and HJC_l
[m,n]=size(x);
for i=1:m
    R_p=[reshape(pp(i,:),3,3)]';
    if length(ThighMrks)>2;
        R_t=[reshape(r_t(i,:),3,3)]';
        HJC_t_p(i,:)=[[R_p]*HJC_t(1:3,1)+d_p(:,i)]';
        HJC_t_t(i,:)=[[R_t]*HJC_t(4:6,1)+d_t(:,i)]';
    else d_t(:,i)=[0;0;0];
    end
    
    time(i,1)=i/sfx;
    
    %Convert the HJC locations into the pelvic anatomical reference frame
    
    for jj=1:3
        
       try
            rasis(jj,1)=[x(i,(RASISnum-1)*3+jj)]';
       catch ME
            keyboard
       end

       lasis(jj,1)=[x(i,(LASISnum-1)*3+jj)]';
        if isempty(SACnum)==1
            RPSISnum=strmatch('RPSIS',mnames);
            LPSISnum=strmatch('LPSIS',mnames);
            sacral(jj,1)=(x(i,(RPSISnum-1)*3+jj)' + x(i,(LPSISnum-1)*3+jj)')/2;
        else
            sacral(jj,1)=[x(i,(SACnum-1)*3+jj)]';
        end
        ric(jj,1)=[x(i, (RICnum-1)*3+jj)]';
%         MidASIS(jj,1)=[x(i, (MidASISnum-1)*3+jj)]';
%         MidIC(jj,1)=[x(i, (MidICnum-1)*3+jj)]';
%         TopASIS(jj,1)=[x(i, (TopASISnum-1)*3+jj)]';
    end
    
%     y = TopASIS - MidASIS;
%     y = y/sqrt(sum(y.*y));
%     z = rasis - MidASIS;
%     z = z/sqrt(sum(z.*z));
%     X = cross(y,z);
%     X = X/sqrt(sum(X.*X));
%     R = [X y z];
    
    midasis=(lasis+rasis)/2;
    y=lasis-rasis;
    y=y/sqrt(sum(y.*y));
    %z = cross((ric-lasis),(rasis-lasis));
    z=cross((sacral-lasis),(rasis-lasis));
    z=z/sqrt(sum(z.*z));
    X=cross(y,z);
    R=[X y z];
    T=[R midasis; 0 0 0 1];

    T = [R midasis; 0 0 0 1];
    
    Rp=[reshape(pp(i,:),3,3)]';
    pRp=[Rp]'*[R];
    p_d_p=[R]'*[d_p(:,i)-midasis];
    pTp=[[pRp]' p_d_p; 0 0 0 1];
    ptp(i,:)=reshape(pTp',1,16);
    
    if length(ThighMrks)>2;
        Rt=[reshape(r_t(i,:),3,3)]';
        tRp=[Rt]'*[R];
        t_d_p=[R]'*[d_t(:,i)-midasis];
        tTp=[[tRp]' t_d_p;0 0 0 1];
        ttp(i,:)=reshape(tTp',1,16);
    end
end

%Average the transformation matrices over the trial
pTp=mean(ptp);
pTp=[reshape(pTp,4,4)]';

if length(ThighMrks)>2;
    tTp=mean(ttp);
    tTp=[reshape(tTp,4,4)]';
    tTp=mean(ttp);
    tTp=[reshape(tTp,4,4)]';
    
    hjc_pt=pTp*[HJC_t(1:3);1];
    hjc_tt=tTp*[HJC_t(4:6);1];
    
    %Estimate the HJC as the mean location between the two HJC estimations
    HJC=[mean([[hjc_pt]';[hjc_tt]'])]';
    
    for i=1:length(x);
        Hdiff(i,1)=sqrt((HJC_t_p(i,1)-HJC_t_t(i,1))^2+(HJC_t_p(i,2)-HJC_t_t(i,2))^2+(HJC_t_p(i,3)-HJC_t_t(i,3))^2);
    end
    Ave=sum((Hdiff(:,1)))/length(x);
    Std_Dev=sqrt((sum((Hdiff(:,1)-Ave).^2))/(length(x)-1));
end
end
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------

function [T,res]=soder(data)
% function [T,res]=soder(data)
%
% Description:	Program calculates the transformation matrix T containing
%		the rotation matrix (3x3) and the translation translation
%		vector d (3x1) for a rigid body segment using a singular
%		value decomposition method (Soederkvist & Wedin 1993).
%
% Input:    data:   columns represent the XYZ positions and the rows
%                   represent time.
% Output:   T:      4x4 Matrix containing the rotation matrix R and the
%                   translation d: T = [R,d; 0 0 0 1]
%           res:    norm of residuals (measure of fit; "rigidity" of body
%
% References:     Soderkvist I. and Wedin P. -A., (1993). Determining the
%                 movements of the skeleton using well-configured markers.
%                 Journal of Biomechanics, 26:1473-1477
%
% Author:	Christoph Reinschmidt, HPL, The University of Calgary
%               (Matlab code adapted from Ron Jacobs, 1993)
% Date:		February, 1995
% Last Changes: December 09, 1996
% Version:      3.1

if (size(data,2)/3)~=fix(size(data,2)/3),
    disp('ERROR: input has to be multiple of 3 (XYZ coordinates)'); return
end



A=[reshape(data(1,:)',3,size(data,2)/3)]';
B=[reshape(data(2,:)',3,size(data,2)/3)]';


% Checking for NaNs and also checking if still 3 pts left and if not
% T=[NaN...];
cut=[0];
qA=isnan(A); qB=isnan(B); qAB=[qA,qB];
qsum=sum(qAB'); cut=find(qsum~=0);
A([cut],:)=[];
B([cut],:)=[];
if size(A,1)<3,
    T=[NaN,NaN,NaN,NaN;NaN,NaN,NaN,NaN;NaN,NaN,NaN,NaN;NaN,NaN,NaN,NaN;]; return;
end


Amean=mean(A)';
Bmean=mean(B)';


for i=1:size(A,1)-size(cut,2),
    Ai(:,i)=[A(i,:)-Amean']';
    Bi(:,i)=[B(i,:)-Bmean']';
end


C=Bi*Ai';
[P,T,Q]=svd(C);
R=P*diag([1 1 det(P*Q')])*Q';
d=Bmean-R*(Amean);

T=[R d;0 0 0 1];

% Calculating the norm of residuals
A=A'; A(4,:)=ones(1,size(A,2));
B=B';
Bcalc=T*A; Bcalc(4,:)=[]; Diff=B-Bcalc; Diffsquare=Diff.^2;
%DOF=3*(number of points)-6 unknowns (Hx,Hy,Hz,alpha,beta,gamma):
DOF=size(B,1)*size(B,2)-6;
res=[sum(Diffsquare(:))/DOF].^0.5;
end
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------

function [A,b]=load_A(r,tx);
% Description: Finds the matrix A and the vector b such that Ax-b=0, so A'b=x
%     Where x is the vector [x y z u v w]', where (x,y,z) is the location of
%     the hip joint center in the pelvis frame, and (u,v,w) is the HJC in
%     the thigh reference frame.
%
% INPUTS: r is a 9 column matrix that is a reshaping of the rotational matrix
%           from the pelvis reference frame to the thigh reference frame at
%           every point in time, organized as:
%           r=[rxx rxy rxz ryx ryy ryz rzx rzy rzz]
%
%        tx is a three column matrix that lists the directional vector from
%           the pelvic reference frame to the thigh reference frame,
%           defined in pelvic coordinates:
%           tx=[tx ty tz]
%
% OUTPUTS: A is a 6 x 6 matrix
%         B is a 6 x 1 vector
%
% Reference: Piazza, S. J., A. Erdemir, et al. (2004). "Assessment of the
%     functional method of hip joint center location subject to reduced range
%     of hip motion." J Biomech 37(3): 349-56.
%
% AUTHOR: Joseph Farron, NMBL, University of Wisconsin-Madison
%
% DATE: October 10, 2005
% VERSION: 1.0
%
% MATLAB Version 7.1


%First find 6x6 matrix A:

A=zeros(6,6);
A(1,1)=length(r); A(2,2)=length(r); A(3,3)=length(r);
A(4,4)=sum(r(:,1).^2)+sum(r(:,4).^2)+sum(r(:,7).^2);
A(4,5)=sum(r(:,2).*r(:,1))+sum(r(:,5).*r(:,4))+sum(r(:,8).*r(:,7));
A(4,6)=sum(r(:,3).*r(:,1))+sum(r(:,6).*r(:,4))+sum(r(:,9).*r(:,7));
A(5,6)=sum(r(:,3).*r(:,2))+sum(r(:,6).*r(:,5))+sum(r(:,9).*r(:,8));
A(5,5)=sum(r(:,2).^2)+sum(r(:,5).^2)+sum(r(:,8).^2);
A(6,6)=sum(r(:,3).^2)+sum(r(:,6).^2)+sum(r(:,9).^2);
A(5,4)=sum(r(:,2).*r(:,1))+sum(r(:,5).*r(:,4))+sum(r(:,8).*r(:,7));
A(6,4)=sum(r(:,3).*r(:,1))+sum(r(:,6).*r(:,4))+sum(r(:,9).*r(:,7));
A(6,5)=sum(r(:,3).*r(:,2))+sum(r(:,6).*r(:,5))+sum(r(:,9).*r(:,8));
for l=1:3
    for m=1:3
        A((l),m+3)=-sum(r(:,((l-1)*3+m)));
        A((m+3),l)=-sum(r(:,((l-1)*3+m)));
    end
end

%Find the vector b:

b=zeros(6,1);
for p=1:3
    b(1,1)=sum(tx(:,1));
    b(2,1)=sum(tx(:,2));
    b(3,1)=sum(tx(:,3));
    b((p+3),1)=-(sum(tx(:,1).*r(:,p))+sum(tx(:,2).*r(:,(p+3)))+sum(tx(:,3).*r(:,(p+6))));
end
end





%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------

function done = writeTRCFile(time,mrkdata,mrknames,directory,file)

if length(time)<2
    time(2)=1;
end
T=time(2)-time(1);

f=1/T;
[mk,nk]=size(mrkdata);
nk=nk/3;
fid = fopen([directory,'/',file,'.trc'],'w');
fprintf(fid,['PathFileType  4	(X/Y/Z) %s\n'],directory);
fprintf(fid,['DataRate	CameraRate	NumFrames	NumMarkers	Units	OrigDataRate	OrigDataStartFrame	OrigNumFrames\n']);
fprintf(fid,['%.1f\t%.1f\t%d\t%d\tmm\t%.1f\t%d\t%d\n'],f,f,mk,nk,f,1,mk);

fprintf(fid,['Frame#	Time ']);
%fprintf(fid,'\t');
for i=1:nk
    temp=strcat(mrknames(i,:));
%     fprintf(fid,'%s\t\t\t',mrknames(i,:));
    fprintf(fid,'%s\t\t\t',temp);
end
fprintf(fid,'\n\t\t');

for i=1:nk
    if (i<10)
        fprintf(fid,'X%1d\tY%1d\tZ%1d\t',i,i,i);
    elseif (i<100)
        fprintf(fid,'X%2d\tY%2d\tZ%2d\t',i,i,i);
    else
        fprintf(fid,'X%3d\tY%3d\tZ%3d\t',i,i,i);
    end
end
fprintf(fid,'\n');
fprintf(fid,'\n');

for i=1:mk
    fprintf(fid,'%d',i);
    fprintf(fid,'\t%.5f',time(i));
    fprintf(fid,'\t%.3f',1000.*mrkdata(i,:));
    fprintf(fid,'\n');
end
fclose(fid);
done=1;
end


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function [Pelvis_Ref,Thigh_Ref]=ref_load(PelvMrks,ThighMrks,dir);
%Description: Reads marker locations averaged from a static calibration
%             trial that are used as a reference.
%Inputs:
%       PelvMrks: A cell list of the pelvis markers to be taken from the
%                 static trial averages.
%       ThighMrks: A cell list of the thigh markers to be taken from the
%                  static trial averages.
%
%Outputs:
%       Pelvis_Ref: A row of the x,y,and z points for each marker in the
%                   pelvis frame based on selected markers.
%       Thigh_Ref:  A row of the x,y, and z point for each marker in the
%                   thigh frame based on selected markers.
%
%('Static_Marker_Reference',dir)
filename='static_ref.trc';
[x2,tx,sfx,nsx,nmrk,mnames,file,inpath]=load_trc(filename,dir);
for i=1:length(PelvMrks)
    Pelv_Ref(1,i)=strmatch(PelvMrks(i),mnames);
    Pelvis_Ref(1,(((i-1)*3+1):((i-1)*3+3)))=x2(1,(((Pelv_Ref(1,i)-1)*3+1):((Pelv_Ref(1,i)-1)*3+3)));
end

for i=1:length(ThighMrks)
    Thi_Ref(1,i)=strmatch(ThighMrks(i),mnames);
    Thigh_Ref(1,(((i-1)*3+1):((i-1)*3+3)))=x2(1,(((Thi_Ref(1,i)-1)*3+1):((Thi_Ref(1,i)-1)*3+3)));
end
status=fclose('all');
end