function [oData] = rotateCoordinateSys(oData,rotation)
% rotateCoordinateSys()
% Rotates xyz data about an axis. Data is assumed to be an nx3 sized array.
% 
% Input - oData - either a nX3 matrix or a strucutre containing matrix variables 
%                 eg oData.LASI = [nx3]
%         rotation - A ordered cell array of axes strings and rotation
%                    values. rotation = [{'z' 90 'x' 90}] is an ordered  
%                    rotation of 90 degrees first about the Z then about the X
%

% Author: Cyril J Donnelly, James Dunne, Thor Besier.  
% Written: Feb 2009     updated: Oct 2014

%%  Determine the coordinate to rotate around   
nRot = length(rotation)/2;

for i = 1 : nRot   
   
    rotAxis = char(rotation(  i*2-1 ));
    Rot     = cell2mat(rotation( i*2 ));
    
    if ischar(Rot)
        a = '''x''';
        b = '''z''';
        error(['input rotation is incorrect; must indicate axis and the amount of rotation i.e. {' a ' 90 ' b ' 90 }']) 
    end
    
    % Create roation matrices according to Rot (degrees)   
    RotAboutX1 = [1,0,0;0,cos(Rot*pi/180),-(sin(Rot*pi/180));0,sin(Rot*pi/180),cos(Rot*pi/180)];
    RotAboutY1 = [cos(Rot*pi/180),0,sin(Rot*pi/180);0,1,0;-(sin(Rot*pi/180)),0,cos(Rot*pi/180)];
    RotAboutZ1 = [cos(Rot*pi/180),-(sin(Rot*pi/180)),0;sin(Rot*pi/180),cos(Rot*pi/180),0;0,0,1];
    % choose which rotation matrix to use based on user input 
    if     strcmpi(rotAxis,'x') 
        rotationMatrix = RotAboutX1;
    elseif strcmpi(rotAxis,'y') 
        rotationMatrix = RotAboutY1;
    elseif strcmpi(rotAxis,'z')
        rotationMatrix = RotAboutZ1;
    end
        
	% Rotate nx3 arrays by the rotation matrix  
    if isstruct(oData)
        fields  = fieldnames(oData);
        nFields = length(fields);
        nData = oData;

        for u = 1:nFields
                % assign the strucutre field to data 
                vectorData = oData.(fields{u});
                % Rotate the data
                rotatedData = [rotationMatrix'*vectorData']';
                % assign the rotated data back to field  
                oData.(fields{u}) = rotatedData;  
        end
    elseif ismatrix(oData)
        % Data is a n*3 matrix
        oData = [rotationMatrix'*oData']';

    end
end
