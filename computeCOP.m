function [structData] = computeCOP(structData)
% Calculates the COP and Free moment from force and moment data 
%   

% Written by Thor Besier, James Dunne, Cyril (Jon) Donnelly, 2008
% Last Modified; James Dunne, September (2014).     

for i = 1 : length(structData.fp_data.GRF_data)
        
        % Place these data into their variable Name. This is slower but is
        % easier to read and make ajustments later.
        
        fx = structData.fp_data.GRF_data(i).F(:,1);
        fy = structData.fp_data.GRF_data(i).F(:,2);
        fz = structData.fp_data.GRF_data(i).F(:,3);
        mx = structData.fp_data.GRF_data(i).M(:,1);
        my = structData.fp_data.GRF_data(i).M(:,2);
        mz = structData.fp_data.GRF_data(i).M(:,3);

        % Get the center of the forceplate
        fpCenter = ...
        [(structData.fp_data.FP_data(i).corners(1,:) +...
        structData.fp_data.FP_data(i).corners(2,:) +...
        structData.fp_data.FP_data(i).corners(3,:) +...
        structData.fp_data.FP_data(i).corners(4,:))]/4;
         
        xo = fpCenter(1);  
        yo = fpCenter(2);
        %zo = fpCenter(3);       
        
        COPx = (-1*(my + fx)./fz) + xo;
        COPy =((mx - fy)./fz) + yo;
        COPz = zeros(length(COPx),1);
        
        % Calculate the free moment (Tz) of the forceplate
        % different meaning of origin among different FP type
        if structData.fp_data.FP_data(i).type == 3   
            % Calculate Torque
            Tz = mz - ((COPx).*fy) + ((COPy).*fx);
        else
            % Calculate Torque
            Tz = mz - ((COPx - xo).*fy) + ((COPy - yo).*fx);
        end

        Tx = zeros(1,length(Tz))';
        Ty = Tx;
        
        % Go from mm to meters
        COPx = COPx/1000;    
        COPy = COPy/1000;
        Tz   = Tz/1000;
        
        % take out any Nans
        nNaN    = find(isnan(COPx));
        COPx(nNaN) = 0;
        COPy(nNaN) = 0;
        COPz(nNaN) = 0;
        Ty(nNaN)   = 0;
        Tz(nNaN)   = 0;
        Tx(nNaN)   = 0;
        
        % save the processed COP to the structure
        structData.fp_data.GRF_data(i).P = [COPx COPy COPz];  
        % save the processed Z torque to the structure
        structData.fp_data.GRF_data(i).M= [Tx Ty Tz];  
         
end
       
end