function [betterissues] = IssuesToCells(Issues)
    % for the issues you have
    % from Jon Stingel
    % 10/13/2020
    betterissues = cell(length(Issues),1);
    for i=1:length(Issues)
        temp = cellstr(strcat(char(Issues(i,1)),'_______',char(Issues(i,2))));
        betterissues(i,1) = temp;
    end
end
