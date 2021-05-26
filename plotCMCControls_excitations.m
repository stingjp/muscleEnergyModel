function plotCMCControls_excitations(filename, outputname)
    import org.opensim.modeling.*
    
    % load the full file
    tempTimeSeriesTable = TimeSeriesTable(filename);
    % grab the time column
    temptime = tempTimeSeriesTable.getIndependentColumn();
    % convert to percent gait cycle
    times = zeros(temptime.size(),1);
    for i=0:temptime.size()-1
        times(i+1) = temptime.get(i);
    end
    timespercent = (times - times(1)) / (times(end) - times(1)) *100;
    timespercent101 = [0:1:100]';
    
    % now for each signal in the model
    numCols = tempTimeSeriesTable.getNumColumns(); % including time
    labels = tempTimeSeriesTable.getColumnLabels();
    
    tempstruct = struct();
    
    % loop through everything and get the excitations.
    for i=0:labels.size()-1
        templab = char(labels.get(i));
        % check if activation 
        if contains(string(templab), '.')
            if string(templab(length(templab)-9:end))=='excitation'
                % now select leg starting with strike or non-sided
                % actuators
                if templab(length(templab)-11)~='r'
                    % we want this one
                    tempcol = tempTimeSeriesTable.getDependentColumn(java.lang.String(templab)).getAsMat();

                    tempcolinterp = interp1(timespercent, tempcol, timespercent101);

                    if ~isfield(tempstruct, templab) % (11:length(templab)-11))
                        % fix the naming
                        tempstruct.(genvarname(templab))=[]; % (11:length(templab)-11))) = [];
                    end
                    % welkexostruct.(genvarname(templab(11:length(templab)-11))) = [welkexostruct.(genvarname(templab(11:length(templab)-11))), tempcolinterp]; 
                    tempstruct.(genvarname(templab)) = [tempstruct.(genvarname(templab)), tempcolinterp];
                end
            end
        end
    end
    
    % make a figure and plot
    
    tempfig = figure('Position',[1,1,1920,1080]);
    % do more stuff
    % averaging and whatnot
    labels = fields(tempstruct);
    for i=1:length(labels)
        subplot(8,9,i);
        templabel = char(labels(i));
        muscleplot1 = tempstruct.(genvarname(templabel));
        % plot(timespercent101, muscleplot1, ':')
        hold on;
        plot(timespercent101, mean(muscleplot1,2), 'k-', 'LineWidth', 2)
        title(templabel)
        xlabel('% gait cycle')
        ylabel('activation')
        grid on;
    end
   print(tempfig, ...
        strcat('C:\Users\JP\code\repos\Stanford\delplab\projects\muscleModel\analysis\', ...
        strcat(outputname, '.png')),...
        '-dpng', '-r500')
    disp('print 1')
end
