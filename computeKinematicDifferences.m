function [] = computeKinematicDifferences(solution)
    import org.opensim.modeling.*
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % solution = MocoTrajectory(solution);
    Time = solution.getTimeMat();
    numColPoints = solution.getNumTimes();

    workingdir = pwd;

    % try storage instead
    solutionstatestable = solution.exportToStatesTable();
    solutionNumRows = solutionstatestable.getNumRows();
    solutionLabels = solutionstatestable.getColumnLabels();
    solutionNumLabels = solutionLabels.size();
    solutionremove = [];
    % find all the labels that are not jointset
    for l = 0:solutionNumLabels-1
        templabel = solutionLabels.get(l);
        if ~contains(string(templabel), 'jointset')
            solutionremove = [solutionremove, templabel];
        end
    end
    % remove all the ones that aren't jointset
    for r = 1:length(solutionremove)
        solutionstatestable.removeColumn(solutionremove(r));
    end
    % overwrite label vector and count
    solutionLabels = solutionstatestable.getColumnLabels();
    solutionNumLabels = solutionLabels.size();
    
    
    % load the tracked or original states for comparison
    trackstatestable = TimeSeriesTable('muscle_statetrack_grfprescribe_tracked_states.sto');
    trackNumRows = trackstatestable.getNumRows();
    trackLabels = trackstatestable.getColumnLabels();
    trackNumLabels = trackLabels.size();
    trackremove = [];
    % find all the labels that are not jointset
    for l = 0:trackNumLabels-1
        templabel = trackLabels.get(l);
        if ~contains(string(templabel), 'jointset')
            trackremove = [trackremove, templabel];
        end
    end
    % remove all the ones that aren't jointset
    for r = 1:length(trackremove)
        trackstatestable.removeColumn(trackremove(r));
    end
    % overwrite label vector and count
    trackLabels = trackstatestable.getColumnLabels();
    trackNumLabels = trackLabels.size();
    
    % have a matching set of jointset values in the tables 
    % compare the two trajectories for each of the values that we care
    % about
    solutiontime = solutionstatestable.getIndependentColumn();
    solutiontime2 = [];
    for t = 0:solutiontime.size()-1
        solutiontime2 = [solutiontime2, solutiontime.get(t)];
    end
    tracktime = trackstatestable.getIndependentColumn();
    tracktime2 = [];
    for t = 0:tracktime.size()-1
        tracktime2 = [tracktime2, tracktime.get(t)];
    end
    
    
    % now want to do the same thing but plot differences between signals
    % have to figure out the trailing ends on the tracked states
    % use the time vectors to see where they overlap
    timesplit1 = length(find(tracktime2 < solutiontime2(1)));
    timesplit2 = find(tracktime2 > solutiontime2(end));
    timesplit2 = timesplit2(1);
    
    if timesplit1 == 0
        timesplit1 = 1;
        temptracktime = tracktime2(timesplit1:timesplit2);
    else
        temptracktime = tracktime2(timesplit1:timesplit2);
    end
%     temptracktime = tracktime2(timesplit1:timesplit2);
    
%     tempfig1 = figure('Position',[1,1,1920,1080]);
%     % go through all the joints
%     for l=0:trackNumLabels-1
%         templabel = trackLabels.get(l);
%         temptrack = trackstatestable.getDependentColumn(templabel);
%         temptrack2 = temptrack.getAsMat();
%         tempsolution = solutionstatestable.getDependentColumn(templabel);
%         tempsolution2 = tempsolution.getAsMat();
%         
%         % okay now what...
%         subplot(6,8,l+1);
%         plot(tracktime2, temptrack2);
%         hold on;
%         plot(solutiontime2, tempsolution2);
%         title(string(templabel));
%         xlabel('time [s]');
%         ylabel('[rad or rad/s]');
%         grid on;
%         
%         
%     end
%     
% %     subplot(6,8,l+2);
%     legend('solution','original','location','eastoutside');
%     workingdirlength = length(workingdir);
%     trial = workingdir(workingdirlength-6:workingdirlength);
%     condition = workingdir(workingdirlength-18:workingdirlength-8);
%     subject = workingdir(workingdirlength-26:workingdirlength-20);
%     
%     print(tempfig1, ...
%             strcat(strcat('C:\Users\JP\code\repos\Stanford\delplab\projects\muscleModel\analysis\', ... 
%             strcat(subject,strcat('\',strcat(string(trial),strcat('_',strcat(condition,'_kinematics_differences.png'))))))),...
%             '-dpng', '-r500')
%     disp('print 1')
%     
%     
%     
%     
%         
%     tempfig2 = figure('Position',[1,1,1920,1080]);
%     % go through all the joints
%     for l=0:trackNumLabels-1
%         templabel = trackLabels.get(l);
%         temptrack = trackstatestable.getDependentColumn(templabel);
%         temptrack2 = temptrack.getAsMat();
%         tempsolution = solutionstatestable.getDependentColumn(templabel);
%         tempsolution2 = tempsolution.getAsMat();
%         
%         % have to figure out where they overlap
%         tempsolution3 = interp1(solutiontime2, tempsolution2, temptracktime);
%         
%         % okay now what...
%         subplot(6,8,l+1);
%         plot(temptracktime, temptrack2(timesplit1:timesplit2));
%         hold on;
%         plot(temptracktime, tempsolution3);
%         title(string(templabel));
%         xlabel('time [s]');
%         ylabel('[rad or rad/s]');
%         grid on;
%     end
%     
% %     subplot(6,8,l+2);
%     legend('solution','original','location','eastoutside');
%     workingdirlength = length(workingdir);
%     trial = workingdir(workingdirlength-6:workingdirlength);
%     condition = workingdir(workingdirlength-18:workingdirlength-8);
%     subject = workingdir(workingdirlength-26:workingdirlength-20);
%     
%     print(tempfig2, ...
%             strcat(strcat('C:\Users\JP\code\repos\Stanford\delplab\projects\muscleModel\analysis\', ... 
%             strcat(subject,strcat('\',strcat(string(trial),strcat('_',strcat(condition,'_kinematics_differences_adj.png'))))))),...
%             '-dpng', '-r500')
%     disp('print 1')
    
    
    % same figure but in degrees
    tempfig3 = figure('Position',[1,1,2800,1080]);
    % go through all the joints
    for l=0:trackNumLabels-1
        templabel = trackLabels.get(l);
        temptrack = trackstatestable.getDependentColumn(templabel);
        temptrack2 = temptrack.getAsMat();
        tempsolution = solutionstatestable.getDependentColumn(templabel);
        tempsolution2 = tempsolution.getAsMat();
        
        % have to figure out where they overlap
        tempsolution3 = interp1(solutiontime2, tempsolution2, temptracktime);
        
        % okay now what...
        subplot(6,8,l+1);
        plot(temptracktime, temptrack2(timesplit1:timesplit2).*180./pi());
        hold on;
        plot(temptracktime, tempsolution3.*180./pi());
        title(string(templabel));
        xlabel('time [s]');
        ylabel('[deg or deg/s]');
        grid on;
    end
    
%     subplot(6,8,l+2);
    legend('solution','original','location','eastoutside');
    workingdirlength = length(workingdir);
    trial = workingdir(workingdirlength-6:workingdirlength);
    if strcmp(workingdir(workingdirlength-15), '\') % welkexo
        condition = workingdir(workingdirlength-14:workingdirlength-8);
        subject = workingdir(workingdirlength-22:workingdirlength-16);
    end
    if strcmp(workingdir(workingdirlength-19),'\') % welknatural
        condition = workingdir(workingdirlength-18:workingdirlength-8);
        subject = workingdir(workingdirlength-26:workingdirlength-20);
    end
    
    
    temptarget = strcat(pwd,'\..\..\..\..\analysis\',subject,'\',string(trial),'_',condition,'_kinematics_differences_adj_deg');
    print(tempfig3, temptarget, '-dpng', '-r500')
    disp('print 1')


    % now a figure that actually takes the differences between them, or the
    % deviation in the tracking problem from the IK angles
    % same figure but in degrees
    tempfig4 = figure('Position',[1,1,2800,1080]);
    % go through all the joints
    for l=0:trackNumLabels-1
        templabel = trackLabels.get(l);
        temptrack = trackstatestable.getDependentColumn(templabel);
        temptrack2 = temptrack.getAsMat();
        tempsolution = solutionstatestable.getDependentColumn(templabel);
        tempsolution2 = tempsolution.getAsMat();
        
        % have to figure out where they overlap
        tempsolution3 = interp1(solutiontime2, tempsolution2, temptracktime);
        
        % okay now what...
        subplot(6,8,l+1);
        plot(temptracktime, ((tempsolution3.*180./pi()) - (temptrack2(timesplit1:timesplit2).*180./pi())')');
        hold on;
%         plot(temptracktime, tempsolution3.*180./pi());
        title(string(templabel));
        xlabel('time [s]');
        ylabel('[deg or deg/s]');
        grid on;
    end
    
%     subplot(6,8,l+2);
    legend('solution-original','location','eastoutside');
    workingdirlength = length(workingdir);
    if strcmp(workingdir(workingdirlength-15), '\') % welkexo
        condition = workingdir(workingdirlength-14:workingdirlength-8);
        subject = workingdir(workingdirlength-22:workingdirlength-16);
    end
    if strcmp(workingdir(workingdirlength-19),'\') % welknatural
        condition = workingdir(workingdirlength-18:workingdirlength-8);
        subject = workingdir(workingdirlength-26:workingdirlength-20);
    end
    
    temptarget2 = strcat(pwd,'\..\..\..\..\analysis\',subject,'\',string(trial),'_',condition,'_kinematics_differences_adj_deg_diff');
    print(tempfig4, temptarget2, '-dpng', '-r500')
    disp('print 1')

    close all
end
