% written by Jon Stingel
% 20211022
import org.opensim.modeling.*
% repodir = 'G:\Shared drives\Exotendon\muscleModel\muscleEnergyModel';
repodir = "C:\Users\jonstingel\code\muscleModel\muscleEnergyModel";
resultsdir = strcat(repodir, '/../results');
cd(resultsdir)
exocolor = '#AB82FF'
natcolor = '#FF7F00'
% conditions
% walsconditions = ['walsslack','walslow','walsmed','walshigh','walsmax']
% jackconditions = ['jackpower1','jackpower2','jackpower3','jackpower4','jackpower5','jackpower6',
%                   'jacktau1','jacktau2','jacktau3','jacktau4','jacktau5']
% dembconditions = ['dembnoloadfree','dembnoloadslow','dembloadedfree','dembloadedmatched']
% sildconditions = ['sildbw0','sildbw5','sildbw10','sild10w0','sild10w5','sild10w10',
%                   'sild20w0','sild20w5','sild20w10','sild30w0','sild30w5','sild30w10',
%                   'sildbwrun0','sild10wrun0','sild20wrun0','sild30wrun0']

%%%%% - remember to only put in the exo conditions that you are looking to see the reductions from
% dembconditions = {'dembnoloadfree', 'dembloadedfree'}; %
% dembsubjects = {'demb010','demb011','demb012','demb014', 'demb005','demb007','demb009'}; %
welkexoconditions = {'welkexo'}; % ,'welkexoexo'}; % ,'welknaturalslow','welknaturalnatural', ...
                  % 'welknaturalexo','welkexonatural','welkexoexo','welkexofast'};
welknaturalconditions = {'welknatural'};% ,'welknaturalnatural'};
welksubjects = {'welk002','welk003','welk005','welk008','welk009','welk010','welk013'};

thingstoplot = {'netjointmoments'}; % 'probes', 'shortening', 'mechanical', 'activation'

load 'G:\Shared drives\Exotendon\muscleModel\muscleEnergyModel\subjectgaitcycles.mat';
load 'G:\Shared drives\Exotendon\muscleModel\muscleEnergyModel\subjectmass.mat';

% loop through each of the things we want to plot
for thing=1:length(thingstoplot)
    tempthing = char(thingstoplot(thing))

    % create stucture for combined subject figures
    welknaturalstruct_combine = struct();
    welkexostruct_combine = struct();

    allsubj_nat = struct();
    allsubj_exo = struct();

    meansall_nat = struct();
    meansall_exo = struct();
    stdsall_nat = struct();
    stdsall_exo = struct();

    % loop through the subjects
    for subj=1:length(welksubjects)
        subject = char(welksubjects(subj));
        subjdir = strcat(resultsdir, strcat('/',subject));
        
        % create the struct for individual figures
        welknaturalstruct = struct();
        welkexostruct = struct();
    

        % loop through conditions - exo first
        for cond=1:length(welkexoconditions)
            condition = char(welkexoconditions(cond));
            conddir = strcat(subjdir, strcat('/',condition));
            trials = fieldnames(subjectgaitcycles.(genvarname(subject)).(genvarname(condition)));
            % loop the trials
            for trial=1:length(trials)
                % what do we actually want to do here
                test = char(trials(trial));
                trialdir = strcat(conddir, strcat('/',test));
                cd(trialdir)
                disp(trialdir)

                % now figure out how to get and plot the signal i want
                % have all the muscle analysis files already
                % do I want to do average or individual?
                tempfile = strcat(trialdir, '/muscle_joint_moment_breakdown', '.sto');
                tempTimeSeriesTable = TimeSeriesTable(tempfile);
                temptime = tempTimeSeriesTable.getIndependentColumn();
                times = zeros(temptime.size(),1);
                for i=0:temptime.size()-1
                    times(i+1) = temptime.get(i);
                end
                timespercent = (times - times(1)) / (times(end) - times(1)) *100;
                timespercent101 = [0:1:100]';
                welkexostruct.time = timespercent101;

                
                % now for each of the things
                numCols = tempTimeSeriesTable.getNumColumns(); % including time
                labels = tempTimeSeriesTable.getColumnLabels();            

                for i=0:labels.size()-1
                    coord = char(labels.get(i));
                        
                    tempcol = tempTimeSeriesTable.getDependentColumn(java.lang.String(coord)).getAsMat();
                    tempcolinterp = interp1(timespercent, tempcol, timespercent101);
                    if ~isfield(welkexostruct, coord)
                        welkexostruct.(genvarname(coord)) = [];
                    end
                    welkexostruct.(genvarname(coord)) = [welkexostruct.(genvarname(coord)), tempcolinterp];
                end
            end
        end
        % done with the exo conditions
        
        % loop through conditions - now for the natural
        for cond=1:length(welknaturalconditions)
           condition = char(welknaturalconditions(cond));
           conddir = strcat(subjdir, strcat('/',condition));
           trials = fieldnames(subjectgaitcycles.(genvarname(subject)).(genvarname(condition)));
            % loop the trials
            for trial=1:length(trials)
                % what do we actually want to do here
                test = char(trials(trial));
                trialdir = strcat(conddir, strcat('/',test));
                cd(trialdir)
                disp(trialdir)

                % now figure out how to get and plot the signal i want
                % have all the muscle analysis files already
                % do I want to do average or individual?
                tempfile = strcat(trialdir, '/muscle_joint_moment_breakdown', '.sto');
                tempTimeSeriesTable = TimeSeriesTable(tempfile);
                temptime = tempTimeSeriesTable.getIndependentColumn();
                times = zeros(temptime.size(),1);
                for i=0:temptime.size()-1
                    times(i+1) = temptime.get(i);
                end
                timespercent = (times - times(1)) / (times(end) - times(1)) *100;
                timespercent101 = [0:1:100]';
                welknaturalstruct.time = timespercent101;

                % now for each of the things
                numCols = tempTimeSeriesTable.getNumColumns(); % including time
                labels = tempTimeSeriesTable.getColumnLabels();            

                
                for i=0:labels.size()-1
                    coord = char(labels.get(i));
                    
                    % we also want the whole body measure
                    tempcol = tempTimeSeriesTable.getDependentColumn(java.lang.String(coord)).getAsMat();
                    tempcolinterp = interp1(timespercent, tempcol, timespercent101);
                    if ~isfield(welknaturalstruct, coord)
                        welknaturalstruct.(genvarname(coord)) = [];
                    end
                    welknaturalstruct.(genvarname(coord)) = [welknaturalstruct.(genvarname(coord)), tempcolinterp];
                end
            end
        end


        newlabels = fields(welkexostruct);
        
        % okay now to plot etc. 
        
        % need to redo the labels
%         tempfig = figure('Position',[1,1,1920,1080]);
%         % do more stuff
%         % averaging and whatnot
%         for i=2:length(newlabels)
%             subplot(5,6,i-1);
%             templabel = newlabels(i);
%             templabel = char(templabel);
%             % plot each of the gait cycles
%             muscleplot_nat = welknaturalstruct.(genvarname(char(templabel)));
%             muscleplot_exo = welkexostruct.(genvarname(char(templabel)));
%             plot(welknaturalstruct.time, muscleplot_nat, 'r:')
%             hold on;
%             plot(welkexostruct.time, muscleplot_exo, 'b:')
%             % plot the subject average
%             plot(welknaturalstruct.time, mean(muscleplot_nat,2), 'r-', 'LineWidth', 1)
%             plot(welkexostruct.time, mean(muscleplot_exo,2), 'b-', 'LineWidth', 1)
%             ylabel('Moment [Nm]');
%             % title(templabel)
%             xlabel('% gait cycle')
%             % select the name out 
%             templabel2 = strrep(templabel,'_',' ');
%             % for activation maintenance rate
%             % templabel2 = templabel()
%             title(templabel2);
%             grid on;
%         end
        
%         print(tempfig, ...
%             strcat(strcat('G:\Shared drives\Exotendon\muscleModel\analysis\', strcat(subject,'\')), strcat(strcat(tempthing, '_acrossconditions'), '.png')),...
%             '-dpng', '-r500')
%         disp('print 1')
        


        % add the subject average to the combined struct?
        welknaturalstruct_combine.(genvarname(subject)) = welknaturalstruct;
        welkexostruct_combine.(genvarname(subject)) = welkexostruct;
        
    end

    keyboard
    % loop through the subjects again?
    markr = {'r:','r--'};
    markb = {'b:','b--'};
    
    % now plot across subjects
    tempfig2 = figure('Position',[1,1,1280,1920]);
    
    % then loop through the muscles inside each subject

    for i=2:length(newlabels)
        subplot(5,6,i-1);
        templabel = newlabels(i);
        templabel = char(templabel);
        holding_nat = [];
        holding_exo = [];
        % loop through the subjects
        for subj=1:length(welksubjects)
            subject = char(welksubjects(subj));
            model_mass = subjectmass.(genvarname(subject));
            % disp(model_mass)
            muscleplot_nat = welknaturalstruct_combine.(genvarname(subject)).(genvarname(char(templabel)));
            muscleplot_exo = welkexostruct_combine.(genvarname(subject)).(genvarname(char(templabel)));
            temp_nat = mean(muscleplot_nat,2);
            temp_nat_norm = temp_nat./model_mass;
            temp_exo = mean(muscleplot_exo,2);
            temp_exo_norm = temp_exo./model_mass;
            holding_nat = [holding_nat, temp_nat_norm];
            holding_exo = [holding_exo, temp_exo_norm];
            % have all of them, want the average plotted for each subject
            hold on;
            plot(welknaturalstruct.time, mean(muscleplot_nat,2), 'r:','LineWidth',1);% char(markr(subj)))
            plot(welkexostruct.time, mean(muscleplot_exo,2), 'b:','LineWidth',1);% char(markb(subj)))
            % plot(welknaturalstruct.time, temp_nat_norm,'r:','LineWidth',.4);
            % plot(welkexostruct.time, temp_exo_norm,'b:','LineWidth',.4);
            
            if ~isfield(allsubj_nat, templabel)
                allsubj_nat.(genvarname(templabel)) = [];
                allsubj_exo.(genvarname(templabel)) = [];
            end
            allsubj_nat.(genvarname(templabel)) = [allsubj_nat.(genvarname(templabel)), temp_nat_norm];
            allsubj_exo.(genvarname(templabel)) = [allsubj_exo.(genvarname(templabel)), temp_exo_norm];

        end
        plot(welknaturalstruct.time, mean(holding_nat,2),'Color', natcolor,'LineWidth',2);
        plot(welkexostruct.time, mean(holding_exo,2),'Color', exocolor,'LineWidth',2);
        % legend(num2str(max(mean(holding_nat,2))),num2str(max(mean(holding_exo,2))));
        % 
        % note that we can do max or min for flexions/extensions
        % disp(strcat('nat: ', templabel))
        % maxnat = min(holding_nat)
        % natmaxavg = mean(maxnat)
        % natsd = std(maxnat)
        % natse = natsd/sqrt(length(maxnat))
        % 
        % disp(strcat('exo: ',templabel))
        % maxexo = min(holding_exo)
        % exomaxavg = mean(maxexo)
        % exosd = std(maxexo)
        % exose = exosd/sqrt(length(maxexo))
        % 
        % doing differences
        % note that we can do max or min for flexions/extensions
        disp(strcat('posdiff: ', templabel))
        posdiff = max(holding_nat) - max(holding_exo)
        posdiffavg = mean(posdiff)
        posdiffsd = std(posdiff)
        posdiffse = posdiffsd/sqrt(length(posdiff))

        disp(strcat('mindiff: ',templabel))
        mindiff = min(holding_nat) - min(holding_exo)
        mindiffavg = mean(mindiff)
        mindiffsd = std(mindiff)
        mindiffse = mindiffsd/sqrt(length(mindiff))


        templabel2 = strrep(templabel,'_',' ');
        title(templabel2)
        
        xlabel('% gait cycle')
        ylabel('Moment [Nm/kg]')
%         grid on;
    end
%     subplot(9,3,i);
%     plot(welkexostruct.time, mean(holding_exo,2),'b');
%     hold on;
%     plot(welknaturalstruct.time, mean(holding_nat,2),'r');
    
    print(tempfig2, ...
        strcat('G:\Shared drives\Exotendon\muscleModel\analysis\', tempthing, '_combined_subjectAVG', '.png'),...
        '-dpng', '-r500')
    disp('print 2')

    % now need to loop through all the labels and get the across subject
    % sample means and standard deviations. 
    for i=2:length(newlabels)
        templabel = newlabels(i);
        templabel = char(templabel);
        if ~isfield(meansall_exo, templabel)
            meansall_nat.(genvarname(templabel)) = [];
            meansall_exo.(genvarname(templabel)) = [];
            stdsall_nat.(genvarname(templabel)) = [];
            stdsall_exo.(genvarname(templabel)) = [];
        end
        meansall_nat.(genvarname(templabel)) = [meansall_nat.(genvarname(templabel)), mean(allsubj_nat.(genvarname(templabel)),2)];
        meansall_exo.(genvarname(templabel)) = [meansall_exo.(genvarname(templabel)), mean(allsubj_exo.(genvarname(templabel)),2)];
        stdsall_nat.(genvarname(templabel)) = [stdsall_nat.(genvarname(templabel)), std(allsubj_nat.(genvarname(templabel)),0,2)];
        stdsall_exo.(genvarname(templabel)) = [stdsall_exo.(genvarname(templabel)), std(allsubj_exo.(genvarname(templabel)),0,2)];
    end
    
    cd("G:\Shared drives\Exotendon\muscleModel\muscleEnergyModel\")

    % save the data from the means and stds
    meansall_nat_table = struct2table(meansall_nat);
    meansall_exo_table = struct2table(meansall_exo);
    stdsall_nat_table = struct2table(stdsall_nat);
    stdsall_exo_table = struct2table(stdsall_exo);

    writetable(meansall_nat_table, 'meansall_moments_nat.csv');
    writetable(meansall_exo_table, 'meansall_moments_exo.csv');
    writetable(stdsall_nat_table, 'stdsall_moments_nat.csv');
    writetable(stdsall_exo_table, 'stdsall_moments_exo.csv');

end
