% written by Jon Stingel
% 20210329
% gather and plot all the muscle activities from simulations. 
import org.opensim.modeling.*
repodir = 'G:\Shared drives\Exotendon\muscleModel\muscleEnergyModel';
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
welkexoconditions = {'welkexo'};%,'welkexoexo'}; % ,'welknaturalslow','welknaturalnatural', ...
                  % 'welknaturalexo','welkexonatural','welkexoexo','welkexofast'};
welknaturalconditions = {'welknatural'};%,'welknaturalnatural'};
% welksubjects = {'welk002','welk003','welk005','welk007','welk008','welk009','welk010','welk013'};
welksubjects = {'welk002','welk003','welk005','welk008','welk009','welk010','welk013'}; % 'welk008
tag = 'muscletrack';
thingstoplot = {'excitation','activation'};

load 'G:\Shared drives\Exotendon\muscleModel\muscleEnergyModel\subjectgaitcycles.mat';



exomeans_excitation = struct();
naturalmeans_excitation = struct();
exomeans_activation = struct();
naturalmeans_activation = struct();

exopeaks_excitation = struct();
exopeaks_activation = struct();
naturalpeaks_excitation = struct();
naturalpeaks_activation = struct();


nat_activitygrouped = struct();
exo_activitygrouped = struct();

nat_maxiso = struct();
exo_maxiso = struct();

 
% loop through the subjects
for subj=1:length(welksubjects)
    subject = char(welksubjects(subj));
    subjdir = strcat(resultsdir, strcat('/',subject));
    
    exomeans_excitation.(genvarname(subject)) = [];
    naturalmeans_excitation.(genvarname(subject)) = [];
    exomeans_activation.(genvarname(subject)) = [];
    naturalmeans_activation.(genvarname(subject)) = [];
    
    exopeaks_excitation.(genvarname(subject)) = [];
    exopeaks_activation.(genvarname(subject)) = [];
    naturalpeaks_excitation.(genvarname(subject)) = [];
    naturalpeaks_activation.(genvarname(subject)) = [];


    %idk 
    nat_maxiso.(genvarname(subject)) = [];
    exo_maxiso.(genvarname(subject)) = [];



    
    % loop through each of the things we want to plot
    for thing=1:length(thingstoplot)
        tempthing = char(thingstoplot(thing))
        welknaturalstruct = struct();
        welkexostruct = struct();

        welknat_isostruct = struct();
        welkexo_isostruct = struct();


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
                % disp(trialdir)

                % now figure out how to get and plot the signal i want
                % have all the muscle analysis files already
                % do I want to do average or individual?
                
                if tempthing == 'activation'
                    % do something for the activations
                    % disp('getting activations...')
                    
                    tempfile = strcat(trialdir, strcat('/',tag,'_states_100con.sto'));
                    % disp(tempfile)
                    % if strcmp(subject,'welk002') || strcmp(subject,'welk003')
                    %     tempfile = strcat(trialdir,'/muscleprescribe_states.sto');
                    % else
                    %     tempfile = strcat(trialdir, strcat('/',tag,'_states.sto'));
                    % end
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
                    
                    % get the model so that we can grab all the max iso forces
                    model = Model(strcat(trialdir,'/simple_model_all_the_probes_adjusted.osim'));
                    modelmuscles = model.getMuscles();

                    % loop through everything and get the activations.
                    for i=0:labels.size()-1
                        templab = char(labels.get(i));
                        % check if it is an activation for the right leg
                        if string(templab(length(templab)-9:end)) == 'activation' && templab(length(templab)-11) == 'r'
                            % we want the activations - doing right leg for
                            % simplicity
                            tempcol = tempTimeSeriesTable.getDependentColumn(java.lang.String(templab)).getAsMat();
                            tempcolinterp = interp1(timespercent, tempcol, timespercent101);
                            
                            % get max iso force as well for muscle averaging.
                            % templab;
                            % get the forcepath name
                            musclename = templab(11:length(templab)-11);
                            tempmusc = modelmuscles.get(java.lang.String(musclename));
                            tempmaxiso = tempmusc.getMaxIsometricForce();

                            % add to the temp struct
                            if ~isfield(welkexostruct, templab(11:length(templab)-11))
                                % fix the naming
                                welkexostruct.(genvarname(templab(11:length(templab)-11))) = [];
                                welkexo_isostruct.(genvarname(templab(11:length(templab)-11))) = [];
                            end
                            welkexostruct.(genvarname(templab(11:length(templab)-11))) = [welkexostruct.(genvarname(templab(11:length(templab)-11))), tempcolinterp]; 
                            welkexo_isostruct.(genvarname(templab(11:length(templab)-11))) = [welkexo_isostruct.(genvarname(templab(11:length(templab)-11))), tempmaxiso];
                        end
                    end
                end

                % now for the case of excitations
                if tempthing == 'excitation'
                    % do something for the excitations
                    % disp('getting excitations...')
                    
                    tempfile = strcat(trialdir, strcat('/',tag,'_controls_100con.sto'));
                    % disp(tempfile)
                    % if strcmp(subject,'welk002') || strcmp(subject,'welk003')
                    %     tempfile = strcat(trialdir,'/muscleprescribe_controls.sto');
                    % else
                    %     tempfile = strcat(trialdir, strcat('/',tag,'_controls.sto'));
                    % end
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
                    
                    % loop through everything and get the activations.
                    for i=0:labels.size()-1
                        templab = char(labels.get(i));
                        % check if it is an activation for the right leg                        
                        if string(templab(length(templab)-1:end)) ~= '_l'
                            % we want the controls - doing right leg for
                            % simplicity
                            tempcol = tempTimeSeriesTable.getDependentColumn(java.lang.String(templab)).getAsMat();
                            tempcolinterp = interp1(timespercent, tempcol, timespercent101);
                            
                            if ~isfield(welkexostruct, templab(11:end))
                                % fix the naming
                                welkexostruct.(genvarname(templab(11:end))) = [];
                            end
                            welkexostruct.(genvarname(templab(11:end))) = [welkexostruct.(genvarname(templab(11:end))), tempcolinterp]; 
                        end
                    end                    
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
                % disp(trialdir)

                % now figure out how to get and plot the signal i want
                % have all the muscle analysis files already
                % do I want to do average or individual?
                
                if tempthing == 'activation'
                    % do something for the activations
                    % disp('getting activations...')
                    
                    tempfile = strcat(trialdir, strcat('/',tag,'_states_100con.sto'));

                    % if strcmp(subject,'welk002') || strcmp(subject,'welk003')
                    %     tempfile = strcat(trialdir,'/muscleprescribe_states.sto');
                    % else
                    %     tempfile = strcat(trialdir, strcat('/',tag,'_states.sto'));
                    % end
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
                    
                    % shouldn't need this but adding anyway
                    % get the model so that we can grab all the max iso forces
                    model = Model(strcat(trialdir,'/simple_model_all_the_probes_adjusted.osim'));
                    modelmuscles = model.getMuscles();


                    % loop through everything and get the activations.
                    for i=0:labels.size()-1
                        templab = char(labels.get(i));
                        % check if it is an activation for the right leg
                        if string(templab(length(templab)-9:end)) == 'activation' && templab(length(templab)-11) == 'r'
                            % we want the activations - doing right leg for
                            % simplicity
                            tempcol = tempTimeSeriesTable.getDependentColumn(java.lang.String(templab)).getAsMat();
                            tempcolinterp = interp1(timespercent, tempcol, timespercent101);
                            
                            % get max iso force as well for muscle averaging.
                            % templab;
                            % get the forcepath name
                            musclename = templab(11:length(templab)-11);
                            tempmusc = modelmuscles.get(java.lang.String(musclename));
                            tempmaxiso = tempmusc.getMaxIsometricForce();

                            if ~isfield(welknaturalstruct, templab(11:length(templab)-11))
                                % fix the naming
                                welknaturalstruct.(genvarname(templab(11:length(templab)-11))) = [];
                                welknat_isostruct.(genvarname(templab(11:length(templab)-11))) = [];
                            end
                            welknaturalstruct.(genvarname(templab(11:length(templab)-11))) = [welknaturalstruct.(genvarname(templab(11:length(templab)-11))), tempcolinterp]; 
                            welknat_isostruct.(genvarname(templab(11:length(templab)-11))) = [welknat_isostruct.(genvarname(templab(11:length(templab)-11))), tempmaxiso];

                        end
                    end
                end

                % now for the case of excitations
                if tempthing == 'excitation'
                    % do something for the excitations
                    % disp('getting excitations...')
                    
                    tempfile = strcat(trialdir, strcat('/',tag,'_controls_100con.sto'));

                    % if strcmp(subject,'welk002') || strcmp(subject,'welk003')
                    %     tempfile = strcat(trialdir,'/muscleprescribe_controls.sto');
                    % else
                    %     tempfile = strcat(trialdir, strcat('/',tag,'_controls.sto'));
                    % end
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
                    
                    % loop through everything and get the activations.
                    for i=0:labels.size()-1
                        templab = char(labels.get(i));
                        % check if it is an activation for the right leg                        
                        if string(templab(length(templab)-1:end)) ~= '_l'
                            % we want the controls - doing right leg for
                            % simplicity
                            tempcol = tempTimeSeriesTable.getDependentColumn(java.lang.String(templab)).getAsMat();
                            tempcolinterp = interp1(timespercent, tempcol, timespercent101);
                            
                            if ~isfield(welknaturalstruct, templab(11:end))
                                % fix the naming
                                welknaturalstruct.(genvarname(templab(11:end))) = [];
                            end
                            welknaturalstruct.(genvarname(templab(11:end))) = [welknaturalstruct.(genvarname(templab(11:end))), tempcolinterp]; 
                        end
                    end                    
                end
            end
        end
        
        % now need to loop through both natural and exo to find the 3 glutes
        disp('can do weighted avgs here')

        if tempthing == 'activation'
            % create maps of all the muscles that we want to grab out 
            quadriceps = {'recfem_r','vasint_r','vaslat_r','vasmed_r'};
            hipflexors = {'grac_r','iliacus_r','psoas_r','sart_r','tfl_r'};
            hipabductors = {'glmed1_r','glmed2_r','glmed3_r','piri_r','glmin1_r','glmin2_r','glmin3_r'};
            hamstrings = {'bflh_r','bfsh_r','semiten_r','semimem_r'};
            hipadductors = {'addbrev_r','addlong_r','addmagDist_r','addmagIsch_r','addmagMid_r','addmagProx_r'};
            hipextensors = {'glmax1_r','glmax2_r','glmax3_r'};
            plantarflexors = {'fdl_r','fhl_r','gaslat_r','gasmed_r','soleus_r','perlong_r','perbrev_r','tibpost_r'};
            dorsiflexors = {'tibant_r','edl_r','ehl_r'};

            % create empty data 
            quadriceps_nat_data = [];
            hipflexors_nat_data = [];
            hipabductors_nat_data = [];
            hamstrings_nat_data = [];
            hipadductors_nat_data = [];
            hipextensors_nat_data = [];
            plantarflexors_nat_data = [];
            dorsiflexors_nat_data = [];

            quadriceps_exo_data = [];
            hipflexors_exo_data = [];
            hipabductors_exo_data = [];
            hamstrings_exo_data = [];
            hipadductors_exo_data = [];
            hipextensors_exo_data = [];
            plantarflexors_exo_data = [];
            dorsiflexors_exo_data = [];

            labels_iso = fields(welknat_isostruct);
            quadcount = 0;
            quadholder_nat = zeros(101,4);
            quadholder_exo = zeros(101,4);
            hipflexcount = 0;
            hipflexholder_nat = zeros(101,4);
            hipflexholder_exo = zeros(101,4);
            hipabdcount = 0;
            hipabdholder_nat = zeros(101,4);
            hipabdholder_exo = zeros(101,4);
            hamstringcount = 0;
            hamstringholder_nat = zeros(101,4);
            hamstringholder_exo = zeros(101,4);
            hipaddcount = 0;
            hipaddholder_nat = zeros(101,4);
            hipaddholder_exo = zeros(101,4);
            hipextcount = 0;
            hipextholder_nat = zeros(101,4);
            hipextholder_exo = zeros(101,4);
            pfcount = 0;
            pfholder_nat = zeros(101,4);
            pfholder_exo = zeros(101,4);
            dfcount = 0;
            dfholder_nat = zeros(101,4);
            dfholder_exo = zeros(101,4);
            
            % can just multiply each muscle activation by the max iso for that subject, 
            for o=1:length(labels_iso)
                templabel_iso = string(labels_iso(o));
                
                
                if any(strcmp(quadriceps, templabel_iso))
                    tempquadact_nat = welknaturalstruct.(genvarname(templabel_iso));
                    tempquadact_exo = welkexostruct.(genvarname(templabel_iso));
                    getisoforce = welknat_isostruct.(genvarname(templabel_iso))(1);
                    % have the curve and the iso - want to multiply and keep a running sum of the multiplier
                    newquad_nat = tempquadact_nat.*getisoforce;
                    quadholder_nat = quadholder_nat + newquad_nat;
                    newquad_exo = tempquadact_exo.*getisoforce;
                    quadholder_exo = quadholder_exo + newquad_exo;

                    quadcount = quadcount + getisoforce;
                end
                if any(strcmp(hipflexors, templabel_iso))
                    temphipflexact_nat = welknaturalstruct.(genvarname(templabel_iso));
                    temphipflexact_exo = welkexostruct.(genvarname(templabel_iso));
                    getisoforce = welknat_isostruct.(genvarname(templabel_iso))(1);
                    % have the curve and the iso - want to multiply and keep a running sum of the multiplier
                    newhipflex_nat = temphipflexact_nat.*getisoforce;
                    hipflexholder_nat = hipflexholder_nat + newhipflex_nat;
                    newhipflex_exo = temphipflexact_exo.*getisoforce;
                    hipflexholder_exo = hipflexholder_exo + newhipflex_exo;
                    
                    hipflexcount = hipflexcount + getisoforce;
                end
                if any(strcmp(hipabductors, templabel_iso))
                    temphipabdact_nat = welknaturalstruct.(genvarname(templabel_iso));
                    temphipabdact_exo = welkexostruct.(genvarname(templabel_iso));
                    getisoforce = welknat_isostruct.(genvarname(templabel_iso))(1);
                    % have the curve and the iso - want to multiply and keep a running sum of the multiplier
                    newhipabd_nat = temphipabdact_nat.*getisoforce;
                    hipabdholder_nat = hipabdholder_nat + newhipabd_nat;
                    newhipabd_exo = temphipabdact_exo.*getisoforce;
                    hipabdholder_exo = hipabdholder_exo + newhipabd_exo;

                    hipabdcount = hipabdcount + getisoforce;
                end
                if any(strcmp(hamstrings, templabel_iso))
                    temphamstringact_nat = welknaturalstruct.(genvarname(templabel_iso));
                    temphamstringact_exo = welkexostruct.(genvarname(templabel_iso));
                    getisoforce = welknat_isostruct.(genvarname(templabel_iso))(1);
                    % have the curve and the iso - want to multiply and keep a running sum of the multiplier
                    newhamstring_nat = temphamstringact_nat.*getisoforce;
                    hamstringholder_nat = hamstringholder_nat + newhamstring_nat;
                    newhamstring_exo = temphamstringact_exo.*getisoforce;
                    hamstringholder_exo = hamstringholder_exo + newhamstring_exo;

                    hamstringcount = hamstringcount + getisoforce;
                end
                if any(strcmp(hipadductors, templabel_iso))
                    temphipaddact_nat = welknaturalstruct.(genvarname(templabel_iso));
                    temphipaddact_exo = welkexostruct.(genvarname(templabel_iso));
                    getisoforce = welknat_isostruct.(genvarname(templabel_iso))(1);
                    % have the curve and the iso - want to multiply and keep a running sum of the multiplier
                    newhipadd_nat = temphipaddact_nat.*getisoforce;
                    hipaddholder_nat = hipaddholder_nat + newhipadd_nat;
                    newhipadd_exo = temphipaddact_exo.*getisoforce;
                    hipaddholder_exo = hipaddholder_exo + newhipadd_exo;
                    
                    hipaddcount = hipaddcount + getisoforce;
                end
                if any(strcmp(hipextensors, templabel_iso))
                    temphipextact_nat = welknaturalstruct.(genvarname(templabel_iso));
                    temphipextact_exo = welkexostruct.(genvarname(templabel_iso));
                    getisoforce = welknat_isostruct.(genvarname(templabel_iso))(1);
                    % have the curve and the iso - want to multiply and keep a running sum of the multiplier
                    newhipext_nat = temphipextact_nat.*getisoforce;
                    hipextholder_nat = hipextholder_nat + newhipext_nat;
                    newhipext_exo = temphipextact_exo.*getisoforce;
                    hipextholder_exo = hipextholder_exo + newhipext_exo;
                    
                    hipextcount = hipextcount + getisoforce;
                end
                if any(strcmp(plantarflexors, templabel_iso))
                    temppfact_nat = welknaturalstruct.(genvarname(templabel_iso));
                    temppfact_exo = welkexostruct.(genvarname(templabel_iso));
                    getisoforce = welknat_isostruct.(genvarname(templabel_iso))(1);
                    % have the curve and the iso - want to multiply and keep a running sum of the multiplier
                    newpf_nat = temppfact_nat.*getisoforce;
                    pfholder_nat = pfholder_nat + newpf_nat;
                    newpf_exo = temppfact_exo.*getisoforce;
                    pfholder_exo = pfholder_exo + newpf_exo;
                    
                    pfcount = pfcount + getisoforce;
                end
                if any(strcmp(dorsiflexors, templabel_iso))
                    tempdfact_nat = welknaturalstruct.(genvarname(templabel_iso));
                    tempdfact_exo = welkexostruct.(genvarname(templabel_iso));
                    getisoforce = welknat_isostruct.(genvarname(templabel_iso))(1);
                    % have the curve and the iso - want to multiply and keep a running sum of the multiplier
                    newdf_nat = tempdfact_nat.*getisoforce;
                    dfholder_nat = dfholder_nat + newdf_nat;
                    newdf_exo = tempdfact_exo.*getisoforce;
                    dfholder_exo = dfholder_exo + newdf_exo;
                    
                    dfcount = dfcount + getisoforce;
                end
            end
            
            % need to divide by count for weighted average
            quadriceps_exo_data = mean(quadholder_exo./quadcount, 2);
            quadriceps_nat_data = mean(quadholder_nat./quadcount, 2);

            hipflexors_nat_data = mean(hipflexholder_nat./hipflexcount, 2);
            hipflexors_exo_data = mean(hipflexholder_exo./hipflexcount, 2);

            hipabductors_nat_data = mean(hipabdholder_nat./hipabdcount, 2);
            hipabductors_exo_data = mean(hipabdholder_exo./hipabdcount, 2);

            hamstrings_nat_data = mean(hamstringholder_nat./hamstringcount, 2);
            hamstrings_exo_data = mean(hamstringholder_exo./hamstringcount, 2);

            hipadductors_nat_data = mean(hipaddholder_nat./hipaddcount, 2);
            hipadductors_exo_data = mean(hipaddholder_exo./hipaddcount, 2);

            hipextensors_nat_data = mean(hipextholder_nat./hipextcount, 2);
            hipextensors_exo_data = mean(hipextholder_exo./hipextcount, 2);

            plantarflexors_nat_data = mean(pfholder_nat./pfcount, 2);
            plantarflexors_exo_data = mean(pfholder_exo./pfcount, 2);

            dorsiflexors_nat_data = mean(dfholder_nat./dfcount, 2);
            dorsiflexors_exo_data = mean(dfholder_exo./dfcount, 2);

            % now add them to the struct so that they are with the rest of
            % the fields
            welknaturalstruct.quadriceps = quadriceps_nat_data;
            welknaturalstruct.hipflexors = hipflexors_nat_data;
            welknaturalstruct.hipabductors = hipabductors_nat_data;
            welknaturalstruct.hamstrings = hamstrings_nat_data;
            welknaturalstruct.hipadductors = hipadductors_nat_data;
            welknaturalstruct.hipextensors = hipextensors_nat_data;
            welknaturalstruct.plantarflexors = plantarflexors_nat_data;
            welknaturalstruct.dorsiflexors = dorsiflexors_nat_data;

            welkexostruct.quadriceps = quadriceps_exo_data;
            welkexostruct.hipflexors = hipflexors_exo_data;
            welkexostruct.hipabductors = hipabductors_exo_data;
            welkexostruct.hamstrings = hamstrings_exo_data;
            welkexostruct.hipadductors = hipadductors_exo_data;
            welkexostruct.hipextensors = hipextensors_exo_data;
            welkexostruct.plantarflexors = plantarflexors_exo_data;
            welkexostruct.dorsiflexors = dorsiflexors_exo_data;
            % and update the fields variable so that they are included.
            testlabels_nat = fields(welknaturalstruct);
            testlabels_exo = fields(welkexostruct);       

        end



        labels_nat = fields(welknaturalstruct);
        glutemax = {'glmax1_r','glmax2_r','glmax3_r'};
        glutemed = {'glmed1_r','glmed2_r','glmed3_r'};
        glutemin = {'glmin1_r','glmin2_r','glmin3_r'};
        
        glutemax_data_nat = [];
        glutemed_data_nat = [];
        glutemin_data_nat = [];
        
        glutemax_data_exo = [];
        glutemed_data_exo= [];
        glutemin_data_exo = [];
        
        % loop the naturals first
        for i=1:length(labels_nat)
            templabel_nat = string(labels_nat(i));
            if any(strcmp(glutemax, templabel_nat))
                tempglute = welknaturalstruct.(genvarname(templabel_nat));
                glutemax_data_nat = [glutemax_data_nat, tempglute];
            end
            if any(strcmp(glutemed, templabel_nat))
                tempglute = welknaturalstruct.(genvarname(templabel_nat));
                glutemed_data_nat = [glutemed_data_nat, tempglute];
            end
            if any(strcmp(glutemin, templabel_nat))
                tempglute = welknaturalstruct.(genvarname(templabel_nat));
                glutemin_data_nat = [glutemin_data_nat, tempglute];
            end
        end
        glutemax_data_nat = mean(glutemax_data_nat, 2);
        glutemed_data_nat = mean(glutemed_data_nat, 2);
        glutemin_data_nat = mean(glutemin_data_nat, 2);
        
        labels_exo = fields(welkexostruct);
        % loop the exos now
        for i=1:length(labels_exo)
            templabel_exo = string(labels_exo(i));
            if any(strcmp(glutemax, templabel_exo))
                tempglute = welkexostruct.(genvarname(templabel_exo));
                glutemax_data_exo = [glutemax_data_exo, tempglute];
            end
            if any(strcmp(glutemed, templabel_exo))
                tempglute = welkexostruct.(genvarname(templabel_exo));
                glutemed_data_exo = [glutemed_data_exo, tempglute];
            end
            if any(strcmp(glutemin, templabel_exo))
                tempglute = welkexostruct.(genvarname(templabel_exo));
                glutemin_data_exo = [glutemin_data_exo, tempglute];
            end
        end
        glutemax_data_exo = mean(glutemax_data_exo, 2);
        glutemed_data_exo = mean(glutemed_data_exo, 2);
        glutemin_data_exo = mean(glutemin_data_exo, 2);

        
        % make sure the new averaged will get into figure
        welknaturalstruct.glmax_avg_r = glutemax_data_nat;
        welknaturalstruct.glmed_avg_r = glutemed_data_nat;
        welknaturalstruct.glmin_avg_r = glutemin_data_nat;
        
        welkexostruct.glmax_avg_r = glutemax_data_exo;
        welkexostruct.glmed_avg_r = glutemed_data_exo;
        welkexostruct.glmin_avg_r = glutemin_data_exo;
        
        % need to get new total labels
        testlabels_nat = fields(welknaturalstruct);
        testlabels_exo = fields(welkexostruct);        
        
        
        % now create a figure 
        % tempfig = figure('Position',[1,1,1920,1080]);
        % for i=2:length(testlabels_nat)
        %     subplot(7,9,i-1);
        %     templabel = char(testlabels_nat(i));
        %     muscleplot1 = welknaturalstruct.(genvarname(templabel));
        %     plot(welknaturalstruct.time, muscleplot1, ':')
        %     hold on;
        %     plot(welknaturalstruct.time, mean(muscleplot1,2), 'k-', 'LineWidth', 2)
        %     title(templabel)
        %     xlabel('% gait cycle')
        %     ylabel(tempthing)
        %     grid on;
        % end
        % print(tempfig, ...
        %     strcat(repodir,'\..\analysis\',subject,'\', tempthing, '_natural', '.png'),...
        %     '-dpng', '-r500')
        % disp('print 1')
        
        
        % tempfig2 = figure('Position',[1,1,1920,1080]);
        % for i=2:length(testlabels_exo)
        %     subplot(7,9,i-1);
        %     templabel = char(testlabels_exo(i));
        %     muscleplot2 = welkexostruct.(genvarname(templabel));
        %     plot(welkexostruct.time, muscleplot2, ':')
        %     hold on;
        %     plot(welkexostruct.time, mean(muscleplot2,2), 'k-', 'LineWidth', 2)
        %     title(templabel)
        %     xlabel('% gait cycle')
        %     ylabel(tempthing)
        %     grid on;
        % end
        % print(tempfig2, ...
        %     strcat(repodir,'\..\analysis\',subject,'\', tempthing, '_exo', '.png'),...
        %     '-dpng', '-r500')
        % disp('print 2')
        
        
        % combined exo and natural fig
        % combineexonaturalfig = figure('Position',[1,1,1920,1080]);
        % title('red=natural, blue=exo');
        labels = fields(welkexostruct);
        for i=2:length(labels)
            % subplot(7,9,i-1);
            templabel = char(labels(i));
            muscleplot2 = welkexostruct.(genvarname(templabel));
            muscleplot1 = welknaturalstruct.(genvarname(templabel));
            % % plot(welkexostruct.time, muscleplot1, ':')
            % hold on;
            % plot(welkexostruct.time, mean(muscleplot2,2), 'b-', 'LineWidth', 2)
            % plot(welknaturalstruct.time, mean(muscleplot1,2), 'r-', 'LineWidth', 2)
            % title(templabel)
            % xlabel('% gait cycle')
            % ylabel(tempthing)
            % grid on;
            % % legend('exo','natural')
            if tempthing == 'excitation'
                exomeans_excitation.(genvarname(subject)) = [exomeans_excitation.(genvarname(subject)), mean(muscleplot2, 2)];
                naturalmeans_excitation.(genvarname(subject)) = [naturalmeans_excitation.(genvarname(subject)), mean(muscleplot1, 2)];
                excitelabels = fields(welkexostruct);
                % exomeans_excitation. = [exomeans_excitation, mean(muscleplot2,2)];
                % naturalmeans_excitation = [naturalmeans_excitation, mean(muscleplot1, 2)];

                exopeaks_excitation.(genvarname(subject)) = [exopeaks_excitation.(genvarname(subject)), max(mean(muscleplot2, 2))];
                naturalpeaks_excitation.(genvarname(subject)) = [naturalpeaks_excitation.(genvarname(subject)), max(mean(muscleplot1, 2))];

            end
            if tempthing == 'activation'
                exomeans_activation.(genvarname(subject)) = [exomeans_activation.(genvarname(subject)), mean(muscleplot2, 2)];
                naturalmeans_activation.(genvarname(subject)) = [naturalmeans_activation.(genvarname(subject)), mean(muscleplot1, 2)];
                activelabels = fields(welkexostruct);
                % exomeans_activation = [exomeans_activation, mean(muscleplot2, 2)];
                % naturalmeans_activation = [naturalmeans_activation, mean(muscleplot1, 2)];

                exopeaks_activation.(genvarname(subject)) = [exopeaks_activation.(genvarname(subject)), max(mean(muscleplot2, 2))];
                naturalpeaks_activation.(genvarname(subject)) = [naturalpeaks_activation.(genvarname(subject)), max(mean(muscleplot1, 2))];

            end
        end
        % print(combineexonaturalfig, ...
        %     strcat(repodir,'\..\analysis\',subject,'\', tempthing, '_combined', '.png'),...
        %     '-dpng', '-r500')
        disp('print combined')
    end
    nat_maxiso.(genvarname(subject)) = [nat_maxiso.(genvarname(subject)), welknat_isostruct];
    exo_maxiso.(genvarname(subject)) = [exo_maxiso.(genvarname(subject)), welkexo_isostruct];


end



% combined exo and natural fig - averages for subjects and mean
subjectcombineexonaturalfig1 = figure('Position',[1,1,1200,1920]);
title('red=natural, blue=exo');
% do more stuff
% averaging and whatnot
    

for i=1:length(excitelabels)-1    
    % make a subplot
    subplot(16,4,i);
    templabel = char(excitelabels(i+1));
    hold on;
    tempsubjavgs1 = [];
    tempsubjavgs2 = [];
    
    % for each muscle go through each subject
    for subj=1:length(welksubjects)
        subject = char(welksubjects(subj));
        plot(welkexostruct.time, exomeans_excitation.(genvarname(subject))(:,i), 'b:', 'LineWidth', 0.4)
        plot(welknaturalstruct.time, naturalmeans_excitation.(genvarname(subject))(:,i), 'r:', 'LineWidth', 0.4)
        
        % add each signal to the temp vector to get means before moving to
        % the next muscle
        tempsubjavgs1 = [tempsubjavgs1, naturalmeans_excitation.(genvarname(subject))(:,i)];
        tempsubjavgs2 = [tempsubjavgs2, exomeans_excitation.(genvarname(subject))(:,i)];
    end
    
    % now plot the means from the tempsubjavgs1/2
    plot(welkexostruct.time, mean(tempsubjavgs2, 2), 'b-', 'LineWidth', 2)
    plot(welknaturalstruct.time, mean(tempsubjavgs1, 2), 'r-', 'LineWidth', 2)
    
    
    title(templabel)
    xlabel('% gait cycle')
    ylabel(tempthing)
%     grid on;
    % legend('exo','natural')
%     legend(strcat('exo peak: ', num2str(max(mean(tempsubjavgs2, 2)))),strcat('nat peak: ',num2str(max(mean(tempsubjavgs1, 2)))))

end

% print(subjectcombineexonaturalfig1, ...
%     strcat(repodir, '\..\analysis\', 'excitation_all_subjects_combined_nolegend.png'),...
%     '-dpng', '-r500')
disp('print combined')


% now for the activations
% combined exo and natural fig - averages for subjects and mean
subjectcombineexonaturalfig2 = figure(100); %'Number',100,'Position',[1,1,1920,1080]);
set(gcf,'WindowStyle','Docked','Position',[1,1,3508,1280,])
% set(gcf,'Position',[1,1,700,4508])

title('red=natural, blue=exo');
% do more stuff
% averaging and whatnot


for i=1:length(activelabels)-1   
    % make a subplot
    subplot(5,11,i);
%     axis('square')
%     axis('tight')
    templabel = char(activelabels(i+1));
    hold on;
    tempsubjavgs1 = [];
    tempsubjavgs2 = [];
    
    % for each muscle go through each subject
    for subj=1:length(welksubjects)
        subject = char(welksubjects(subj));
%         plot(welkexostruct.time, exomeans_activation.(genvarname(subject))(:,i), 'b:', 'LineWidth', 0.4)
%         plot(welknaturalstruct.time, naturalmeans_activation.(genvarname(subject))(:,i), 'r:', 'LineWidth', 0.4)
        
        % add each signal to the temp vector to get means before moving to
        % the next muscle
        tempsubjavgs1 = [tempsubjavgs1, naturalmeans_activation.(genvarname(subject))(:,i)];
        tempsubjavgs2 = [tempsubjavgs2, exomeans_activation.(genvarname(subject))(:,i)];
    end
    
    disp(templabel)
    % spm tests for simulated activities
    addpath(genpath('G:\Shared drives\Exotendon\Matlab\common\')) ;
    spm = spm1d.stats.nonparam.ttest_paired(squeeze(tempsubjavgs1)', squeeze(tempsubjavgs2)');
    spmi   = spm.inference(0.05, 'two_tailed',true, 'interp',true);
    disp(spmi)
%     tempfig = figure;
%     spmi.plot();
%     spmi.plot_threshold_label();
%     spmi.plot_p_values();
%     title(templabel)
%     print(tempfig, strcat(repodir, '\..\analysis\activitySPM\',templabel,'.png'), ...
%         '-dpng', '-r500')
    % now plot the means from the tempsubjavgs1/2
    figure(100);
    plot(welkexostruct.time, mean(tempsubjavgs2, 2), 'b-', 'LineWidth', 2)
    plot(welknaturalstruct.time, mean(tempsubjavgs1, 2), 'r-', 'LineWidth', 2)
    
    
    title(templabel)
    xlabel('% gait cycle')
    ylabel(tempthing)
%     grid on;
    % legend('exo','natural')
%     legend(strcat('exo peak: ', num2str(max(mean(tempsubjavgs2, 2)))),strcat('nat peak: ',num2str(max(mean(tempsubjavgs1, 2)))))
end

% print(subjectcombineexonaturalfig2, ...
%     strcat(repodir, '\..\analysis\', 'testgroups_activation_all_subjects_combined_nolegend.png'),...
%     '-dpng', '-r500')
disp('print combined')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% now need to get the EMG stuff in here and grouped as well to add to the same figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
% EMG_analysis.m
% This program loads processed EMG signals and determines portions of the
% gait cycle when they are significantly different 
% Written by: Cara Welker, updated Jon Stingel
% 3/9/2023

% close all; clc; clear
% set(0,'DefaultFigureWindowStyle','docked')

if ismac
    basedir = '/Volumes/GoogleDrive/Shared drives/Exotendon/DATA';
    addpath(genpath('/Volumes/GoogleDrive/Shared drives/Exotendon/Matlab/common'));
    addpath('G:\Shared drives\Exotendon\Matlab\plotting');
else
    basedir = 'G:\Shared drives\Exotendon\DATA\' ;
    addpath(genpath('G:\Shared drives\Exotendon\Matlab\common\')) ;
    addpath('G:\Shared drives\Exotendon\Matlab\plotting');
end

% load('EMG_muscles_removed_withS4.mat');
% load('EMG_muscles_removed.mat');
load('EMG_all_updated.mat')
[nSubs, nTrials] = size(EMG.data);
nMuscl = length(EMG.muscleNames);
timesteps = 101;

%% idk what this is
EMG.avg_avg = zeros(timesteps, nMuscl);
EMG.avg = EMG.data;
% Loop through the subjects
for sub = 1:nSubs

        %for statistic testing, make the number of gait cycles tested
        %the same for exo and non-exo cases
        
        %ToDO; rewrite so not hard-coded in!
        
        EMG_exo = EMG.data{sub}; %exo
        EMG_nat = EMG.data{sub+nSubs}; %natural
        nTrials = min(length(squeeze(EMG_exo(1,1,:))),...
            length(squeeze(EMG_nat(1,1,:))));
        EMG_exo = EMG_exo(:,:,length(squeeze(EMG_exo(1,1,:)))-nTrials+1:end);
        EMG_nat = EMG_nat(:,:,length(squeeze(EMG_nat(1,1,:)))-nTrials+1:end);
        EMG.data{sub} = EMG_exo;
        EMG.data{sub+nSubs} = EMG_nat;

        % take averages and add to struct
        EMG.avg{sub} = mean(EMG_exo,3);
        EMG.avg{sub+nSubs} = mean(EMG_nat,3);
    
end % sub

%% statistical analysis and plotting for all subjects
%TODO: change so generalizable for different trials!
EMG_plotting = zeros(timesteps, nMuscl, nSubs, 2);
%EMG_plotting_avg = zeros(timesteps, nMuscl, 2);

Trials = [1,2]; % 1,2 bc natural used to be the 5th trial
for trial = 1:length(Trials)
    nTrial = Trials(trial);
    temp = zeros(timesteps,nSubs,nMuscl);
    for sub = 1:nSubs
        for muscl = 1:nMuscl
            EMG_plotting(:,muscl,sub,trial) = EMG.avg{sub,nTrial}(:,muscl);
            %temp(:,sub,muscl) = EMG.avg{sub,nTrial}(:,muscl);
        end
    end
  %  EMG_plotting_avg(:,:,trial) = squeeze(nanmean(temp,2));
end

params = [];
params.fh = 1;
data = [];
data.time = 0:100;
data.colheaders = EMG.muscleNames;

% EMG_plotting is 101, 14, 7, 2
quadEMG = {'VL','VM','RF'};
hipflexEMG = {'PS'};
hipabdEMG = {'GMED'};
hamEMG = {'BF','ST'};
hipaddEMG = {'ADD'};
hipextEMG = {'GMAX'};
pfEMG = {'SOL','LG','MG','PER'};
dfEMG = {'TA'};


secondtime = false;
for i = 1:2
    if i==2
        secondtime = true;
    end
    data.data = squeeze(EMG_plotting(:,:,:,i));
    if i==2
        %comparison data for statistical test
        comp_data = squeeze(EMG_plotting(:,:,:,1)); 
        params.color = [0.5, 0, 0];
    else
        comp_data = squeeze(EMG_plotting(:,:,:,2));
        params.color = [0, 0, 0.5];
    end
    
    % need a way to get the iso forces to multiply
    for subj=1:length(welksubjects)
        subject = char(welksubjects(subj));
        % get data
        tempPS = data.data(:,1,subj);
        tempADD = data.data(:,2,subj);
        tempRF = data.data(:,3,subj);
        tempVL = data.data(:,4,subj);
        tempVM = data.data(:,5,subj);
        tempGMAX = data.data(:,6,subj);
        tempGMED = data.data(:,7,subj);
        tempTA = data.data(:,8,subj);
        tempBF = data.data(:,9,subj);
        tempST = data.data(:,10,subj);
        tempLG = data.data(:,11,subj);
        tempMG = data.data(:,12,subj);
        tempSOL = data.data(:,13,subj);
        tempPER = data.data(:,14,subj);
        
        % get the max iso for each muscle for each subject
        subjisos = nat_maxiso.(genvarname(subject));

        PSiso = subjisos.psoas_r(1);
        ADDiso = subjisos.addlong_r(1);
        RFiso = subjisos.recfem_r(1);
        VLiso = subjisos.vaslat_r(1);
        VMiso = subjisos.vasmed_r(1);
        GMAXiso = mean([subjisos.glmax1_r(1), subjisos.glmax2_r(1), subjisos.glmax3_r(1)]);
        GMEDiso = mean([subjisos.glmed1_r(1), subjisos.glmed2_r(1), subjisos.glmed3_r(1)]);
        TAiso = subjisos.tibant_r(1);
        BFiso = subjisos.bflh_r(1);
        STiso = subjisos.semiten_r(1);
        LGiso = subjisos.gaslat_r(1);
        MGiso = subjisos.gasmed_r(1);
        SOLiso = subjisos.soleus_r(1);
        PERiso = subjisos.perlong_r(1);

        % get averages by hand easily
        quadsEMGavg = ((tempVL.*VLiso) + (tempVM.*VMiso) + (tempRF.*RFiso))/(VLiso + VMiso + RFiso);
        hipflexEMGavg = tempPS;
        hipabdEMGavg = tempGMED;
        hamEMGavg = ((tempBF.*BFiso) + (tempST.*STiso))/(BFiso + STiso);
        hipaddEMGavg = tempADD;
        hipextEMGavg = tempGMAX;
        pfEMGavg = ((tempSOL.*SOLiso) + (tempPER.*PERiso) + (tempLG.*LGiso) + (tempMG.*MGiso))/(SOLiso + PERiso + LGiso + MGiso);
        dfEMGavg = tempTA;

        % add to struct
        if i==1
            % we have the exo stuff so add to the exo struct
            exomeans_activation.(genvarname(subject)) = [exomeans_activation.(genvarname(subject)), ...
                                quadsEMGavg, hipflexEMGavg, hipabdEMGavg, hamEMGavg, hipaddEMGavg, hipextEMGavg, pfEMGavg, dfEMGavg];
        else
            naturalmeans_activation.(genvarname(subject)) = [naturalmeans_activation.(genvarname(subject)), ...
                                quadsEMGavg, hipflexEMGavg, hipabdEMGavg, hamEMGavg, hipaddEMGavg, hipextEMGavg, pfEMGavg, dfEMGavg];
        end
    end


    % plot_results(data, comp_data, params, secondtime); 
end
activelabels = [activelabels; {'quadEMG'};{'hipflexEMG'};{'hipabdEMG'};{'hamEMG'};{'hipaddEMG'};{'hipextEMG'};{'pfEMG'};{'dfEMG'}];


% print(strcat('G:\Shared drives\Exotendon\muscleModel\analysis\', 'emg_updated_all_w_spm.png'),...
    % '-dpng', '-r500')


keyboard
%% finally plot them both together. 

% now for the activations
% combined exo and natural fig - averages for subjects and mean
subjectcombineexonaturalfig3 = figure(200); %'Number',100,'Position',[1,1,1920,1080]);
set(gcf,'WindowStyle','Docked','Position',[1,1,1500,1280,])
% set(gcf,'Position',[1,1,700,4508])

title('red=natural, blue=exo');
% do more stuff
% averaging and whatnot
% 'hipflexors'
% 'hipflexEMG'
stuffwewant = {'quadriceps','hipabductors','hamstrings','hipadductors','hipextensors','plantarflexors','dorsiflexors'};
emgwewant = {'quadEMG','hipabdEMG','hamEMG','hipaddEMG','hipextEMG','pfEMG','dfEMG'};

plotcount = 0;
emgplotcount = 0
for i=1:length(activelabels)-1   
    templabel = char(activelabels(i+1));
    if any(strcmp(templabel, stuffwewant))
        
        plotcount = plotcount+1;
        % make a subplot
        subplot(3,4,plotcount);
    %     axis('square')
    %     axis('tight')
        hold on;
        tempsubjavgs1 = [];
        tempsubjavgs2 = [];
        
        % for each muscle go through each subject
        for subj=1:length(welksubjects)
            subject = char(welksubjects(subj));
    %         plot(welkexostruct.time, exomeans_activation.(genvarname(subject))(:,i), 'b:', 'LineWidth', 0.4)
    %         plot(welknaturalstruct.time, naturalmeans_activation.(genvarname(subject))(:,i), 'r:', 'LineWidth', 0.4)
            
            % add each signal to the temp vector to get means before moving to
            % the next muscle
            tempsubjavgs1 = [tempsubjavgs1, naturalmeans_activation.(genvarname(subject))(:,i)];
            tempsubjavgs2 = [tempsubjavgs2, exomeans_activation.(genvarname(subject))(:,i)];
        end
        
        % disp(templabel)
        % spm tests for simulated activities
        % addpath(genpath('G:\Shared drives\Exotendon\Matlab\common\')) ;
        % spm = spm1d.stats.nonparam.ttest_paired(squeeze(tempsubjavgs1)', squeeze(tempsubjavgs2)');
        % spmi   = spm.inference(0.05, 'two_tailed',true, 'interp',true);
        % disp(spmi)
    %     tempfig = figure;
    %     spmi.plot();
    %     spmi.plot_threshold_label();
    %     spmi.plot_p_values();
    %     title(templabel)
    %     print(tempfig, strcat(repodir, '\..\analysis\activitySPM\',templabel,'.png'), ...
    %         '-dpng', '-r500')
        % now plot the means from the tempsubjavgs1/2
        figure(200);
        % if any(strcmp(stuffwewant, templabel))
            disp(templabel)
            hold on;
            plot(welkexostruct.time, mean(tempsubjavgs2, 2), 'Color',exocolor,'LineStyle',':', 'LineWidth', 2)
            plot(welknaturalstruct.time, mean(tempsubjavgs1, 2), 'Color',natcolor,'LineStyle',':', 'LineWidth', 2)
        
        
            title(templabel)
            xlabel('% gait cycle')
            ylabel(tempthing)
        %     grid on;
            % legend('exo','natural')
        %     legend(strcat('exo peak: ', num2str(max(mean(tempsubjavgs2, 2)))),strcat('nat peak: ',num2str(max(mean(tempsubjavgs1, 2)))))
        % end
    end
    if any(strcmp(templabel, emgwewant))
        
        emgplotcount = emgplotcount+1;
        % make a subplot
        subplot(3,4,emgplotcount);
    %     axis('square')
    %     axis('tight')
        hold on;
        tempsubjavgs1 = [];
        tempsubjavgs2 = [];
        
        % for each muscle go through each subject
        for subj=1:length(welksubjects)
            subject = char(welksubjects(subj));
    %         plot(welkexostruct.time, exomeans_activation.(genvarname(subject))(:,i), 'b:', 'LineWidth', 0.4)
    %         plot(welknaturalstruct.time, naturalmeans_activation.(genvarname(subject))(:,i), 'r:', 'LineWidth', 0.4)
            
            % add each signal to the temp vector to get means before moving to
            % the next muscle
            tempsubjavgs1 = [tempsubjavgs1, naturalmeans_activation.(genvarname(subject))(:,i)];
            tempsubjavgs2 = [tempsubjavgs2, exomeans_activation.(genvarname(subject))(:,i)];
        end
        
        % disp(templabel)
        % spm tests for simulated activities
        % addpath(genpath('G:\Shared drives\Exotendon\Matlab\common\')) ;
        % spm = spm1d.stats.nonparam.ttest_paired(squeeze(tempsubjavgs1)', squeeze(tempsubjavgs2)');
        % spmi   = spm.inference(0.05, 'two_tailed',true, 'interp',true);
        % disp(spmi)
    %     tempfig = figure;
    %     spmi.plot();
    %     spmi.plot_threshold_label();
    %     spmi.plot_p_values();
    %     title(templabel)
    %     print(tempfig, strcat(repodir, '\..\analysis\activitySPM\',templabel,'.png'), ...
    %         '-dpng', '-r500')
        % now plot the means from the tempsubjavgs1/2
        figure(200);
        % if any(strcmp(emgwewant, templabel))
            disp(templabel)
            hold on
            plot(welkexostruct.time, mean(tempsubjavgs2, 2), 'LineStyle','-','Color',exocolor,'LineWidth', 2)
            plot(welknaturalstruct.time, mean(tempsubjavgs1, 2), 'LineStyle','-','Color', natcolor, 'LineWidth', 2)
        
        
%             title(templabel)
            xlabel('% gait cycle')
            ylabel(tempthing)
        %     grid on;
            % legend('exo','natural')
        %     legend(strcat('exo peak: ', num2str(max(mean(tempsubjavgs2, 2)))),strcat('nat peak: ',num2str(max(mean(tempsubjavgs1, 2)))))
        % end
    end
end

wantlegend = true;
if wantlegend
    templabel = 'quadriceps';
    if any(strcmp(templabel, stuffwewant))

    plotcount = plotcount+1;
    % make a subplot
    subplot(3,4,plotcount);
    %     axis('square')
    %     axis('tight')
    hold on;
    tempsubjavgs1 = [];
    tempsubjavgs2 = [];

    % for each muscle go through each subject
    for subj=1:length(welksubjects)
        subject = char(welksubjects(subj));
    %         plot(welkexostruct.time, exomeans_activation.(genvarname(subject))(:,i), 'b:', 'LineWidth', 0.4)
    %         plot(welknaturalstruct.time, naturalmeans_activation.(genvarname(subject))(:,i), 'r:', 'LineWidth', 0.4)
        
        % add each signal to the temp vector to get means before moving to
        % the next muscle
        tempsubjavgs1 = [tempsubjavgs1, naturalmeans_activation.(genvarname(subject))(:,i)];
        tempsubjavgs2 = [tempsubjavgs2, exomeans_activation.(genvarname(subject))(:,i)];
    end

    % disp(templabel)
    % spm tests for simulated activities
    % addpath(genpath('G:\Shared drives\Exotendon\Matlab\common\')) ;
    % spm = spm1d.stats.nonparam.ttest_paired(squeeze(tempsubjavgs1)', squeeze(tempsubjavgs2)');
    % spmi   = spm.inference(0.05, 'two_tailed',true, 'interp',true);
    % disp(spmi)
    %     tempfig = figure;
    %     spmi.plot();
    %     spmi.plot_threshold_label();
    %     spmi.plot_p_values();
    %     title(templabel)
    %     print(tempfig, strcat(repodir, '\..\analysis\activitySPM\',templabel,'.png'), ...
    %         '-dpng', '-r500')
    % now plot the means from the tempsubjavgs1/2
    figure(200);
    % if any(strcmp(stuffwewant, templabel))
        disp(templabel)

        plot(welkexostruct.time, mean(tempsubjavgs2, 2), 'Color',exocolor,'LineStyle',':', 'LineWidth', 2)
        plot(welknaturalstruct.time, mean(tempsubjavgs1, 2), 'Color',natcolor,'LineStyle',':', 'LineWidth', 2)


        title(templabel)
        xlabel('% gait cycle')
        ylabel(tempthing)
    %     grid on;
        % legend('exo','natural')
    %     legend(strcat('exo peak: ', num2str(max(mean(tempsubjavgs2, 2)))),strcat('nat peak: ',num2str(max(mean(tempsubjavgs1, 2)))))
    % end
    end
    templabel = 'quadEMG';
    if any(strcmp(templabel, emgwewant))

    emgplotcount = emgplotcount+1;
    % make a subplot
    subplot(3,4,emgplotcount);
    %     axis('square')
    %     axis('tight')
    hold on;
    tempsubjavgs1 = [];
    tempsubjavgs2 = [];

    % for each muscle go through each subject
    for subj=1:length(welksubjects)
        subject = char(welksubjects(subj));
    %         plot(welkexostruct.time, exomeans_activation.(genvarname(subject))(:,i), 'b:', 'LineWidth', 0.4)
    %         plot(welknaturalstruct.time, naturalmeans_activation.(genvarname(subject))(:,i), 'r:', 'LineWidth', 0.4)
        
        % add each signal to the temp vector to get means before moving to
        % the next muscle
        tempsubjavgs1 = [tempsubjavgs1, naturalmeans_activation.(genvarname(subject))(:,i)];
        tempsubjavgs2 = [tempsubjavgs2, exomeans_activation.(genvarname(subject))(:,i)];
    end

    % disp(templabel)
    % spm tests for simulated activities
    % addpath(genpath('G:\Shared drives\Exotendon\Matlab\common\')) ;
    % spm = spm1d.stats.nonparam.ttest_paired(squeeze(tempsubjavgs1)', squeeze(tempsubjavgs2)');
    % spmi   = spm.inference(0.05, 'two_tailed',true, 'interp',true);
    % disp(spmi)
    %     tempfig = figure;
    %     spmi.plot();
    %     spmi.plot_threshold_label();
    %     spmi.plot_p_values();
    %     title(templabel)
    %     print(tempfig, strcat(repodir, '\..\analysis\activitySPM\',templabel,'.png'), ...
    %         '-dpng', '-r500')
    % now plot the means from the tempsubjavgs1/2
    figure(200);
    % if any(strcmp(emgwewant, templabel))
        disp(templabel)

        plot(welkexostruct.time, mean(tempsubjavgs2, 2), 'Color',exocolor,'LineStyle','-', 'LineWidth', 2)
        plot(welknaturalstruct.time, mean(tempsubjavgs1, 2), 'Color',natcolor,'LineStyle','-', 'LineWidth', 2)


    %             title(templabel)
        xlabel('% gait cycle')
        ylabel(tempthing)
    %     grid on;
        % legend('exo','natural')
    %     legend(strcat('exo peak: ', num2str(max(mean(tempsubjavgs2, 2)))),strcat('nat peak: ',num2str(max(mean(tempsubjavgs1, 2)))))
    % end
    end
    
end
legend('Simulation exotendon','Simulation natural','EMG exotendon','EMG natural','Location','eastoutside')

print(subjectcombineexonaturalfig3, ...
    strcat(repodir, '\..\analysis\', 'MUSCLEGROUPS_activationstuff.png'),...
    '-dpng', '-r500')
disp('print combined')



