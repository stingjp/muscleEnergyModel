%% Repeated Measures Analysis of Variance  

%% 
% Load the sample data. 
% load fisheriris 

% for i=1:length(subjects)
%     fprintf(subjects(i))
%     fprintf('\n')
% end




%%
% The column vector |species| consists of iris flowers of three different
% species: setosa, versicolor, virginica. The double matrix |meas| consists
% of four types of measurements on the flowers: the length and width of
% sepals and petals in centimeters, respectively.

%% 
% Store the data in a table array. 
% t = table(species,meas(:,1),meas(:,2),meas(:,3),meas(:,4),...
% 'VariableNames',{'species','meas1','meas2','meas3','meas4'});
t = table(subjects,conditions(:,1),conditions(:,2),conditions(:,3),...
    conditions(:,4),conditions(:,5),'VariableNames',...
    {'subjects','slack','low','med','high','max'});

% Meas = table([1 2 3 4]','VariableNames',{'Measurements'});  
Meas = table([1 2 3 4 5]','VariableNames',{'Measurements'});
 
% %% 
% % Fit a repeated measures model, where the measurements are the responses
% % and the species is the predictor variable. 
% rm = fitrm(t,'meas1-meas4~species','WithinDesign',Meas);  
rm = fitrm(t,'slack-max~subjects','WithinDesign',Meas);

% 
% %% 
% % Perform repeated measures analysis of variance. 
% ranovatbl = ranova(rm) 
ranovatbl = ranova(rm)
oneanova = anova1(conditions)

%%
% There are four measurements, three types of species, and 150
% observations. So, degrees of freedom for measurements is (4&#8211;1) = 3,
% for species-measurements interaction it is (4&#8211;1)*(3&#8211;1) = 6,
% and for error it is (150&#8211;4)*(3&#8211;1) = 441. |ranova| computes
% the last three $p$-values using Greenhouse-Geisser, Huynh-Feldt, and
% Lower bound corrections, respectively. You can check the compound
% symmetry (sphericity) assumption using the |mauchly| method, and display
% the epsilon corrections using the |epsilon| method.
%%
% fprintf('\ntangent modulus for 20 percent stretch (stretch = 1.2):\n')
% [pval_tanmod12,~,stats_tanmod12] = anova1([tanmod12_01, tanmod12_4, tanmod12_160], category_peakstress_01);
% fprintf(['we see a significant difference in the tangent modulus when stretch = 1.2\n'...
%     'between groups, as p=0.1881. Additionally, we can see that there is a larger\n'...
%     'variation in the 4 percent strain rate group once again.'])
% multcompare_tanmod12 = multcompare(stats_tanmod12);
% fprintf(['\nfrom the plot we can see that both the 0.1 percent rate and the 160\n'...
%     'percent rate have different means than the 4 percent rate.\n'])
[idk,~,statsidk] = anova1(conditions(:,1:5))
multcompare_idk = multcompare(statsidk)