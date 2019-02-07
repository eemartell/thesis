clc
clear
close all

% Load nparam and priors
load('options.mat','F')
nparam = F.nparam;
priors = F.priors;
ndraws = 50000;
ntasks = 50;

% Set Geweke diagnostic options and subsamples
%   % of sample to burn off from beginning
perburn = 0;
%   Tapering percentages to remove noise from moments 
taper_pers = [2,4,8,15];
geweke_interval = [.2,.5];
nburn = round(perburn(1)*ndraws);
nhalf = round(geweke_interval(1)*(ndraws-nburn));
n1 = nburn + 1;
n2 = nburn + round(geweke_interval(1)*(ndraws-nburn));
n3 = nburn + round((1-geweke_interval(2))*(ndraws-nburn))+1;
n4 = ndraws;
disp('Sample 1 indices')
disp([n1,n2])
disp('Sample 2 indices')
disp([n3,n4])

for taskid = 1:ntasks
    % Load draws and posterior log-likelihoods
    filename = ['Estimation-global-me5-zlb1-mh1/' num2str(taskid) '-MHdraws1.txt'];
    MHdraws = textread(filename,'%f');
    MHdraws = reshape(MHdraws,[numel(MHdraws)/nparam,nparam]);

    % Dynare defaults for Geweke test
    for iparam = 1:nparam
        [results_vec1, results_struct1(iparam)] = geweke_mom(MHdraws(n1:n2,iparam),taper_pers);
        [results_vec2, results_struct2(iparam)] = geweke_mom(MHdraws(n3:n4,iparam),taper_pers);
        results_struct3(iparam) = geweke_test(results_vec1,results_vec2,[],taper_pers);
        results(iparam,:,taskid) = results_struct3(iparam).prob_chi2_test;
    end
end

%% Results
format bank;

% Column headers
temp = [1,taper_pers];
for itemp = 1:length(temp)
    VariableNames{itemp} = ['taper',num2str(temp(itemp))];
end

% p-values for particular task
taskid = 1;

tab = array2table(results(:,:,taskid),'VariableNames',VariableNames);
tab.Properties.RowNames = F.priors(:,1);

disp('Geweke Chi-squared test (p values)')
disp('  Column header indicates tapering % of NSE (smooths draws to remove noise)')
disp('  If p > confidence level, e.g., 0.05, then test fails to reject equal means.')
disp(tab)

% percent of samples where test fails to reject equal means at 5% level
conflevs = [0.01,0.05];
for ilev = 1:numel(conflevs)
    conflev = conflevs(ilev);
    indresults = 100*sum(all(results > conflev,1),3)/ntasks;
    tab2 = array2table(indresults,'VariableNames',VariableNames);
    disp(['Percent of samples where Chi-squared test fails to reject null at ' num2str(100*conflev) '% level'])
    disp('  Column header indicates tapering % of NSE (smooths draws to remove noise)')
    disp(tab2)
end