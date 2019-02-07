clc
clear
close all

% Load nparam and priors and Geweke test results
load('options.mat','F')
nparam = F.nparam;
priors = F.priors;
ndraws = 50000;
ntasks = 50;
load('geweke_results')

%% Results
format bank;

% Column headers
temp = [1,taper_pers];
for itemp = 1:length(temp)
    VariableNames{itemp} = ['taper',num2str(temp(itemp))];
end

% percent of samples where test fails to reject equal means at 5% level
conflevs = [0.01,0.05];
for ilev = 1:numel(conflevs)
    conflev = conflevs(ilev);
    tempindresults = all(results > conflev,1);
    indresults = reshape(tempindresults,[5,300,3]);
    indresults = squeeze(100*mean(indresults,2,'omitnan'))';
    tab2 = array2table(indresults,'VariableNames',VariableNames);
    tab2.Properties.RowNames = methods;
    disp(['Percent of samples where Chi-squared test fails to reject null at ' num2str(100*conflev) '% level'])
    disp('  Column header indicates tapering % of NSE (smooths draws to remove noise)')
    disp(tab2)
end