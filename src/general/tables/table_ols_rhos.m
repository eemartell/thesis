%  Creates a figure of the percent of actuals outside of forecast credible sets
% from the different model %  specifications from the periods just prior to the ZLB binding.

clear
clc
close all

npers = 120;
nsimsall = 100;
nsims = 50;
nZLBpers = 6;
nqs = [0,6,12,18,24,30];
ndraws = 1000;

% Load parameter names
load('../options.mat','V','P','F')

% Load actual data
simscap = textread('../../data-artificial/sims.txt','%f');
simscap = reshape(simscap,[npers,V.nvar-V.nfore,nsimsall,nZLBpers]);

Vcap= V;

%Load in V from misspecified model
load('../../Results-nocap/options.mat','V','P','F')
Vnstate = V.nvar-V.nfore;

%Create artificial data matrix (true model) that aligns with variable order
%from the true model. Note that yg is named cg in misspec model
sims = nan(npers,Vnstate,nsims,nZLBpers);
for ivar = 1:(Vnstate)
    if strcmp(V.names{ivar},'cg')
        sims(:,ivar,:,:) = simscap(:,Vcap.yg,1:nsims,:);
    else
        eval(['sims(:,ivar,:,:) = simscap(:,Vcap.',V.names{ivar},',1:nsims,:);']);
    end
end

s = squeeze(sims(2:120,V.s,:,:));
lags = squeeze(sims(1:119,V.s,:,:));

for izlb = 1:6
    for itask = 1:50
        b = regress(s(:,itask,izlb),[ones(119,1),lags(:,itask,izlb)]);
        rhos(itask,izlb) = b(2);
    end
end
median(rhos,1)