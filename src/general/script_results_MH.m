% RW-MH results
clear all
close all
clc
plotops

% Method and Task ID
% method = 'Estimation-global';
% method = 'Estimation-pwlinear';
method = 'Estimation-global-ukf';
iZLBper = 1;
taskid = 5;

% Results options
nmodesearchdraws = 5000;

% Load parameters
load('options','P','F');

% Load priors, data, draws, and posterior log-likelihoods
sfthin = ['../' method '/' num2str(iZLBper) '/thin/'];
sfms = ['../' method '/' num2str(iZLBper) '/ms/'];
filename = [sfthin num2str(taskid) '-MHdraws.txt'];
draws = textread(filename,'%f');
draws = reshape(draws,[numel(draws)/F.nparam,F.nparam]);
filename = [sfthin num2str(taskid) '-MHpostlogliks.txt'];
postlogliks = textread(filename,'%f');
ndiscarddraws = textread([sfms num2str(taskid) '-ndiscarddraws.txt'],'%f');
idraw = find(sum(abs(draws),2) == 0,1,'first') - 1;
if isempty(idraw); idraw = size(draws,1); end

% Adjust posterior density for truncated prior density
ntotalmodearchdraws = ndiscarddraws + nmodesearchdraws;
postlogliks = postlogliks + log(1 - ndiscarddraws/ntotalmodearchdraws);

% Date modified
file = dir(filename);
disp(['Loaded: ' filename])
disp(['Last Modified: ' file.date ])

% John Geweke's Modified Harmonic Mean
mhm = geweke_mhm(draws,postlogliks);
disp(['Data density: ' num2str(mean(mhm))]);

% Posterior Distribution
probs = [.05,.95];
drawpostmean = mean(draws,1);
quantMHdraws = [drawpostmean' quantile(draws,probs)'];
tab = table(quantMHdraws(:,1),quantMHdraws(:,2),quantMHdraws(:,3));
tab.Properties.VariableNames = {'Mean','p5','p95'};
tab.Properties.RowNames = F.priors(:,1);
format shortg; disp(tab)


% Figure padding and dimensions
subpad.bot = .175; % Increase if xaxis lables are cut off
subpad.top = .175; % Increase if subplot title is cut off
subpad.left = .15; % Increase if ylabel is cut off
subpad.right = .05; % Decrease if white-space on RHS
subpad.legend = 0; % Increase if legend overlaps subplot titles
plotdim = [3,3];
ndraws = size(draws,1);
xaxis = 0:ndraws-1;
titlestr = strcat(F.paramstrs,' ($', F.paramtex, '$)');

% Chains
figure(1)
set(gcf,'Position',[3,.5,10,7.5])
for iparam = 1:F.nparam
    subplot(plotdim(1),plotdim(2),iparam); grid on; hold on; box on;
    plot(xaxis,draws(:,iparam))
    title(titlestr{iparam},'interpreter','latex')
    axis('tight')
end
for iparam = 1:F.nparam
    [col,row] = ind2sub(plotdim([2 1]),iparam);
    left = (col-1+subpad.left)/plotdim(2);
    bottom = (1-(row-subpad.bot)/plotdim(1))/(1+subpad.legend);
    width = (1-(subpad.left+subpad.right))/plotdim(2);
    height = (1-subpad.bot-subpad.top)/(plotdim(1)*(1+subpad.legend));
    subplot(plotdim(1),plotdim(2),iparam);
    set(gca,'Position',[left bottom width height])
end


% Prior vs. posterior
npts = 1000;
[priorxs,priorpdfs] = pdfs(F.priors,npts);
figure(2)
subpad.top = .175; % Increase if subplot title is cut off
subpad.bot = .175; % Increase if xaxis lables are cut off
set(gcf,'Position',[1,1,10,7.5])
for iparam = 1:F.nparam
    subplot(plotdim(1),plotdim(2),iparam); grid on; hold on; box on;
    title(titlestr{iparam},'interpreter','latex')
    % Prior
    plot(priorxs(iparam,:),priorpdfs(iparam,:),'k');
    % Posterior
    [f,xi] = ksdensity(draws(:,iparam));
    plot(xi,f,'--b')
    axis('tight')
    ylims = get(gca,'ylim');
end
for iparam = 1:F.nparam
    [col,row] = ind2sub(plotdim([2 1]),iparam);
    left = (col-1+subpad.left)/plotdim(2);
    bottom = (1-(row-subpad.bot)/plotdim(1))/(1+subpad.legend);
    width = (1-(subpad.left+subpad.right))/plotdim(2);
    height = (1-subpad.bot-subpad.top)/(plotdim(1)*(1+subpad.legend));
    subplot(plotdim(1),plotdim(2),iparam);
    set(gca,'Position',[left bottom width height])
end

% Log-likelihoods
figure(3); hold on; grid on
set(gcf,'Position',[9.75,.5,6.5,3.25])
set(gca,'position',[.05,.1,.9,.85])
xmax = min(5000,ndraws-1);
xaxis = 0:xmax;
plot(xaxis,postlogliks(end-xmax:end))
axis('tight')