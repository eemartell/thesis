% Plot distributions of the state variables and interest rate
clear all
close all
plotops
set(0,'DefaultLineColor','black')
set(0,'DefaultLineLineWidth',2)

% Saving 'on' or 'off'
saving = 'off';

%% ------------------------------------------------------------------------
% Run simulations
%--------------------------------------------------------------------------
if exist('Rules\results.mat','file')
    % Load baseline solution
    load('Rules\results.mat')
    load([baselinename '.mat'])
    V = variables;
    if ~exist('paths','var')
        runsims = 'on';
    else
        runsims = 'off';
    end
else
    disp('Please run script_TL_onetime.m first to get baseline solution.')
    break
end
    
% Simulation parameters
nburn = 10000;
nplot = 500000;
npers = nburn+nplot;

if strcmp(runsims,'on')
    %   Production function
    y = G.z_gr.*G.k_gr.^P.alpha.*pf.n.^(1-P.alpha); 
    %   Interest rate rule     
    r = max(1,S.r*(pf.pi/P.pi).^P.phipi.*(y/S.y).^P.phiy);
    zlbbinds = length(find(r == 1))/G.nodes;
    %--------------------------------------------------------------------------
    % Simulate
    %--------------------------------------------------------------------------
    % Draw shocks
    % Set random number seed
    mtstream = RandStream('mt19937ar','seed',2);
    RandStream.setGlobalStream(mtstream);
    e = zeros(npers,1);
    u = zeros(npers,1);
    e(nburn+1:end) = P.sige*randn(nplot,1);
    u(nburn+1:end) = P.sigu*randn(nplot,1);
    % Simulate
    paths = simnonlin(pf,O,P,S,G,V,npers,e,u,0);

    % Truncate paths
    paths = paths(nburn+1:end,:);
    
    % Save paths
    save('Rules\results.mat','paths','-append')
end

%% ------------------------------------------------------------------------
% State Histograms (Unconditional vs. Conditional)
%--------------------------------------------------------------------------
% Probability of hitting ZLB
binds = paths(:,V.r) == 1;
nbinds = sum(binds(:));

% Longest ZLB event
events = 0;
ievent = 1; 
temp = 0;
for j = 2:size(paths,1)
    if binds(j) == 1 && binds(j-1) == 0
        b = j;
    elseif binds(j) == 0 && binds(j-1) == 1
        events(ievent) = j - b;
        ievent = ievent + 1;
    end
end

% Plot data
nbins = 20;
%   Capital
%       Unconditional
ks_uncond = paths(:,V.k);
[ku_bins,ku_binloc] = hist(ks_uncond(:),nbins);
%       Conditional
ks_cond = paths(binds,V.k);
[k_bins,k_binloc] = hist(ks_cond(:),nbins);
%   Productivity
%       Unconditional
zs_uncond = linspace(O.zbound(1),O.zbound(2),1000);
z_pdf = normpdf(zs_uncond,P.zbar,sqrt(P.sige^2/(1-P.rhoz^2)));
[zu_bins,zu_binloc] = hist(paths(:,V.z),nbins);
%       Conditional
zs_cond = paths(binds,V.z);
[z_bins,z_binloc] = hist(zs_cond(:),nbins);
%   Disount factor
%       Unconditional
betas_uncond = linspace(O.betabound(1),O.betabound(2),1000);
beta_pdf = normpdf(betas_uncond,P.beta,sqrt(P.sige^2/(1-P.rhobeta^2)));
[betau_bins,betau_binloc] = hist(paths(:,V.beta),nbins);
%       Conditional
betas_cond = paths(binds,V.beta);
[beta_bins,beta_binloc] = hist(betas_cond(:),nbins);
%   Interest rate
%       Unconditional
[ru_bins,ru_binloc] = hist(paths(:,V.r),nbins);

% Display key statistics
disp('Maximum net nominal interest rate')
maxr = 100*(max(paths(:,V.r))-1);
disp([num2str(maxr) '%'])
rub = ceil(100*(max(paths(:,V.r))-1))/nbins;
disp(['% Interest rates between 0 and ' num2str(rub*100) ' basis points'])
disp([num2str(100*numel(find(paths(:,V.r) <= 1+rub/100))/nplot) '%'])
disp('Minimum/Maximum values of capital (% of steady state)')
disp([min(paths(:,V.k))/S.k max(paths(:,V.k))/S.k])
disp('State space boundary (% of steady state)')
disp(O.kbound/S.k)

% Conditional histograms
figure(1)
set(gcf,'Position',[0 0 3 5.33])

subplot(3,1,1); hold on; box on;
bar(100*(k_binloc-S.k)/S.k,100*k_bins/nplot);
ylims = get(gca,'ylim');
plot(100*(O.kbound([1 1])-S.k)/S.k,ylims,'--k');
plot(100*(O.kbound([2 2])-S.k)/S.k,ylims,'--k');
xlabel('Capital')
set(gca,'xlim',[100*(O.kbound(1)-S.k)/S.k - 1 100*(O.kbound(2)-S.k)/S.k + 1]);
set(gca,'ylim',[0 0.025],'ytick',0:.01:.02,'xtick',-5:2.5:5)

subplot(3,1,2); hold on; box on;
bar(100*(z_binloc-P.zbar)/P.zbar,100*z_bins/nplot);
adj = max(100*z_bins/nplot)/max(z_pdf);
plot(100*(zs_uncond-P.zbar)/P.zbar,adj*z_pdf,'linewidth',2);
ylims = get(gca,'ylim');
plot(100*(O.zbound([1 1])-P.zbar)/P.zbar,ylims,'--k');
plot(100*(O.zbound([2 2])-P.zbar)/P.zbar,ylims,'--k');
xlabel('Technology')
set(gca,'xlim',[100*(O.zbound(1)-P.zbar)/P.zbar - 1 100*(O.zbound(2)-P.zbar)/P.zbar + 1]);
set(gca,'ylim',[0 0.03],'ytick',0:.01:.03,'xtick',-3:1:3)

subplot(3,1,3); hold on; box on;
bar(100*(beta_binloc-P.beta)/P.beta,100*beta_bins/nplot);
adj = max(100*beta_bins/nplot)/max(beta_pdf);
plot(100*(betas_uncond-P.beta)/P.beta,adj*beta_pdf,'linewidth',2);
ymax = ceil(max(adj*beta_pdf)*10)/10; set(gca,'ylim',[0 ymax]); set(gca,'ytick',0:.01:ymax)
plot(100*(O.betabound([1 1])-P.beta)/P.beta,[0 ymax],'--k');
plot(100*(O.betabound([2 2])-P.beta)/P.beta,[0 ymax],'--k');
xlabel('Discount Factor')
set(gca,'xlim',[100*(O.betabound(1)-P.beta)/P.beta - 1 100*(O.betabound(2)-P.beta)/P.beta + 1]);
set(gca,'ylim',[0 0.03],'ytick',0:.01:.03,'xtick',-2:1:2)

% Unconditional histograms
figure(2)
set(gcf,'Position',[0 0 3 5.33])

subplot(4,1,1); hold on; box on;
bar(100*(ku_binloc-S.k)/S.k,100*ku_bins/nplot);
ylims = get(gca,'ylim');
plot(100*(O.kbound([1 1])-S.k)/S.k,ylims,'--k');
plot(100*(O.kbound([2 2])-S.k)/S.k,ylims,'--k');
xlabel('Capital')
set(gca,'xlim',[100*(O.kbound(1)-S.k)/S.k - 1 100*(O.kbound(2)-S.k)/S.k + 1]);
set(gca,'ylim',[0 17.5],'ytick',0:5:15,'xtick',-5:2.5:5)

subplot(4,1,2); hold on; box on;
bar(100*(zu_binloc-P.zbar)/P.zbar,100*zu_bins/nplot);
ylims = get(gca,'ylim');
plot(100*(O.zbound([1 1])-P.zbar)/P.zbar,ylims,'--k');
plot(100*(O.zbound([2 2])-P.zbar)/P.zbar,ylims,'--k');
xlabel('Technology')
set(gca,'xlim',[100*(O.zbound(1)-P.zbar)/P.zbar - 1 100*(O.zbound(2)-P.zbar)/P.zbar + 1]);
set(gca,'ylim',[0 20],'ytick',0:10:20,'xtick',-3:1:3)

subplot(4,1,3); hold on; box on;
bar(100*(betau_binloc-P.beta)/P.beta,100*betau_bins/nplot);
ylims = get(gca,'ylim');
plot(100*(O.betabound([1 1])-P.beta)/P.beta,ylims,'--k');
plot(100*(O.betabound([2 2])-P.beta)/P.beta,ylims,'--k');
xlabel('Discount Factor')
set(gca,'xlim',[100*(O.betabound(1)-P.beta)/P.beta - 1 100*(O.betabound(2)-P.beta)/P.beta + 1]);
set(gca,'ylim',[0 20],'ytick',0:10:20)

subplot(4,1,4); hold on; box on;
bar(100*(ru_binloc-1),100*ru_bins/nplot);
ylims = get(gca,'ylim');
xlabel('Interest Rate')
axis('tight')
set(gca,'ylim',[0 17.5],'ytick',0:5:15)

% Save
if strcmp(saving,'on')
    figure(1); print(gcf,'-depsc2','-painters','Figs\Model3b_DistCond.eps')
    figure(2); print(gcf,'-depsc2','-painters','Figs\Model3b_DistUncond.eps')
end
