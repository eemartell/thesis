% Simulation comparison of normal and ZLB times
clear
clc
close all

% Saving 'on' or 'off'
saving = 'on';

% Simulation ZLB duration: 1SDx4, 1SDx7, 1p5SDx4
dur = '1p5SDx4';

% Load V from capital model
load('../options.mat','V')
V.nstate = V.nvar - V.nfore;
Vcap = V;

% Load structures for no capital model
load('../../Results-nocap/options.mat','V','P','F')
V.nstate = V.nvar - V.nfore;

%--------------------------------------------------------------------------
% True simulation
%--------------------------------------------------------------------------
T = 21;
simstemp = textread(['../../data-artificial/simulation_',dur,'.txt'],'%f');
simstemp = reshape(simstemp,[T,Vcap.nstate]);
sims{1} = simstemp(:,[Vcap.c,Vcap.n,Vcap.y,Vcap.yf,Vcap.yg,Vcap.w,...
    Vcap.pi,Vcap.i,Vcap.in,Vcap.lam,Vcap.g,Vcap.s,Vcap.mp]);

%--------------------------------------------------------------------------
% NL-PF-5% simulations
%--------------------------------------------------------------------------
nsims = 100;
nZLBpers = 6;
simstemp = textread(['../../Estimation-global/me5/simulation_',dur,'.txt'],'%f');
simstemp = reshape(simstemp,[T,V.nstate,nsims,nZLBpers]);
simstemp(:,:,51:end,:) = [];
simstemp = reshape(simstemp,[T,V.nstate,(nsims/2)*nZLBpers]);
sims{2}(:,:,1) = mean(simstemp,3);
sims{2}(:,:,2:3) = quantile(simstemp,[.05,.95],3);

%--------------------------------------------------------------------------
% PW-IF-0% simulations
%--------------------------------------------------------------------------
simstemp = textread(['../../Estimation-pwlinear/simulation_',dur,'.txt'],'%f');
simstemp = reshape(simstemp,[T,V.nstate,nsims,nZLBpers]);
simstemp(:,:,51:end,:) = [];
simstemp = reshape(simstemp,[T,V.nstate,(nsims/2)*nZLBpers]);
sims{3}(:,:,1) = mean(simstemp,3);
sims{3}(:,:,2:3) = quantile(simstemp,[.05,.95],3);

%% ------------------------------------------------------------------------
% Plot the results
%--------------------------------------------------------------------------
% close all
plotops
set(0,'DefaultLineColor','black')
set(0,'DefaultLineLineWidth',1)

%   Subplot padding
figpos = [1 1 6.5 6];
pad.topmarg = .08; % Increase if subplot title is cut off
pad.leftmarg = .08; % Increase if ylabel is cut off
pad.vspace = 0.06;
pad.hspace = .1;
pad.rightmarg = .025; % Decrease if white-space on RHS
pad.botmarg = .05; % Increase if xaxis lables are cut off
legadj = [0,.06,0,0];

ylims = {[-10,5],[-4,2],[-6,4]};
yticks = [7,7,6];

figure(1)
set(gcf,'Position',figpos)
plotvars = [V.yg,V.pi,V.in];
nplotvars = numel(plotvars);
plotdim = [nplotvars,2];
x=0:T-1;
X=[x,fliplr(x)];
titlestrs = {'NL-PF-5\%','PW-IF-0\%'};
legstrs = {'True Simulation','Mean Estimated Simulation'};
for isubplot = 1:prod(plotdim)
    [icol,irow] = ind2sub(plotdim([2 1]),isubplot);
    subplot(plotdim(1),plotdim(2),isubplot); hold on; box on; grid on;
    
    % Data to plot
    y0 = 400*(sims{1}(:,plotvars(irow))-1);
    y1 = 400*(sims{icol+1}(:,plotvars(irow),1)'-1);
    y2 = 400*(sims{icol+1}(:,plotvars(irow),2)'-1);
    y3 = 400*(sims{icol+1}(:,plotvars(irow),3)'-1);
    
    % Shade credible sets
    Y=[y3,fliplr(y2)];             
    fill(X,Y,[0.8 0.8 0.8],'linestyle','none');
    set(gca,'Layer','top')
    
    % Plot true response and median estimate
    p1(1) = plot(x,y0,'k-');
    p1(2) = plot(x,y1,'b--');
    
    % Labels and axes
    ylabel(V.desc{plotvars(irow)})
    set(gca,'ylim',ylims{irow},...
        'ytick',linspace(ylims{irow}(1),ylims{irow}(2),yticks(irow)))
    set(gca,'xlim',[0,T-1])
    if irow == 1
        title(titlestrs{icol},'interpreter','latex')
    end
    xticks(0:2:20)
end
% Resize subplots
subplotpad(gcf,pad)
% Add Data Legend
leg1 = legend(p1,legstrs,'Orientation','horizontal','Interpreter','latex','fontsize',10);
set(leg1,'Position',[.25 1-pad.topmarg 0.5 0]);
op = get(leg1,'Position');
set(leg1,'Position',op + legadj);

%% Save
if strcmp(saving,'on') 
    savename = ['simulation_' dur];
    print(gcf,'-depsc2','-painters',[savename '.eps'])
    saveas(gcf,[savename '.fig'])
    saveas(gcf,[savename '.pdf'])
end