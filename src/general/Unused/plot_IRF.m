% IRF comparison of normal and ZLB times
clear all
close all
clc

% Saving 'on' or 'off'
saving = 'off';

%% ------------------------------------------------------------------------
% Calculate stochastic steady state
%--------------------------------------------------------------------------
% Set random number seed
mtstream = RandStream('mt19937ar','seed',0);
RandStream.setGlobalStream(mtstream);

nburn = 10000;
% Set shocks to zero
e = zeros(nburn,1);
u = zeros(nburn,1);
% Load baseline solution
if exist('Rules\results.mat','file')
    load('Rules\results.mat');
    load([baselinename '.mat']);
else
    disp('Please run script_TL_onetime.m first to get baseline solution.')
    break
end
V = variables;

% Production function
y = G.z_gr.*G.k_gr.^P.alpha.*pf.n.^(1-P.alpha); 
% Interest rate rule
r = max(1,S.r*(pf.pi/P.pi).^P.phipi.*(y./S.y).^P.phiy);
% Fix beta cross-section
betaCS = find(r(O.k_pts,ceil(O.z_pts/2),:)==1,1,'first');
beta_zlb = G.beta_grid(betaCS);

% Find stochastic steady state
paths = simnonlin(pf,O,P,S,G,V,nburn,e,u,1);
ss = paths(nburn,:);

%% ------------------------------------------------------------------------
% Impulse Responses
%--------------------------------------------------------------------------
% Number of periods
npers = 21;

% Shocks
e = zeros(npers,1);
e(2) = 0.01;
unorm = zeros(npers,1);
uZLB = zeros(npers,1);
u_fixbeta = log((beta_zlb/P.beta)^(1-P.rhobeta));
uZLB(:) = u_fixbeta;

path = zeros(npers,V.nplotvar,2);
% Fixed beta shock ss
path_zlb = simnonlin(pf,O,P,S,G,V,nburn,zeros(nburn,1),u_fixbeta*ones(nburn,1),1,ss);
ss_zlb = path_zlb(end,:);
%   Normal Times
path(:,:,1) = simnonlin(pf,O,P,S,G,V,npers,e,unorm,1,ss);
%   ZLB Times
path(:,:,2) = simnonlin(pf,O,P,S,G,V,npers,e,uZLB,1,ss_zlb);
% Store steady states
ss_all(:,1) = ss;
ss_all(:,2) = ss_zlb;

%% ------------------------------------------------------------------------
% Plot the results
%--------------------------------------------------------------------------
% Plot variables
plotvars = [V.ytil   V.n 
            V.realr  V.w 
            V.pi     V.k
            V.r      V.rk];
plottypes = [1 1;
             2 1;
             2 1
             2 1];
plotdim = size(plotvars); 
varstr(1) = strcat(V.desc(V.ytil),' ($\tilde{y}$)');
varstr(2) = strcat(V.desc(V.realr),' ($r/E[\pi]$)');
varstr(3) = strcat(V.desc(V.pi),' ($\pi$)');
varstr(4) = strcat(V.desc(V.r),' ($r$)');
varstr(5) = strcat(V.desc(V.n),' ($n$)');
varstr(6) = strcat(V.desc(V.w),' ($w$)');
varstr(7) = strcat(V.desc(V.k),' ($k$)');
varstr(8) = strcat(V.desc(V.rk),' ($r^k$)');
path_plotvar = path(:,plotvars,:);
ss_plotvar = ss_all(plotvars,:);

% Plot data
plotdata = 0*path_plotvar;
for ivar = 1:numel(plotvars)
    if plottypes(ivar) == 1
        ss_rep = permute(ss_plotvar(ivar,:,ones(npers,1)),[3 1 2]);
        plotdata(:,ivar,:) = 100*(path_plotvar(:,ivar,:)-ss_rep)./ss_rep;              
    elseif plottypes(ivar) == 2
        plotdata(:,ivar,:) = 100*(path_plotvar(:,ivar,:) - 1);
        ss_plotvar(ivar,:) = 100*(ss_plotvar(ivar,:)-1);
    end
end

%  Adjusted output
idx=find(plotvars==V.ytil);
ylims{idx} = [-1 1.25];
yticks{idx} = -1:0.5:1;
%  Real interest rate
idx=find(plotvars==V.realr);
ylims{idx} = [0.15 0.55];
yticks{idx} = 0.2:.1:0.5;
%  Inflation
idx=find(plotvars==V.pi);
ylims{idx} = [-1 .75];
yticks{idx} = -1:.5:.5;
%  Nominal interest rate
idx=find(plotvars==V.r);
ylims{idx} = [-0.05 1.15];
yticks{idx} = 0:.25:1;
%  Labor hours
idx=find(plotvars==V.n);
ylims{idx} = [-2.3 0.5];
yticks{idx} = -2:1:0;
%  Real wage rate
idx=find(plotvars==V.w);
ylims{idx} = [-1.1 .9];
yticks{idx} = -0.8:.4:.8;
%  Capital
idx=find(plotvars==V.k);
ylims{idx} = [-0.3 0.6];
yticks{idx} = -0.25:.25:0.5;
%  Real rental rate
idx=find(plotvars==V.rk);
ylims{idx} = [-3.25 1.5];
yticks{idx} = -3:1.5:1.5;

% Plot options
plotposition = [0 0 6.5 plotdim(1)*1.5];
subpad.bot = .175; % Increase if xaxis lables are cut off
subpad.top = .20; % Increase if subplot title is cut off
subpad.left = .175; % Increase if ylabel is cut off
subpad.right = .05; % Decrease if white-space on RHS
subpad.legend = .06; % Increase if legend overlaps subplot titles
fontsize = 10;

% Plot selected variables
figure(1)
set(gcf,'Position',plotposition);
xaxis = 0:npers-1;
colors = {[0 0 0] [0 0 1]};
type_line = {'-','--'};
p1 = zeros(1,2);
for ivar = 1:numel(plotvars)
    [col,row] = ind2sub(plotdim([2 1]),ivar);
    left = (col-1+subpad.left)/plotdim(2);
    bottom = (1-(row-subpad.bot)/plotdim(1))/(1+subpad.legend);
    width = (1-(subpad.left+subpad.right))/plotdim(2);
    height = (1-subpad.bot-subpad.top)/(plotdim(1)*(1+subpad.legend));
    subplot(plotdim(1),plotdim(2),ivar,'Position',[left bottom width height]); grid on; hold on; box on;
    for ipath = 1:2
        temp = plotdata(:,ivar,ipath);
        if plottypes(ivar) == 2
            plot(xaxis([1 end]),[ss_plotvar(ivar,ipath) ss_plotvar(ivar,ipath)],':','linewidth',2,'color',colors{ipath});
        end
        p1(ipath) = plot(xaxis,temp,type_line{ipath},'color',colors{ipath});  
    end
    set(gca,'ylim',ylims{ivar},'ytick',yticks{ivar});
    title(varstr{ivar},'interpreter','latex');
end
% Add Data Legend
legstr = {'Baseline','ZLB'};
leg1 = legend(p1,legstr,'Orientation','horizontal','Position',[.25 (1-subpad.legend/2) 0.5 0]);
obj1 = findobj(leg1,'type','text');
set(obj1,'FontSize',12)

% Save
if strcmp(saving,'on')
    print(gcf,'-depsc2','-painters','Figs\Model2_IRF.eps')
    saveas(gcf,'Figs\Model2_IRF.fig')
end