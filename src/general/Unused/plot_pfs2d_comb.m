clear all
close all
clc
plotops
set(0,'DefaultLineColor','black')
set(0,'DefaultLineLineWidth',2)

% Saving 'on' or 'off'
saving = 'on';

% Plot type
%   M: Manuscript
%   P: Presentation
plottype = 'M';

%% ------------------------------------------------------------------------
% Load/calulate decision rules
%--------------------------------------------------------------------------
% Load decision rules
load('Rules\pf_baseline.mat');
V = variables;

%   Production function
y = G.z_gr.*G.k_gr.^P.alpha.*pf.n.^(1-P.alpha); 
%   Interest rate rule
r = max(1,S.r*(pf.pi/P.pi).^P.phipi.*(y/S.y).^P.phiy);
% Tobin's q
q = 1+P.nu*(pf.i./G.k_gr-P.delta);  
% Adjusted output
ytil = y.*(1-(P.varphi*(pf.pi/P.pi-1).^2)/2);
% Capital adjustment costs
kac = P.nu*(pf.i./G.k_gr-P.delta).^2/2;
% Aggregate resource constraint
c = ytil - pf.i - kac.*G.k_gr;   
% FOC Labor
w = S.chi*pf.n.^P.eta.*c.^P.sigma;
% Firm FOC labor
psi = w.*pf.n./((1-P.alpha)*y);
% Firm FOC capital
rk = P.alpha*psi.*y./G.k_gr;
% Expected future variables
realr = 0*r;
Ec = 0*r;
Er = 0*r;
Erk = 0*r;
for i = 1:numel(realr)
    [Einvpi,Ec(i),Er(i),Erk(i)] = expvar(G.k_gr(i),G.z_gr(i),G.beta_gr(i),pf.i(i),...
                                         P.alpha,P.pi,P.eta,P.varphi,P.phipi,P.phiy,...  
                                         P.delta,P.nu,P.zbar,P.rhoz,P.beta,P.rhobeta,...
                                         S.chi,S.y,S.r,...       
                                         G.e_weight,G.e_nodes,G.u_weight,G.u_nodes,...
                                         G.k_grid,G.z_grid,G.beta_grid,...
                                         pf.n,pf.pi,pf.i);
    realr(i) = r(i)*Einvpi;
end
% Investment
kp = pf.i + (1-P.delta)*G.k_gr; 

%% ------------------------------------------------------------------------
% Plot decision rules (compare steady state and diagonal cross sections)
%--------------------------------------------------------------------------
% Fix capital cross section
%  steady state
betaCS = find(r(ceil(O.k_pts/2),ceil(O.z_pts/2),:)==1,1,'first');
kCS(:,1) = sub2ind(G.griddim,round(O.k_pts/2)*ones(1,O.k_pts),1:O.z_pts,betaCS*ones(1,O.beta_pts));

betaCS(1) = O.z_pts;
for i=2:O.k_pts
    betaCS(i) = find(r(i,ceil(O.z_pts/2),:)==1,1,'first');
end
kCS(:,2) = sub2ind(G.griddim,1:O.k_pts,1:O.z_pts,betaCS);

% Plot data
idx = 0;
%  Adjusted output
idx = idx+1;
titles{idx} = strcat(V.desc(V.ytil),' ($\hat{y}^{adj}$)');
plotme1{idx} = 100*(ytil(kCS(:,1))-S.y)/S.y;
plotme2{idx} = 100*(ytil(kCS(:,2))-S.y)/S.y;
ylims{idx} = [-3 2];
yticks{idx} = -3:1:2;    
%  Real interest rate
idx = idx+1;
titles{idx} = strcat(V.desc(V.realr),' ($\widetilde{r/E[\pi]}$)');
plotme1{idx} = 100*(realr(kCS(:,1))-1);
plotme2{idx} = 100*(realr(kCS(:,2))-1);
ylims{idx} = [0 1.25];
yticks{idx} = 0:0.25:1.25;    
%  Inflation
idx = idx+1;
titles{idx} = strcat(V.desc(V.pi),' ($\tilde{\pi}$)');
plotme1{idx} = 100*(pf.pi(kCS(:,1))-1);
plotme2{idx} = 100*(pf.pi(kCS(:,2))-1);
ylims{idx} = [-2 0.25];
yticks{idx} = -2:1:0;  
%  Nominal interest rate
idx = idx+1;
titles{idx} = strcat(V.desc(V.r),' ($\tilde{r}$)');
plotme1{idx} = 100*(r(kCS(:,1))-1);
plotme2{idx} = 100*(r(kCS(:,2))-1);
ylims{idx} = [0 0.4];
yticks{idx} = 0:.1:.4;  
%  Consumption
idx = idx+1;
titles{idx} = strcat(V.desc(V.c),' ($\hat{c}$)');
plotme1{idx} = 100*(c(kCS(:,1))-S.c)/S.c;
plotme2{idx} = 100*(c(kCS(:,2))-S.c)/S.c;
ylims{idx} = [-10 0];
yticks{idx} = -10:2:0; 
%  Investment
idx = idx+1;
titles{idx} = strcat(V.desc(V.i),' ($\hat{i}$)');
plotme1{idx} = 100*(pf.i(kCS(:,1))-S.i)/S.i;
plotme2{idx} = 100*(pf.i(kCS(:,2))-S.i)/S.i;
ylims{idx} = [-6 18];
yticks{idx} = -6:6:18;
%  Labor hours
idx = idx+1;
titles{idx} = strcat(V.desc(V.n),' ($\hat{n}$)');
plotme1{idx} = 100*(pf.n(kCS(:,1))-P.n)/P.n;
plotme2{idx} = 100*(pf.n(kCS(:,2))-P.n)/P.n;
ylims{idx} = [-6 3];
yticks{idx} = -6:3:3;     
%  Real rental rate
idx = idx+1;
titles{idx} = strcat(V.desc(V.rk),' ($r^k$)');
plotme1{idx} = 100*rk(kCS(:,1));
plotme2{idx} = 100*rk(kCS(:,2));
ylims{idx} = [2.5 3.10];
yticks{idx} = 2.5:0.2:3.1;  

% %  Real wage rate
% idx = idx+1;
% titles{idx} = strcat(V.desc(V.w),' ($w$)');
% plotme1{idx} = 100*(w(kCS(:,1))-S.w)/S.w;
% plotme2{idx} = 100*(w(kCS(:,2))-S.w)/S.w;
% ylims{idx} = [-13 4];
% yticks{idx} = -12:4:4;   
% %  Capital
% idx = idx+1;
% titles{idx} = strcat(V.desc(V.k),' ($k$)');
% plotme1{idx} = 100*(kp(kCS(:,1))-S.k)/S.k;
% plotme2{idx} = 100*(kp(kCS(:,2))-S.k)/S.k;
% ylims{idx} = [-6 4.5];
% yticks{idx} = -6:3:3; 

% Plot
%   Axis and bounds
zper = 100*(G.z_grid-P.zbar)/P.zbar;
zbound = [min(zper) max(zper)];
%   Find ZLB region
ZLB_ss = find(r(kCS(:,1)) == 1,1,'first');
shadeint_ss = [zper(ZLB_ss) zper(end)];
ZLB_diag = find(r(kCS(:,2)) == 1,1,'first');
shadeint_diag = [zper(ZLB_diag) zper(ZLB_ss)];

%   Position and size

set(0,'DefaultAxesFontSize',12)
legfont = 12;
plotdim = [4 2];
%   Subplot padding
subpadbot = .275; % Increase if xaxis labels are cut off
subpadtop = .15; % Increase if subplot title is cut off
subpadleft = .1; % Increase if ylabel is cut off
subpadright = .05; % Decrease if white-space on RHS
subpadlegend = .07; % Increase if legend overlaps subplot titles

% Plot
figure(1); hold on;
set(gcf,'Position',[0 0 6.5 8])
for iplot = 1:numel(plotme1)
    [col,row] = ind2sub(plotdim([2 1]),iplot);
    left = (col-1+subpadleft)/plotdim(2);
    bottom = (1-(row-subpadbot)/plotdim(1))/(1+subpadlegend);
    width = (1-(subpadleft+subpadright))/plotdim(2);
    height = (1-subpadbot-subpadtop)/(plotdim(1)*(1+subpadlegend));
    subplot(plotdim(1),plotdim(2),iplot,'Position',[left bottom width height]); grid on; hold on; box on 
    set(gca,'xlim',zbound,'xtick',-2:1:2);
    set(gca,'ylim',ylims{iplot},'ytick',yticks{iplot});
    area(shadeint_ss,ones(2,1)*max(ylims{iplot}),min(ylims{iplot}),...
        'FaceColor',[0.5 0.5 0.5],'EdgeColor','none');
    area(shadeint_diag,ones(2,1)*max(ylims{iplot}),min(ylims{iplot}),...
        'FaceColor',[0.8 0.8 0.8],'EdgeColor','none');
    set(gca,'Layer','top')
    p(1) = plot(zper,plotme1{iplot},'-k');
    p(2) = plot(zper,plotme2{iplot},'--b');
    title(titles{iplot},'interpreter','latex');
    xlabel(strcat(V.desc(V.z),' ($\hat{z}_{-1}$)'),'interpreter','latex');
end
% Add Legend
leg = legend(p,{'$\hat{k}_{-1} = 0$','$\hat{k}_{-1} = \hat{k}_{diag}$'},'Orientation','horizontal','Position',[.25 (1-subpadlegend/2) 0.5 0],'Interpreter','latex'); 
obj = findobj(leg,'type','text');
set(obj,'FontSize',12)

% Save
if strcmp(saving,'on')
    figname = 'Figs\Model2_pfs_tech';
    print(gcf,'-depsc2','-painters',[figname '.eps'])
    saveas(gcf,[figname '.fig'])
end
    