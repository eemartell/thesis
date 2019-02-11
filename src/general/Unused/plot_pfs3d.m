clear all
close all
clc
plotops
set(0,'DefaultLineColor','black')
set(0,'DefaultLineLineWidth',2)

% Saving 'on' or 'off'
saving = 'off';

% Plot type
%   M: Manuscript
%   P: Presentation
plottype = 'M';

%% ------------------------------------------------------------------------
% Load/calulate decision rules
%--------------------------------------------------------------------------
% Load decision rules
load('Rules\pf_baseline_phiy25_phipi150.mat');
V = variables;

%   Production function
y = G.z_gr.*G.k_gr.^P.alpha.*pf.n.^(1-P.alpha); 
%   Interest rate rule
r = max(1,S.r*(pf.pi/P.pi).^P.phipi.*(y/S.y).^P.phiy);
disp(['ZLB bind (% nodes)' num2str(100*sum(r(:)==1)/G.nodes)]) 
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
rk = (1-P.alpha)*psi.*y./G.k_gr;
% Expected future variables
realr = 0*r;
Ec = 0*r;
Er = 0*r;
Erk = 0*r;
for i = 1:numel(realr)
    [Einvpi,Ec(i),Er(i),Erk(i)] = expvar(...
        G.k_gr(i),G.z_gr(i),G.beta_gr(i),pf.i(i),...
        P.alpha,P.pi,P.eta,P.varphi,P.phipi,P.phiy,...  
        P.delta,P.nu,P.zbar,P.rhoz,P.beta,P.rhobeta,...
        S.chi,S.r,S.y,...       
        G.e_weight,G.e_nodes,G.u_weight,G.u_nodes,...
        G.k_grid,G.z_grid,G.beta_grid,...
        pf.n,pf.pi,pf.i);
    realr(i) = r(i)*Einvpi;
end
% Investment
kp = pf.i + (1-P.delta)*G.k_gr; 

%% ------------------------------------------------------------------------
% Plot decision rules
%--------------------------------------------------------------------------
% Axes
kper = 100*(G.k_gr-S.k)/S.k;
zper = 100*(G.z_gr-P.zbar)/P.zbar;
betaper = 100*(G.beta_gr-P.beta)/P.beta;
betab = [min(betaper(:)) max(betaper(:))];
kmid = round(O.k_pts/2);
zmid = round(O.z_pts/2);
betaZLB = find(r(kmid,zmid,:)==1,1,'first');


% Plot data
idx = 0;
vars = {'r','ytil'};
%  Nominal interest rate
idx = idx+1;
titles{idx} = strcat(V.desc(V.r),' ($\tilde{r}$)');
plotme{idx} = 100*(r(:,:,betaZLB)-1);
%  Adjusted output
idx = idx+1;
titles{idx} = strcat(V.desc(V.ytil),' ($\hat{y}^{adj}$)');
plotme{idx} = 100*(ytil(:,:,betaZLB)-S.y)/S.y;
%  Consumption
idx = idx+1;
titles{idx} = strcat(V.desc(V.c),' ($\hat{c}$)');
plotme{idx} = 100*(c(:,:,betaZLB)-S.c)/S.c;
%  Investment
idx = idx+1;
titles{idx} = strcat(V.desc(V.i),' ($\hat{i}$)');
plotme{idx} = 100*(pf.i(:,:,betaZLB)-S.i)/S.i;

% Locate ZLB region
kper_zlb =  100*(G.k_grid-S.k)/S.k;
zper_zlb = 0*kper_zlb;
rper = 100*(r-1);
dellocs = [];
j = 1;
for ik = 1:numel(kper_zlb)
    loc = find(abs(rper(ik,:,betaZLB)-0)<P.tol,1,'first');
    if ~isempty(loc)
        zper_zlb(ik) = zper(ik,loc,betaZLB);
    else
        dellocs(j) = ik;
        j = j+1;
    end
end
kper_zlb(dellocs) = [];
zper_zlb(dellocs) = [];

%   Locate cross sections
kCS(:,1) = sub2ind(G.griddim,kmid*ones(1,O.k_pts),1:O.z_pts,betaZLB*ones(1,O.beta_pts));
kCS(:,2) = sub2ind(G.griddim,1:O.k_pts,1:O.z_pts,betaZLB*ones(1,O.beta_pts));

%   Position and size
if strcmp(plottype,'M')
    set(0,'DefaultAxesFontSize',12)
    legfont = 12;
    plotdim = [2 2];
    plotposition = [0 0 6.5 3];
    % subplot padding
    subpadbot = .175; % Increase if xaxis labels are cut off
    subpadtop = .15; % Increase if subplot title is cut off
    subpadleft = 0.21; % Increase if ylabel is cut off
    subpadright = .05; % Decrease if white-space on RHS
    subpadlegend = .05; % Increase if legend overlaps subplot titles

    % Plot
    figure(1); hold on;
    set(gcf,'Position',[0 0 6.5 6])
    for iplot = 1:numel(plotme)        
        [col,row] = ind2sub(plotdim([2 1]),iplot);
        left = (col-1+subpadleft)/plotdim(2);
        bottom = (1-(row-subpadbot)/plotdim(1))/(1+subpadlegend);
        width = (1-(subpadleft+subpadright))/plotdim(2);
        height = (1-subpadbot-subpadtop)/(plotdim(1)*(1+subpadlegend));
        subplot(plotdim(1),plotdim(2),iplot,'Position',[left bottom width height]); grid on; hold on; box on 
        area(kper_zlb,zper_zlb,max(zper(:)),'FaceColor',0.8*ones(1,3),'EdgeColor','none');
        set(gca,'Layer','top')
        [cl,hl] = contour(kper(:,:,betaZLB),zper(:,:,betaZLB),plotme{iplot}); clabel(cl,hl);
        title(titles{iplot},'interpreter','latex'); 
        xlabel(strcat(V.desc(V.k),' ($\hat{k}_{-1}$)'),'interpreter','latex')
        ylabel(strcat(V.desc(V.z),' ($\hat{z}_{-1}$)'),'interpreter','latex')
        axis('tight'); set(gca,'ytick',-2:2,'xtick',-5:2.5:5)
        p(1) = plot(kper(kCS([1 end],1)),zper([1 end]),'k-','linewidth',2);
        p(2) = plot(kper(kCS([1 end],2)),zper([1 end]),'b--','linewidth',2);
    end
    % Add Legend
    leg = legend(p,{'$\hat{k}_{-1} = 0$','$\hat{k}_{-1} = \hat{k}_{diag}$'},...
        'Orientation','horizontal','Interpreter','latex');
    set(leg,'Position',[.25 1-subpadlegend 0.5 subpadlegend])
    set(leg,'units','pixels')
    op = get(leg,'outerposition');
    set(leg,'outerposition',[op(1)+30,op(2),op(3)-60,.7*op(4)]); 
    obj = findobj(leg,'type','text');
    set(obj,'FontSize',legfont)

    % Save
    if strcmp(saving,'on')
        figname = 'Figs\Model2_pfs3d';
        print(gcf,'-depsc2','-painters',[figname '.eps'])
        saveas(gcf,[figname '.fig'])
    end

% Individual subplots
elseif strcmp(plottype,'P')
    ncslines = 0;
    p = 0;
    % Set local defaults
    set(0,'DefaultAxesFontSize',8)
    legfont = 8;
    plotdim = [1 2];
    plotposition = 0.7*[0 0 6.5 3.75];
    % Subplot padding
    subpadbot = .18; % Increase if xaxis labels are cut off
    subpadtop = .12; % Increase if subplot title is cut off
    subpadleft = 0.25; % Increase if ylabel is cut off
    subpadright = .05; % Decrease if white-space on RHS
    subpadlegend = .14; % Increase if legend overlaps subplot titles
    
    % Plot
    for ifig=1:2     
        figure(ifig+1); hold on; grid on; box on
        set(gcf,'Position',plotposition)
        for iplot = 1*(2*ifig-1):ifig*numel(plotme)/2        
            [col,row] = ind2sub(plotdim([2 1]),floor(iplot/ifig));
            left = (col-1+subpadleft)/plotdim(2);
            bottom = (1-(row-subpadbot)/plotdim(1))/(1+subpadlegend);
            width = (1-(subpadleft+subpadright))/plotdim(2);
            height = (1-subpadbot-subpadtop)/(plotdim(1)*(1+subpadlegend));
            subplot(plotdim(1),plotdim(2),floor(iplot/ifig),'Position',[left bottom width height]); grid on; hold on; box on 
            area(kper_zlb,zper_zlb,max(zper(:)),'FaceColor',0.8*ones(1,3),'EdgeColor','none');
            set(gca,'Layer','top')
            [cl,hl] = contour(kper(:,:,betaZLB),zper(:,:,betaZLB),plotme{iplot}); clabel(cl,hl,'FontSize',8);
            title(titles{iplot},'interpreter','latex'); 
            xlabel(strcat(V.desc(V.k),' ($\hat{k}_{-1}$)'),'interpreter','latex')
            ylabel(strcat(V.desc(V.z),' ($\hat{z}_{-1}$)'),'interpreter','latex')
            axis('tight'); set(gca,'ytick',-2:2,'xtick',-5:2.5:5)
            if ncslines == 1
                p = plot(kper(kCS([1 end],1,betaZLB)),zper([1 end]),'k-','linewidth',2);
                legstr = {'$\hat{k}_{-1} = 0$'};
            end
            if ncslines == 2
                p(1) = plot(kper(kCS([1 end],1)),zper([1 end]),'k-','linewidth',2);
                p(2) = plot(kper(kCS([1 end],2)),zper([1 end]),'b--','linewidth',2);
                legstr = {'$\hat{k}_{-1} = 0$','$\hat{k}_{-1} = \hat{k}_{diag}$'};
            end
        end
        % Add Legend
        if ncslines > 0
            leg = legend(p,legstr,'Orientation','horizontal','Position',[.25 (1-subpadlegend/2) 0.5 0],'Interpreter','latex'); 
            obj = findobj(leg,'type','text');
            set(obj,'FontSize',legfont)        
        end
        % Save
        if strcmp(saving,'on')
            figname = ['Figs\Model2_pfs3d_tech' num2str(ifig) '_ncslines=' num2str(ncslines) '_pres'];
            print(gcf,'-depsc2','-painters',[figname '.eps'])
        end    
    end
end    
