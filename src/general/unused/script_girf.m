% GIRF comparison of normal and ZLB times
clear
clear mex
clc

% % Load GIRFs, otherwise run
% loadflag = 0;

% Set random number seed
mtstream = RandStream('mt19937ar','seed',0);
RandStream.setGlobalStream(mtstream);

% Number of periods and simulations
nplot = 21;     % i.e., 5 years
nmcsims = 1000;  % Number of Monte Carlo simulations for GIRF
shocksize = -2; % 2SD negative shock
ninitstate = 2;     % Initial state vectors
npers = 2000;  % Periods for stochastic SS and ergodic distribution

% Artificial dataset parameters
nZLBper = 6;
nsims = 3;

% Load V from capital model
load('../Results/options.mat','V','S')
V.nstate = V.nvar - V.nfore;
Vcap = V;
Scap = S;

% Load true parameters and solution
load('options.mat')
load('solution_test.mat')
V = variables;

% if ~loadflag
    %--------------------------------------------------------------------------
    % Set initial states using true solution
    %--------------------------------------------------------------------------
%     initstate = zeros(V.nplotvar,ninitstate);
%     %   Stochastic steady state (zero shocks)
%     epsg = zeros(npers,1);
%     epss = zeros(npers,1);
%     epsmp = zeros(npers,1);
%     sims = simulation(pf,P,S,G,V,epsg,epss,epsmp);
% %     sims = Fsimulation_test(...
% %         zeros(V.nplotvar,1),...
% %         P.eta,P.varphi,P.phipi,P.phiy,...
% %         P.h,P.rhog,P.rhoa,P.rhoi,P.rhomp,P.sige,P.sigu,P.sigv,P.g,P.pi,...
% %         S.c,S.chi,S.i,...
% %         V.n,V.w,V.i,V.in,V.c,V.g,V.a,V.y,V.pi,V.cg,V.mp,...
% %         G.g_grid,G.a_grid,G.mp_grid,G.in_grid,G.c_grid,...
% %         pf.c,pf.pigap,...
% %         e,u,v,P.zlbflag);
%     initstate(:,1) = sims(end,:);
%     
%     %   State after 4x1.5SD epss shocks
%     epsg = zeros(npers,1);
%     epss = zeros(npers,1);
%     epss(1:5) = 1.5;
%     epsmp = zeros(npers,1);
%     sims = simulation(pf,P,S,G,V,epsg,epss,epsmp);
% %     sims = Fsimulation_test(...
% %         zeros(V.nplotvar,1),...
% %         P.eta,P.varphi,P.phipi,P.phiy,...
% %         P.h,P.rhog,P.rhoa,P.rhoi,P.rhomp,P.sige,P.sigu,P.sigv,P.g,P.pi,...
% %         S.c,S.chi,S.i,...
% %         V.n,V.w,V.i,V.in,V.c,V.g,V.a,V.y,V.pi,V.cg,V.mp,...
% %         G.g_grid,G.a_grid,G.mp_grid,G.in_grid,G.c_grid,...
% %         pf.c,pf.pigap,...
% %         e,u,v,P.zlbflag);
%     initstate(:,2) = sims(5,:);
    
    initstate_dgp = textread('../data-artificial/girf_initstate_nlpf.txt','%f');
    initstate_dgp = reshape(initstate_dgp,[23,2]);
    initstate(V.c,:) = initstate_dgp(Vcap.c,:);
	initstate(V.n,:) = initstate_dgp(Vcap.n,:);
	initstate(V.y,:) = initstate_dgp(Vcap.y,:);
	initstate(V.yf,:) = initstate_dgp(Vcap.yf,:);
	initstate(V.yg,:) = initstate_dgp(Vcap.yg,:);
	initstate(V.w,:) = initstate_dgp(Vcap.w,:);
	initstate(V.pi,:) = initstate_dgp(Vcap.pi,:);
	initstate(V.i,:) = initstate_dgp(Vcap.i,:);
	initstate(V.in,:) = initstate_dgp(Vcap.in,:);
	initstate(V.lam,:) = initstate_dgp(Vcap.lam,:);
	initstate(V.g,:) = initstate_dgp(Vcap.g,:);
	initstate(V.s,:) = initstate_dgp(Vcap.s,:);
	initstate(V.mp,:) = initstate_dgp(Vcap.mp,:);
    
    %--------------------------------------------------------------------------
    % Compute true GIRF
    %--------------------------------------------------------------------------
    baseline = randn(nplot,nmcsims,O.nshocks);
    paths_temp = zeros(nplot,V.nplotvar,nmcsims,2);
    meanpaths_true = zeros(nplot,V.nplotvar,2,ninitstate,O.nshocks);
    for ishock = 1:O.nshocks
        for istate = 1:ninitstate
            % Baseline
            paths_temp(:,:,:,1) = simulation(...
                pf,P,S,G,V,baseline(:,:,1),baseline(:,:,2),baseline(:,:,3),initstate(:,istate));
%             paths_temp(:,:,:,istate,ishock,1) = Fsimulation_test(...
%                 initstate(:,istate),...
%                 P.eta,P.varphi,P.phipi,P.phiy,...
%                 P.h,P.rhog,P.rhoa,P.rhoi,P.rhomp,P.sige,P.sigu,P.sigv,P.g,P.pi,...
%                 S.c,S.chi,S.i,...
%                 V.n,V.w,V.i,V.in,V.c,V.g,V.a,V.y,V.pi,V.cg,V.mp,...
%                 G.g_grid,G.a_grid,G.mp_grid,G.in_grid,G.c_grid,...
%                 pf.c,pf.pigap,...
%                 baseline(:,:,1),baseline(:,:,2),baseline(:,:,3),P.zlbflag);
            % Impulse
            impulse = baseline;
            impulse(2,:,ishock) = shocksize;
            paths_temp(:,:,:,2) = simulation(...
                pf,P,S,G,V,impulse(:,:,1),impulse(:,:,2),impulse(:,:,3),initstate(:,istate));
%             paths_temp(:,:,:,istate,ishock,2) = Fsimulation_test(...
%                 initstate(:,istate),...
%                 P.eta,P.varphi,P.phipi,P.phiy,...
%                 P.h,P.rhog,P.rhoa,P.rhoi,P.rhomp,P.sige,P.sigu,P.sigv,P.g,P.pi,...
%                 S.c,S.chi,S.i,...
%                 V.n,V.w,V.i,V.in,V.c,V.g,V.a,V.y,V.pi,V.cg,V.mp,...
%                 G.g_grid,G.a_grid,G.mp_grid,G.in_grid,G.c_grid,...
%                 pf.c,pf.pigap,...
%                 impulse(:,:,1),impulse(:,:,2),impulse(:,:,3),P.zlbflag);

            % Take average across sims
            meanpaths_true(:,:,:,istate,ishock) = squeeze(mean(paths_temp,3)); 
        end
    end
    
    %--------------------------------------------------------------------------
    % Posterior means
    %--------------------------------------------------------------------------
    postmeans = zeros(F.nparam,nsims,nZLBper);
    for iZLBper = 1:nZLBper
        for taskid = 1:nsims
            fname = ...
                ['../Estimation-global/me5/zlb' num2str(iZLBper) '/mean/' ...
                num2str(taskid) '-postmean.txt'];
            postmeans(:,taskid,iZLBper) = textread(fname,'%f');
        end
    end

    %--------------------------------------------------------------------------
    % Compute IRFs/GIRFs
    %--------------------------------------------------------------------------
    meanpaths = zeros(nplot,V.nplotvar,2,ninitstate,O.nshocks,nsims,nZLBper);
    for iZLBper = 1:nZLBper
        disp(['iZLBper: ' num2str(iZLBper)])
        for taskid = 1:nsims
            disp(['taskid: ' num2str(taskid)])
            
            postmean = postmeans(:,taskid,iZLBper);
            % Solve model at posterior mean
            %   Parameters, steady state, grids, and initial conjectures
            for iparam = 1:F.nparam
                eval(['P.' F.priors{iparam,1} ' = postmean(iparam);']);
            end
            S = steadystate(P);
            G = grids(O,P,S);
            pf = guess(P,S,G);
            [pf.c(:),pf.pigap(:)] = Fsolution_test( ...
                O.nstates,O.npfs,O.nshocks,...
                P.tol,P.beta,P.thetap,P.eta,P.g,P.pi,P.s,...
                P.varphip,P.h,P.rhos,P.rhoi,P.phipi,P.phiy,...
                S.chi,S.i,...
                G.g_grid,G.s_grid,G.mp_grid,G.in_grid,G.c_grid,...
                G.e_weight,G.e_nodes,...
                G.u_weight,G.u_nodes,...
                G.v_weight,G.v_nodes,...
                pf.c,pf.pigap,...
                G.g_gr,G.mp_gr,G.in_gr,G.c_gr);
            paths_temp = zeros(nplot,V.nplotvar,nmcsims,2);
            for ishock = 1:O.nshocks
                impulse = baseline;
                impulse(2,:,ishock) = shocksize;
                for istate = 1:ninitstate
                    % Baseline
                    paths_temp(:,:,:,1) = simulation(...
                        pf,P,S,G,V,baseline(:,:,1),baseline(:,:,2),baseline(:,:,3),initstate(:,istate));
%                     paths_temp(:,:,:,istate,ishock,1) = Fsimulation_test(...
%                         initstate(:,istate),...
%                         P.eta,P.varphi,P.phipi,P.phiy,...
%                         P.h,P.rhog,P.rhoa,P.rhoi,P.rhomp,P.sige,P.sigu,P.sigv,P.g,P.pi,...
%                         S.c,S.chi,S.i,...
%                         V.n,V.w,V.i,V.in,V.c,V.g,V.a,V.y,V.pi,V.cg,V.mp,...
%                         G.g_grid,G.a_grid,G.mp_grid,G.in_grid,G.c_grid,...
%                         pf.c,pf.pigap,...
%                         baseline(:,:,1),baseline(:,:,2),baseline(:,:,3),P.zlbflag);
                    % Impulse
                    paths_temp(:,:,:,2) = simulation(...
                        pf,P,S,G,V,impulse(:,:,1),impulse(:,:,2),impulse(:,:,3),initstate(:,istate));
%                     paths_temp(:,:,:,istate,ishock,2) = Fsimulation_test(...
%                         initstate(:,istate),...
%                         P.eta,P.varphi,P.phipi,P.phiy,...
%                         P.h,P.rhog,P.rhoa,P.rhoi,P.rhomp,P.sige,P.sigu,P.sigv,P.g,P.pi,...
%                         S.c,S.chi,S.i,...
%                         V.n,V.w,V.i,V.in,V.c,V.g,V.a,V.y,V.pi,V.cg,V.mp,...
%                         G.g_grid,G.a_grid,G.mp_grid,G.in_grid,G.c_grid,...
%                         pf.c,pf.pigap,...
%                         impulse(:,:,1),impulse(:,:,2),impulse(:,:,3),P.zlbflag);

                    % Take average across sims
                    meanpaths(:,:,:,istate,ishock,taskid,iZLBper) = ...
                        squeeze(mean(paths_temp,3)); 
                end
            end
        end
    end
%     save('script_girf3.mat',...
%         'initstate','truth','linkf0','nlpf5','nlpf5b','nlpf5b_keep',...
%         'idxzlb15','nzlb15')
% else
%     load('script_girf3.mat')
% end

% Compute GIRF
girf_true = 400*squeeze(diff(meanpaths_true,1,3));
girfs = reshape(400*squeeze(diff(meanpaths,1,3)),...
        [nplot,V.nplotvar,ninitstate,O.nshocks,nsims*nZLBper]);

%% ------------------------------------------------------------------------
% Plot the results
%--------------------------------------------------------------------------
% close all
plotops
set(0,'DefaultLineColor','black')
set(0,'DefaultLineLineWidth',1)

% Saving 'on' or 'off'
saving = 'on';

% GIRFs (percentage point difference from mean baseline)
varidx = [V.yg,V.pi,V.in];

% Plot options
plotdim = [3,2];
legstrs = {'True Response','Estimated Median Response'};
% titlestrs = {'Steady State','Recession State'};
titlestrs = {'Growth Shock','Interest Rate Shock'};
ishocks = [1,3];
statestrs = {[ 'Steady State ($i_0^n =' num2str(round(1e3*(initstate(V.in,1)-1))/10) '\%$)'],...
          [ 'ZLB State ($i_0^n =' num2str(round(1e3*(initstate(V.in,2)-1))/10) '\%$)']};
%   Adjust ylims and yticks
ylims{1} = [-1.5,.5];
yticks{1} = ylims{1}(1):.5:ylims{1}(2);
ylims{2} = [-2,.5];
yticks{2} = ylims{2}(1):.5:ylims{2}(2);
%   Subplot padding
plotposition = [1,1,6.5,6.75];
pad.topmarg = .08; % Increase if subplot title is cut off
pad.leftmarg = .08; % Increase if ylabel is cut off
pad.vspace = 0.06;
pad.hspace = .1;
pad.rightmarg = .025; % Decrease if white-space on RHS
pad.botmarg = .05; % Increase if xaxis lables are cut off
legadj = [0,.06,0,0];
legfont = 10;

% Quantiles
probs = [.5,.05,.95];
girfs_quants = quantile(girfs,probs,5);

figure(2)
set(gcf,'Position',plotposition);
xaxis = 0:(17-1);
colors = {[0,0,0],[0,0,1],[1,0,0]};
type_line = {'-','--','-.'};
p1 = zeros(2,1);
istate = 2; % Recession state
% ishock = 2; % Preference shock
for isubplot = 1:prod(plotdim)
    [icol,irow] = ind2sub(plotdim([2,1]),isubplot);
    ivar = irow;
    ishock = ishocks(icol);
    subplot(plotdim(1),plotdim(2),isubplot); grid on; hold on; box on;
    
    % Estimated Credible sets
    y1 = girfs_quants(xaxis+1,varidx(ivar),istate,ishock,2)';
    y2 = girfs_quants(xaxis+1,varidx(ivar),istate,ishock,3)';
    X = [xaxis,fliplr(xaxis)];
    Y = [y1,fliplr(y2)];
    fill(X,Y,.82*ones(1,3),'LineStyle','none');
    
    % Truth
    temp = girf_true(:,varidx(ivar),istate,ishock);
    p1(1) = plot(xaxis,temp(xaxis+1),type_line{1},'color',colors{1});  
    % Estimated  Median
    temp = girfs_quants(:,varidx(ivar),istate,ishock,1);
    p1(2) = plot(xaxis,temp(xaxis+1),type_line{2},'color',colors{2});

    set(gca,'xlim',xaxis([1,end]),'xtick',xaxis(1:4:end));
%     set(gca,'ylim',ylims{istate},'ytick',yticks{istate});
    ylabel(V.desc(varidx(irow)),'interpreter','latex','fontsize',10);
    title(titlestrs{icol},'Interpreter','latex','fontsize',10);
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
    savename = 'Figs\girf3';
    print(gcf,'-depsc2','-painters',[savename '.eps'])
    saveas(gcf,[savename '.fig'])
end