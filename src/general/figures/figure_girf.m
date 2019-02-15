% Simulation comparison of normal and ZLB times
clear
clc
close all

% Number of periods and simulations
nplot = 21;         % i.e., 5 years
ninitstate = 2;     % Initial state vectors

% Load V from capital model
load('../options.mat','V','S')
V.nstate = V.nvar - V.nfore;
Vcap = V;
Scap = S;

% Load structures for no capital model
load('../../Results-nocap/options.mat','V','P','S')
V.nstate = V.nvar - V.nfore;

%--------------------------------------------------------------------------
% True GIRF paths
%--------------------------------------------------------------------------
initstate = textread('../../data-artificial/girf_initstate_nlpf.txt','%f');
initstate_nlpf = reshape(initstate,[Vcap.nstate,2]);
initstate = textread('../../data-artificial/girf_initstate_pwlin.txt','%f');
initstate_pwif = reshape(initstate,[Vcap.nstate,2]);
meanpaths_true = textread('../../data-artificial/girf_meanpaths.txt','%f');
meanpaths_true = reshape(meanpaths_true,[nplot,Vcap.nstate,2,ninitstate,V.nshock]);

%--------------------------------------------------------------------------
% NL-PF-5% GIRF paths
%--------------------------------------------------------------------------
nsimstemp = 100;
nsims = 50;
% nsims = 10;
nZLBpers = 6;
meanpathstemp = textread('../../Estimation-global/me5/girf_meanpaths.txt','%f');
meanpathstemp = reshape(meanpathstemp,[nplot,V.nstate,2,ninitstate,V.nshock,nsimstemp,nZLBpers]);
meanpathstemp(:,:,:,:,:,nsims+1:end,:) = [];
meanpaths{1} = reshape(meanpathstemp,[nplot,V.nstate,2,ninitstate,V.nshock,nsims*nZLBpers]);
% meanpaths{1} = reshape(meanpathstemp(:,:,:,:,:,:,[1,2]),[nplot,V.nstate,2,ninitstate,V.nshock,nsims*2]);

%--------------------------------------------------------------------------
% PW-IF-0% GIRF paths
%--------------------------------------------------------------------------
meanpathstemp = textread('../../Estimation-pwlinear/girf_meanpaths.txt','%f');
meanpathstemp = reshape(meanpathstemp,[nplot,V.nstate,2,ninitstate,V.nshock,nsimstemp,nZLBpers]);
meanpathstemp(:,:,:,:,:,nsims+1:end,:) = [];
meanpaths{2} = reshape(meanpathstemp,[nplot,V.nstate,2,ninitstate,V.nshock,nsims*nZLBpers]);
% temp = load('../../Estimation-pwlinear-MATLAB/girf_meanpaths.mat');
% meanpaths{2} = reshape(temp.meanpaths,[nplot,V.nstate,2,ninitstate,V.nshock,nsims*nZLBpers]);
% meanpaths{2} = reshape(meanpathstemp(:,:,:,:,:,:,[1,2]),[nplot,V.nstate,2,ninitstate,V.nshock,nsims*2]);

% Compute GIRF
girf_true = 400*squeeze(diff(meanpaths_true,1,3));
girfs{1} = 400*squeeze(diff(meanpaths{1},1,3));
girfs{2} = 400*squeeze(diff(meanpaths{2},1,3));

%% ------------------------------------------------------------------------
% Plot the results
%--------------------------------------------------------------------------
% close all
plotops
set(0,'DefaultLineColor','black')
set(0,'DefaultLineLineWidth',1)

% Saving 'on' or 'off'
saving = 'on';
%   What to plot
istate = 2;         % 1) Steady state, 2) Recession/ZLB state
ishocks = [1,2,3];    % 1) Growth, 2) Risk Premium, 3) Interest Rate shocks

% GIRFs to plot
varidx_true = [Vcap.yg,Vcap.pi,Vcap.in];
varidx = [V.yg,V.pi,V.in];
probs = [.05,.95];

% Plot options
plotdim = [numel(varidx),2];
legstrs = {'True Response','Mean Estimated Response'};
shockstrs = {'epsg','epss','epsi'};
methodstrs = {'NL-PF-$5\%$','PW-IF-$0\%$'};
colors = {[0,0,0],[0,0,1],[1,0,0]};
type_line = {'-','--','-.'};
%   Adjust lims and ticks
ylims{1} = {[-1.75,.25],[-.1,.35],[-.2,.1]};
yspace{1} = [.25,.1,.05];
ylims{2} = {[-.5,2],[0,.7],[-1.5,.1]};
yspace{2} = [.5,.1,.25];
xaxis = 0:(17-1);
%   Subplot padding
plotposition = [1,1,6.5,6];
pad.topmarg = .08; % Increase if subplot title is cut off
pad.leftmarg = .08; % Increase if ylabel is cut off
pad.vspace = 0.06;
pad.hspace = .1;
pad.rightmarg = .025; % Decrease if white-space on RHS
pad.botmarg = .05; % Increase if xaxis lables are cut off
legadj = [0,.06,0,0];
legfont = 10;

for ifig = 1:numel(ishocks)
    ishock = ishocks(ifig);
    figure(ifig)
    set(gcf,'Position',plotposition);
    p1 = zeros(2,1);
    for isubplot = 1:prod(plotdim)
        [icol,irow] = ind2sub(plotdim([2,1]),isubplot);
        ivar = irow;
        imethod = icol;
        subplot(plotdim(1),plotdim(2),isubplot); grid on; hold on; box on;

        % Estimated Credible sets
        girfs_quants(:,:,:,:,1) = mean(girfs{imethod},5);
        girfs_quants(:,:,:,:,2:3) = quantile(girfs{imethod},probs,5);
        y1 = girfs_quants(xaxis+1,varidx(ivar),istate,ishock,2)';
        y2 = girfs_quants(xaxis+1,varidx(ivar),istate,ishock,3)';
        X = [xaxis,fliplr(xaxis)];
        Y = [y1,fliplr(y2)];
        fill(X,Y,.82*ones(1,3),'LineStyle','none');

        % Truth
        temp = girf_true(:,varidx_true(ivar),istate,ishock);
        p1(1) = plot(xaxis,temp(xaxis+1),type_line{1},'color',colors{1});  
        % Estimated  Median
        temp = girfs_quants(:,varidx(ivar),istate,ishock,1);
        p1(2) = plot(xaxis,temp(xaxis+1),type_line{2},'color',colors{2});

        set(gca,'xlim',xaxis([1,end]),'xtick',xaxis(1:4:end));
        ylabel(V.desc(varidx(ivar)),'interpreter','latex','fontsize',10);
%         set(gca,'ylim',ylims{ifig}{ivar},...
%             'ytick',ylims{ifig}{ivar}(1):yspace{ifig}(ivar):ylims{ifig}{ivar}(2))
        if irow == 1
            title(methodstrs{imethod},'Interpreter','latex','fontsize',10);
        end
    end
    % Resize subplots
    subplotpad(gcf,pad)
    % Legend
    leg1 = legend(p1,legstrs,'Orientation','horizontal','Interpreter','latex','fontsize',10);
    set(leg1,'Position',[.25 1-pad.topmarg 0.5 0]);
    op = get(leg1,'Position');
    set(leg1,'Position',op + legadj);
    
    % Save
    if strcmp(saving,'on')
        savename = ['girf_' num2str(shockstrs{ishock})];
        print(gcf,'-depsc2','-painters',[savename '.eps'])
        saveas(gcf,[savename '.fig'])
    end
end