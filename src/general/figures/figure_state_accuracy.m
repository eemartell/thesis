%  Creates a figure of the RMSE of RMSE of filtered state estimates across
%  ZLB bins and specification. Bottom of the code plots various simulations

clear
clc
close all

%Specifications to include: path, method, filter, measurement error,%label
specs = {...
    '../../Estimation-global/me5/','global','pf','me5','Global, Particle Filter, ME 5$\%$';  
    '../../Estimation-pwlinear/','pwlinear','if','me0','Piecewise Linear, Inversion Filter, ME 0$\%$';...
    '../../Estimation-linear/me0/','linear','kf','me0','Level Linear, Kalman Filter, ME 0$\%$';...
    };

npers = 120;
nsimsall = 100;
nsims = 50;
nZLBpers = 6;
nqs = [0,6,12,18,24,30];
nspecs = size(specs,1);
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

states = nan(npers,Vnstate,nsims,nZLBpers,nspecs);
%Loop through specifications
for ifolder = 1:size(specs,1)
    disp(['Specification: ' specs{ifolder,1}])
    %Load state estimate for given specification
    filename = [specs{ifolder,1},'/states/states.txt'];
    if exist(filename, 'file') == 2
        statestemp = textread(filename,'%f');
        statestemp = reshape(statestemp,[npers,Vnstate,nsims,nZLBpers]);
        % For global specs which may be missing, replace zeros with nans.
        if strcmp(specs{ifolder,2},'global')
            statestemp(statestemp==0) = NaN;
        end
        states(:,:,:,:,ifolder) = statestemp;
    end
end   

actual = repmat(sims,1,1,1,1,nspecs);

% Convert interest rates rate to net annualized percent and shocks to annualized
% percent
actual(:,[V.in,V.i],:,:,:) = (actual(:,[V.in,V.i],:,:,:) - 1)*400; 
states(:,[V.in,V.i],:,:,:) = (states(:,[V.in,V.i],:,:,:) - 1)*400; 
shocks = [V.s, V.g, V.mp];
actual(:,[V.s, V.g, V.mp],:,:,:) = actual(:,[V.s, V.g, V.mp],:,:,:)*400;
states(:,[V.s, V.g, V.mp],:,:,:) = states(:,[V.s, V.g, V.mp],:,:,:)*400;
V.plotnames{V.in} = [V.plotnames{V.in},', net annualized percent'];
V.plotnames{V.i} = [V.plotnames{V.i},', net annualized percent'];
V.plotnames{V.s} = [V.plotnames{V.s},', annualized percent'];
V.plotnames{V.g} = [V.plotnames{V.g},', annualized percent'];
V.plotnames{V.mp} = [V.plotnames{V.mp},', annualized percent'];

%Dummy for if ZLB is binding o rnot
zlb = repmat(actual(:,V.i,:,:,:)==0,1,Vnstate,1,1,1);
%Alternate actuals that are NaN if the ZLB is binding (or not binding)
actualzlb = actual;
actualzlb(~zlb) = NaN;
actualnzlb = actual;
actualnzlb(zlb) = NaN;

%RMSE across datasets (dim 2)
rmse = squeeze(mean(mean((states-actual).^2,1,'omitnan'),3,'omitnan').^0.5);
rmsezlb = squeeze(mean(mean((states-actualzlb).^2,1,'omitnan'),3,'omitnan').^0.5);
rmsenzlb = squeeze(mean(mean((states-actualnzlb).^2,1,'omitnan'),3,'omitnan').^0.5);

%Plot accuracy for the three categories of data
plotops
ploterrs(rmsezlb,[V.in,V.s,V.g,V.mp],V,'RMSE of state esimates for ZLB periods','state-accuracy-zlb.pdf');
ploterrs(rmsenzlb,[V.i,V.s,V.g,V.mp],V,'RMSE of state esimates for non-ZLB periods','state-accuracy-nzlb.pdf');
ploterrs(rmse,[V.in,V.s,V.g,V.mp],V,'RMSE of state esimates for all periods','state-accuracy.pdf');

%% ------------------------------------------------------------------------
% Plot lines for a given variable, simulation and ZLB bin. Change three
% values below to change what is plotted.
% -----------------------------------------------------------------------
close all
rmsetemp = squeeze(mean((states(:,:,:,:,1)-states(:,:,:,:,2)).^2,1,'omitnan').^0.5);
[~,sortidx] = sort(rmsetemp(V.in,:,6),'descend')

isim = 22;
inq = 6;

legstrs = {'Actual','NL-PF-5\%','PW-IF-0\%','Lin-KF-0\%'};
figure
set(gcf,'Position',[1 1 6.5 5.5]);

subplot(3,1,1); box on; grid on; hold on;
var = V.g;
p2 = plot([actual(:,var,isim,inq,1),states(:,var,isim,inq,1),states(:,var,isim,inq,2),states(:,var,isim,inq,3)]);
title([V.desc{var},' (',V.plotnames{var},')'],'Interpreter','latex','fontsize',10);
p2(1).Color = 'k';p2(1).LineWidth = 1.5;
p2(2).Color= 'r';p2(2).LineStyle = '--';
p2(3).Color = 'b';p2(3).LineStyle = '--';
p2(4).Color = 'g';p2(4).LineStyle = '--';
leg1 = legend(p2,legstrs,'Orientation','horizontal','Interpreter','latex','fontsize',10,'Location','South');


subplot(3,1,2); box on; grid on; hold on;
var = V.s;
p2 = plot([actual(:,var,isim,inq,1),states(:,var,isim,inq,1),states(:,var,isim,inq,2),states(:,var,isim,inq,3)]);
title([V.desc{var},' (',V.plotnames{var},')'],'Interpreter','latex','fontsize',10);
p2(1).Color = 'k';p2(1).LineWidth = 1.5;
p2(2).Color= 'r';p2(2).LineStyle = '--';
p2(3).Color = 'b';p2(3).LineStyle = '--';
p2(4).Color = 'g';p2(4).LineStyle = '--';
% leg1 = legend(p2,legstrs,'Orientation','horizontal','Interpreter','latex','fontsize',10);

subplot(3,1,3); box on; grid on; hold on;
var = V.in;
p2 = plot([actual(:,var,isim,inq,1),states(:,var,isim,inq,1),states(:,var,isim,inq,2),states(:,var,isim,inq,3)]);
title([V.desc{var},' (',V.plotnames{var},')'],'Interpreter','latex','fontsize',10);
p2(1).Color = 'k';p2(1).LineWidth = 1.5;
p2(2).Color= 'r';p2(2).LineStyle = '--';
p2(3).Color = 'b';p2(3).LineStyle = '--';
p2(4).Color = 'g';p2(4).LineStyle = '--';
% leg1 = legend(p2,legstrs,'Orientation','horizontal','Interpreter','latex','fontsize',10);


%% ------------------------------------------------------------------------
function [] = ploterrs(rmse,statesplot,V,windowname,savename)
% Plot the results
%--------------------------------------------------------------------------
% close all
    % Spec strings for legend
    legstrs = {'NL-PF-5\%','PW-IF-0\%','Lin-KF-0\%'};
    % Cat strings for x axis
    list = {'ZLB = 0q','ZLB = 6q','ZLB = 12q','ZLB = 18q','ZLB = 24q','ZLB = 30q'};
    c = categorical(list);
    c = reordercats(c,list);

    % Plot options
    plotdim = [size(statesplot,2),1];
    %   Subplot padding
    plotposition = [1 1 6.5 5.5];
    pad.topmarg = .12; % Increase if subplot title is cut off
    pad.leftmarg = .08; % Increase if ylabel is cut off
    pad.vspace = 0.06;
    pad.hspace = .10;
    pad.rightmarg = .025; % Decrease if white-space on RHS
    pad.botmarg = .05; % Increase if xaxis lables are cut off
    legfont = 10;

    figure('Name',windowname,'NumberTitle','off')
    set(gcf,'Position',plotposition);

    for isubplot = 1:prod(plotdim)
        subplot(plotdim(1),plotdim(2),isubplot); grid on; hold on; box on;
        istate = statesplot(isubplot);

        p1 = bar(c,squeeze(rmse(istate,:,:)));

        %Only label bottom xaxis
        if isubplot < prod(plotdim)
            set(gca,'xticklabel',{[]});
        end
        title([V.desc{istate},' (',V.plotnames{istate},')'],'Interpreter','latex','fontsize',10);
        set(gca,'ylim',[0,6],'ytick',0:2:6)
    end
    % Resize subplots
    subplotpad(gcf,pad)
    % Add Data Legend
    leg1 = legend(p1,legstrs,'Orientation','horizontal','Interpreter','latex','fontsize',10);
    set(leg1,'Position',[.25 1-pad.topmarg 0.5 0]);
    op = get(leg1,'Position');
    set(leg1,'Position',op + [0,.09,0,0]);

    saveas(gcf,savename)
end
