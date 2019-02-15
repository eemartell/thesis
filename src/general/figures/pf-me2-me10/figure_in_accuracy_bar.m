%  Creates a figure of the RMSE of RMSE of filtered state estimates across
%  ZLB bins and specification. Bottom of the code plots various simulations

clear
clc
close all
savename = 'in_accuracy_bar'; 
%Specifications to include: path, method, filter, measurement error,%label
specs = {...
    '../../../Estimation-global/me2/','global','pf','me2','Global, Particle Filter, ME 2$\%$';  ...
    '../../../Estimation-global/me5/','global','pf','me5','Global, Particle Filter, ME 5$\%$';  ...
    '../../../Estimation-global/me10/','global','pf','me10','Global, Particle Filter, ME 10$\%$';  ...
    '../../../Estimation-pwlinear/','pwlinear','if','me0','Piecewise Linear, Inversion Filter, ME 0$\%$';...
    };

npers = 120;
nsimsall = 100;
nsims = 50;
nZLBpers = 6;
nqs = [0,6,12,18,24,30];
nspecs = size(specs,1);
ndraws = 1000;

% Load parameter names
load('../../options.mat','V','P','F')

% Load actual data
simscap = textread('../../../data-artificial/sims.txt','%f');
simscap = reshape(simscap,[npers,V.nvar-V.nfore,nsimsall,nZLBpers]);

Vcap= V;

%Load in V from misspecified model
load('../../../Results-nocap/options.mat','V','P','F')
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
ploterrs(rmsezlb,[V.in],V,'RMSE of state esimates for ZLB periods',savename);

%% ------------------------------------------------------------------------
% Plot lines for a given variable, simulation and ZLB bin. Change three
% values below to change what is plotted.
% -----------------------------------------------------------------------

function [] = ploterrs(rmse,statesplot,V,windowname,savename)
%% ------------------------------------------------------------------------
% Plot the results
%--------------------------------------------------------------------------
% close all
    % Spec strings for legend
    legstrs = {'NL-PF-2\%','NL-PF-5\%','NL-PF-10\%','PW-IF-0\%'};
    %legstrs = {'NL-PF-5\%','PW-IF-0\%'};
    % Cat strings for x axis
    list = {'6Q','12Q','18Q','24Q','30Q'};
    c = categorical(list);
    c = reordercats(c,list);

    % Plot options
    plotdim = [size(statesplot,2),1];
    %   Subplot padding
    plotposition = [1 1 6.5 2.85];
    pad.topmarg = .05; % Increase if subplot title is cut off
    pad.leftmarg = .08; % Increase if ylabel is cut off
    pad.vspace = 0.06;
    pad.hspace = .10;
    pad.rightmarg = .025; % Decrease if white-space on RHS
    pad.botmarg = .1; % Increase if xaxis lables are cut off
    legfont = 10;

    figure('Name',windowname,'NumberTitle','off')
    set(gcf,'Position',plotposition);
    
    for isubplot = 1:prod(plotdim)
        subplot(plotdim(1),plotdim(2),isubplot); grid on; hold on; box on;
        istate = statesplot(isubplot);

        p1 = bar(squeeze(rmse(istate,2:end,:)),'group','FaceColor','flat');
        colors = {[0,0,0],[1,1,1]};
        %for ip1 = 1:numel(p1)
        %   p1(ip1).CData = colors{ip1};
        %end
        colormap(gray)
        
        %Only label bottom xaxis
        if isubplot < prod(plotdim)
            set(gca,'xticklabel',{[]});
        end
%         title(['RMSE of ',V.desc{istate},' at ZLB (',V.plotnames{istate},')'],'Interpreter','latex','fontsize',10);
        set(gca,'ylim',[0,1.75],'ytick',0:0.25:1.75)
        
        set(gca, 'XTickLabel',list, 'XTick',1:numel(list))
        %offset(1) = p1(1).XOffset
        %offset(2) = p1(2).XOffset
        offset(1) = -0.2727;
        offset(2) = -0.0909;
        offset(3) = 0.0909;
        offset(4) = 0.2727;
        
        hT = [];
        for i=1:length(p1)  % iterate over number of bar objects
            hT=[hT text(p1(i).XData+offset(i),p1(i).YData,num2str(p1(i).YData.','%.2f'), ...
                          'VerticalAlignment','bottom','horizontalalign','center','fontsize',7)];
        end
        %hatchfill2(p1(1),'single','HatchAngle',0,'HatchDensity',30,'HatchColor','w','HatchLineWidth',1.5); 
        %hatchfill2(p1(2),'single','HatchAngle',135,'HatchDensity',30,'HatchColor','w','HatchLineWidth',1.5)
    end
    % Resize subplots
    subplotpad(gcf,pad)
    % Add Data Legend
    leg1 = legend(p1,legstrs,'Orientation','horizontal',...
        'Location','Northwest','Interpreter','latex','fontsize',10);
    op = get(leg1,'Position');
    set(leg1,'Position',op + [0,.03,0,0]);

    print(gcf,'-depsc2','-painters',[savename '.eps'])
    saveas(gcf,[savename '.pdf'])
end
