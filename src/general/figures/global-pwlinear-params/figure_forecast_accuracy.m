%  Creates a figure of the RMSE forecasts from the different model
%  specifications from the periods just prior to the ZLB binding.

clear
clc
close all
savename = 'forecast_accuracy'; 
%Specifications to include: path, method, filter, measurement error,%label
specs = {...
    %'../../Estimation-global/me2/','global','pf','me2','Global, Particle Filter, ME 2$\%$';  ...
    '../../../data-artificial/Ezlb_duration/','dgp','na','na','DGP';  ...
    '../../../Estimation-global/me5/states/','global','pf','me5','Global, Particle Filter, ME 5$\%$';  ...
    '../../../Estimation-pwlinear/states/','pwlinear','if','me0','Piecewise Linear, Inversion Filter, ME 0$\%$';...
    '../../../Estimation-linear/me0/states/','linear','kf','me0','Lin-KF-$0\%$';...
    '../../../global-pwlinear-params/me5/states/','global','pf','me5','NL-PF using PW-IF estimated parameters';...
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

%Load in forecasts from the specifications
temp = textread([specs{1,1},'forecast.txt'],'%f');
forecast_dgp = reshape(temp,[npers,3,4,nsimsall,nZLBpers]);

forecast = nan(npers,3,Vnstate,nsims,nZLBpers,nspecs);
forecast(:,:,V.yg,:,:,1) = forecast_dgp(:,:,1,1:nsims,:);
forecast(:,:,V.pi,:,:,1) = forecast_dgp(:,:,2,1:nsims,:);
forecast(:,:,V.i,:,:,1) = forecast_dgp(:,:,3,1:nsims,:);
forecast(:,:,V.in,:,:,1) = forecast_dgp(:,:,4,1:nsims,:);
%Loop through specifications
for ifolder = 2:size(specs,1)
    disp(['Specification: ' specs{ifolder,1}])
    %Load state estimate for given specification
    filename = [specs{ifolder,1},'forecast.txt'];
    if exist(filename, 'file') == 2
        temp = textread(filename,'%f');
        temp = reshape(temp,[npers,3,Vnstate,nsims,nZLBpers]);
        % For global specs which may be missing, replace zeros with nans.
        temp(temp==0) = NaN;
        forecast(:,:,:,:,:,ifolder) = temp;
    end
end   

%Find periods where ZLB binds
zlb = repmat(sims(:,V.i,:,:,:)==1,1,Vnstate,1,1);
% Find periods where the ZLB doesn't bind, but does the next period
prezlb = false(size(zlb));
prezlb(1:end-1,:,:,:) = ~zlb(1:end-1,:,:,:) & zlb(2:end,:,:,:); 
%prezlb(1:end-2,:,:,:) = ~zlb(1:end-2,:,:,:) & ~zlb(2:end-1,:,:,:) & zlb(3:end,:,:,:); 
%prezlb(1:end-4,:,:,:) = ~zlb(1:end-4,:,:,:) & ~zlb(2:end-3,:,:,:)  & ~zlb(3:end-2,:,:,:) & ~zlb(4:end-1,:,:,:) & zlb(5:end,:,:,:); 
%prezlb(2:end,:,:,:) = ~zlb(1:end-1,:,:,:) & zlb(2:end,:,:,:); 
%prezlb(1:end-4,:,:,:) = ~zlb(1:end-4,:,:,:)  & zlb(5:end,:,:,:); 

%Replace actual data with NaN for all other periods so it doesn't enter
%into RMSE calculations
sims(~prezlb(:)) = NaN;
%sims(~zlb(:)) = NaN;

%Assign actuals for forecasts 1,4 and 8 periods ahead
actual = NaN(npers,3,Vnstate,nsims,nZLBpers);
actual(1:end-1,1,:,:,:) = sims(2:end,:,:,:);
actual(1:end-4,2,:,:,:) = sims(5:end,:,:,:);
actual(1:end-8,3,:,:,:) = sims(9:end,:,:,:);
actual = repmat(actual,1,1,1,1,1,nspecs);

% Convert growth and interest rates to net annualized percent 
vars = [V.yg V.pi V.in,V.i];
actual(:,:,vars,:,:,:) = (actual(:,:,vars,:,:,:).^4 - 1)*100; 
forecast(:,:,vars,:,:,:) = (forecast(:,:,vars,:,:,:).^4 - 1)*100; 

V.plotnames{V.in} = [V.plotnames{V.in},', net annualized percent'];
V.plotnames{V.i} = [V.plotnames{V.i},', net annualized percent'];
V.plotnames{V.yg} = [V.plotnames{V.yg},', net annualized percent'];
V.plotnames{V.pi} = [V.plotnames{V.pi},', net annualized percent'];

%RMSE across periods (dim1) and datasets (dim 4)
rmse = squeeze(mean(mean((forecast-actual).^2,1,'omitnan'),4,'omitnan').^0.5);

%Plot accuracy for the three forecast horizons
plotops
ploterrs(squeeze(rmse(1,:,:,:)),[V.yg V.pi V.i],V,'RMSE of 1Q ahead forecasts',[savename,'1Q']);
ploterrs(squeeze(rmse(2,:,:,:)),[V.yg V.pi V.i],V,'RMSE of 4Q ahead forecasts',[savename,'4Q']);
ploterrs(squeeze(rmse(3,:,:,:)),[V.yg V.pi V.i],V,'RMSE of 8Q ahead forecasts',[savename,'8Q']);


function [] = ploterrs(rmse,statesplot,V,windowname,savename)
%% ------------------------------------------------------------------------
% Plot the results
%--------------------------------------------------------------------------
% close all
    % Spec strings for legend
    %legstrs = {'NL-PF-2\%','NL-PF-5\%','PW-IF-0\%'};
    %legstrs = {'NL-PF-5\%','PW-IF-0\%','Lin-KF-0\%'};
    legstrs = {'DGP','NL-PF-5\%','PW-IF-0\%','Lin-KF-0\%','NL-PW-params'};
    % Cat strings for x axis
    list = {'6Q','12Q','18Q','24Q','30Q'};
    c = categorical(list);
    c = reordercats(c,list);

    % Plot options
    plotdim = [size(statesplot,2),1];
    %   Subplot padding
    plotposition = [1 1 6.5 5];
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

        p1 = bar(c,squeeze(rmse(istate,2:end,:)),'FaceColor','flat');
        %for ip1 = 1:numel(p1)
        %    p1(ip1).CData = ip1;
        %end
        
        %Only label bottom xaxis
        if isubplot < prod(plotdim)
            set(gca,'xticklabel',{[]});
        end
        title([V.desc{istate},' (',V.plotnames{istate},')'],'Interpreter','latex','fontsize',10);
        %set(gca,'ylim',[0,1.75],'ytick',0:0.25:1.75)
    end
    % Resize subplots
    subplotpad(gcf,pad)
    % Add Data Legend
    leg1 = legend(p1,legstrs,'Orientation','horizontal',...
        'Location','South','Interpreter','latex','fontsize',10);
    op = get(leg1,'Position');
    set(leg1,'Position',op - [0,.120,0,0]);

    print(gcf,'-depsc2','-painters',[savename '.eps'])
    saveas(gcf,[savename '.pdf'])
end
