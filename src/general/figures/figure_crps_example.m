%  Creates a figure of the percent of actuals outside of forecast credible sets
% from the different model %  specifications from the periods just prior to the ZLB binding.

clear
clc
close all
plotops

savename = 'crps_example'; 

npers = 120;
nsimsall = 100;
nsims = 50;
nZLBpers = 6;
nqs = [0,6,12,18,24,30];

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

%Load forecast distribution
forecast_dist = textread('../../Estimation-global/me5/states/forecast_dist.txt','%f');
forecast_dist = 100*(forecast_dist.^4-1);
actual = 100*(sims(89+8,V.i,2,2)^4-1);
disp(['mean forecast =    ',num2str(mean(forecast_dist))])
disp(['true realization = ',num2str(actual)])
pts = 0:0.01:8;
%%
%   Subplot padding
plotposition = [1 1 6.5 2.75];
pad.topmarg = .12; % Increase if subplot title is cut off
pad.leftmarg = .05; % Increase if ylabel is cut off
pad.vspace = 0.04;
pad.hspace = .08;
pad.rightmarg = .025; % Decrease if white-space on RHS
pad.botmarg = .15; % Increase if xaxis lables are cut off
legfont = 10;
plotdim = [1,2];

figure(1);
set(gcf,'Position',plotposition);

% PDF
subplot(1,2,1); hold on; box on; grid on;
[f,xi] = ksdensity(forecast_dist,pts,'Bandwidth',0.2); 
plot(xi,f,'k-');
%   True Realization
ylim = get(gca,'ylim');
plot([actual,actual],ylim,'k--');
set(gca,'xlim',[0,8],'ylim',ylim)
set(gca,'Layer','top')

%   Labels
% text(actual+.1,.1*diff(ylim),'True Realization')
title('PDF of 8-Quarter Ahead Forecasts','Interpreter','latex','fontsize',10);
xlabel('Nominal Interest Rate')

% CDF
subplot(1,2,2); hold on; box on; grid on;
[f,xi] = ksdensity(forecast_dist,pts,'Function','cdf','Bandwidth',0.2);
plot(xi,f,'k-');
%   Shade
y1 = f;
y1(xi<=actual) = 0;
y2 = f;
y2(xi>=actual) = 1;
X=[xi,fliplr(xi)];               
Y=[y1,fliplr(y2)];             
fill(X,Y,[0.8 0.8 0.8],'linestyle','none');
plot(xi,f,'k-');
set(gca,'Layer','top')
%   True Realization
ylim = get(gca,'ylim');
plot([actual,actual],ylim,'k--');
set(gca,'xlim',[0,8],'ylim',ylim)

% 	Labels
% text(actual+.1,.1*diff(ylim),'True Realization')
title('CDF of 8-Quarter Ahead Forecasts','Interpreter','latex','fontsize',10);
xlabel('Nominal Interest Rate')

% Resize subplots
subplotpad(gcf,pad)

saveas(gcf,[savename '.pdf'])

print(gcf,'-depsc2','-painters',[savename '.eps'])
saveas(gcf,[savename '.fig'])