% Plot histogram of Euler equation errors
clear
clc
close all

plotops
saving = 'on';
plottype = 'P';

% Iteration
%   ti: time iteration
%   fp: fixed point
O.it = 'fp';

% Figure name
figname = ['eeerrors' O.it 'bothalgs'];

% Load solution at posterior mean
load('../solutions/eeerrors_simfpART_5')
EE_ART(:,1) = R.EE1(:);
EE_ART(:,2) = R.EE2(:);
EE_ART(:,3) = R.EE3(:);
EE_ART(:,4) = R.EE4(:);
RMeanEEART(:,1) = R.meanEE;
RMaxEEART(:,1) = R.maxEE;

load('../solutions/eeerrors_simfpGust_5')
EE_Gust(:,1) = R.EE1(:);
EE_Gust(:,2) = R.EE3(:);
EE_Gust(:,3) = R.EE4(:);
EE_Gust(:,4) = R.EE5(:);

RMeanEEGust(:,1) = R.meanEE(:,1);
RMeanEEGust(:,2) = R.meanEE(:,3);
RMeanEEGust(:,3) = R.meanEE(:,4);
RMeanEEGust(:,4) = R.meanEE(:,5);

RMaxEEGust(:,1) = R.maxEE(:,1);
RMaxEEGust(:,2) = R.maxEE(:,3);
RMaxEEGust(:,3) = R.maxEE(:,4);
RMaxEEGust(:,4) = R.maxEE(:,5);

bins{1} = -6:.125:-1.25;
bins{2} = -7.25:.125:-2.25;
bins{3} = -4.5:.125:-.5;
bins{4} = -4.5:.125:-.5;

nEE = size(EE_ART,1);

titles = {'HH FOC bond','HH FOC capital','HH FOC investment','Price Phillips Curve','HH FOC Bond','HH FOC capital','HH FOC investment', 'Price Phillips Curve'};

% Plot options
if strcmp(plottype,'M')
    savename = [figname];
    figbox = [1,1,6.5,2.5];
    subpad.bot = .175; % Increase if xaxis lables are cut off
    subpad.top = .125; % Increase if subplot title is cut off
    subpad.left = .075; % Increase if ylabel is cut off
    subpad.right = .025; % Decrease if white-space on RHS 
    subpad.legend = 0; % Increase if legend overlaps subplot titles
    fontsize = 10; 
elseif strcmp(plottype,'P')
    savename = [figname '_pres'];
    set(0,'DefaultAxesFontSize',8)
    figbox = [1,1,4.5,2.5];
    subpad.bot = .175; % Increase if xaxis lables are cut off
    subpad.top = .125; % Increase if subplot title is cut off
    subpad.left = .15; % Increase if ylabel is cut off
    subpad.right = .05; % Decrease if white-space on RHS
    subpad.legend = 0; % Increase if legend overlaps subplot titles
    fontsize = 8;
end

plotdim = [2,4];

figure;
set(gcf,'position',figbox);
iEE1 = 1;
iEE2 = 1;
for isubplot = 1:prod(plotdim)
    [col,row] = ind2sub(plotdim([2 1]),isubplot);
    left = (col-1+subpad.left)/plotdim(2);
    bottom = (1-(row-subpad.bot)/plotdim(1))/(1+subpad.legend);
    width = (1-(subpad.left+subpad.right))/plotdim(2);
    height = (1-subpad.bot-subpad.top)/(plotdim(1)*(1+subpad.legend));
    subplot(plotdim(1),plotdim(2),isubplot,'Position',[left bottom width height]); grid on; box on; hold on;
    %[N,X] = hist(EE(:,iEE),nbars);
    %bar(X,100*N/nEE)
    subplot(plotdim(1),plotdim(2),isubplot,'Position',[left bottom width height]); grid on; box on; hold on;
    if isubplot < 5
        N1 = hist(EE_ART(:,iEE1),bins{iEE1}); %ART
        bar(bins{iEE1},100*N1/length(EE_ART(:,iEE1)),'grouped')
        text(-7.5,16.5,['Mean: ', num2str(RMeanEEART(iEE1))],'Fontsize',8);
        text(-7.5,13.5,['Max: ', num2str(RMaxEEART(iEE1))],'Fontsize',8);
        iEE1 = iEE1 + 1;    
        ylabel('ART','interpreter','latex','fontsize',fontsize)
        ylim([0, 20])
    else
        N2 = hist(EE_Gust(:,iEE2),bins{iEE2}); %GustEtAl
        bar(bins{iEE2},100*N2/length(EE_Gust(:,iEE2)),'grouped')
        text(-7.5,11,['Mean: ', num2str(RMeanEEGust(iEE2))],'Fontsize',8);        
        text(-7.5,8.5,['Max: ', num2str(RMaxEEGust(iEE2))],'Fontsize',8);
        iEE2 = iEE2 + 1;
        ylabel('GustEtAl','interpreter','latex','fontsize',fontsize)
        ylim([0, 15])
    end
    %title('Euler Equation','interpreter','latex','fontsize',fontsize)
    xlabel('Errors ($\log_{10}$)','interpreter','latex','fontsize',fontsize)
    set(gca,'xlim',[-8,0],'xtick',-8:2:0);
    title(titles{isubplot})
    %lgd = legend('ART', 'GustEtAl');
    %lgd.FontSize = 6;
    %lgd.Location = 'northwest';
end

%disp('Mean Euler Equation Error')
%disp(['c: ', num2str(RMeanEE(1,:))])
%disp(['pigap: ', num2str(RMeanEE(2,:))])
%disp('Max Euler Equation Error')
%disp(['c: ', num2str(RMaxEE(1,:))])
%disp(['pigap: ', num2str(RMaxEE(2,:))])
% disp('Integral Euler Equation Error')
% disp(R.intEE)

%% Save figure
if strcmp(saving,'on')
    savefig([savename '.fig'])
end