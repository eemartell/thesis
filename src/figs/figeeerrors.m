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
O.alg = 'Gust';

% Figure name
figname = ['eeerrors' O.it O.alg];

% Load solution at posterior mean
sim = true;
if sim == true
    if strcmp(O.it,'fp') && strcmp(O.alg, 'ART')
        load('../solutions/eeerrors_simfpART')
    elseif strcmp(O.it,'fp') && strcmp(O.alg, 'Gust')
        load('../solutions/eeerrors_simfpGust')
    end
else
    load('../solutions/eeerrors')
end

EE(:,1) = R.EE1(:);
EE(:,2) = R.EE2(:);
if strcmp(O.alg, 'Gust')
    EE(:,3) = R.EE3(:);
end
nEE = size(EE,1);

EE_noZLB(:,1) = R.EE1(R.notZLBlocs);
EE_noZLB(:,2) = R.EE2(R.notZLBlocs);
nEE_noZLB = size(EE_noZLB,1);
EE_ZLB(:,1) = R.EE1(R.ZLBlocs);
EE_ZLB(:,2) = R.EE2(R.ZLBlocs);
if strcmp(O.alg, 'Gust')
    EE_noZLB(:,3) = R.EE3(R.notZLBlocs);
    EE_ZLB(:,3) = R.EE3(R.ZLBlocs);
end
nEE_ZLB = size(EE_ZLB,1);
nbars = 10;%5;
if strcmp(O.alg,'ART')
    titles = {'Consumption Euler','Firm Pricing'};
elseif strcmp(O.alg,'Gust')
     titles = {'Consumption Euler non-ZLB','Consumption Euler ZLB','Firm Pricing'};
end
% Plot options
if strcmp(plottype,'M')
    savename = ['Figs/' figname];
    figbox = [1,1,6.5,2.5];
    subpad.bot = .175; % Increase if xaxis lables are cut off
    subpad.top = .125; % Increase if subplot title is cut off
    subpad.left = .075; % Increase if ylabel is cut off
    subpad.right = .025; % Decrease if white-space on RHS 
    subpad.legend = 0; % Increase if legend overlaps subplot titles
    fontsize = 10; 
elseif strcmp(plottype,'P')
    savename = ['Figs/' figname '_pres'];
    set(0,'DefaultAxesFontSize',8)
    figbox = [1,1,4.5,2.5];
    subpad.bot = .175; % Increase if xaxis lables are cut off
    subpad.top = .125; % Increase if subplot title is cut off
    subpad.left = .15; % Increase if ylabel is cut off
    subpad.right = .05; % Decrease if white-space on RHS
    subpad.legend = 0; % Increase if legend overlaps subplot titles
    fontsize = 8;
end
if strcmp(O.alg,'ART')
    num = 2;
elseif strcmp(O.alg,'Gust')
    num = 3;
end
plotdim = [1,num];

figure;
set(gcf,'position',figbox);
for iEE = 1:num
    [col,row] = ind2sub(plotdim([2 1]),iEE);
    left = (col-1+subpad.left)/plotdim(2);
    bottom = (1-(row-subpad.bot)/plotdim(1))/(1+subpad.legend);
    width = (1-(subpad.left+subpad.right))/plotdim(2);
    height = (1-subpad.bot-subpad.top)/(plotdim(1)*(1+subpad.legend));
    subplot(plotdim(1),plotdim(2),iEE,'Position',[left bottom width height]); grid on; box on; hold on;
    %[N,X] = hist(EE(:,iEE),nbars);
    %bar(X,100*N/nEE)
    [N1,X1] = hist(EE_noZLB(:,iEE),nbars); %no ZLB
    bar(X1,100*N1/nEE_noZLB,'grouped') 
    [N2,X2] = hist(EE_ZLB(:,iEE),nbars); %ZLB
    bar(X2,100*N2/nEE_ZLB,'grouped')
    title('Euler Equation','interpreter','latex','fontsize',fontsize)
    ylabel('Frequency (\%)','interpreter','latex','fontsize',fontsize)
    xlabel('Errors ($\log_{10}$)','interpreter','latex','fontsize',fontsize)
    set(gca,'xlim',[-8,0],'xtick',-8:2:0);
    title(titles{iEE})
    lgd = legend('ZLB not binding', 'ZLB binding');
    lgd.FontSize = 6;
    lgd.Location = 'northwest';
end

disp('Mean Euler Equation Error')
disp(R.meanEE)
disp('Max Euler Equation Error')
disp(R.maxEE)
% disp('Integral Euler Equation Error')
% disp(R.intEE)

%% Save figure
if strcmp(saving,'on')
    print(gcf,'-depsc2','-painters',[savename '.eps'])
    saveas(gcf,[savename '.fig'])
end