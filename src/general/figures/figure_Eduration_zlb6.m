%  Creates a figure of the Expectated duration, actual and estimated
%  scatter pot

clear
clc
close all

savename = 'Eduration_zlb6';

%Specifications to include: path, method, filter, measurement error,%label
specs = {...
    '../../Estimation-pwlinear/','pwlinear','if','me0','Piecewise Linear, Inversion Filter, ME 0$\%$';...
    '../../Estimation-global/me5/','global','pf','me5','Global, Particle Filter, ME 5$\%$';  
    '../../Estimation-global/me2/','global','pf','me2','Global, Particle Filter, ME 2$\%$';  
    };

npers = 120;
nsimsall = 100;
nsims = 50;
nZLBpers = 6;
nqs = [0,6,12,18,24,30];
nspecs = size(specs,1);
ndraws = 1000;

% Load actual expected duration from DGP
Eduration = textread('../../data-artificial/Ezlb_duration/Eduration.txt','%f');
Eduration = reshape(Eduration,[npers,nsimsall,nZLBpers]);

%Load in V from misspecified model
load('../../Results-nocap/options.mat','V','P','F')

Eduration_est = nan(npers,nsims,nZLBpers,nspecs);
%Loop through specifications
for ifolder = 1:size(specs,1)
    disp(['Specification: ' specs{ifolder,1}])
    %Load state estimate for given specification
    filename = [specs{ifolder,1},'/states/Eduration.txt'];
    if exist(filename, 'file') == 2
        Eduration_temp = textread(filename,'%f');
        Eduration_temp = reshape(Eduration_temp,[npers,nsims,nZLBpers]);
        %Global data may be missing, replace zeros with nan
        %if strcmp(specs{ifolder,2},'global')
        %   Eduration_temp(Eduration_temp==0) = NaN;
        %end
        %Error in global and linear expectation calculation means Eduration
        %needs to be reduced by 1
        %if (strcmp(specs{ifolder,2},'global') || strcmp(specs{ifolder,2},'linear'))
        %    Eduration_temp = Eduration_temp - 1;
        %end
        Eduration_est(:,:,:,ifolder) = Eduration_temp;
    end
end   

filename = ['../../Estimation-pwlinear/states/zlblasts.txt'];
zlblasts = textread(filename,'%f');
zlblasts = reshape(zlblasts,[npers,nsims,nZLBpers]);

Eduration_est(:,50,5,2) = NaN;

%Replace zeros with NaNs
Eduration = Eduration(:,1:nsims,:);
zlblasts(Eduration==0) = NaN;
Eduration(Eduration==0) = NaN;
Eduration(:,:,1:5) = NaN;

Eduration_pw = Eduration_est(:,:,:,1);
%test = ((Eduration_pw == 0) & (zlblasts > 0));
%Eduration_test = Eduration(test);
%in = squeeze(sims(:,V.in,:,:));
%in= in(test);

%Plot accuracy for the three categories of data
plotops

% %% ------------------------------------------------------------------------
% % Plot the results
% %--------------------------------------------------------------------------
% % close all
% % Spec strings for legend
% legstrs = {'Global-PF-5\%','PwLin-IF-0\%','Linear-KF-0\%'};
% 
% %   Subplot padding
% plotposition = [1 1 6.5 5.5];
% pad.topmarg = .12; % Increase if subplot title is cut off
% pad.leftmarg = .08; % Increase if ylabel is cut off
% pad.vspace = 0.06;
% pad.hspace = .10;
% pad.rightmarg = .025; % Decrease if white-space on RHS
% pad.botmarg = .05; % Increase if xaxis lables are cut off
% legfont = 10;
% 
% figure(1); hold on;
% set(gcf,'Position',plotposition);
% 
% p1(1) = scatter(Eduration(:),reshape(Eduration_est(:,:,:,2),npers*nsims*nZLBpers,1),'MarkerEdgeColor','r');
% p1(2) = scatter(Eduration(:),reshape(Eduration_est(:,:,:,1),npers*nsims*nZLBpers,1),'MarkerEdgeColor','b');
% p1(3) = scatter(Eduration(:),reshape(Eduration_est(:,:,:,3),npers*nsims*nZLBpers,1),'MarkerEdgeColor','g');
% plot([0,15],[0,15],'k')
% 
% title('Expected duration of ZLB (number of quarters)','Interpreter','latex','fontsize',10);
% xlabel('Actual')
% ylabel('Estimated')
% 
% % Add Data Legend
% leg1 = legend(p1,legstrs,'Orientation','horizontal','Interpreter','latex','fontsize',10);
% set(leg1,'Position',[.25 1-pad.topmarg 0.5 0]);
% 
% saveas(gcf,savename)

% %% ------------------------------------------------------------------------
% % Plot the results
% %--------------------------------------------------------------------------
% % close all
% % Spec strings for legend
% legstrs = {'Global-PF-5\%','PwLin-IF-0\%','Linear-KF-0\%'};
% 
% %   Subplot padding
% plotposition = [1 1 6.5 5.5];
% pad.topmarg = .12; % Increase if subplot title is cut off
% pad.leftmarg = .08; % Increase if ylabel is cut off
% pad.vspace = 0.06;
% pad.hspace = .10;
% pad.rightmarg = .025; % Decrease if white-space on RHS
% pad.botmarg = .05; % Increase if xaxis lables are cut off
% legfont = 10;
% 
% figure(2); hold on;
% set(gcf,'Position',plotposition);
% 
% % Compute mean/credible set for each actual duration
% ZLBdurs = 1:20;
% for iZLBdur = ZLBdurs
%     idxs = Eduration(:) > iZLBdur - .5 & Eduration(:) <= iZLBdur + .5;
%     for imethod = 1:3
%         temp = Eduration_est(:,:,:,imethod);
% %         meanEduration_est(iZLBdur,:,imethod) = mean(temp(idxs),'omitnan');
%         temp = temp(idxs(:));
%         temp(isnan(temp)) = [];
%         qEduration_est(iZLBdur,:,imethod) = quantile(temp,[.5,.05,.95]);
%     end
% end
% 
% % p1(1) = scatter(ZLBdurs,meanEduration_est(:,2),'MarkerEdgeColor','r');
% % p1(2) = scatter(ZLBdurs,meanEduration_est(:,1),'MarkerEdgeColor','b');
% % p1(3) = scatter(ZLBdurs,meanEduration_est(:,3),'MarkerEdgeColor','g');
% p1(1) = plot(ZLBdurs,qEduration_est(:,1,2),'rd-');
% p1(2) = plot(ZLBdurs,qEduration_est(:,1,1),'bd-');
% p1(3) = plot(ZLBdurs,qEduration_est(:,1,3),'gd-');
% plot(ZLBdurs,qEduration_est(:,2,2),'rx--');
% plot(ZLBdurs,qEduration_est(:,2,1),'bx--');
% plot(ZLBdurs,qEduration_est(:,2,3),'gx--');
% plot(ZLBdurs,qEduration_est(:,3,2),'rx--');
% plot(ZLBdurs,qEduration_est(:,3,1),'bx--');
% plot(ZLBdurs,qEduration_est(:,3,3),'gx--');
% plot([0,15],[0,15],'k')
% 
% title('Expected duration of ZLB (number of quarters)','Interpreter','latex','fontsize',10);
% xlabel('Actual')
% ylabel('Estimated')
% 
% % Add Data Legend
% leg1 = legend(p1,legstrs,'Orientation','horizontal','Interpreter','latex','fontsize',10);
% set(leg1,'Position',[.25 1-pad.topmarg 0.5 0]);
% 
% saveas(gcf,savename)

%% ------------------------------------------------------------------------
% Plot the results
%--------------------------------------------------------------------------
% close all
% Spec strings for legend
titles = {'PW-IF-0\%','NL-PF-5\%','NL-PF-2\%'};
order = [2,1,3];

%   Subplot padding
plotposition = [1 1 6.5 2.25];
pad.topmarg = .12; % Increase if subplot title is cut off
pad.leftmarg = .06; % Increase if ylabel is cut off
pad.vspace = 0.04;
pad.hspace = .08;
pad.rightmarg = .025; % Decrease if white-space on RHS
pad.botmarg = .15; % Increase if xaxis lables are cut off
legfont = 10;
plotdim = [1,3];

% Compute mean/credible set for each actual duration
ZLBdurs = 1:20;
for iZLBdur = ZLBdurs
    idxs = Eduration(:) > iZLBdur - .5 & Eduration(:) <= iZLBdur + .5;
    for imethod = 1:3
        temp = Eduration_est(:,:,:,imethod);
        temp = temp(idxs(:));
        qEduration_est(iZLBdur,:,imethod) = quantile(temp,[.5,.05,.95]);
    end
end

figure(2);
set(gcf,'Position',plotposition);

for isubplot = 1:3
    imethod = order(isubplot);
    subplot(plotdim(1),plotdim(2),isubplot); hold on; box on; grid on;
    
    % Shade       
    x=ZLBdurs(1:12);
    y1 = qEduration_est(1:12,1,imethod)';
    y2 = qEduration_est(1:12,2,imethod)';
    y3 = qEduration_est(1:12,3,imethod)';
    X=[x,fliplr(x)];               
    Y=[y3,fliplr(y2)];             
    fill(X,Y,[0.8 0.8 0.8],'linestyle','none');
    set(gca,'Layer','top')
    
    % Plot median
    plot(x,y1,'k-');
    plot(0:1:15,0:1:15,'k--','MarkerSize',2)
    
    % Labels
    title(titles{imethod},'Interpreter','latex','fontsize',10);
    xlabel('Actual')
    ylabel('Estimated')
    set(gca,'xlim',[1,12],'ylim',[1,12])
end

% Set subplot dimensions
for j = 1:3
    [col,row] = ind2sub(plotdim([2 1]),j);
    width = (1-pad.leftmarg-pad.hspace*(plotdim(2)-1)-pad.rightmarg)/plotdim(2);
    height = (1-pad.topmarg-pad.vspace*(plotdim(1)-1)-pad.botmarg)/plotdim(1);
    left = pad.leftmarg + (col-1)*(width+pad.hspace);
    bottom = 1-pad.topmarg-row*height-(row-1)*pad.vspace;
    subplot(plotdim(1),plotdim(2),j,'Position',[left bottom width height]);
end

saveas(gcf,[savename '.pdf'])

print(gcf,'-depsc2','-painters',[savename '.eps'])
saveas(gcf,[savename '.fig'])
