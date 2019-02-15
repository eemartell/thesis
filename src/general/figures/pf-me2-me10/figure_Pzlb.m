%  Creates a figure of the Expectated duration, actual and estimated
%  scatter pot

clear
clc
close all

savename = 'Pzlb';

%Specifications to include: path, method, filter, measurement error,%label
specs = {...
    '../../../Estimation-global/me2/','global','pf','me2','Global, Particle Filter, ME 2$\%$';,...  
    '../../../Estimation-global/me5/','global','pf','me5','Global, Particle Filter, ME 5$\%$';,...  
    '../../../Estimation-global/me10/','global','pf','me10','Global, Particle Filter, ME 10$\%$';,...  
    '../../../Estimation-pwlinear/','pwlinear','if','me0','Piecewise Linear, Inversion Filter, ME 0$\%$';...
    '../../../Estimation-linear/me0/','linear','kf','me0','Level Linear, Kalman Filter, ME 0$\%$';...
    };

npers = 120;
nsimsall = 100;
nsims = 50;
nZLBpers = 6;
nqs = [0,6,12,18,24,30];
nspecs = size(specs,1);
ndraws = 1000;

% Load actual expected duration from DGP
Pzlb = textread('../../../data-artificial/Ezlb_duration/Pzlb.txt','%f');
Pzlb = reshape(Pzlb,[npers,nsimsall,nZLBpers]);

% Load parameter names
load('../../options.mat','V','P','F')

% Load actual data
simscap = textread('../../../data-artificial/sims.txt','%f');
simscap = reshape(simscap,[npers,V.nvar-V.nfore,nsimsall,nZLBpers]);

Vcap= V;

%Load in V from misspecified model
load('../../../Results-nocap/options.mat','V','P','F')
Vnstate = V.nvar-V.nfore

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

Pzlb_est = nan(npers,nsims,nZLBpers,nspecs);
%Loop through specifications
for ifolder = 1:size(specs,1)
    disp(['Specification: ' specs{ifolder,1}])
    %Load state estimate for given specification
    filename = [specs{ifolder,1},'/states/Pzlb.txt'];
    if exist(filename, 'file') == 2
        Pzlb_temp = textread(filename,'%f');
        Pzlb_temp = reshape(Pzlb_temp,[npers,nsims,nZLBpers]);
        %Global data may be missing, replace zeros with nan
        %if strcmp(specs{ifolder,2},'global')
        %   Eduration_temp(Eduration_temp==0) = NaN;
        %end
        %Error in global and linear expectation calculation means Eduration
        %needs to be reduced by 1
        %if (strcmp(specs{ifolder,2},'global') || strcmp(specs{ifolder,2},'linear'))
        %    Eduration_temp = Eduration_temp - 1;
        %end
        Pzlb_temp(Pzlb_temp>1 | Pzlb_temp<0) = NaN;
        Pzlb_est(:,:,:,ifolder) = Pzlb_temp;
    end
end   

Pzlb_est(:,50,5,2) = NaN;

in = squeeze(sims(:,V.in,:,:));
%Replace zeros with NaNs
Pzlb = Pzlb(:,1:nsims,:);
Pzlb(in<=1) = NaN;
%Pzlb(in>(1.25/100+1)^0.25) = NaN;

Pzlb_error = Pzlb_est(:,:,:,2) - Pzlb;
i_error = squeeze(states(:,V.i,:,:,2) - sims(:,V.i,:,:));


%test = ((Eduration_pw == 0) & (zlblasts > 0));
%Eduration_test = Eduration(test);
%in = squeeze(sims(:,V.in,:,:));
%in= in(test);


plotops
%% ------------------------------------------------------------------------
% Plot the results
%--------------------------------------------------------------------------
% close all
% Spec strings for legend
titles = {'NL-PF-2\%','NL-PF-5\%','NL-PF-10\%','PW-IF-0\%','Lin-KF-0\%'};
order = [1,2,3,4,5];

%   Subplot padding
plotposition = [1 1 6.5 2.25];
pad.topmarg = .12; % Increase if subplot title is cut off
pad.leftmarg = .07; % Increase if ylabel is cut off
pad.vspace = 0.04;
pad.hspace = .08;
pad.rightmarg = .025; % Decrease if white-space on RHS
pad.botmarg = .15; % Increase if xaxis lables are cut off
legfont = 10;
plotdim = [1,4];

% Compute mean/credible set for each actual duration
Pzlbbins = 1:15;
for iPzlb = Pzlbbins
    idxs = Pzlb(:) >= (iPzlb-1)/20 & Pzlb(:) < iPzlb/20;
    for imethod = 1:4
        temp = Pzlb_est(:,:,:,imethod);
        temp = temp(idxs(:));
        qPzlb_est(iPzlb,:,imethod) = [mean(temp,'omitnan'),quantile(temp,[.05,.95])];
    end
end

figure(3);
set(gcf,'Position',plotposition);

for isubplot = 1:4
    imethod = order(isubplot);
    subplot(plotdim(1),plotdim(2),isubplot); hold on; box on; grid on;
    
    % Shade       
    x=(Pzlbbins-0.5)/20;
    y1 = qPzlb_est(:,1,imethod)';
    y2 = qPzlb_est(:,2,imethod)';
    y3 = qPzlb_est(:,3,imethod)';
    X=[x,fliplr(x)];               
    Y=[y3,fliplr(y2)];             
    fill(X,Y,[0.8 0.8 0.8],'linestyle','none');
    %fill(X,Y,'k');
    set(gca,'Layer','top')
    
    % Plot median
    plot(x,y1,'k-');
    plot(0:0.05:1,0:0.05:1,'k--','MarkerSize',2)
    
    % Labels
    title(titles{imethod},'Interpreter','latex','fontsize',10);
    xlabel('Actual')
    ylabel('Estimated')
    set(gca,'xlim',[0,0.5],'ylim',[0,0.5])
end

% Set subplot dimensions
for j = 1:4
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

