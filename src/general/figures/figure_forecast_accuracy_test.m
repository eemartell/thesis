%  Creates a figure of the RMSE forecasts from the different model
%  specifications from the periods just prior to the ZLB binding.

clear
clc
close all
savename = 'forecast_accuracy'; 
%Specifications to include: path, method, filter, measurement error,%label
specs = {...
    %'../../Estimation-global/me2/','global','pf','me2','Global, Particle Filter, ME 2$\%$';  ...
    '../../data-artificial/Ezlb_duration/','dgp','na','na','DGP';  ...
    '../../Estimation-global/me5/states/','global','pf','me5','Global, Particle Filter, ME 5$\%$';  ...
    '../../Estimation-pwlinear/states/','pwlinear','if','me0','Piecewise Linear, Inversion Filter, ME 0$\%$';...
    %'../../Estimation-linear/me0/','linear','kf','me0','Lin-KF-$0\%$';...
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

temp = textread([specs{1,1},'forecast.txt'],'%f');
% forecast_dgp = reshape(temp,[npers,3,V.nvar-V.nfore,nsimsall,nZLBpers]);

Vcap= V;

%Load in V from misspecified model
load('../../Results-nocap/options.mat','V','P','F')
Vnstate = V.nvar-V.nfore;

%Create artificial data matrix (true model) that aligns with variable order
%from the true model. Note that yg is named cg in misspec model
sims = nan(npers,Vnstate,nsims,nZLBpers);
forecast = nan(npers,3,Vnstate,nsims,nZLBpers,nspecs);
for ivar = 1:(Vnstate)
    if strcmp(V.names{ivar},'cg')
        sims(:,ivar,:,:) = simscap(:,Vcap.yg,1:nsims,:);
%         forecast(:,:,ivar,:,:,1) = forecast_dgp(:,:,Vcap.yg,1:nsims,:);
    else
        eval(['sims(:,ivar,:,:) = simscap(:,Vcap.',V.names{ivar},',1:nsims,:);']);
%         eval(['forecast(:,:,ivar,:,:,1) = forecast_dgp(:,:,Vcap.',V.names{ivar},',1:nsims,:);']);
    end
end

%Load in forecast cred sets from the specifications
temp = textread([specs{1,1},'p05.txt'],'%f');
p05_dgp = reshape(temp,[npers,3,4,nsimsall,nZLBpers]);
temp = textread([specs{1,1},'p95.txt'],'%f');
p95_dgp = reshape(temp,[npers,3,4,nsimsall,nZLBpers]);

p05 = nan(npers,3,Vnstate,nsims,nZLBpers,nspecs);
p05(:,:,V.yg,:,:,1) = p05_dgp(:,:,1,1:nsims,:);
p05(:,:,V.pi,:,:,1) = p05_dgp(:,:,2,1:nsims,:);
p05(:,:,V.i,:,:,1) = p05_dgp(:,:,3,1:nsims,:);
p05(:,:,V.in,:,:,1) = p05_dgp(:,:,4,1:nsims,:);

p95 = nan(npers,3,Vnstate,nsims,nZLBpers,nspecs);
p95(:,:,V.yg,:,:,1) = p95_dgp(:,:,1,1:nsims,:);
p95(:,:,V.pi,:,:,1) = p95_dgp(:,:,2,1:nsims,:);
p95(:,:,V.i,:,:,1) = p95_dgp(:,:,3,1:nsims,:);
p95(:,:,V.in,:,:,1) = p95_dgp(:,:,4,1:nsims,:);

%Load in forecasts from the specifications
%Loop through specifications
for ifolder = 2:size(specs,1)
    disp(['Specification: ' specs{ifolder,1}])
    %Load state estimate for given specification
    filename = [specs{ifolder,1},'p05.txt'];
    if exist(filename, 'file') == 2
        temp = textread(filename,'%f');
        temp = reshape(temp,[npers,3,Vnstate,nsims,nZLBpers]);
        % For global specs which may be missing, replace zeros with nans.
        temp(temp==0) = NaN;
        p05(:,:,:,:,:,ifolder) = temp;
        
        filename = [specs{ifolder,1},'p95.txt'];
        temp = textread(filename,'%f');
        temp = reshape(temp,[npers,3,Vnstate,nsims,nZLBpers]);
        % For global specs which may be missing, replace zeros with nans.
        temp(temp==0) = NaN;
        p95(:,:,:,:,:,ifolder) = temp;
        
        filename = [specs{ifolder,1},'forecast.txt'];
        temp = textread(filename,'%f');
        temp = reshape(temp,[npers,3,Vnstate,nsims,nZLBpers]);
        % For global specs which may be missing, replace zeros with nans.
        temp(temp==0) = NaN;
        forecast(:,:,:,:,:,ifolder) = temp;
    end
end   


%%
clf

isim = 12;
iZLBper = 6;
sim_true = sims(:,:,isim,iZLBper);
tZLBstart = find(sim_true(:,V.in) <= 1,1,'first') - 1;

% forecast = nan(npers,3,Vnstate,nsims,nZLBpers,nspecs);
sim_nlpf = forecast(:,3,V.in,isim,iZLBper,2);

figure(1)
plotvars = [V.yg,V.pi,V.in];
for isubplot = 1:numel(plotvars)
    subplot(3,1,isubplot); hold on; grid on; box on;
    plot(sim_true(tZLBstart:tZLBstart+12,plotvars(isubplot)))
%     plot(9,forecast(tZLBstart,3,V.in,isim,iZLBper,1),'d');
    plot([2,5,9],forecast(tZLBstart,:,plotvars(isubplot),isim,iZLBper,2),'x');
    plot([2,5,9],forecast(tZLBstart,:,plotvars(isubplot),isim,iZLBper,3),'o');
    
    plot([2,5,9],p05(tZLBstart,:,plotvars(isubplot),isim,iZLBper,2),'x');
    plot([2,5,9],p05(tZLBstart,:,plotvars(isubplot),isim,iZLBper,3),'o');
    
    plot([2,5,9],p95(tZLBstart,:,plotvars(isubplot),isim,iZLBper,2),'x');
    plot([2,5,9],p95(tZLBstart,:,plotvars(isubplot),isim,iZLBper,3),'o');
end
