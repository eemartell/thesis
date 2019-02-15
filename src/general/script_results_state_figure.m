%  This file creates a mat file with summary statistics of the accuracy
%  parameter estimates of different specifications.
%  The output is tables/tables.mat

clear
clc
close all

%Baseline specification
% baseline = 'lin-kf-me0'; 
%  baseline = 'pw-if-me0';
 
%save the output matlab table
save_flag = 1;

% Load parameter names
load('options.mat','V','P','F')

% Load artificial data
npers = 120;
nsimsall = 100;
nsims = 50;
nZLBpers = 6;
 nqs = [0,6,12,18,24,30];
%nqs = [0,6,12,18,24];
simstemp = textread('../data-artificial/sims.txt','%f');
simscap = reshape(simstemp,[npers,V.nvar-V.nfore,nsimsall,nZLBpers]);


Vcap= V;
load('options_nocap.mat','V','P','F')
sims = nan(npers,V.nvar-V.nfore,nsims,nZLBpers);
for ivar = 1:(V.nvar-V.nfore)
    if strcmp(V.names{ivar},'cg')
        sims(:,ivar,:,:) = simscap(:,Vcap.yg,1:nsims,:);
    else
        eval(['sims(:,ivar,:,:) = simscap(:,Vcap.',V.names{ivar},',1:nsims,:);']);
    end
end

%Specifications to include: path, method, filter, measurement error
save_name = 'tables/tables.mat';
specs = {...
    '../Estimation-pwlinear/','pwlinear','if','me0';...
    '../Estimation-global/','global','pf','me5';         
    };

states = nan(npers,V.nvar-V.nfore,nsims,nZLBpers,size(specs,1));
errors = nan(npers,V.nvar-V.nfore,nsims,nZLBpers,size(specs,1));


%Loop through specifications
for ifolder = 1:size(specs,1)
    disp(['Specification: ' specs{ifolder,1}])
    % Loop through ZLB durations
    for inq = 1:numel(nqs)
        % ZLB duration
        nq = nqs(inq);
        disp(['ZLB duration = ' num2str(nq)])
        for isim = 1:nsims
            filename = [specs{ifolder,1},num2str(inq),'/states/',num2str(isim),'-states.txt'];
            if exist(filename, 'file') == 2
                statestemp = textread(filename,'%f');
                states(:,:,isim,nZLBpers,ifolder) = reshape(statestemp,[npers,V.nvar-V.nfore]);
            end
        end
    end
    errors(:,:,:,:,ifolder) = states(:,:,:,:,ifolder) - sims;
end   

plot(states(:,V.n,1,1,1))
%plot(sims(:,V.c,1,1))

%if save_flag 
%    save(save_name,'total')
%end