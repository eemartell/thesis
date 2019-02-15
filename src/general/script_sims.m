clear
clc

% Set random number seed
mtstream = RandStream('mt19937ar','seed',123456);
RandStream.setGlobalStream(mtstream);

% Load options
load('options')

% Load Fortran solution
T = 120;
nsims = 100;
ZLBpers = [0,.05,.1,.15,.20,.25];
nZLBpers = numel(ZLBpers);
sims_tmp = textread('../data-artificial/sims.txt','%f');
sims = reshape(sims_tmp,[T,V.nvar-V.nfore,nsims,nZLBpers]);

%% ZLB duration
clc
iZLBpers = 6;

disp('Growth Shock')
disp([mean(sum(sims(:,V.g,:,iZLBpers)<min(G.g_grid)))...
      mean(sum(sims(:,V.g,:,iZLBpers)>max(G.g_grid)))])

disp('Risk Premium Shock')
disp([mean(sum(sims(:,V.s,:,iZLBpers)<min(G.s_grid)))...
      mean(sum(sims(:,V.s,:,iZLBpers)>max(G.s_grid)))])

disp('Monetary Policy Shock')
disp([mean(sum(sims(:,V.mp,:,iZLBpers)<min(G.mp_grid)))...
      mean(sum(sims(:,V.mp,:,iZLBpers)>max(G.mp_grid)))])

disp('Notional Interest Rate')
temp = squeeze(sims(:,V.in,:,iZLBpers));
disp([mean(sum(temp<min(G.in_grid)))...
      mean(sum(temp>max(G.in_grid)))])
disp(100*(quantile(temp(:),[.025,.975])./S.in-1));

disp('Consumption')
temp = squeeze(sims(:,V.c,:,iZLBpers));
disp([mean(sum(temp<min(G.c_grid)))...
      mean(sum(temp>max(G.c_grid)))])
disp(100*(quantile(temp(:),[.025,.975])./S.c-1));

disp('Capital')
temp = squeeze(sims(:,V.k,:,iZLBpers));
disp([mean(sum(temp<min(G.k_grid)))...
      mean(sum(temp>max(G.k_grid)))])
disp(100*(quantile(temp(:),[.025,.975])./S.k-1));

disp('Investment')
temp = squeeze(sims(:,V.x,:,iZLBpers));
disp([mean(sum(temp<min(G.x_grid)))...
      mean(sum(temp>max(G.x_grid)))])
disp(100*(quantile(temp(:),[.025,.975])./S.x-1));
  
disp('Real wage')
temp = squeeze(sims(:,V.w,:,iZLBpers));
disp([mean(sum(temp<min(G.w_grid)))...
      mean(sum(temp>max(G.w_grid)))])
disp(100*(quantile(temp(:),[.025,.975])./S.w-1));