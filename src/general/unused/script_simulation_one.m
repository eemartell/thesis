% Counterfactual simulation with ZLB event of about 6 quarters
clear
clc

% Load options
load('options')
V.nstate = V.nvar-V.nfore;

% Load solution
load('solution_test.mat')

%% Do a short simulation
T = 22;
nsims = 1;

% Shocks
epsg = zeros(T,nsims);
epss = zeros(T,nsims);
epss(1:7,:) = 1;
epsi = zeros(T,nsims);
sims = simulation(pf,P,S,G,V,epsg,epss,epsi);
sims = sims(2:end,:);

%%
xaxis = 0:T-2;
figure(1)
subplot(3,1,1); box on; grid on;
plot(xaxis,sims(:,V.yg))
title(V.desc{V.yg})
axis('tight')
subplot(3,1,2); box on; grid on;
plot(xaxis,sims(:,V.pi))
title(V.desc{V.pi})
axis('tight')
subplot(3,1,3); box on; grid on;
plot(xaxis,sims(:,V.in))
title(V.desc{V.in})
axis('tight')