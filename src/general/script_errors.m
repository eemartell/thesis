% Counterfactual simulation with ZLB event of about 6 quarters
clear
clc
rng(0)

% Load options
load('options')
V.nstate = V.nvar-V.nfore;

% Load solution
load('solution_test.mat')

%% Ergodic distribution

nergdraws = 2500;
nergpers = 500;
epsg = randn(nergpers,nergdraws);
epss = randn(nergpers,nergdraws);
epsmp = randn(nergpers,nergdraws);
ergsims = simulation(pf,P,S,G,V,epsg,epss,epsmp);
ergdraws = squeeze(ergsims(end,:,:));

%% Euler equation errors

nreps = 20000;
bondexp = zeros(nreps,1);
pcexp = zeros(nreps,1);
bonderr = zeros(nergdraws,1);
pcerr = zeros(nergdraws,1);
for idraw = 1:nergdraws
    if mod(idraw,100) == 1
        disp(idraw)
    end
    % Time-t variables in expectation operator
    g = ergdraws(V.g,idraw);
    lam = ergdraws(V.lam,idraw);
    pi = ergdraws(V.pi,idraw);
    s = ergdraws(V.s,idraw);
    i = ergdraws(V.i,idraw);
    y = ergdraws(V.y,idraw);
    w = ergdraws(V.w,idraw);    
    % Shocks
    epsg = randn(2,nreps);
    epss = randn(2,nreps);
    epsmp = randn(2,nreps);
    % Time-t+1 variables in expectation operator
    varp = simulation(pf,P,S,G,V,epsg,epss,epsmp,ergdraws(:,idraw));
    gp = squeeze(varp(2,V.g,:));
    lamp = squeeze(varp(2,V.lam,:));
    pip = squeeze(varp(2,V.pi,:));
    yp = squeeze(varp(2,V.y,:));
    bondexp = (lam./lamp).*(s*i./(gp.*pip));
    pcexp = (lam./lamp).*(pip/P.pi-1).*(pip/P.pi).*(yp/y);    
    % Euler equation residuals
    bonderr(idraw) = 1-P.beta*mean(bondexp);
    pcerr(idraw) = P.varphip*(pi/P.pi-1)*pi/P.pi-(1-P.thetap)-P.thetap*w-P.beta*P.varphip*mean(pcexp);
end
R.ergdraws = ergdraws;
R.bonderr = log10(abs(bonderr));
R.pcerr = log10(abs(pcerr));
R.meanbonderr = mean(log10(abs(bonderr)));
R.meanpcerr = mean(log10(abs(pcerr)));
R.maxbonderr = max(log10(abs(bonderr)));
R.maxpcerr = max(log10(abs(pcerr)));
save('simulation_test.mat','R','V');
