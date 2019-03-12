clear all
close all
clc
%simulates 10,000 periods and stores Euler equation errors
%creates histogram of simulation and Euler equation errors

% Save rule {'on','off'}
saving = 'on';

O.it = 'fp';
O.alg = 'ART';

if strcmp(O.it,'fp') && strcmp(O.alg, 'ART')
    load('solutions/solutionfpART_5.mat')
elseif strcmp(O.it,'fp') && strcmp(O.alg, 'Gust')
    load('solutions/solutionfpGust_5.mat')    
end

O.it = 'fp';
O.alg = 'ART';
V = variables;
% Numerical pdf of state variables
%   Simulation parameters
npers = 10000;
nobs = npers;
%V = variables_gustetal;
mtstream = RandStream('mt19937ar','seed',2);
RandStream.setGlobalStream(mtstream);
shockse = randn(npers,1);
shocksu = randn(npers,1);
shocksv = randn(npers,1);
%   Growth rate shocks
e = shockse;
%   Preference shocks
u = shocksu;
%   Interest rate shocks
v = shocksv;
%   Simulate for densities
if strcmp(O.alg, 'ART')
    sims = simulation(pf,P,S,G,V,e,u,v);
elseif strcmp(O.alg, 'Gust')
    V.n_zlb = V.nplotvar+1;    
    sims = simulation_gustetal(pf,P,S,G,V,e,u,v);
end

% Shocks
GH.shockpts = 11;
[e_nodes,e_weight] = ghquad(GH.shockpts);
GH.e_nodes = P.sigg*e_nodes;
e_weight = e_weight/sqrt(pi);
[u_nodes,u_weight] = ghquad(GH.shockpts);
GH.u_nodes = P.sigs*u_nodes;
u_weight = u_weight/sqrt(pi);
[v_nodes,v_weight] = ghquad(GH.shockpts);
GH.v_nodes = P.sigmp*v_nodes;
v_weight = v_weight/sqrt(pi);

% Create array of weights for integration
e_weightArr3 = e_weight(:,ones(GH.shockpts,1),ones(GH.shockpts,1));
u_weightArr3 = permute(u_weight(:,ones(GH.shockpts,1),ones(GH.shockpts,1)), [2,1,3]);
v_weightArr3 = permute(v_weight(:,ones(GH.shockpts,1),ones(GH.shockpts,1)),[2,3,1]);
weightArr3 = e_weightArr3.*u_weightArr3.*v_weightArr3;

% Exogenous processes
gpArr3 = e_nodes(:,ones(GH.shockpts,1),ones(GH.shockpts,1));
spArr3 = permute(u_nodes(:,ones(GH.shockpts,1),ones(GH.shockpts,1)),[2,1,3]);
mpArr3 = permute(v_nodes(:,ones(GH.shockpts,1),ones(GH.shockpts,1)),[2,3,1]); %???

EE1 = zeros(npers,1);
EE2 = zeros(npers,1);
EE3 = zeros(npers,1);
EE4 = zeros(npers,1);
if strcmp(O.alg,'Gust')
    EE3 = zeros(npers,1);
end
for time = 2:npers
    state = [sims(time,V.g),sims(time,V.s),sims(time,V.mp),sims(time,V.in),sims(time,V.c),sims(time,V.k),sims(time,V.x)];
        if strcmp(O.alg, 'ART')
            start = [sims(time,V.pi)/P.pi,sims(time,V.n),sims(time,V.q),sims(time,V.mc)]';
            % Approximate solution
            EE_temp = eqm(start,state,O,P,S,G,pf,gpArr3,weightArr3,GH);      
        elseif strcmp(O.alg, 'Gust')
            %%%cs and pigap here
            start = [sims(time,V.pi)/P.pi,sims(time,V.n),sims(time,V.n_zlb,:),sims(time,V.q),sims(time,V.mc)]';
             % Approximate solution           
             EE_temp = eqm_gustetal(start,state,O,P,S,G,pf,gpArr3,mpArr3,weightArr3,GH);    
        end
         % Store Euler Equation errors
        EE1(time) = abs(EE_temp(1));
        EE2(time) = abs(EE_temp(2));
        EE3(time) = abs(EE_temp(3));
        EE4(time) = abs(EE_temp(4));
        if strcmp(O.alg,'Gust')
            EE5(time) = abs(EE_temp(5));
        end
end
% Find where ZLB binds
%inp = sims(1:end-1,V.in).^P.rhoi.*(S.i*sims(2:end,V.pi).^P.phipi).^(1-P.rhoi).*exp(sims(2:end,V.mp));

%R.ZLBlocs = find(inp <= 1);
%R.notZLBlocs = find(inp > 1);
%   Percent nodes binding
%R.perbind = 100*numel(R.ZLBlocs)/npers;

R.EE1 = log10(EE1(2:end));
R.EE2 = log10(EE2(2:end));
R.EE3 = log10(EE3(2:end));
R.EE4 = log10(EE4(2:end));
if strcmp(O.alg,'Gust')
    R.EE5 = log10(EE5(2:end));
    R.meanEE = [mean(R.EE1),mean(R.EE2),mean(R.EE3),mean(R.EE4),mean(R.EE5)];
    R.maxEE = [max(R.EE1),max(R.EE2),max(R.EE3),max(R.EE4),max(R.EE5)];
elseif strcmp(O.alg,'ART')
    R.meanEE = [mean(R.EE1),mean(R.EE2),mean(R.EE3),mean(R.EE4)];
    R.maxEE = [max(R.EE1),max(R.EE2),max(R.EE3),max(R.EE4)];
end
%% Save results
if strcmp(saving,'on')
    fname = ['eeerrors_sim' O.it O.alg '_5'];
    save(['solutions/' fname],'R');    
end