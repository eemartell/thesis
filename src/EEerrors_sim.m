clear all
close all
clc
%simulates 10,000 periods and stores Euler equation errors
%creates histogram of simulation and Euler equation errors

% Save rule {'on','off'}
saving = 'on';

% Iteration
%   ti: time iteration
%   fp: fixed point
O.it = 'fp';

if strcmp(O.it,'fp')
    disp('Using fixed point solution')
    load('solutions/solutionfp1.mat')
elseif strcmp(O.it,'ti')
    disp('Using time iteration solution') 
    load('solutions/solutionti1.mat')
end     
% Numerical pdf of state variables
%   Simulation parameters
npers = 2000;
nobs = npers;
V = variables;
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
sims = simulation_test(pf,P,S,G,V,e,u,v);

% Shocks
GH.shockpts = 11;
[e_nodes,e_weight] = ghquad(GH.shockpts);
GH.e_nodes = P.sige*e_nodes;
e_weight = e_weight/sqrt(pi);
[u_nodes,u_weight] = ghquad(GH.shockpts);
GH.u_nodes = P.sigu*u_nodes;
u_weight = u_weight/sqrt(pi);
[v_nodes,v_weight] = ghquad(GH.shockpts);
GH.v_nodes = P.sigv*v_nodes;
v_weight = v_weight/sqrt(pi);

% Create array of weights for integration
e_weightArr3 = e_weight(:,ones(GH.shockpts,1),ones(GH.shockpts,1));
u_weightArr3 = permute(u_weight(:,ones(GH.shockpts,1),ones(GH.shockpts,1)), [2,1,3]);
v_weightArr3 = permute(v_weight(:,ones(GH.shockpts,1),ones(GH.shockpts,1)),[2,3,1]);
weightArr3 = e_weightArr3.*u_weightArr3.*v_weightArr3;

% Exogenous processes
gpArr3 = e_nodes(:,ones(GH.shockpts,1),ones(GH.shockpts,1));
apArr3 = permute(u_nodes(:,ones(GH.shockpts,1),ones(GH.shockpts,1)),[2,1,3]);

EE1 = zeros(npers,1);
EE2 = zeros(npers,1);
for time = 2:npers
        start = [sims(time,V.c),sims(time,V.pi)/P.pi]';
        state = [sims(time,V.g),sims(time,V.s),sims(time,V.mp),sims(time-1,V.in)];
        % Approximate solution
        EE_temp = eqm(start,state,O,P,S,G,pf,gpArr3,apArr3,weightArr3,GH);      
         % Store Euler Equation errors
        EE1(time) = abs(EE_temp(1));
        EE2(time) = abs(EE_temp(2));
end
% Find where ZLB binds
inp = sims(1:end-1,V.in).^P.rhoi.*(S.i*sims(2:end,V.pi).^P.phipi).^(1-P.rhoi).*exp(sims(2:end,V.mp));
R.ZLBlocs = find(inp <= 1);
R.notZLBlocs = find(inp > 1);
%   Percent nodes binding
R.perbind = 100*numel(R.ZLBlocs)/npers;

R.EE1 = log10(EE1(2:end));
R.EE2 = log10(EE2(2:end));
R.meanEE = [mean(R.EE1),mean(R.EE2)];
R.maxEE = [max(R.EE1),max(R.EE2)];

%% Save results
if strcmp(saving,'on')
    fname = ['eeerrors_sim' O.it];
    save(['solutions/' fname],'R');    
end