clear
clc

% from Plante et al\Data\DSGE\Results\script_solution_postmean.m
load('solutions\solution_test.mat')

% Grids
evalpts = 21;
eval_g_grid = linspace(G.g_grid(1),G.g_grid(end),evalpts);
eval_a_grid = linspace(G.a_grid(1),G.a_grid(end),evalpts);
eval_mp_grid = linspace(G.mp_grid(1),G.mp_grid(end),evalpts);
[Geval.g_gr,Geval.a_gr,Geval.mp_gr] = ndgrid(...
    eval_g_grid,eval_a_grid,eval_mp_grid);
Geval.nodes = numel(Geval.g_gr);
Geval.griddim = size(Geval.g_gr);
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

pfeval.c = zeros(Geval.griddim);
pfeval.pigap = zeros(Geval.griddim);
[pfeval.c(:),pfeval.pigap(:)] = Fallterp32c_F(...
                G.g_grid,G.a_grid,G.mp_grid,...
                Geval.g_gr(:),Geval.a_gr(:),Geval.mp_gr(:),...
                pf.c,pf.pigap);

EE1 = zeros(Geval.nodes,1);
EE2 = zeros(Geval.nodes,1);
for inode = 1:Geval.nodes
        start = [pfeval.c(inode),pfeval.pigap(inode)]';
        state = [Geval.g_gr(inode),Geval.a_gr(inode),Geval.mp_gr(inode)];
        % Approximate solution
        EE_temp = eqm(start,state,O,P,S,G,pf,gpArr3,apArr3,weightArr3,GH);
        % Store Euler Equation errors
        EE1(inode) = abs(EE_temp(1));
        EE2(inode) = abs(EE_temp(2));
end
R.EE1 = log10(EE1);
R.EE2 = log10(EE2);
R.meanEE = [mean(R.EE1),mean(R.EE2)];
R.maxEE = [max(R.EE1),max(R.EE2)];
save('solutions\eeerrors.mat','R')

