    function G = grids(O,P,S)

% G = grids(O,P)
%   Constructs an evenly-spaced discretized state space
% Inputs:
%	O :   Structure of user-specified options
%   P :   Structure of parameters
%   S :   Structure of steady state values
% Output:
%   G :   Structure of grid points

%--------------------------------------------------------------------------
% Set endogenous state space bounds
%--------------------------------------------------------------------------
G.inmin = (1-O.in_boundper)*S.i;
G.inmax = (1+O.in_boundper)*S.i;
G.cmin = (1-O.c_boundper)*S.c;
G.cmax = (1+O.c_boundper)*S.c;
G.kmin = (1-O.k_boundper)*S.k;
G.kmax = (1+O.k_boundper)*S.k;
G.xmin = (1-O.x_boundper)*S.x;
G.xmax = (1+O.x_boundper)*S.x;
G.wmin = (1-O.w_boundper)*S.w;
G.wmax = (1+O.w_boundper)*S.w;
%--------------------------------------------------------------------------
% Nodes and weights for numerical integration
%--------------------------------------------------------------------------
% Rouwenhorst for growth shock
[G.epsg_nodes,G.epsg_weight] = rouwenhorst(P.g,0,P.sigg,O.g_pts);
% Rouwenhorst for risk premium shock
[G.epss_nodes,G.epss_weight] = rouwenhorst(P.s,P.rhos,P.sigs,O.s_pts);
% Rouwenhorst for monetary policy shock
[G.epsmp_nodes,G.epsmp_weight] = rouwenhorst(0,0,P.sigmp,O.mp_pts);
%--------------------------------------------------------------------------
% Set up grids by discretizing the state space
%--------------------------------------------------------------------------
% Define grid points for growth
G.g_grid = G.epsg_nodes';
% Define grid points for the preference shock
G.s_grid = G.epss_nodes';
% Define grid points for the monetary policy shock
G.mp_grid = G.epsmp_nodes';
% Define grid points for notional interest rate
G.in_grid = linspace(G.inmin,G.inmax,O.in_pts);
% Define grid points for consumption
G.c_grid = linspace(G.cmin,G.cmax,O.c_pts);
% Define grid points for capital
G.k_grid = linspace(G.kmin,G.kmax,O.k_pts);
% Define grid points for investment
G.x_grid = linspace(G.xmin,G.xmax,O.x_pts);
% Define grid points for real wage
G.w_grid = linspace(G.wmin,G.wmax,O.w_pts);
%--------------------------------------------------------------------------
% Construct grid arrays
%--------------------------------------------------------------------------
% Full nonlinear model
[G.g_gr,G.s_gr,G.mp_gr,G.in_gr,G.c_gr,G.k_gr,G.x_gr,G.w_gr] = ndgrid(...
    G.g_grid,G.s_grid,G.mp_grid,...
	G.in_grid,G.c_grid,G.k_grid,G.x_grid,G.w_grid);
G.nodes = numel(G.g_gr);
G.griddim = size(G.g_gr);