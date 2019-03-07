function G = grids(O,P)

% G = grids(O,P)
%   Constructs an evenly-spaced discretized state space
% Inputs:
%	O :   Structure of user-specified options
%   P :   Structure of parameters
% Output:
%   G :   Structure of grid points

%--------------------------------------------------------------------------
% Nodes and weights for numerical integration
%--------------------------------------------------------------------------
% Rouwenhorst for growth shock
[G.e_nodes,G.e_weight] = rouwenhorst(P.g,P.rhog,P.sige,O.g_pts);
% Rouwenhorst for discount factor shock
[G.u_nodes,G.u_weight] = rouwenhorst(P.beta,P.rhobeta,P.sigu,O.beta_pts);
% Rouwenhorst for monetary policy shock
[G.v_nodes,G.v_weight] = rouwenhorst(0,0,P.sigv,O.mp_pts);
%--------------------------------------------------------------------------
% Set up grids by discretizing the state space
%--------------------------------------------------------------------------
% Define grid points for consumption
G.c_grid = linspace(O.cbound(1),O.cbound(2),O.c_pts);
% Define grid points for capital
G.k_grid = linspace(O.kbound(1),O.kbound(2),O.k_pts);
% Define grid points for investment
G.i_grid = linspace(O.ibound(1),O.ibound(2),O.i_pts);
% Define grid points for notional interest rate
G.rn_grid = linspace(O.rnbound(1),O.rnbound(2),O.rn_pts);
% Define grid points for growth
G.g_grid = G.e_nodes';
% Define grid points for the discount factor
G.beta_grid = G.u_nodes';
% Define grid points for the monetary policy shock
G.mp_grid = G.v_nodes';
%--------------------------------------------------------------------------
% Construct grid arrays
%--------------------------------------------------------------------------
% Full nonlinear model
[G.c_gr,G.k_gr,G.i_gr,G.rn_gr,G.g_gr,G.beta_gr,G.mp_gr] = ...
    ndgrid(G.c_grid,G.k_grid,G.i_grid,G.rn_grid,...
           G.g_grid,G.beta_grid,G.mp_grid);
G.nodes = numel(G.c_gr);
G.griddim = size(G.c_gr);