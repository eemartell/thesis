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
%--------------------------------------------------------------------------
% Nodes and weights for numerical integration
%--------------------------------------------------------------------------
% Rouwenhorst for productivity growth process
[G.e_nodes,G.e_weight] = rouwenhorst(P.g,0,P.sigg,O.g_pts);
% Rouwenhorst for risk premium process
[G.u_nodes,G.u_weight] = rouwenhorst(P.s,P.rhos,P.sigs,O.s_pts);
% Rouwenhorst for monetary policy shock
[G.v_nodes,G.v_weight] = rouwenhorst(0,0,P.sigmp,O.mp_pts);
%--------------------------------------------------------------------------
% Set up grids by discretizing the state space
%--------------------------------------------------------------------------
% Define grid points for technology growth
G.g_grid = G.e_nodes';
% Define grid points for the risk premium
G.s_grid = G.u_nodes';
% Define grid points for the monetary policy shock
G.mp_grid = G.v_nodes';
% Define grid points for notional interest rate
G.in_grid = linspace(G.inmin,G.inmax,O.in_pts);
% Define grid points for consumption
G.c_grid = linspace(G.cmin,G.cmax,O.c_pts);
%--------------------------------------------------------------------------
% Construct grid arrays
%--------------------------------------------------------------------------
[G.g_gr,G.s_gr,G.mp_gr,G.in_gr,G.c_gr] = ...
    ndgrid(G.g_grid,G.s_grid,G.mp_grid,G.in_grid,G.c_grid);
G.nodes = numel(G.g_gr);
G.griddim = size(G.g_gr);