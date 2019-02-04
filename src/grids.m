function G = grids(O,P,S)

% G = grids(O,P)
%   Constructs an evenly-spaced discretized state space
% Inputs:
%	O :   Structure of user-specified options
%   P :   Structure of parameters
%   P :   Structure of steady state values
% Output:
%   G :   Structure of grid points

%--------------------------------------------------------------------------
% Set endogenous state space bounds
%--------------------------------------------------------------------------
G.inmin = (1-O.in_boundper)*S.i;
G.inmax = (1+O.in_boundper)*S.i;
%--------------------------------------------------------------------------
% Nodes and weights for numerical integration
%--------------------------------------------------------------------------
% Rouwenhorst for productivity growth process
[G.e_nodes,G.e_weight] = rouwenhorst(P.g,P.rhog,P.sige,O.g_pts);
% Rouwenhorst for preference shock
[G.u_nodes,G.u_weight] = rouwenhorst(1,P.rhoa,P.sigu,O.a_pts);
% Rouwenhorst for monetary policy shock
[G.v_nodes,G.v_weight] = rouwenhorst(0,P.rhomp,P.sigv,O.mp_pts);
%--------------------------------------------------------------------------
% Set up grids by discretizing the state space
%--------------------------------------------------------------------------
% Define grid points for prodcutivity growth
G.g_grid = G.e_nodes';
% Define grid points for the risk premium
G.a_grid = G.u_nodes';
% Define grid points for the monetary policy shock
G.mp_grid = G.v_nodes';
% Define grid points for notional interest rate
G.in_grid = linspace(G.inmin,G.inmax,O.in_pts);
%--------------------------------------------------------------------------
% Construct grid arrays
%--------------------------------------------------------------------------
[G.g_gr,G.a_gr,G.mp_gr,G.in_gr] = ...
    ndgrid(G.g_grid,G.a_grid,G.mp_grid,G.in_grid);
G.nodes = numel(G.g_gr);
G.griddim = size(G.g_gr);
