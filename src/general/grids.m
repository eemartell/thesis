function G = grids(O,P)

% G = grids(O,P)
%   Constructs an evenly-spaced discretized state space
% Inputs:
%	O :   Structure of user-specified options
%   P :   Structure of parameters
% Output:
%   G :   Structure of grid points

%--------------------------------------------------------------------------
% Set up grids by discretizing the state space
%--------------------------------------------------------------------------
% Define grid points for capital
G.k_grid = linspace(O.kbound(1),O.kbound(2),O.k_pts);
% Define grid points for technology shock
G.z_grid = linspace(O.zbound(1),O.zbound(2),O.z_pts);
% Define grid points for discount factor shock
G.beta_grid = linspace(O.betabound(1),O.betabound(2),O.beta_pts);
%--------------------------------------------------------------------------
% Weights for numerical integration from truncated normal
%--------------------------------------------------------------------------
[e_nodes,G.e_weight] = ghquad(O.e_pts);
G.e_nodes = (2^.5) * P.sige * e_nodes;
[u_nodes,G.u_weight] = ghquad(O.u_pts);
G.u_nodes = (2^.5) * P.sigu * u_nodes;
%--------------------------------------------------------------------------
% Construct grid arrays
%--------------------------------------------------------------------------
[G.k_gr,G.z_gr,G.beta_gr] = ndgrid(G.k_grid,G.z_grid,G.beta_grid);
G.nodes = numel(G.k_gr);
G.griddim = size(G.k_gr);