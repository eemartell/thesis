clear; clc;

% Problem size
O.nstates = 7;
O.npfs = 4;
O.nshocks = 3;

% Load parameters, steady state and grids
V = variables;
P = parameters;
S = steadystate(P);

% Specify grid options
%   Bounds
O.in_boundper = 0.025;
O.c_boundper = 0.06;
O.k_boundper = 0.08;
O.x_boundper = 0.15;

%   Density
O.g_pts = 7;            % Growth shock state
O.s_pts = 13;           % Risk premium shock state
O.mp_pts = 7;           % Monetary policy shock state
O.in_pts = 7;           % Notional interest rate
O.c_pts = 7;            % Consumption
O.k_pts = 7;            % Capital
O.x_pts = 11;            % Investment

O.epsg_pts = O.g_pts;	% Growth shock
O.epss_pts = O.s_pts;	% Risk premium shock
O.epsmp_pts = O.mp_pts;	% Monetary Policy shock

% Grids
G = grids(O,P,S);

save('options.mat','O','V','P','S','G')