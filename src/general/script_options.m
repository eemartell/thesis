clear; clc;

% Problem size
O.nstates = 8;
O.npfs = 5;
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
%O.w_boundper = 0.04;

%   Density
O.g_pts = 3;            % Growth shock state
O.s_pts = 5;%13;           % Risk premium shock state
O.mp_pts = 3;           % Monetary policy shock state
O.in_pts = 5;           % Notional interest rate
O.c_pts = 5;            % Consumption
O.k_pts = 5;            % Capital
O.x_pts = 5;%11;            % Investment
%O.w_pts = 7;            % Real wage state

O.epsg_pts = O.g_pts;	% Growth shock
O.epss_pts = O.s_pts;	% Risk premium shock
O.epsmp_pts = O.mp_pts;	% Monetary Policy shock

% Grids
G = grids(O,P,S);

save('options.mat','O','V','P','S','G')