clear; clc;

% Problem size
O.nstates = 5;
O.npfs = 2;
O.nshocks = 3;

% Load parameters, steady state and grids
V = variables;
P = parameters;
S = steadystate(P);

% Specify grid options
%   Bounds
O.in_boundper = 0.025;
O.c_boundper = 0.06;

%   Density
O.g_pts = 7;            % Growth shock state
O.s_pts = 8;           % Risk premium shock state
O.mp_pts = 7;           % Monetary policy shock state
O.in_pts = 8;           % Notional interest rate
O.c_pts = 7;            % Consumption

O.epsg_pts = O.g_pts;	% Growth shock
O.epss_pts = O.s_pts;	% Risk premium shock
O.epsmp_pts = O.mp_pts;	% Monetary Policy shock

% Grids
G = grids(O,P,S);

% Estimation options
F.nparam = 9;
F.paramstrs = {...
    'Rotemberg Price Adj. Cost',...
    'Habit Persistence','Risk Premium Persistence','Interest Rate Persistence',...
    'Growth Shock SD','Preference Shock SD','MP Shock SD',...
    'Inflation Gap Response','Output Gap Response'};
F.paramtex = {...
    '\varphi_p','h','\rho_s','\rho_i',...
    '\sigma_g','\sigma_s','\sigma_i',...
    '\phi_\pi','\phi_y',};
F.priors = {...  
    'varphip',  'norm',   100,      25,      0.01,   1000;   % Price adjustment cost
    'h',        'beta',   0.8,      0.1,     0,      0.99;   % Habit parameter
    'rhos',     'beta',   0.8,      0.1,     0,      0.99;   % Preference shock persistence
    'rhoi',     'beta',   0.8,      0.1,     0,      0.99;   % Notional rate persistence
    'sigg',     'invgam', 0.0050, 	0.0050,	 1e-5,   1;      % Growth shock std
    'sigs',     'invgam', 0.0050,   0.0050,  1e-5,   1;      % Preference shock std
    'sigmp',    'invgam', 0.0020,   0.0020,  1e-5,   1;      % Interest rate shock std
    'phipi',    'norm',   2.0,      0.25,    1.01,   20;     % Inflation response
    'phiy',     'norm',   0.5,      0.25,    0,      5};     % Output response

save('options.mat','O','V','P','S','G','F')