% Data and options structure
clear
clc

% Model solution variables and options
V = variables;

% Solution space
%   Problem size
O.nstates = 4;
O.npfs = 2;
O.nshocks = 3;
%   Endogenous variable bounds
O.in_boundper = 0.025;
%   Density
O.g_pts = 7;     	% Growth state
O.s_pts = 7;     	% Risk premium shock state
O.mp_pts = 7;    	% Policy shock state
O.in_pts = 7;       	% Notional interest rate
O.e_pts = O.g_pts; 	% Growth shock
O.u_pts = O.s_pts;	% Risk premium shock
O.v_pts = O.mp_pts;	% MP shock

% Parameters, steady state, and grids
P = parameters;
P.zlbflag = 1;
S = steadystate(P);
G = grids(O,P,S);


% Save
save('solutions/options.mat','V','O','P','S','G')
