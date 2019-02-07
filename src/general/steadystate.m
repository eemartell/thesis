function S = steadystate(P)

% S = steadystate(P) 
%   Computes the deterministic steady state
% Input:
%   P : Structure of parameters
% Output:
%   S : structure of steady state values

% FOC bond 
S.i = P.g*P.pi/(P.beta*P.s);
S.in = S.i;
% Price Phillips curve
S.w = (P.thetap-1)/P.thetap;
% Production function
S.yf = P.n;
% Output definition
S.y = S.yf;
% Aggregate resouce constraint
S.c = S.y;
% Marginal utility of consumption
S.lam = (1-P.h/P.g)*S.c;
% HH FOC labor
S.chi = S.w/(P.n^P.eta*S.lam);