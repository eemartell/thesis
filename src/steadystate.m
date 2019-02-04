function S = steadystate(P)

% S = steadystate(P) 
%   Computes the deterministic steady state
% Input:
%   P : Structure of parameters
% Output:
%   S : structure of steady state values

% Firm pricing
S.w = (P.theta-1)/P.theta;
% Production function
S.y = P.n;
% Aggregate resouce constraint
S.c = S.y;
% FOC labor
S.chi = S.w/(P.n^P.eta*S.c);
% FOC bond 
S.i = P.g*P.pi/(P.beta*(1+P.s));