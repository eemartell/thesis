function [pf,eu] = guess(P,S,G)

% pf = guess(P,S,G) 
%   Sets the initial policy functions
% Inputs:
%   P : structure of parameters
%   S : structure of steady state values 
%   G : structure of grids
% Outputs:
%  pf : structure of policy functions

%----------------------------------------------------------------------
% log-linear solution - ZLB not imposed
%----------------------------------------------------------------------
V = variables;
[T,~,eu] = linmodel(P,S,V);

% Transform discretized state space to level deviation from steady state
g_gr_per = G.g_gr - P.g;
s_gr_per = G.s_gr - P.s;
mp_gr_per = G.mp_gr;
in_gr_per = G.in_gr - S.i;
c_gr_per = G.c_gr - S.c;

% Calculate linear policy functions on discretized state space    
linpf_c = zeros(G.griddim);
linpf_pi = zeros(G.griddim);
state = [g_gr_per(:),s_gr_per(:),mp_gr_per(:),in_gr_per(:),c_gr_per(:)]';
linpf_c(:) = T(V.c,[V.g,V.s,V.mp,V.in,V.c])*state;
linpf_pi(:) = T(V.pi,[V.g,V.s,V.mp,V.in,V.c])*state;

% Convert back level deviations to levels
pf.c = S.c + linpf_c;
pf.pigap = (P.pi + linpf_pi)/P.pi;