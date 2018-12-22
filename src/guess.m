function [pf,eu] = guess(P,S,G,O)

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
if strcmp(O.alg,'ART')
    [T,M,eu] = linmodel(P,S,V);
elseif strcmp(O.alg,'Gust')
    [T,M,eu] = linmodel_gustetal(P,S,V);
end    
% Transform discretized state space to percent deviation from steady state
g_gr_per = G.g_gr./P.g - 1;
s_gr_per = G.s_gr/P.s - 1;
mp_gr_per = G.mp_gr;
in_gr_per = G.in_gr./S.i - 1;

% Calculate linear policy functions on discretized state space    
linpf_hh = zeros(G.griddim);
linpf_firm = zeros(G.griddim);
state = [g_gr_per(:),s_gr_per(:),mp_gr_per(:),in_gr_per(:)]';
linpf_hh(:) = T(V.c,[V.g,V.s,V.mp,V.in])*state;
linpf_firm(:) = T(V.pi,[V.g,V.s,V.mp,V.in])*state;

% Convert back to levels
if strcmp(O.alg,'ART')
    pf.hh = S.c*(1 + linpf_hh); 
    pf.firm = 1 + linpf_firm;
elseif strcmp(O.alg,'Gust')
    pf.hh = S.c*(1 + linpf_hh); %??? 
    pf.firm = 1 + linpf_firm; %??? CHECK THIS
end
