function pf = guess(V,P,S,G)

% pf = guess(V,P,S,G) 
%   Sets the initial policy functions
% Inputs:
%   V : structure of variable indices
%   P : structure of parameters
%   S : structure of steady state values 
%   G : structure of grids
% Outputs:
%  pf : structure of policy functions

%----------------------------------------------------------------------
% log-linear solution - ZLB not imposed
%----------------------------------------------------------------------
[T,~,eu] = linmodel(P,S,V);
disp(['Linear Solution: eu=' mat2str(eu)])

% Transform discretized state space to level deviation from steady state
g_gr_per = G.g_gr - P.g;
s_gr_per = G.s_gr - P.s;
mp_gr_per = G.mp_gr;
in_gr_per = G.in_gr - S.i;
c_gr_per = G.c_gr - S.c;
k_gr_per = G.k_gr - S.k;
x_gr_per = G.x_gr - S.x;
%w_gr_per = G.w_gr - S.w;

% Calculate linear policy functions on discretized state space 
linpf_c = zeros(G.griddim);
linpf_pi = zeros(G.griddim);
linpf_n = zeros(G.griddim);
linpf_q = zeros(G.griddim);
%linpf_ups = zeros(G.griddim);
state = [g_gr_per(:),s_gr_per(:),mp_gr_per(:),in_gr_per(:),c_gr_per(:),k_gr_per(:),x_gr_per(:)]';
linpf_c(:) = T(V.c,[V.g,V.s,V.mp,V.in,V.c,V.k,V.x])*state;
linpf_pi(:) = T(V.pi,[V.g,V.s,V.mp,V.in,V.c,V.k,V.x])*state;
linpf_n(:) = T(V.n,[V.g,V.s,V.mp,V.in,V.c,V.k,V.x])*state;
linpf_q(:) = T(V.q,[V.g,V.s,V.mp,V.in,V.c,V.k,V.x])*state;
%linpf_ups(:) = T(V.ups,[V.g,V.s,V.mp,V.in,V.c,V.k,V.x])*state;

% Convert back to levels
pf.c = S.c + linpf_c;
pf.pigap = (P.pi + linpf_pi)/P.pi;
pf.n = P.n + linpf_n;
pf.q = 1 + linpf_q;
%pf.ups = 1 + linpf_ups;  