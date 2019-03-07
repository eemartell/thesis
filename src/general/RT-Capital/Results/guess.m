function pf = guess(P,S,G)

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
[T,M,eu] = linmodel(P,S,V);
disp(['Linear Solution: eu=' mat2str(eu)])

% Transform discretized state space to percent deviation from steady state
c_gr_per = G.c_gr./S.c - 1;
k_gr_per = G.k_gr./S.k - 1;
i_gr_per = G.i_gr./S.i - 1;
rn_gr_per = G.rn_gr./S.r - 1;
g_gr_per = (G.g_gr - P.g)./P.g;
beta_gr_per = (G.beta_gr - P.beta)/P.beta;
mp_gr_per = G.mp_gr;

% Calculate linear policy functions on discretized state space    
linpf_n = zeros(G.griddim);
linpf_q = zeros(G.griddim);
linpf_pi = zeros(G.griddim);
linpf_psi = zeros(G.griddim);
state = [c_gr_per(:),k_gr_per(:),i_gr_per(:),rn_gr_per(:),g_gr_per(:),beta_gr_per(:)]';
linpf_n(:) = T(V.n,[V.c,V.k,V.i,V.rn,V.g,V.beta])*state + M(V.n,V.epsmp)*mp_gr_per(:)';
linpf_q(:) = T(V.q,[V.c,V.k,V.i,V.rn,V.g,V.beta])*state + M(V.q,V.epsmp)*mp_gr_per(:)';
linpf_pi(:) = T(V.pi,[V.c,V.k,V.i,V.rn,V.g,V.beta])*state + M(V.pi,V.epsmp)*mp_gr_per(:)';
linpf_psi(:) = T(V.psi,[V.c,V.k,V.i,V.rn,V.g,V.beta])*state + M(V.psi,V.epsmp)*mp_gr_per(:)';

% Convert back to levels
pf.n = P.n*(1 + linpf_n);
pf.q = S.q*(1 + linpf_q);
pf.pi = P.pi*(1 + linpf_pi);
pf.psi = S.psi*(1 + linpf_psi);  