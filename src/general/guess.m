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
[T,~,eu] = linmodel(P,S,V);
disp(['Linear Solution: eu=' mat2str(eu)])

% Transform discretized state space to percent deviation from steady state     
k_gr_per = (G.k_gr - S.k)./S.k;   
z_gr_per = (G.z_gr - P.zbar)./P.zbar;
beta_gr_per = (G.beta_gr - P.beta)./P.beta;

% Calculate linear policy functions on discretized state space    
linpf_n = zeros(G.griddim);
linpf_pi = zeros(G.griddim);
linpf_i = zeros(G.griddim);
state = [k_gr_per(:) z_gr_per(:) beta_gr_per(:)]';
linpf_n(:) = T(V.n,[V.k V.z V.beta])*state;
linpf_pi(:) = T(V.pi,[V.k V.z V.beta])*state;
linpf_i(:) = T(V.i,[V.k V.z V.beta])*state;

% Convert back to levels
pf.n = P.n*(1 + linpf_n);
pf.pi = P.pi*(1 + linpf_pi);
pf.i = S.i*(1 + linpf_i);