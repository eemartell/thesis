function R = eqm(x,state,O,P,S,G,pf,gpArr3,spArr3,weightArr3,varargin)

% Get original grids and GH nodes
if ~isempty(varargin)
    EEflag = 1;
    GH = varargin{1};
else
    EEflag = 0;
end

% Preallocate function output
R = zeros(size(x));

% State Values
g = state(1);       %Growth state current period
s = state(2);       %Preference state current period
mp = state(3);      %Monetary policy state current period
in = state(4);  	%Notional interest rate last period %%%no in

% Policy Function Guesses
cp = x(1);     %Consumption policy current period 
pigap = x(2);	%Inflation policy current period     
%----------------------------------------------------------------------
% Solve for variables
%----------------------------------------------------------------------
% Aggregate resource constraint (1)    
y = cp/(1-P.varphi*(pigap-1)^2/2);
% Interest rate rule (2,3)
inp = in^P.rhoi*(S.i*pigap^P.phipi)^(1-P.rhoi)*exp(mp); %%%make like Model.pdf
if P.zlbflag
    i = max(1,inp);
else
    i = inp;
end
% FOC Labor (4,5)
lam = cp;
w = S.chi*y^P.eta*lam;
%----------------------------------------------------------------------
% Linear interpolation of the policy variables
%----------------------------------------------------------------------
if ~EEflag  %%%no interpolation
    [cppArr3,pigappArr3] = Fallterp423_R(...
        O.g_pts,O.s_pts,O.mp_pts,O.in_pts,...
        G.in_grid,...
        inp,...
        pf.c,pf.pigap);
else
    % Growth rate (9)
    gpVec = ((1-P.rhog)*P.g + P.rhog*g + GH.e_nodes);
    % Risk premium (10)
    spVec = ((1-P.rhos)*P.s + P.rhos*s + GH.u_nodes);
    [cppArr3,pigappArr3] = allterp423(... %%%no interpolation
                            G.g_grid,G.s_grid,G.mp_grid,G.in_grid,...
                            gpVec,spVec,GH.v_nodes,inp,...
                            pf.c,pf.pigap);
    gpArr3 = gpVec(:,ones(GH.shockpts,1),ones(GH.shockpts,1));
end
%----------------------------------------------------------------------        
% Solve for variables inside expectations
%----------------------------------------------------------------------    
% Aggregate resource constraint (1)  
ypArr3 = cppArr3./(1-P.varphi*(pigappArr3-1).^2/2);
% Stochastic discount factor (4,5)
lampArr3 = cppArr3;
sdfArr3 = P.beta*lam./lampArr3;
%----------------------------------------------------------------------
% Numerical integration
%----------------------------------------------------------------------
% Compute all combinations of shocks
EbondArr3 = weightArr3.*sdfArr3./(gpArr3.*pigappArr3);
EfpArr3 = weightArr3.*sdfArr3.*(pigappArr3-1).*pigappArr3.*ypArr3;
% Integrate
Ebond = sum(EbondArr3(:));
Efp = sum(EfpArr3(:));
%----------------------------------------------------------------------
% First-order conditions
%----------------------------------------------------------------------
% Consumption Euler equation (6)
R(1) = 1 - s*i*Ebond/P.pi;
% Firm Pricing (Philips Curve, 7)
R(2) = P.varphi*(pigap-1)*pigap-(1-P.theta)-P.theta*w-P.varphi*Efp/y;
end
