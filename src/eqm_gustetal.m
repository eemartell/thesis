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
Rdim = size(R,2);

% State Values
g = state(1);       %Growth state current period
s = state(2);       %Preference state current period
mp = state(3);      %Monetary policy state current period
in = state(4);  	%Notional interest rate last period

for icol = 1:Rdim   
    % Policy Function Guesses
% Put V functions here instead
    Vlamp = x(1,icol); %Vlam policy current period
    Vpip = x(2,icol);     % Vpi policy current period
    %cp = x(1,icol);     %Consumption policy current period 
    %pigap = x(2,icol);	%Inflation policy current period     
    %----------------------------------------------------------------------
    % Solve for variables
    %----------------------------------------------------------------------
% solve for other time t variables, in particular c and pigap, given state and the V functions
% with goal of updating the state variables
    lam = 1/Vlamp;
    c = lam;
    pigap = (1+sqrt((P.varphi + 4*Vpip)/P.varphi))/2;
    % Aggregate resource constraint (6)    
    y = c/(1-P.varphi*(pigap-1)^2/2);
    % Interest rate rule (8)
    inp = in^P.rhoi*(S.i*pigap^P.phipi)^(1-P.rhoi)*exp(mp);
    if P.zlbflag
        i = max(1,inp);
    else
        i = inp;
    end
    % FOC Labor (1,2)
    w = S.chi*y^P.eta*lam;
    %----------------------------------------------------------------------
    % Linear interpolation of the policy variables
    %----------------------------------------------------------------------
% interpolate V functions (instead of c and pigap) at updated state
    if ~EEflag         
        %[cppArr3,pigappArr3] = Fallterp423_R(...
        [VlampArr3,VpipArr3] = Fallterp423_R(...
            O.g_pts,O.s_pts,O.mp_pts,O.in_pts,...
            G.in_grid,...
            inp,...
            pf.Vlam,pf.Vpi);
            %pf.c,pf.pigap);
    else
        % Growth rate
        gpVec = ((1-P.rhog)*P.g + P.rhog*g + GH.e_nodes);
        % Preferences 
        spVec = ((1-P.rhos)*P.s + P.rhos*s + GH.u_nodes);
        %[cppArr3,pigappArr3] = allterp423(...
        [VlampArr3,VpipArr3] = allterp423(...
                                G.g_grid,G.s_grid,G.mp_grid,G.in_grid,...
                                gpVec,spVec,GH.v_nodes,inp,...
                                pf.Vlam,pf.Vpi);
                                %pf.c,pf.pigap);
        gpArr3 = gpVec(:,ones(GH.shockpts,1),ones(GH.shockpts,1));
    end

    %----------------------------------------------------------------------        
    % Solve for variables inside expectations
    %---------------------------------------------------------------------- 
% replicate the update at time t for t+1 variables instead only using the
% equations you need to get variables to evaulate the expectation operators
    lampArr3 = VlampArr3;
    cppArr3 = lampArr3;
    pigappArr3 = (1+sqrt((P.varphi + 4*VpipArr3)/P.varphi))/2;
    % Aggregate resource constraint  
    ypArr3 = cppArr3./(1-P.varphi*(pigappArr3-1).^2/2);
    % Stochastic discount factor
    sdfArr3 = P.beta*lam./lampArr3;
    %----------------------------------------------------------------------
    % Numerical integration
    %----------------------------------------------------------------------
% Integration probably won't change because there are no Vs at this point
    % Compute all combinations of shocks
    EbondArr3 = weightArr3.*sdfArr3./(gpArr3.*pigappArr3);
    EfpArr3 = weightArr3.*sdfArr3.*(pigappArr3-1).*pigappArr3.*ypArr3;
    % Integrate
    Ebond = sum(EbondArr3(:));
    Efp = sum(EfpArr3(:));
    %----------------------------------------------------------------------
    % First-order conditions
    %----------------------------------------------------------------------
% probably doesn't change unless they are writtin in terms of V somehow
    % Consumption Euler equation (3)
    R(1,icol) = 1 - s*i*Ebond/P.pi;
    % Firm Pricing (Philips Curve, 5)
    R(2,icol) = P.varphi*(pigap-1)*pigap-(1-P.theta)-P.theta*w-P.varphi*Efp/y;
end
