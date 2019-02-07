clear
clc
% close all

postmeanflag = 0;

% Initialize simulation
nplot = 21;
shocksize = -2;

% Simulate

% Parameters, Steadystate, and Variables
% load('options');
load('options','V','P','F');
postmean = [...
          166.136743422136
          0.63609169716914
         0.797968531695273
         0.760553680680901
       0.00556322247731273
        0.0047217628003337
       0.00202764440559228
          1.95998009550232
         0.389219799038252];         
for iparam = 1:F.nparam
    eval(['P.' F.priors{iparam,1} ' = postmean(iparam);']);
end
S = steadystate(P);
        
% Solution
V = variables;
[T,M,eu] = linmodel(P,S,V);
disp(eu)

% Compute IRF
for ishock = 1:V.nshock
    eps = zeros(V.nshock,nplot);
    eps(ishock,2) = shocksize;
    y = zeros(V.nvar,nplot);
    for t = 2:nplot
        y(:,t) = T*y(:,t-1) + M*eps(:,t);
    end
    % Store solution
    irfs(ishock).y = y;
end

% Plot simulations
x = 0:(nplot-1);
figure(4);
clf
varidx = [V.yg,V.pi,V.i];
plotdim = [numel(varidx),V.nshock];
for isubplot = 1:prod(plotdim)
    subplot(plotdim(1),plotdim(2),isubplot); hold on; box on; grid on;
    [col,row] = ind2sub(plotdim([2 1]),isubplot);
    ishock = col; ivar = row;
    plot(x,400*irfs(ishock).y(varidx(ivar),:),'k-');
    if col == 1
        ylabel(V.desc(varidx(ivar)))
    end
    if row == 1
        title(V.shocktypes(ishock))
    end
end