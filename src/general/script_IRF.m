clear
clc
% close all

postmeanflag = 0;

% Initialize simulation
nplot = 21;
shocksize = -2;

% Simulate

% Parameters, Steadystate, and Variables
load('options');
% P.alpha = 0.001;
% P.thetaw = 1e8;
S = steadystate(P);

if postmeanflag
    % Results options
    MHstage = 1;
    nburn = 0;
    thin = 1;

    filename = ['estimation_MHtemp' num2str(MHstage) '.mat'];
    load(filename)
    
    % Remove burn and thin
    disp(['Total draws:     ' num2str(idraw)]);
    draws = draws(nburn+1:thin:idraw,:); 
    [ndraws,nparam] = size(draws);
    drawpostmean = mean(draws,1);
    
    % Update parameters and steady state
    for iparam = 1:nparam
        eval(['P.' F.priors{iparam,1} ' = drawpostmean(iparam);']);
    end
    S = steadystate(P);
end

% Solution
thisdir = cd;
%   Capital
disp('Capital solution')
[cap.T,cap.M,eu] = linmodel(P,S,V);
disp(eu)
%   No capital
disp('No capital solution')
cd('../Results-nocap')
nocap = load('options.mat');
[nocap.T,nocap.M,eu] = linmodel(P,S,nocap.V);	
disp(eu)
cd(thisdir)

% Compute IRF
for ishock = 1:V.nshock
    eps = zeros(V.nshock,nplot);
    eps(ishock,2) = shocksize;
    ycap = zeros(V.nvar,nplot);
    ynocap = zeros(nocap.V.nvar,nplot);
    for t = 2:nplot
        ycap(:,t) = cap.T*ycap(:,t-1) + cap.M*eps(:,t);
        ynocap(:,t) = nocap.T*ynocap(:,t-1) + nocap.M*eps(1:nocap.V.nshock,t);
    end
    % Store solution
    irfs(ishock).ycap = ycap;
    irfs(ishock).ynocap = ynocap;
end

% Plot simulations
x = 0:(nplot-1);
figure(3);
% clf
varidx = [V.yg,V.pi,V.i,V.x,V.u];
varidx2 = [nocap.V.yg,nocap.V.pi,nocap.V.i];
plotdim = [numel(varidx),V.nshock];
for isubplot = 1:prod(plotdim)
    subplot(plotdim(1),plotdim(2),isubplot); hold on; box on; grid on;
    [col,row] = ind2sub(plotdim([2 1]),isubplot);
    ishock = col; ivar = row;
    plot(x,400*irfs(ishock).ycap(varidx(ivar),:),'r-.');
%     if ivar <= numel(varidx2)
%         plot(x,400*irfs(ishock).ynocap(varidx2(ivar),:),'b--');
%     end
    if col == 1
        ylabel(V.desc(varidx(ivar)))
    end
    if row == 1
        title(V.shocktypes(ishock))
    end
end
