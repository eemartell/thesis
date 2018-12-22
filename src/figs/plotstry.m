clear all
close all
clc

%%%calculate the nominal interest rate policy (41)
load('solutions/solution_test.mat')
inp = S.i*pf.pigap.^P.phipi;
r = max(1, inp);

%%%pick a cross section of the state space. look at the preferences and
%%%set the other shock processes to steady state (67)
gslice = ceil(O.g_pts/2);
mpslice = ceil(O.mp_pts/2);
findlastZLB = find(r(gslice,:,mpslice)==1,1,'last');

%%%find the ZLB region (143-144)
aper = 100*(G.a_grid-1);
shadeint = [aper(1) aper(findlastZLB)];

set(0, 'DefaultAxesFontSize', 12)
x = 100*(G.a_grid-1);
y{1} = 100*(pf.c(gslice,:,mpslice)-S.c)/S.c;
y{2} = 100*(pf.pigap(gslice,:,mpslice)-1);
y{3} = 100*(r(gslice,:,mpslice)-S.i)/S.i;
ytks{1} = -4:1:1;
ytks{2} = -1.5:.5:.5;
ytks{3} = -.75:.5:1;
plotline{1} = [0 0];
plotline{2} = [0 0];
plotline{3} = [0 0];
%txtbounds_a{1} = [-1.8, .35];
%txtbounds_a{2} = [-1.8, .16];
%txtbounds_a{3} = [-1.8, .3]; 
%txt{1} = '$\bar{c}$ = 0';
%txt{2} = '$\bar{\pi}_{gap}$ = 0';
%txt{3} = '$\bar{r}$ = 0';
%txtbounds_b{1} = [-1.8, .18];
%txtbounds_b{2} = [-1.8, .08];
%txtbounds_b{3} = [-1.8, .18]; 
ylabtxt{1} = '\textbf{Consumption ($c$)}';
ylabtxt{2} = '\textbf{Inflation gap ($\pi_{gap}$)}';
ylabtxt{3} = '\textbf{Interest rate ($r$)}';
%   subplot padding
subpadbot = .295; % Increase if xaxis labels are cut off
subpadtop = .175; % Increase if subplot title is cut off
subpadleft = .1; % Increase if ylabel is cut off
subpadright = .05; % Decrease if white-space on RHS 

figure('Color', 'white'); hold on
plotdim = [3 1];
   
for i = 1:3
[col,row] = ind2sub(plotdim([2 1]),i);
left = (col-1+subpadleft)/plotdim(2);
bottom = 1-(row-subpadbot)/plotdim(1);
width = (1-(subpadleft+subpadright))/plotdim(2);
height = (1-subpadbot-subpadtop)/plotdim(1);
subplot(3,1,i, 'Position', [left bottom width height]); hold on; box on
xlm = [x(1) x(end)];
xlim(xlm)
ylm = [y{i}(1) y{i}(end)];
ylim(ylm)
%%% shade the zlb area (179)
area(shadeint, ones(2,1)*max(ylm), min(ylm),...
    'FaceColor',[0.8 0.8 0.8],'EdgeColor','none',...
    'HandleVisibility','off')
plot(x, y{i},... 
    'Color','k',...
    'Linewidth', 1.8)
xlabel('Preferences ($a$)', 'Interpreter', 'latex',...
        'FontName','Times New Roman')
title(ylabtxt{i}, 'Interpreter', 'latex',...
      'FontName', 'Times New Roman')
%plot([(G.a_grid(1)-1)*100 (G.a_grid(end)-1)*100] ,plotline{i}, 'Color','k', 'Linewidth', 1)
xtks = ((0.97:0.01:1.03)-1)*100;
xticks(xtks)
yticks(ytks{i})
for xt = xtks
    plot([xt xt], ylm,...
        'LineStyle', ':',...
        'Color', 'k')
end
for yt = ytks{i}
    plot(xlim, [yt yt],...
        'LineStyle', ':',...
        'Color', 'k')
end
%text(txtbounds_a{i}(1), txtbounds_a{i}(2), txt{i}, 'Interpreter', 'Latex', 'Fontsize', 10)
%text(txtbounds_b{i}(1), txtbounds_b{i}(2), '\downarrow', 'Fontsize', 10)

legtext = ['(\rho_a, \sigma_e, \sigma_u, \sigma_v) = (' ,num2str(P.rhoa),', ', num2str(P.sige),', ', num2str(P.sigu), ', ' num2str(P.sigv), ')'];
set(legend(legtext),'Interpreter','tex',...
    'Location', 'southeast'); 
end