load('../solutions/solutionfpART.mat')
c_ART = pf.hh;
pigap_ART = pf.firm;
load('../solutions/solutionfpGust.mat')
c_Gust = pf.hh;
c_zlb_Gust = pf.hh_zlb;
pigap_Gust = pf.firm;

gslice = ceil(O.g_pts/2);
mpslice = ceil(O.mp_pts/2);
inslice = 3; %in = 1.0001

set(0, 'DefaultAxesFontSize', 12)
x = 100*(G.s_grid-1);
y{1} = 100*(c_ART(gslice,:,mpslice,inslice)-S.c)/S.c;
z{1} = 100*(c_Gust(gslice,:,mpslice,inslice)-S.c)/S.c;
y{2} = 100*(pigap_ART(gslice,:,mpslice,inslice)-1);
z{2} = 100*(pigap_Gust(gslice,:,mpslice,inslice)-1);

ylabtxt{1} = '\textbf{Consumption ($c$)}';
ylabtxt{2} = '\textbf{Inflation gap ($\pi_{gap}$)}';

%   subplot padding
subpadbot = .295; % Increase if xaxis labels are cut off
subpadtop = .175; % Increase if subplot title is cut off
subpadleft = .1; % Increase if ylabel is cut off
subpadright = .05; % Decrease if white-space on RHS 

figure('Color', 'white'); hold on
set(gcf,'position',figbox);
plotdim = [2 1];

for i = 1:2
[col,row] = ind2sub(plotdim([2 1]),i);
left = (col-1+subpadleft)/plotdim(2);
bottom = 1-(row-subpadbot)/plotdim(1);
width = (1-(subpadleft+subpadright))/plotdim(2);
height = (1-subpadbot-subpadtop)/plotdim(1);
subplot(2,1,i, 'Position', [left bottom width height]); hold on; box on
xlm = [x(1) x(end)];
xlim(xlm)
ylm = [y{i}(end) y{i}(1)];
ylim(ylm)
plot(x, y{i},... 
    'Color','b')
    %'Linewidth', 1.8)
plot(x, z{i},... 
    'Color','r')
xlabel('Risk Premium ($s$)', 'Interpreter', 'latex',...
        'FontName','Times New Roman')
title(ylabtxt{i}, 'Interpreter', 'latex',...
      'FontName', 'Times New Roman')

legend('ART','GustEtAl')
%legtext = ['(\rho_s, \sigma_e, \sigma_u, \sigma_v) = (' ,num2str(P.rhos),', ', num2str(P.sige),', ', num2str(P.sigu), ', ' num2str(P.sigv), ')'];
%set(legend(legtext),'Interpreter','tex',...
%    'Location', 'southeast');
end
print(gcf,'-depsc2','-painters','figpfs.eps')
saveas(gcf,[savename 'figpfs.fig'])
