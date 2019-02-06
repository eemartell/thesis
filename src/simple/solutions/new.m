% Smoothness tests
clear
clc
close all

% Load solutions
load('solutionfpART.mat')
y_ART = pf.c(:);
c_ART = pf.c;

load('solutionfpGust.mat')
y_Gust = pf.c(:);
c_Gust = pf.c;
y_Gust_zlb = pf.c_zlb(:);

g = G.g_gr(:);
s = G.s_gr(:);
mp = G.mp_gr(:);
in = G.in_gr(:);
ds_Gust = dataset(y_Gust,g,s,mp,in);
ds_Gust_zlb = dataset(y_Gust_zlb,g,s,mp,in);
ds_ART = dataset(y_ART,g,s,mp,in);
% Linear regressions using fitlm
lm_Gust_lin = fitlm(ds_Gust,'y_Gust~g+s+mp+in');
lm_Gust_lin_zlb = fitlm(ds_Gust_zlb,'y_Gust_zlb~g+s+mp+in');
lm_ART_lin = fitlm(ds_ART,'y_ART~g+s+mp+in');
disp('RMSE (residual standard error) in linear models')
disp(['Gust, c: ', num2str(lm_Gust_lin.RMSE)])
disp(['Gust, c_zlb: ',num2str(lm_Gust_lin_zlb.RMSE)])
disp(['ART, c: ', num2str(lm_ART_lin.RMSE)])

ZLB_boundary = 3;
z{1} =  squeeze(c_ART(4,:,4,ZLB_boundary:end));
z{2} = squeeze(c_Gust(4,:,4,ZLB_boundary:end));
label{1} = 'ART solution method';
label{2} = 'GustEtAl solution method';
RMSE{1} = lm_ART_lin.RMSE;
RMSE{2} = lm_Gust_lin.RMSE;

% corresponds to ZLB points
in_small = squeeze(G.in_gr(4,:,4,1:3)); 
g_small = squeeze(G.s_gr(4,:,4,1:3));
c_small{1} = squeeze(c_ART(4,:,4,1:3));
c_small{2} = squeeze(c_Gust(4,:,4,1:3));

figure('Renderer', 'painters', 'Position', [100 100 900 450])

subplot(1,2,1)
surf(in_small, g_small, c_small{1})
colormap(gca,'gray'); 
colormap gray
freezeColors
xlim([.98, 1.04])
ylim([.98, 1.04])
zlim([.28, .38])
hold on 
surf(squeeze(G.in_gr(4,:,4,ZLB_boundary:end)),squeeze(G.s_gr(4,:,4,ZLB_boundary:end)), z{1})
colormap(gca,'jet'); 
colormap jet
freezeColors 
%hcb1=colorbar;
%hbc1.Label.String = 'non-ZLB locations';
hold off
xlabel('Interest rate')
ylabel('Risk premium')
zlabel('Consumption')
title(label{1})
text(1.03,1.04,.36,['RMSE: ',num2str(RMSE{1})])
text(1.03,1.04,.35,'Grayscale: ZLB')
[az, el] = view;
view(az-90,el-10)

subplot(1,2,2)
surf(squeeze(G.in_gr(4,:,4,ZLB_boundary:end)),squeeze(G.s_gr(4,:,4,ZLB_boundary:end)), z{2})
colormap(gca,'jet'); 
colormap jet
freezeColors
xlim([.98, 1.04])
ylim([.98, 1.04])
zlim([.28, .38])
hold on 
surf(in_small, g_small, c_small{2})
colormap(gca,'gray');
%hcb2=colorbar;
%hbc2.Label.String = 'ZLB locations';
hold off
xlabel('Interest rate')
ylabel('Risk premium')
zlabel('Consumption')
title(label{2})
text(1.03,1.04,.36,['RMSE: ',num2str(RMSE{2})])
text(1.03,1.04,.35,'Grayscale: ZLB')
[az, el] = view;
view(az-90,el-10)


% load clown.mat
% image(X)
% load penny.mat
% image(P)
% figure(1)
% hFig=gcf;
% set(hFig, 'Position', [50 50 300 500]);
% subplot(2,1,1)
% image(X);
% colormap(gca,'gray')
% hcb1=colorbar;
% subplot(2,1,2)
% image(P);
% colormap(gca,'jet')
% hcb2=colorbar;

saving = 'on';
savename = '../figs/Figs/pfs3D';
%% Save figure
if strcmp(saving,'on')
    print(gcf,'-depsc2','-painters',[savename '.eps'])
    saveas(gcf,[savename '.fig'])
end
