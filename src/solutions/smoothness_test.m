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

z{1} =  squeeze(c_ART(4,:,4,:));
z{2} = squeeze(c_Gust(4,:,4,:));
label{1} = 'ART solution method';
label{2} = 'GustEtAl solution method';
RMSE{1} = lm_ART_lin.RMSE;
RMSE{2} = lm_Gust_lin.RMSE;

figure('Renderer', 'painters', 'Position', [100 100 900 450])

for i = 1:2
subplot(1,2,i)
surf(squeeze(G.in_gr(4,:,4,:)),squeeze(G.s_gr(4,:,4,:)), z{i})
xlabel('Interest rate')
ylabel('Risk premium')
zlabel('Consumption')
xlim([.98, 1.04])
ylim([.98, 1.04])
zlim([.28, .38])
title(label{i})
text(.985,1.03,.36,['RMSE: ',num2str(RMSE{i})])
end

% % Linear regressions using regress
% [b_Gust_lin,~,resid_Gust,~,stats_Gust_lin] = regress(y_Gust, [ones(size(G.g_gr(:))) g s mp in]);
% [b_Gust_zlb_lin,~,resid_Gust_zlb,~,stats_Gust_zlb_lin] = regress(y_Gust, [ones(size(G.g_gr(:))) g s mp in]);
% [b_ART_lin,~,resid_ART,~,stats_ART_lin] = regress(y_ART, [ones(size(G.g_gr(:))) g s mp in]);
% % Linear regressions with quadratic terms using regress
% [b_Gust,~,~,~,stats_Gust] = regress(y_Gust, [ones(size(G.g_gr(:))) g s mp in s.^2 in.^2]);
% [b_Gust_zlb,~,~,~,stats_Gust_zlb] = regress(y_Gust_zlb, [ones(size(G.g_gr(:))) g s mp in s.^2 in.^2]);
% [b_ART,~,~,~,stats_ART] = regress(y_ART, [ones(size(G.g_gr(:))) g s mp in s.^2 in.^2]);

% ds_Gust = dataset(y_Gust(:),g,s,mp,in);
% ds_Gust_zlb = dataset(y_Gust_zlb(:),g,s,mp,in);
% ds_ART = dataset(y_ART(:),g,s,mp,in);
% % Linear regressions using fitlm
% lm_Gust_lin = fitlm(ds_Gust,'y_Gust~g+s+mp+in');
% lm_Gust_lin_zlb = fitlm(ds_Gust_zlb,'y_Gust_zlb~g+s+mp+in');
% lm_ART_lin = fitlm(ds_ART,'y_ART~g+s+mp+in');
% disp('RMSE (residual standard error) in linear models')
% disp(['Gust, c: ', num2str(lm_Gust_lin.RMSE)])
% disp(['Gust, c_zlb: ',num2str(lm_Gust_lin_zlb.RMSE)])
% disp(['ART, c: ', num2str(lm_ART_lin.RMSE)])
% %Linear regressions with quadratic terms using fitlm
% lm_Gust = fitlm(ds_Gust,'y_Gust~g+s+mp+in+s^2+in^2+mp^2');
% lm_Gust_zlb = fitlm(ds_Gust_zlb,'y_Gust_zlb~g+s+mp+in+s^2+in^2+mp^2');
% lm_ART = fitlm(ds_ART,'y_ART~g+s+mp+in+s^2+in^2+mp^2');

