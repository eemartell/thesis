% Smoothness tests
clear
clc
close all

% Load solutions
load('solution_test.mat')
y_ART = pf.c(:);
c_ART = pf.c;
inp_ART = G.in_gr.^P.rhoi.*(S.i*pf.pigap.^P.phipi).^(1-P.rhoi).*exp(G.mp_gr);

g = G.g_gr(:);
s = G.s_gr(:);
mp = G.mp_gr(:);
in = G.in_gr(:);
c = G.c_gr(:);
ds_ART = dataset(y_ART,g,s,mp,in,c);
% Linear regressions using fitlm
lm_ART_lin = fitlm(ds_ART,'y_ART~g+s+mp+in+c');
disp('RMSE (residual standard error) in linear models')
disp(['ART, c: ', num2str(lm_ART_lin.RMSE)])

z{1} =  squeeze(c_ART(4,:,4,:,4));
z_zlb{1} = z{1};
inp{1} = squeeze(inp_ART(4,:,4,:,4));
z_zlb{1}(inp{1}>1) = nan;

label{1} = 'ART, risk premium shock';

RMSE{1} = lm_ART_lin.RMSE;

figure('Renderer', 'painters', 'Position', [100 100 825 525])

in_vec = squeeze(G.in_gr(4,:,4,:,4));
in_vec = in_vec(:);
s_vec = squeeze(G.s_gr(4,:,4,:,4));
s_vec = s_vec(:);

for i = 1:1
%hold on
surf(squeeze(G.in_gr(4,:,4,:,4)),squeeze(G.s_gr(4,:,4,:,4)), z{i})
colormap winter
freezeColors %from MATLAB file exchange
%xlim([.98, 1.04])
%ylim([.98, 1.04])
%zlim([.28, .39])
hold on
scatter3(in_vec,s_vec, z_zlb{i}(:),15,'MarkerEdgeColor','k',...
        'MarkerFaceColor','k')
hold off
xlabel('Interest rate')
ylabel('Risk premium')
zlabel('Consumption')
title(label{i})
[az, el] = view;
view(az-90,el-10)
text(1.03,1.04,.37,['RMSE: ',num2str(RMSE{i})])
text(1.03,1.04,.36,'Black dots: ZLB')
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
saving = 'on';
savename = 'pfs3D';
%% Save figure
if strcmp(saving,'on')
    print(gcf,'-depsc2','-painters',[savename '.eps'])
    saveas(gcf,[savename '.fig'])
end
