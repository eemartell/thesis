% Smoothness tests
clear
clc
close all

% Load solutions
load('solutionfpART.mat')
y_ART = pf.c(:);
c_ART = pf.c;
inp_ART = G.in_gr.^P.rhoi.*(S.i*pf.pigap.^P.phipi).^(1-P.rhoi).*exp(G.mp_gr);

load('solutionfpGust.mat')
y_Gust = pf.c(:);
c_Gust = pf.c;
inp_Gust = G.in_gr.^P.rhoi.*(S.i*pf.pigap.^P.phipi).^(1-P.rhoi).*exp(G.mp_gr);
y_Gust_zlb = pf.c_zlb(:);
c_Gust_zlb = pf.c_zlb;
ZLBlocs = find(inp_Gust <= 1);
y_Gust_agg = y_Gust;
y_Gust_agg(ZLBlocs) = y_Gust_zlb(ZLBlocs);
c_Gust_agg = y_Gust;
c_Gust_agg(ZLBlocs) = y_Gust_zlb(ZLBlocs);
c_Gust_agg = reshape(c_Gust_agg,size(c_Gust));
g = G.g_gr(:);
s = G.s_gr(:);
mp = G.mp_gr(:);
in = G.in_gr(:);
ds_Gust = dataset(y_Gust,g,s,mp,in);
ds_Gust_zlb = dataset(y_Gust_zlb,g,s,mp,in);
ds_Gust_agg = dataset(y_Gust_agg,g,s,mp,in);
ds_ART = dataset(y_ART,g,s,mp,in);
% Linear regressions using fitlm
lm_Gust_lin = fitlm(ds_Gust,'y_Gust~g+s+mp+in');
lm_Gust_lin_zlb = fitlm(ds_Gust_zlb,'y_Gust_zlb~g+s+mp+in');
lm_Gust_lin_agg = fitlm(ds_Gust_agg,'y_Gust_agg~g+s+mp+in');
lm_ART_lin = fitlm(ds_ART,'y_ART~g+s+mp+in');

z{1} =  squeeze(c_ART(4,:,4,:));
z_zlb{1} = z{1};
inp{1} = squeeze(inp_ART(4,:,4,:));
z_zlb{1}(inp{1}>1) = nan;

z{2} = squeeze(c_Gust(4,:,4,:));
z_zlb{2} = z{2};
inp{2} = squeeze(inp_Gust(4,:,4,:));
z_zlb{2}(inp{2}>1) = nan;

z{3} = squeeze(c_Gust_zlb(4,:,4,:));
z_zlb{3} = z{3};
z_zlb{3}(inp{2}>1) = nan;

z{4} = squeeze(c_Gust_agg(4,:,4,:));
z_zlb{4} = z{4};
z_zlb{4}(inp{2}>1) = nan;

label{1} = 'ART policy function';
label{2} = 'GustEtAl non-ZLB policy function';
label{3} = 'GustEtAl ZLB policy function';
label{4} = 'GustEtAl aggregated policy function';

RMSE{1} = lm_ART_lin.RMSE;
RMSE{2} = lm_Gust_lin.RMSE;
RMSE{3} = lm_Gust_lin_zlb.RMSE;
RMSE{3} = lm_Gust_lin_agg.RMSE;
FIT{1} = lm_ART_lin.Fitted;
FIT{2} = lm_Gust_lin.Fitted;
FIT{3} = lm_Gust_lin_zlb.Fitted;
FIT{4} = lm_Gust_lin_agg.Fitted;
perc_errART = abs(y_ART - FIT{1})./S.y*100;
meanperc_err{1} = mean(perc_errART(:));
perc_errGust = abs(y_Gust - FIT{2})./S.y*100;
meanperc_err{2} = mean(perc_errGust(:));
perc_errGust_zlb = abs(y_Gust_zlb - FIT{3})./S.y*100;
meanperc_err{3} = mean(perc_errGust_zlb(:));
perc_errGust_agg = abs(y_Gust_agg - FIT{4})./S.y*100;
meanperc_err{4} = mean(perc_errGust_agg(:));

% disp('RMSE (residual standard error) in linear models')
% disp(['ART, c: ', num2str(lm_ART_lin.RMSE), ' consumption units'])
% disp(['Gust, c: ', num2str(lm_Gust_lin.RMSE), ' consumption units'])
% disp(['Gust, c_zlb: ',num2str(lm_Gust_lin_zlb.RMSE), ' consumption units'])

disp('Average percent error of consumption policy function from linear fitted model')
disp(['ART, c: ', num2str(meanperc_err{1}), '%'])
disp(['Gust, c: ',num2str(meanperc_err{2}), '%'])
disp(['Gust, c_zlb: ', num2str(meanperc_err{3}), '%'])
disp(['Gust, c_agg: ', num2str(meanperc_err{4}), '%'])

figure('Renderer', 'painters', 'Position', [100 100 825 525])

in_vec = squeeze(G.in_gr(4,:,4,:));
in_vec = in_vec(:);
s_vec = squeeze(G.s_gr(4,:,4,:));
s_vec = s_vec(:);

for i = 1:3
    if i < 3
    subplot(2,2,i)
    else
    subplot(2,2,i+1)
    end
%hold on
surf(squeeze(G.in_gr(4,:,4,:)),squeeze(G.s_gr(4,:,4,:)), z{i})
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
text(1.013,1.04,.37,['avg % err: ',num2str(meanperc_err{i})])
text(1.013,1.04,.36,'Black dots: ZLB')
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
savename = '../figs/Figs/pfs3D';
%% Save figure
if strcmp(saving,'on')
    print(gcf,'-depsc2','-painters',[savename '.eps'])
    saveas(gcf,[savename '.fig'])
end
