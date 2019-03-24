% Smoothness tests
clear
clc
close all

% Load solutions
load('solutionfpART.mat')
y_ART = pf.n(:);
n_ART = pf.n;
% Production function (2)
y = (G.k_gr./G.g_gr).^P.alpha.*pf.n.^(1-P.alpha); 
% Real gdp
rgdp = G.c_gr + G.x_gr;
rgdpp = (1-P.varphi.*(pf.pigap-1).^2/2).*y;
% Output growth
rgdpg = G.g_gr.*rgdpp./(P.g.*rgdp); 
 %   Interest rate rule
inp_ART = G.in_gr.^P.rhoi.*(S.i*pf.pigap.^P.phipi.*rgdpg.^P.phiy).^(1-P.rhoi).*exp(G.mp_gr);   
%inp_ART = G.in_gr.^P.rhoi.*(S.i*pf.pigap.^P.phipi).^(1-P.rhoi).*exp(G.mp_gr);

load('solutionfpGust.mat')
y_Gust = pf.n(:);
n_Gust = pf.n;
%inp_Gust = G.in_gr.^P.rhoi.*(S.i*pf.pigap.^P.phipi).^(1-P.rhoi).*exp(G.mp_gr);
y_Gust_zlb = pf.n_zlb(:);
n_Gust_zlb = pf.n_zlb;
% Production function (2)
y = (G.k_gr./G.g_gr).^P.alpha.*pf.n.^(1-P.alpha); 
% Real gdp
rgdp = G.c_gr + G.x_gr;
rgdpp = (1-P.varphi.*(pf.pigap-1).^2/2).*y;
% Output growth
rgdpg = G.g_gr.*rgdpp./(P.g.*rgdp); 
 %   Interest rate rule
inp_Gust = G.in_gr.^P.rhoi.*(S.i*pf.pigap.^P.phipi.*rgdpg.^P.phiy).^(1-P.rhoi).*exp(G.mp_gr);  
ZLBlocs = find(inp_Gust<=1);
y_Gust_agg = y_Gust;
y_Gust_agg(ZLBlocs) = y_Gust_zlb(ZLBlocs);
n_Gust_agg = reshape(y_Gust_agg,size(n_Gust));

g = G.g_gr(:);
s = G.s_gr(:);
mp = G.mp_gr(:);
in = G.in_gr(:);
c = G.c_gr(:);
k = G.k_gr(:);
x = G.x_gr(:);
ds_Gust = dataset(y_Gust,g,s,mp,in,c,k,x);
ds_Gust_zlb = dataset(y_Gust_zlb,g,s,mp,in,c,k,x);
ds_Gust_agg = dataset(y_Gust_agg,g,s,mp,in,c,k,x);
ds_ART = dataset(y_ART,g,s,mp,in,c,k,x);
% Linear regressions using fitlm
lm_Gust_lin = fitlm(ds_Gust,'y_Gust~g+s+mp+in+c+k+x');
lm_Gust_lin_zlb = fitlm(ds_Gust_zlb,'y_Gust_zlb~g+s+mp+in+c+k+x');
lm_Gust_lin_agg = fitlm(ds_Gust_agg,'y_Gust_agg~g+s+mp+in+c+k+x');
lm_ART_lin = fitlm(ds_ART,'y_ART~g+s+mp+in+c+k+x');
disp('RMSE (residual standard error) in linear models')
disp(['Gust, n: ', num2str(lm_Gust_lin.RMSE)])
disp(['Gust, n_zlb: ',num2str(lm_Gust_lin_zlb.RMSE)])
disp(['Gust, n_agg: ',num2str(lm_Gust_lin_agg.RMSE)])
disp(['ART, n: ', num2str(lm_ART_lin.RMSE)])

z{1} =  squeeze(n_ART(2,:,2,:,3,3,3));
z_zlb{1} = z{1};
inp{1} = squeeze(inp_ART(2,:,2,:,3,3,3));
z_zlb{1}(inp{1}>1) = nan;

z{2} = squeeze(n_Gust(2,:,2,:,3,3,3));
z_zlb{2} = z{2};
inp{2} = squeeze(inp_Gust(2,:,2,:,3,3,3));
z_zlb{2}(inp{2}>1) = nan;

z{3} = squeeze(n_Gust_zlb(2,:,2,:,3,3,3));
z_zlb{3} = z{3};
z_zlb{3}(inp{2}>1) = nan;

z{4} = squeeze(n_Gust_agg(2,:,2,:,3,3,3));
z_zlb{4} = z{4};
z_zlb{4}(inp{2}>1) = nan;

label{1} = 'ART policy function';
label{2} = 'GustEtAl non-ZLB policy function';
label{3} = 'GustEtAl ZLB policy function';
label{4} = 'GustEtAl combined policy function';

RMSE{1} = lm_ART_lin.RMSE;
RMSE{2} = lm_Gust_lin.RMSE;
RMSE{3} = lm_Gust_lin_zlb.RMSE;
RMSE{4} = lm_Gust_lin_agg.RMSE;
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

disp('RMSE (residual standard error) in linear models')
disp(['ART, n: ', num2str(lm_ART_lin.RMSE), ' consumption units'])
disp(['Gust, n: ', num2str(lm_Gust_lin.RMSE), ' consumption units'])
disp(['Gust, n_zlb: ',num2str(lm_Gust_lin_zlb.RMSE), ' consumption units'])
disp(['Gust, n_agg: ',num2str(lm_Gust_lin_agg.RMSE), ' consumption units'])

disp('Average percent error of consumption policy function from linear fitted model')
disp(['ART, n: ', num2str(meanperc_err{1}), '%'])
disp(['Gust, n: ',num2str(meanperc_err{2}), '%'])
disp(['Gust, n_zlb: ', num2str(meanperc_err{3}), '%'])
disp(['Gust, n_agg: ', num2str(meanperc_err{4}), '%'])

figure('Renderer', 'painters', 'Position', [100 100 825 525])

in_vec = squeeze(G.in_gr(2,:,2,:,3,3,3));
in_vec = in_vec(:);
s_vec = squeeze(G.s_gr(2,:,2,:,3,3,3));
s_vec = s_vec(:);



for i = 1:4
subplot(2,2,i)
%     if i ~= 3
%     subplot(2,2,i)
%     else
%     subplot(2,2,i+1)
%     end
%hold on
surf(squeeze(G.in_gr(2,:,2,:,3,3,3)),squeeze(G.s_gr(2,:,2,:,3,3,3)), z{i})
zlim([.32,.37])
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
zlabel('Labor')
title(label{i})
[az, el] = view;
view(az+90-25,el-20)
% view(az-90,el-10)
%view(-37.5-90,30-10)
text(0.999,1.01,.365,['mean \% err: ',num2str(meanperc_err{i},'%.6f'),'\%'])
text(0.999,1.01,.36,'Black dots: ZLB')
end

saving = 'on';
savename = '../figs/pfs3D';
%% Save figure
if strcmp(saving,'on')
    savefig([savename '.fig'])
end
