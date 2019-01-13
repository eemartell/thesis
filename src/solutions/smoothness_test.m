% Smoothness tests
clear
clc
close all

% Load solutions
load('solutionfpART.mat')
y_ART = pf.hh(:);
load('solutionfpGust.mat')
y_Gust = pf.hh(:);
y_Gust_zlb = pf.hh_zlb(:);

g = G.g_gr(:);
s = G.s_gr(:);
mp = G.mp_gr(:);
in = G.in_gr(:);

% Linear regressions using regress
[b_Gust_lin,~,~,~,stats_Gust_lin] = regress(y_Gust, [ones(size(G.g_gr(:))) g s mp in]);
[b_Gust_zlb_lin,~,~,~,stats_Gust_zlb_lin] = regress(y_Gust, [ones(size(G.g_gr(:))) g s mp in]);
[b_ART_lin,~,~,~,stats_ART_lin] = regress(y_ART, [ones(size(G.g_gr(:))) g s mp in]);
% Linear regressions with quadratic terms using regress
[b_Gust,~,~,~,stats_Gust] = regress(y_Gust, [ones(size(G.g_gr(:))) g s mp in s.^2 in.^2]);
[b_Gust_zlb,~,~,~,stats_Gust_zlb] = regress(y_Gust_zlb, [ones(size(G.g_gr(:))) g s mp in s.^2 in.^2]);
[b_ART,~,~,~,stats_ART] = regress(y_ART, [ones(size(G.g_gr(:))) g s mp in s.^2 in.^2]);

ds_Gust = dataset(y_Gust,g,s,mp,in);
ds_Gust_zlb = dataset(y_Gust_zlb,g,s,mp,in);
ds_ART = dataset(y_ART,g,s,mp,in);

% Linear regressions using fitlm
lm_Gust_lin = fitlm(ds_Gust,'y_Gust~g+s+mp+in');
lm_Gust_lin_zlb = fitlm(ds_Gust_zlb,'y_Gust_zlb~g+s+mp+in');
lm_ART_lin = fitlm(ds_ART,'y_ART~g+s+mp+in');
%Linear regressions with quadratic terms using fitlm
lm_Gust = fitlm(ds_Gust,'y_Gust~g+s+mp+in+s^2+in^2+mp^2');
lm_Gust_zlb = fitlm(ds_Gust_zlb,'y_Gust_zlb~g+s+mp+in+s^2+in^2+mp^2');
lm_ART = fitlm(ds_ART,'y_ART~g+s+mp+in+s^2+in^2+mp^2');

