% Smoothness tests
clear
clc
close all

% Load solution at posterior mean
load('solutionfpART.mat')
pf_ART.c = pf.hh;
load('solutionfpGust.mat')
pf_Gust.c = pf.hh;
pf_Gust.c_zlb = pf.hh_zlb;

% Regression models with slice of data
% 
% gslice = ceil(O.g_pts/2);
% mpslice = ceil(O.mp_pts/2);
% inslice = 3; %in = 1.0001
% 
% x = 100*(G.s_grid-1)';
% y_Gust_sl = 100*(pf_Gust.c(gslice,:,mpslice,inslice)'-S.c)/S.c;
% y_ART_sl = 100*(pf_ART.c(gslice,:,mpslice,inslice)'-S.c)/S.c;
% 
% [b_Gust_sl,~,~,~,stats_Gust_sl] = regress(y_Gust_sl,[x x.^2 ones(size(x))]);
% [b_ART_sl,~,~,~,stats_ART_sl] = regress(y_ART_sl,[x x.^2 ones(size(x))]);
% 
% ds_Gust_sl = dataset(y_Gust_sl,x);
% lm_Gust_sl = fitlm(ds_Gust_sl,'y_Gust_sl~x+x^2');
% ds_ART_sl = dataset(y_ART_sl,x);
% lm_ART_sl = fitlm(ds_ART_sl,'y_ART_sl~x+x^2');

%Regression models with full data
g = G.g_gr(:);
s = G.s_gr(:);
mp = G.mp_gr(:);
in = G.in_gr(:);
y_Gust = pf_Gust.c(:);
y_ART = pf_ART.c(:);

[b_Gust_lin,~,~,~,stats_Gust_lin] = regress(y_Gust, [ones(size(G.g_gr(:))) g s mp in]);
[b_ART_lin,~,~,~,stats_ART_lin] = regress(y_ART, [ones(size(G.g_gr(:))) g s mp in]);

[b_Gust,~,~,~,stats_Gust] = regress(y_Gust, [ones(size(G.g_gr(:))) g s mp in s.^2 in.^2]);
[b_ART,~,~,~,stats_ART] = regress(y_ART, [ones(size(G.g_gr(:))) g s mp in s.^2 in.^2]);

ds_Gust = dataset(y_Gust,g,s,mp,in);
lm_Gust_lin = fitlm(ds_Gust,'y_Gust~g+s+mp+in');
lm_Gust = fitlm(ds_Gust,'y_Gust~g+s+mp+in+s^2+in^2');
ds_ART = dataset(y_ART,g,s,mp,in);
lm_ART_lin = fitlm(ds_ART,'y_ART~g+s+mp+in');
lm_ART = fitlm(ds_ART,'y_ART~g+s+mp+in+s^2+in^2');