%  This file creates a mat file with summary statistics of the accuracy
%  parameter estimates of different specifications.
%  The output is tables/tables.mat

clear
clc
close all

%Baseline specification
% baseline = 'lin-kf-me0'; 
%  baseline = 'pw-if-me0';
 
%save the output matlab table
save_flag = 1;

% Load parameter names
load('options.mat','V','P','F')

% Load artificial data
npers = 120;
nsimsall = 100;
nZLBpers = 6;
% nqs = [0,6,12,18,24,30];
nqs = [0,30];
simstemp = textread('../data-artificial/sims.txt','%f');
sims = reshape(simstemp,[npers,V.nvar-V.nfore,nsimsall,nZLBpers]);

%Specifications to include: path, method, filter, measurement error
nsims = 50;
save_name = 'tables/tables.mat';
specs = {...
    '../Estimation-pwlinear/','pwlinear','if','me0';...
    '../Estimation-global/','global','pf','me5';         
    };

% Table headers and probabilities
headerstrs = {...
    'Mean','p005','p025','p05','p95','p975','p995',...
    'NRMSE','NRMSE_005','NRMSE_025','NRMSE_05','NRMSE_95','NRMSE_975','NRMSE_995',...
    'mdd','mdd_005','mdd_025','mdd_05','mdd_95','mdd_975','mdd_995',...
    'reject0595',...
    'postvar','postvar_005','postvar_025','postvar_05','postvar_95','postvar_975','postvar_995',...
    };
p = [.5,0.005,0.025,0.05,0.95,0.975,0.995];
np = numel(p);

total =[];
% Loop through ZLB durations
for inq = 1:numel(nqs)
    % ZLB duration
    nq = nqs(inq);
    disp(['ZLB duration = ' num2str(nq)])
    %Loop through specifications
    for ifolder = 1:size(specs,1)
        disp(['Specification: ' specs{ifolder,1}])
        %Get accuracy stats
        [draws_all,nzlb,NRMSE,MDD,reject0595,postvar,MSE,bias2,isimdone] = ...
            postdraws_conditional(...
                V,P,F,nsims,specs{ifolder,1},sims(:,:,:,inq),0,inq);
        % Marginal Data Densities
        MDDs{ifolder} = MDD;
        %Root-mean qquare error of posteriors
        RMSEs{ifolder} = NRMSE;
        %Share of datasets where truth is outside of credible set
        reject0595s{ifolder} = reject0595;
        %Posterior variances
        postvars{ifolder} = postvar;
        %MSE
        MSEs{ifolder} = MSE;
        %Bias squared
        bias2s{ifolder} = bias2;

        % Quantiles
        means = NaN(F.nparam+1,1);
        post_quant = NaN(np-1,F.nparam+1);
        % Quantiles of pooled (across datasets) quantiles for each
        % parameter
        for j = 1:F.nparam
            meandraws = squeeze(mean(draws_all(:,j,isimdone==1),1));
            means(j) = mean(meandraws);
            post_quant(:,j)=quantile(meandraws,p(2:end));
        end
        %Quantiles across datasets 
        NRMSE_quant = quantile(NRMSE,p,1);
        MDD_quant = repmat(quantile(MDD',p,1),[1 size(NRMSE,2)]);
        meanreject0595 = mean(reject0595,1);
        postvar_quant = quantile(postvar,p,1);

        % Create table
        temp = table(...
            means,post_quant(1,:)',post_quant(2,:)',post_quant(3,:)',post_quant(4,:)',post_quant(5,:)',post_quant(6,:)',...
            NRMSE_quant(1,:)',NRMSE_quant(2,:)',NRMSE_quant(3,:)',NRMSE_quant(4,:)',NRMSE_quant(5,:)',NRMSE_quant(6,:)',NRMSE_quant(7,:)',...
            MDD_quant(1,:)',MDD_quant(2,:)',MDD_quant(3,:)',MDD_quant(4,:)',MDD_quant(5,:)',MDD_quant(6,:)',MDD_quant(7,:)',...
            meanreject0595', ...
            postvar_quant(1,:)',postvar_quant(2,:)',postvar_quant(3,:)',postvar_quant(4,:)',postvar_quant(5,:)',postvar_quant(6,:)',postvar_quant(7,:)');
        
        %Name columns
        temp.Properties.VariableNames = headerstrs;
        temp.var = [F.priors(:,1);'Sum'];
        
        % Add specification labels to table
        temp.method(:,1) = string(specs(ifolder,2));
        temp.filter(:,1) = string(specs(ifolder,3));
        temp.me(:,1) = string(specs(ifolder,4)); 
        temp.nq(:,1) = string(nq); 
        % Combine
        total = [total;temp];
    end
end   

if save_flag 
    save(save_name,'total')
end