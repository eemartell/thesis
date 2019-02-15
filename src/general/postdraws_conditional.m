function [draws_all,nzlb,NRMSE,MDD,reject0595,postvar,MSE,bias2,isimdone] = ...
    postdraws_conditional(V,P,F,nsims,folder,sims,isimdone,inq)
    % Calculates accuracy statistics for a single dataset
    % Inputs:
    %   V:        Structure of variable indexes
    %   P:        Structure of model parameters
    %   F:        Structure of estimation information (e.g., priors)
    %   nsims:    Number of datasets, integer
    %   folder:   File location of posterior estimate files, string
    %   sims:     Artificial data
    %   subset:   Subset of datasets to consider, nsims x 1 logical
    %   inq:      ZLB duration index in sims
    %
    % Outputs:
    %   draws_all: posterior draws, 1000 x nparam x nsims
    %   nzlb:      number of zlb events in the datasets, nsims x 1 integer
    %   NRMSE:     Normalized root mean squared error of the posterior
    %                draws, nsims x nparam+1
    %   MDD:       Marginal data density, nsims x 1
    %   reject0595: Share of datasets that truth is outside of credible
    %                 set, nparam+1 x 1
    %   postvar:   Posterior variance, nparam+1 x 1
    %   MSE:   Mean-squared error, nparam+1 x 1
    %   bias2:   Square of the bias of the mean, nparam+1 x 1

    % Finished datasets
    if isimdone == 0
        isimdone = ones(nsims,1);
    end
    draws_all = NaN(1000,F.nparam,nsims);
    nmodesearchdraws = 5000;

    % True parameters
    for j = 1:F.nparam
       eval(['thetatilde(j) = P.' F.priors{j,1} ';']);
    end

    %Loop through simulations
    itemp = 1;
    nzlb = zeros(nsims,1);
    for isim = 1:nsims
        % Check if posterior files exist
        filename = [folder,'zlb',num2str(inq),'/thin/' num2str(isim) '-MHdraws.txt'];
        if (exist(filename,'file') && isimdone(isim)==1)
            % Load draws
            draws = textread(filename,'%f');
            draws = reshape(draws,[numel(draws)/F.nparam,F.nparam]);
            draws_all(:,:,isim) = draws;
            % Load data
            data = sims(:,:,isim);
            % Load discarded draws and log likelihoods for MDD
            ndiscarddraws = textread([folder,num2str(inq),'/ms/',num2str(isim),'-ndiscarddraws.txt'],'%f'); 
            postlogliks = textread([folder,num2str(inq),'/thin/',num2str(isim),'-MHpostlogliks.txt'],'%f');        

            % Check ZLB condition
            nzlb(itemp,1) = sum(data(:,V.i)==1);
            % NRMSE for parameters
            for j = 1:F.nparam
                % Error
                D = draws(:,j)-thetatilde(j);
                % Root mean squared error, normalized by true value
                NRMSE(itemp,j) = sqrt(mean(D.^2))/thetatilde(j);
                % Check if true value outside of cred set
                reject0595(itemp,j) = double((thetatilde(j) < quantile(draws(:,j),0.05)) || (thetatilde(j) > quantile(draws(:,j),0.95)));
                % Posterior variance
                postvar(itemp,j) = var(draws(:,j));
                % mean squared error, normalized by true value
                MSE(itemp,j) = mean(D.^2);
                % Bias squared
                bias2(itemp,j) = (mean(draws(:,j))-thetatilde(j)).^2;
            end
            % Sum or average across parameters
            NRMSE(itemp,F.nparam+1) = sum(NRMSE(itemp,1:F.nparam));
            reject0595(itemp,F.nparam+1) = mean(reject0595(itemp,1:F.nparam));
            postvar(itemp,F.nparam+1) = NaN;
            MSE(itemp,F.nparam+1) = NaN;
            bias2(itemp,F.nparam+1) = NaN;

            % Adjust posterior density for truncated prior density
            ntotalmodearchdraws = ndiscarddraws + nmodesearchdraws;
            postlogliks = postlogliks + log(1 - ndiscarddraws/ntotalmodearchdraws);
            % John Geweke's Modified Harmonic Mean
            MDD(itemp) = mean(geweke_mhm(draws,postlogliks));

            itemp = itemp + 1;
        else
            isimdone(isim) = 0;
        end
    end