%  This file creates a table with the median, credible sets and RMSE of
%  mean posterior estimates for different ZLB bins and specifications.

clear
clc
close all

savename = 'table_param_sum_neg_in.tex';

%Specifications to include: path, method, filter, measurement error, table
%label
specs = {...
    '../../Estimation-global/me5/','global','pf','me5','Global, Particle Filter, ME 5$\%$';...         
    '../../Estimation-pwlinear/','pwlinear','if','me0','Piecewise Linear, Inversion Filter, ME 0$\%$';...
    };

npers = 120;
nsimsall = 100;
nsims = 50;
nZLBpers = 6;
nqs = [0,6,12,18,24,30];
nspecs = size(specs,1);
ndraws = 1000;

% Load parameter names
load('../options.mat','V','P','F')

% Load actual data
simscap = textread('../../data-artificial/sims.txt','%f');
simscap = reshape(simscap,[npers,V.nvar-V.nfore,nsimsall,nZLBpers]);

Vcap= V;

%Load in V from misspecified model
load('../../Results-nocap/options.mat','V','P','F')
Vnstate = V.nvar-V.nfore

%Create artificial data matrix (true model) that aligns with variable order
%from the true model. Note that yg is named cg in misspec model
sims = nan(npers,Vnstate,nsims,nZLBpers);
for ivar = 1:(Vnstate)
    if strcmp(V.names{ivar},'cg')
        sims(:,ivar,:,:) = simscap(:,Vcap.yg,1:nsims,:);
    else
        eval(['sims(:,ivar,:,:) = simscap(:,Vcap.',V.names{ivar},',1:nsims,:);']);
    end
end

in = squeeze(sims(:,V.in,:,:));
neg_in = in;
neg_in(in>1) = NaN;
sumin = squeeze(sum(1-neg_in,1,'omitnan'));

Qsumin = reshape(discretize(sumin(:),[min(sumin(:)),quantile(sumin(:),5),max(sumin(:))]),nsims,nZLBpers);

means = nan(F.nparam,nsims,nZLBpers,nspecs);

%Loop through specifications
for ifolder = 1:size(specs,1)
    simnum = [1,1,1,1,1,1];
    disp(['Specification: ' specs{ifolder,1}])
    % Loop through ZLB durations
    for inq = 1:numel(nqs)
        % ZLB duration
        nq = nqs(inq);
        disp(['ZLB duration = ' num2str(nq)])
        for isim = 1:nsims
            %Load posterior distributions they exist, otherwise keep as
            %NaNs
            newbin = Qsumin(isim,inq);
            newsim = simnum(newbin);
            simnum(newbin) = simnum(newbin) + 1;
            filename = [specs{ifolder,1},'zlb',num2str(inq),'/mean/',num2str(isim),'-postmean.txt'];
            if exist(filename, 'file') == 2
                %Calculate mean across posterior distribution
                means(:,newsim,newbin,ifolder) = textread(filename,'%f');
            end
        end
    end
end   

%Median and credible sets across datasets (dim 2)
p50 = squeeze(quantile(means,0.5,2));
p05 = squeeze(quantile(means,0.05,2));
p95 = squeeze(quantile(means,0.95,2));
%RMSE across datasets (dim 2), ignore nans
actual = repmat(cell2mat(F.priors(:,3)),1,nsims,nZLBpers,nspecs);
rmse = squeeze(mean((means-actual).^2,2,'omitnan').^0.5);
actual = repmat(cell2mat(F.priors(:,3)),1,nZLBpers,nspecs);
%Sum of normalized RMSE
nrmse = rmse./actual;
sigma_nrmse = squeeze(sum(nrmse,1));
%Number of nonmissing datasets
N = nsims - squeeze(sum(isnan(means(1,:,:,:)),2));

%number format (# of decimals) for parameter estimates and true value
nformat = {'%.1f','%.3f','%.3f','%.3f','%.4f','%.4f','%.4f','%.3f','%.3f'};
nformat_true = {'%.0f','%.1f','%.1f','%.1f','%.4f','%.4f','%.4f','%.1f','%.1f'};

%Header of table
Tc = cell(5,8);
Tc{1,1} = '\begin{table}[!htb]\footnotesize';
Tc{2,1} = '\sisetup{group-digits=false}';
Tc{3,1} = '\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}*{8}{c}}\toprule';
Tc{5,1} = [' Ptr & Truth & Bin 1',...
        ' & Bin 2  ', ...
        ' & Bin 3 ', ...
        ' & Bin 4 ', ...
        ' & Bin 5 ', ...
        ' & Bin 6 \\'];

    %Loop through specifications
for ispec = 1:nspecs
    %Add specification label row
    Tc0 = cell(1,8);
    Tc0{1,1} = ['\midrule \multicolumn{8}{c}{',specs{ispec,5},'} \\ \midrule'];
    Tc = [Tc;Tc0];
    %Loop through parameters         
    for ipar = 1:F.nparam

        var = F.priors{ipar,1};

        %Add mean, cred set and RMSE for each zlb bin
        Tc0 = cell(3,nZLBpers);
        for iZLBper = 1:nZLBpers
            Tc0{1,iZLBper} = ['$',num2str(p50(ipar,iZLBper,ispec),nformat{ipar}),'$'];
            Tc0{2,iZLBper} = ['\scs$(',num2str(p05(ipar,iZLBper,ispec),nformat{ipar}),',',num2str(p95(ipar,iZLBper,ispec),nformat{ipar}),')$'];
            Tc0{3,iZLBper} = ['\scs$[',num2str(nrmse(ipar,iZLBper,ispec),'%.3f'),']$'];
        end

        %Label rows in the first two columns, parameter name and true value
        rowlab = cell(3,2);
        rowlab{1,1} = ['$',F.paramtex{ipar},'$'];
        rowlab{1,2} = ['$',num2str(F.priors{ipar,3},nformat_true{ipar}),'$'];
        Tc0 = [rowlab,Tc0];

        %Add & between cells and \\ at end or rows
        Tc0(:,1:end-1) = strcat(Tc0(:,1:end-1),' & ');
        Tc0(:,end) = strcat(Tc0(:,end),'\\ ');
        Tc0(1:2,end) = strcat(Tc0(1:2,end),'[-4pt] ');

        %Apend this parameters row to the full table
        Tc = [Tc;Tc0];
    end
    % Sum of NRMSE and number of datasets
    rowlab = cell(1,2);
    rowlab{1,1} = ' $\Sigma$';
    Tc0 = cell(1,nZLBpers);
    for iZLBper = 1:nZLBpers
        Tc0{1,iZLBper} = ['\scs$[',num2str(sigma_nrmse(iZLBper,ispec),'%.3f'),']$'];
    end
    Tc0 = [rowlab,Tc0];

    %Add & between cells and \\ at end or rows
    Tc0(:,1:end-1) = strcat(Tc0(:,1:end-1),' & ');
    Tc0(:,end) = strcat(Tc0(:,end),'\\ ');
    
    %Combine parameter rows
    Tc = [Tc;Tc0];
end

%Table footer
TcEnd = cell(4,size(Tc0,2));
TcEnd{1,1} = '\bottomrule \end{tabular*}';
TcEnd{2,1} = '\caption{Median, $(5\%,95\%)$ credible sets and $[NRMSE]$ of the mean posterior estimated parameters. Column number refers to quantile of the summed notional rate below 0.}';
TcEnd{3,1} = '\label{tab:estimates}';
TcEnd{4,1} = '\end{table}';

%add footnote
Tc = [Tc;TcEnd];

%Save table as tex file 
fid = fopen(savename, 'w');
[nrows,ncols] = size(Tc);
format = [repmat('%s ',1,ncols),' \n'];
for row = 1:nrows
    fprintf(fid,format,Tc{row,:});
end
fclose(fid);

% p50 = squeeze(quantile(means,0.5,2));
% figure
% scatter(sum_neg_in(:,1),means(9,:,1,1))
% hold on
% scatter(sum_neg_in(:,2),means(9,:,2,1))
% x = [min(sum_neg_in(:,2)),max(sum_neg_in(:,2))];
% y = [p50(9,1,1),p50(9,2,1)];
% plot(x,y)
% scatter(sum_neg_in(:,3),means(9,:,3,1))
% scatter(sum_neg_in(:,4),means(9,:,4,1))
% scatter(sum_neg_in(:,5),means(9,:,5,1))
% scatter(sum_neg_in(:,6),means(9,:,6,1))
% set(gcf,'Position',[1 1 6.5 3]);
% 
% x = [min(sum_neg_in(:,1)),max(sum_neg_in(:,1))];
% y = [p50(9,1,1),p50(9,1,1)];
% plot(x,y) 