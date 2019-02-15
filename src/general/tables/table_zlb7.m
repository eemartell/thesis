%  This file creates a table with the median, credible sets and RMSE of
%  mean posterior estimates for different ZLB bins and specifications.

clear
clc
close all

savename = 'table_zlb7.tex';
%Specifications to include: path, method, filter, measurement error, table
%label
specs = {...
    '../../Estimation-global/me5/','global','pf','me5','NL-PF-5$\%$'; ...   
    '../../Estimation-pwlinear/','pwlinear','if','me0','PW-IF-0$\%$';...   
    };

npers = 120;
nsimsall = 100;
nsims = 50;
nZLBpers = 7;
nqs = [0,6,12,18,24,30,30];
nspecs = size(specs,1);
ndraws = 1000;

%Load in V from misspecified model
load('../../Results-nocap/options.mat','V','P','F')

means = nan(F.nparam,nsims,nZLBpers,nspecs);

%Loop through specifications
for ifolder = 1:size(specs,1)
    disp(['Specification: ' specs{ifolder,1}])
    % Loop through ZLB durations
    for inq = [1:numel(nqs)]
        % ZLB duration
        nq = nqs(inq);
        disp(['ZLB duration = ' num2str(nq)])
        for isim = 1:nsims
            %Load posterior distributions they exist, otherwise keep as
            %NaNs
            filename = [specs{ifolder,1},'zlb',num2str(inq),'/mean/',num2str(isim),'-postmean.txt'];
            if exist(filename, 'file') == 2
                %Calculate mean across posterior distribution
                means(:,isim,inq,ifolder) = textread(filename,'%f');
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
nformat = {'%.1f','%.2f','%.2f','%.2f','%.4f','%.4f','%.4f','%.2f','%.2f'};
nformat_true = {'%.0f','%.1f','%.1f','%.1f','%.3f','%.3f','%.3f','%.1f','%.1f'};

%Header of table
Tc = cell(5,8);
Tc{1,1} = '\begin{table}[!htb]\footnotesize';
Tc{2,1} = '\sisetup{group-digits=false}';
Tc{3,1} = '\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}l*{7}{c}}\toprule';
Tc{4,1} = [' Ptr & Truth & \multicolumn{3}{c}{', specs{1,5},...
        '} &  \multicolumn{3}{c}{', specs{2,5},'}  \\ \midrule'];
Tc{5,1} = ['  & & ', num2str(nqs(1)), 'Q &  ', num2str(nqs(6)), 'Q & 5x6Q & ', ...
            num2str(nqs(1)), 'Q &  ', num2str(nqs(6)), 'Q & 5x6Q  \\ \midrule'];

%Loop through parameters         
for ipar = 1:F.nparam

    var = F.priors{ipar,1};

    %Label rows in the first two columns, parameter name and true value
    Tc0 = cell(3,2);
    Tc0{1,1} = ['$',F.paramtex{ipar},'$'];
    Tc0{1,2} = ['$',num2str(F.priors{ipar,3},nformat_true{ipar}),'$'];

    %Loop through specifications
    for ispec = 1:nspecs
        %Add mean, cred set and RMSE for each zlb bin
        Tc1 = cell(3,2);
        itab=1;
        for iZLBper = [1,nZLBpers-1,nZLBpers]
            Tc1{1,itab} = ['$',num2str(p50(ipar,iZLBper,ispec),nformat{ipar}),'$'];
            Tc1{2,itab} = ['\scs$(',num2str(p05(ipar,iZLBper,ispec),nformat{ipar}),',',num2str(p95(ipar,iZLBper,ispec),nformat{ipar}),')$'];
            Tc1{3,itab} = ['\scs$[',num2str(nrmse(ipar,iZLBper,ispec),'%.2f'),']$'];
            itab=itab+1;
        end
        Tc0 = [Tc0,Tc1];
    end

    %Add & between cells and \\ at end or rows
    Tc0(:,1:end-1) = strcat(Tc0(:,1:end-1),' & ');
    Tc0(:,end) = strcat(Tc0(:,end),'\\ ');
    Tc0(1:2,end) = strcat(Tc0(1:2,end),'[-4pt] ');

    %Apend this parameters row to the full table
    Tc = [Tc;Tc0];
end
 % Sum of NRMSE and number of datasets
 Tc0 = cell(1,2);
 Tc0{1,1} = '\midrule $\Sigma$';
 %Loop through specifications
 for ispec = 1:nspecs
     Tc1 = cell(1,2);
     itab = 1;
     for iZLBper = [1,nZLBpers-1,nZLBpers]
         Tc1{1,itab} = ['\scs$[',num2str(sigma_nrmse(iZLBper,ispec),'%.2f'),']$'];
         itab = itab+1;
     end
     Tc0 = [Tc0,Tc1];
 end

% %Add & between cells and \\ at end or rows
 Tc0(:,1:end-1) = strcat(Tc0(:,1:end-1),' & ');
 Tc0(:,end) = strcat(Tc0(:,end),'\\ ');
% 
% %Combine parameter rows
 Tc = [Tc;Tc0];

%Table footer
TcEnd = cell(4,size(Tc0,2));
TcEnd{1,1} = '\bottomrule \end{tabular*}';
TcEnd{2,1} = '\caption{Median, $(5\%,95\%)$ credible sets and $[\operatorname{NRMSE}]$ of the posterior mean parameter estimates.}';
TcEnd{3,1} = '\label{tab:Mestimates}';
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