%  This file creates a figure with the RMSE of filtered state estimates. 
%  The output is tbd

clear
clc
close all

savename = 'table4b.tex';
%Specifications to include: path, method, filter, measurement error
specs = {...
    '../Estimation-pwlinear/','pwlinear','if','me0','Piecewise Linear, Inversion Filter, ME 0$\%$';...
    '../Estimation-global/me5/','global','pf','me5','Global, Particle Filter, ME 5$\%$';         
    };

% First load in the simulated date from the true model
npers = 120;
nsimsall = 100;
nsims = 50;
nZLBpers = 6;
nqs = [0,6,12,18,24,30];
nspecs = size(specs,1);
ndraws = 1000;

% Load parameter names - true model
load('options.mat','V','P','F')
%Rename V to indicate it is from true model
Vcap= V;
Pcap = P;
Fcap = F;

%Load in V from misspecified model
load('../Results-nocap/options.mat','V','P','F')

means = nan(F.nparam,nsims,nZLBpers,nspecs);

%Loop through specifications
for ifolder = 1:size(specs,1)
    disp(['Specification: ' specs{ifolder,1}])
    % Loop through ZLB durations
    for inq = 1:numel(nqs)
        % ZLB duration
        nq = nqs(inq);
        disp(['ZLB duration = ' num2str(nq)])
        for isim = 1:nsims
            %Load filtered state estimates if they exist, otherwise keep as
            %NaNs
            filename = [specs{ifolder,1},'zlb',num2str(inq),'/thin/',num2str(isim),'-MHdraws.txt'];
            if exist(filename, 'file') == 2
                draws = textread(filename,'%f');
                means(:,isim,inq,ifolder) = mean(reshape(draws,[ndraws,F.nparam]),1);
            end
        end
    end
end   

p50 = squeeze(quantile(means,0.5,2));
p05 = squeeze(quantile(means,0.05,2));
p95 = squeeze(quantile(means,0.95,2));
actual = repmat(cell2mat(F.priors(:,3)),1,nsims,nZLBpers,nspecs);
rmse = squeeze(mean((means-actual).^2,2,'omitnan').^0.5);
actual = repmat(cell2mat(F.priors(:,3)),1,nZLBpers,nspecs);
nrmse = rmse./actual;
nrmse = squeeze(sum(nrmse,1));
N = nsims - squeeze(sum(isnan(means(1,:,:,:)),2));

%number format (# of decimals) for parameter estimates and true value
nformat = {'%.1f','%.3f','%.3f','%.3f','%.4f','%.4f','%.4f','%.3f','%.3f'};
nformat_true = {'%.0f','%.1f','%.1f','%.1f','%.4f','%.4f','%.4f','%.1f','%.1f'};

%Header of table
Tc = cell(5,8);
Tc{1,1} = '\begin{table}[!htb]\footnotesize';
Tc{2,1} = '\sisetup{group-digits=false}';
Tc{3,1} = '\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}*{8}{c}}\toprule';
Tc{5,1} = [' Ptr & Truth & ZLB duration = ', nqs(1), 'Q & ZLB duration = ', nqs(2), ...
        'Q & ZLB duration = ', nqs(3), 'Q & ZLB duration = ', nqs(4),...
        'Q & ZLB duration = ', nqs(5), 'Q & ZLB duration = ', nqs(6) 'Q  \\ \midrule'];

for ispec = 1:nspecs
    Tc0 = cell(1,8);
    Tc0{1,1} = ['\midrule \multicolumn{8}{c}{',specs{ispec,5},'} \\ \midrule'];
    %Combine parameter rows
    Tc = [Tc;Tc0];
    %Loop through parameters         
    for ipar = 1:F.nparam

        var = F.priors{ipar,1};

        Tc0 = cell(3,nZLBpers);
        for iZLBper = 1:nZLBpers
            Tc0{1,iZLBper} = ['$',num2str(p50(ipar,iZLBper,ispec),nformat{ipar}),'$'];
            Tc0{2,iZLBper} = ['\scs$(',num2str(p05(ipar,iZLBper,ispec),nformat{ipar}),',',num2str(p95(ipar,iZLBper,ispec),nformat{ipar}),')$'];
            Tc0{3,iZLBper} = ['$\[',num2str(rmse(ipar,iZLBper,ispec),nformat{ipar}),'\]$'];
        end

        %Label rows in the first two columns, parameter name and true value
        rowlab = cell(3,2);
        rowlab{1,1} = ['$',F.paramtex{ipar},'$'];
        rowlab{1,2} = ['$',num2str(F.priors{ipar,3},nformat_true{ipar}),'$'];
        Tc0 = [rowlab,Tc0];

        %Add & between cells and \\ at end or rows
        Tc0(:,1:end-1) = strcat(Tc0(:,1:end-1),' & ');
        Tc0(:,end) = strcat(Tc0(:,end),'\\ ');
        Tc0(1,end) = strcat(Tc0(1,end),'[-4pt] ');

        %Apend this parameters row to the full table
        Tc = [Tc;Tc0];
    end
    % MDD
    % Create subtable for two ZLB condtionals with median MDD
    %   Label rows in the first two columns, parameter name and true value
    rowlab = cell(2,2);
    rowlab{1,1} = ' $\Sigma$';
    rowlab{2,1} = ' $N$';   
    Tc0 = cell(2,nZLBpers);
    for iZLBper = 1:nZLBpers
        Tc0{1,iZLBper} = ['$\[',num2str(nrmse(iZLBper,ispec),'%.3f'),'\]$'];
        Tc0{2,iZLBper} = ['$',num2str(N(iZLBper,ispec)),'$'];
    end
    %   Combine subtables
    Tc0 = [rowlab,Tc0];

    %Add & between cells and \\ at end or rows
    Tc0(:,1:end-1) = strcat(Tc0(:,1:end-1),' & ');
    Tc0(:,end) = strcat(Tc0(:,end),'\\ ');
    Tc0(1:3:end,end) = strcat(Tc0(1:3:end,end),'[-4pt] ');

    %Combine parameter rows
    Tc = [Tc;Tc0];
end

%Table footer
TcEnd = cell(4,size(Tc0,2));
TcEnd{1,1} = '\bottomrule \end{tabular*}';
TcEnd{2,1} = '\caption{Median, $(5\%,95\%)$ credible sets and $\[RMSE\]$ of the mean posterior estimated parameters for $N$ datasets. $\Sigma$ is the sum of the normalized RMSE.}';
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