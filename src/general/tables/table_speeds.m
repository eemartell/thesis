%  This file creates a table with the median, credible sets of estimation speed.

clear
clc
close all

savename = 'table_speeds.tex';
%Specifications to include: path, method, filter, measurement error, table
%label
specs = {...
    '../../Estimation-linear/me0/','linear','kf','me0','Level Linear, Kalman Filter, ME 0$\%$ (1 core)';...
    '../../Estimation-pwlinear/','pwlinear','if','me0','Piecewise Linear, Inversion Filter, ME 0$\%$ (1 core)';...
    '../../Estimation-global/me5/','global','pf','me5','Global, Particle Filter, ME 5$\%$ (16 cores)';         
    };

npers = 120;
nsimsall = 100;
nsims = 50;
nZLBpers = 6;
nqs = [0,6,12,18,24,30];
nspecs = size(specs,1);
ndraws = 1000;

%Load in V from misspecified model
load('../../Results-nocap/options.mat','V','P','F')

speed = nan(nsims,nZLBpers,nspecs);

%Loop through specifications
for ifolder = 1:size(specs,1)
    disp(['Specification: ' specs{ifolder,1}])
    % Loop through ZLB durations
    filename = [specs{ifolder,1},'states/speed.txt'];
    speed(:,:,ifolder) = reshape(textread(filename,'%f'),[nsims,nZLBpers]);
end   
s = speed;
h = s*80000/(60^2);
%Median and credible sets across datasets (dim 2)
s_m = squeeze(mean(s,1));
s_p05 = squeeze(quantile(s,0.05,1));
s_p95 = squeeze(quantile(s,0.95,1));
h_m = squeeze(mean(h,1));
h_p05 = squeeze(quantile(h,0.05,1));
h_p95 = squeeze(quantile(h,0.95,1));

%number format (# of decimals) for parameter estimates and true value
nformat = {'%.3f','%.3f','%.1f'};

%Header of table
Tc = cell(5,7);
Tc{1,1} = '\begin{table}[!htb]\footnotesize';
Tc{2,1} = '\sisetup{group-digits=false}';
Tc{3,1} = '\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}*{7}{c}}\toprule';
Tc{5,1} = ['  & ', num2str(nqs(1)),...
        'Q &  ', num2str(nqs(2)), ...
        'Q &  ', num2str(nqs(3)), ...
        'Q &  ', num2str(nqs(4)), ...
        'Q &  ', num2str(nqs(5)), ...
        'Q &  ', num2str(nqs(6)), 'Q  \\'];

    %Loop through specifications
for ispec = 1:nspecs
    %Add specification label row
    Tc0 = cell(1,7);
    Tc0{1,1} = ['\midrule \multicolumn{7}{c}{',specs{ispec,5},'} \\ \midrule'];
    Tc = [Tc;Tc0];

    %Add mean, cred set and RMSE for each zlb bin
    Tc0 = cell(2,7);
    Tc0{1,1} = 'Seconds per draw';
    Tc0{3,1} = 'Hours per dataset';
    for iZLBper = 1:nZLBpers
        Tc0{1,iZLBper+1} = ['$',num2str(s_m(iZLBper,ispec),nformat{ispec}),'$'];
        Tc0{2,iZLBper+1} = ['\scs$(',num2str(s_p05(iZLBper,ispec),nformat{ispec}),',',num2str(s_p95(iZLBper,ispec),nformat{ispec}),')$'];
        Tc0{3,iZLBper+1} = ['$',num2str(h_m(iZLBper,ispec),nformat{ispec}),'$'];
        Tc0{4,iZLBper+1} = ['\scs$(',num2str(h_p05(iZLBper,ispec),nformat{ispec}),',',num2str(h_p95(iZLBper,ispec),nformat{ispec}),')$'];
    end

    %Add & between cells and \\ at end or rows
    Tc0(:,1:end-1) = strcat(Tc0(:,1:end-1),' & ');
    Tc0(:,end) = strcat(Tc0(:,end),'\\ ');
    Tc0(1,end) = strcat(Tc0(1,end),'[-4pt] ');
    Tc0(3,end) = strcat(Tc0(3,end),'[-4pt] ');
    Tc = [Tc;Tc0];
end

%Table footer
TcEnd = cell(4,size(Tc0,2));
TcEnd{1,1} = '\bottomrule \end{tabular*}';
TcEnd{2,1} = '\caption{Mean, $(5\%,95\%)$ quantiles of estimation times.}';
TcEnd{3,1} = '\label{tab:speed}';
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