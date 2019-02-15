% Creates table 4, as latex fragment table4.tex, which displays the median
% and 5-95% credible set of posterior estimates of the parameters for the
% specifications that have results for 640 datasets.
% Input to create the table is param_acc_640.mat, which is created by
% script_postdraws_conditional.m.
clear
clc

% Table options
nqs = {'0','24'};
spec = {'pwlinear-if-me0';'global-pf-me5'};
savename = 'table4.tex';

% Load options and tables
load('../options.mat','V','P','F')
load('tables.mat');
nparam = F.nparam;

%number format (# of decimals) for parameter estimates and true value
nformat = {'%.1f','%.3f','%.3f','%.3f','%.4f','%.4f','%.4f','%.3f','%.3f'};
nformat_true = {'%.0f','%.1f','%.1f','%.1f','%.4f','%.4f','%.4f','%.1f','%.1f'};

%Header of table
Tc = cell(5,6);
Tc{1,1} = '\begin{table}[!htb]\footnotesize';
Tc{2,1} = '\sisetup{group-digits=false}';
Tc{3,1} = '\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}*{6}{c}}\toprule';
Tc{5,1} = ' Ptr & Truth & PwLin-IF-$0\%$ & Global-PF-$5\%$ & PwLin-IF-$0\%$ & Global-PF-$5\%$ \\ \midrule';

%Loop through parameters         
for ipar = 1:nparam
    
    var = F.priors{ipar,1};
        
    %Create subtable for different ZLB condtionals with median and cred sets
    T1 = subtable(total,nqs{1},var,nformat{ipar},'med',spec);
    T2 = subtable(total,nqs{2},var,nformat{ipar},'med',spec);

    %Combine subtables
    Tc0 = [T1,T2];
    
    %Label rows in the first two columns, parameter name and true value
    rowlab = cell(2,2);
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
rowlab{1,1} = '\midrule $\ell$';
T1 = subtable(total,nqs{1},'Sum','%.1f','mdd',spec);
T2 = subtable(total,nqs{2},'Sum','%.1f','mdd',spec);
%   Combine subtables
Tc0 = [rowlab,T1,T2];

%Add & between cells and \\ at end or rows
Tc0(:,1:end-1) = strcat(Tc0(:,1:end-1),' & ');
Tc0(:,end) = strcat(Tc0(:,end),'\\ ');
Tc0(1:2:end,end) = strcat(Tc0(1:2:end,end),'[-4pt] ');

%Combine parameter rows
Tc = [Tc;Tc0];

%Table footer
TcEnd = cell(4,size(Tc0,2));
TcEnd{1,1} = '\bottomrule \end{tabular*}';
TcEnd{2,1} = '\caption{Posterior medians and $(5\%,95\%)$ credible sets of the estimated parameters and marginal data densities. An asterisk indicates the differences in the marginal data densities (values different from 0) are significant at a $10\%$ level.}';
TcEnd{3,1} = '\label{tab:estimates}';
TcEnd{4,1} = '\end{table}';

%Add row to header with sampel size
Tc{4,1} = ['  &       & \multicolumn{2}{c}{ZLB duration = ' nqs{1} 'Q}  & \multicolumn{2}{c}{ZLB duration = ' nqs{2} 'Q}\\'];

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


function Tc0 = subtable(total,nq,var,format,stat,spec)
    %Function to create a subtable of estimation results
    % Inputs:
    %    total:  accuracy results, matlab table
    %    nq:   ZLB duration in quarters
    %    var:    parameter name, string
    %    format: number format for results, string
    %    stat:   accurarcy statistic, string
    %    spec: specifications to report
    %
    % Output:
    %    Tc0: formated strings for table, cell array
    %    samp: number of datasets included in subsample, ingeger
    
    % Blank cell array to store output
    Tc0 = cell(2,size(spec,1));
    %Only keep observations that match condition and variable
    T = total(strcmp(total.nq, nq) & strcmp(total.var, var),:);
    %Combine method, filter and me
    T.spec = cellstr(strcat(T.method,'-',T.filter,'-',T.me));

    %loop through specifications
    for ispec = 1:size(spec,1)
        % Median and cred sets
        if strcmp(stat,'med')                           
            Tc0{1,ispec} = ['$',num2str(T.Mean(strcmp(T.spec,spec{ispec})),format),'$'];
            Tc0{2,ispec} = ['\scs$(',num2str(T.p05(strcmp(T.spec,spec{ispec})),format),',',num2str(T.p95(strcmp(T.spec,spec{ispec})),format),')$'];
        % Marginal data density
        elseif strcmp(stat,'mdd')
            Tc0{1,ispec} = ['$',num2str(T.mdd(strcmp(T.spec,spec{ispec})),format),'$'];
            Tc0{2,ispec} = ['\scs$(',num2str(T.mdd_05(strcmp(T.spec,spec{ispec})),format),',',num2str(T.mdd_95(strcmp(T.spec,spec{ispec})),format),')$'];
        end

    end   
end