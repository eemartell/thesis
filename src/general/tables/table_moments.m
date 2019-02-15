%  This file creates a table with the median, credible sets and RMSE of
%  mean posterior estimates for different ZLB bins and specifications.

clear
clc
close all

savename = 'table_moments.tex';
%Specifications to include: path, method, filter, measurement error, table
%label
specs = {...
    '../../Estimation-global/me5/','NL_PF_5','pf','me5','NL-PF-5$\%$';  
    '../../Estimation-pwlinear/','PW_IF_0','if','me0','PW-IF-0$\%$';...  
    '../../Estimation-linear/me0/','Lin_KF_0','kf','me2','Lin-KF-0$\%$';... 
    };

npers = 120;
nsimsall = 100;
nsims = 50;
nZLBpers = 6;
nqs = [0,6,12,18,24,30];
nspecs = size(specs,1);
nmoments = 4;

% Load parameter names
load('../options.mat','V','P','F')

% Load actual data
momentscap = textread('../../data-artificial/Ezlb_duration/moments.txt','%f');
momentscap = reshape(momentscap,[nmoments,V.nvar-V.nfore]);

Vcap= V;

%Load in V from misspecified model
load('../../Results-nocap/options.mat','V','P','F')
Vnstate = V.nvar-V.nfore

%Create artificial data matrix (true model) that aligns with variable order
%from the true model. Note that yg is named cg in misspec model
moments = nan(nmoments,Vnstate);
for ivar = 1:(Vnstate)
    if strcmp(V.names{ivar},'cg')
        moments(:,ivar) = momentscap(:,Vcap.yg);
    else
        eval(['moments(:,ivar) = momentscap(:,Vcap.',V.names{ivar},');']);
    end
end

moments_est = nan(nmoments,Vnstate,nsims,nZLBpers,nspecs);
%Loop through specifications
for ifolder = 1:size(specs,1)
    disp(['Specification: ' specs{ifolder,1}])
    %Load state estimate for given specification
    filename = [specs{ifolder,1},'/states/moments.txt'];
    if exist(filename, 'file') == 2
        temp = textread(filename,'%f');
        temp = reshape(temp,[nmoments,Vnstate,nsims,nZLBpers]);
        moments_est(:,:,:,:,ifolder) = temp;
    end
end   

%Median and credible sets across datasets (dim 3)
p50 = squeeze(quantile(moments_est,0.5,3));
p05 = squeeze(quantile(moments_est,0.05,3));
p95 = squeeze(quantile(moments_est,0.95,3));

actual = repmat(moments,1,1,nsims,nZLBpers,nspecs);
rmse = squeeze(mean((moments_est-actual).^2,3,'omitnan').^0.5);

moment_labels = {'mean($y^g$)',1,V.yg;...
              'mean($\pi$)',1,V.pi;...
              'mean($i$)',1,V.i;...
              'std($y^g$)',2,V.yg;...
              'std($\pi$)',2,V.pi;...
              'std($i$)',2,V.i;...
              'skew($y^g$)',3,V.yg;...
              'skew($\pi$)',3,V.pi;...
              'skew($i$)',3,V.i;...
              'AC($y^g$)',4,V.yg;...
              'AC($\pi$)',4,V.pi;...
              'AC($i$)',4,V.i;};


%number format (# of decimals) for parameter estimates and true value
%nformat = {'%.1f','%.3f','%.3f','%.3f','%.4f','%.4f','%.4f','%.3f','%.3f'};
%nformat_true = {'%.0f','%.1f','%.1f','%.1f','%.4f','%.4f','%.4f','%.1f','%.1f'};
nformat = '%.3f';



%Header of table
Tc = cell(5,8);
Tc{1,1} = '\begin{table}[!htb]\footnotesize';
Tc{2,1} = '\sisetup{group-digits=false}';
Tc{3,1} = '\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}*{8}{c}}\toprule';
Tc{4,1} = [' Ptr & Truth & \multicolumn{2}{c}{', specs{1,5},...
        '} &  \multicolumn{2}{c}{', specs{2,5},...
        '} &  \multicolumn{2}{c}{', specs{3,5},'}  \\ \midrule'];
Tc{5,1} = ['  & & ', num2str(nqs(1)),...
        'Q &  ', num2str(nqs(6)), ...
        'Q &  ', num2str(nqs(1)), ...
        'Q &  ', num2str(nqs(6)), ...
        'Q &  ', num2str(nqs(1)), ...
        'Q &  ', num2str(nqs(6)), 'Q  \\ \midrule'];

%Loop through parameters   
for imoment = 1:size(moment_labels,1)

    ivar = moment_labels{imoment,3};
    imoment2 = moment_labels{imoment,2};

    %Label rows in the first two columns, parameter name and true value
    Tc0 = cell(3,2);
    Tc0{1,1} = moment_labels{imoment,1};
    Tc0{1,2} = ['$',num2str(moments(imoment2,ivar),nformat),'$'];

    %Loop through specifications
    for ispec = 1:nspecs
        %Add mean, cred set and RMSE for each zlb bin
        Tc1 = cell(3,2);
        itab=1;
        for iZLBper = [1,nZLBpers]
            Tc1{1,itab} = ['$',num2str(p50(imoment2,ivar,iZLBper,ispec),nformat),'$'];
            Tc1{2,itab} = ['\scs$(',num2str(p05(imoment2,ivar,iZLBper,ispec),nformat),',',num2str(p95(imoment2,ivar,iZLBper,ispec),nformat),')$'];
            Tc1{3,itab} = ['\scs$[',num2str(rmse(imoment2,ivar,iZLBper,ispec),'%.3f'),']$'];
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
%Tc0 = cell(1,2);
%Tc0{1,1} = ' $\Sigma$';
%Tc0{2,1} = ' $N$';   
%Loop through specifications
% for ispec = 1:nspecs
%     Tc1 = cell(1,2);
%     itab = 1;
%     for iZLBper = [1,nZLBpers]
%         Tc1{1,itab} = ['\scs$[',num2str(sigma_nrmse(iZLBper,ispec),'%.3f'),']$'];
%         %Tc1{2,itab} = ['\scs$',num2str(N(iZLBper,ispec)),'$'];
%         itab = itab+1;
%     end
%     Tc0 = [Tc0,Tc1];
% end

%Add & between cells and \\ at end or rows
%Tc0(:,1:end-1) = strcat(Tc0(:,1:end-1),' & ');
%Tc0(:,end) = strcat(Tc0(:,end),'\\ ');

%Combine parameter rows
%Tc = [Tc;Tc0];

%Table footer
TcEnd = cell(4,size(Tc0,2));
TcEnd{1,1} = '\bottomrule \end{tabular*}';
    TcEnd{2,1} = '\caption{Median, $(5\%,95\%)$ credible sets and $[NRMSE]$ of moments.}';
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
%   
% %Loop through specifications
% for ispec = 1:nspecs
% %Header of table
%     Tc = cell(5,8);
%     Tc{1,1} = '\begin{table}[!htb]\footnotesize';
%     Tc{2,1} = '\sisetup{group-digits=false}';
%     Tc{3,1} = '\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}*{8}{c}}\toprule';
%     Tc{5,1} = [' Ptr & Truth & ', num2str(nqs(1)),...
%             'Q &  ', num2str(nqs(2)), ...
%             'Q &  ', num2str(nqs(3)), ...
%             'Q &  ', num2str(nqs(4)), ...
%             'Q &  ', num2str(nqs(5)), ...
%             'Q &  ', num2str(nqs(6)), 'Q  \\'];
% 
%     %Add specification label row
%     Tc0 = cell(1,8);
%     Tc0{1,1} = ['\midrule \multicolumn{8}{c}{',specs{ispec,5},'} \\ \midrule'];
%     Tc = [Tc;Tc0];
%     %Loop through parameters         
%     for imoment = 1:size(moment_labels,1)
% 
%         ivar = moment_labels{imoment,3};
%         imoment2 = moment_labels{imoment,2};
% 
%         %Add mean, cred set and RMSE for each zlb bin
%         Tc0 = cell(3,nZLBpers);
%         for iZLBper = 1:nZLBpers
%             Tc0{1,iZLBper} = ['$',num2str(p50(imoment2,ivar,iZLBper,ispec),nformat),'$'];
%             Tc0{2,iZLBper} = ['\scs$(',num2str(p05(imoment2,ivar,iZLBper,ispec),nformat),',',num2str(p95(imoment2,ivar,iZLBper,ispec),nformat),')$'];
%             Tc0{3,iZLBper} = ['\scs$[',num2str(rmse(imoment2,ivar,iZLBper,ispec),'%.3f'),']$'];
%         end
% 
%         %Label rows in the first two columns, parameter name and true value
%         rowlab = cell(3,2);
%         rowlab{1,1} = moment_labels{imoment,1};
%         rowlab{1,2} = ['$',num2str(moments(imoment2,ivar),nformat),'$'];
%         Tc0 = [rowlab,Tc0];
% 
%         %Add & between cells and \\ at end or rows
%         Tc0(:,1:end-1) = strcat(Tc0(:,1:end-1),' & ');
%         Tc0(:,end) = strcat(Tc0(:,end),'\\ ');
%         Tc0(1:2,end) = strcat(Tc0(1:2,end),'[-4pt] ');
% 
%         %Apend this parameters row to the full table
%         Tc = [Tc;Tc0];
%     end
% %     % Sum of NRMSE and number of datasets
% %     rowlab = cell(1,2);
% %     rowlab{1,1} = ' $\Sigma$';
% %     Tc0 = cell(1,nZLBpers);
% %     for iZLBper = 1:nZLBpers
% %         Tc0{1,iZLBper} = ['\scs$[',num2str(nrmse(iZLBper,ispec),'%.3f'),']$'];
% %     end
% %     Tc0 = [rowlab,Tc0];
% % 
% %     %Add & between cells and \\ at end or rows
% %     Tc0(:,1:end-1) = strcat(Tc0(:,1:end-1),' & ');
% %     Tc0(:,end) = strcat(Tc0(:,end),'\\ ');
% %     
% %     %Combine parameter rows
% %     Tc = [Tc;Tc0];
% 
%     %Table footer
%     TcEnd = cell(4,size(Tc0,2));
%     TcEnd{1,1} = '\bottomrule \end{tabular*}';
%     TcEnd{2,1} = '\caption{Median, $(5\%,95\%)$ credible sets and $[NRMSE]$ of moments.}';
%     TcEnd{3,1} = '\label{tab:estimates}';
%     TcEnd{4,1} = '\end{table}';
% 
%     %add footnote
%     Tc = [Tc;TcEnd];
% 
%     %Save table as tex file 
%     fid = fopen([savename,'_',specs{ispec,2},'.tex'], 'w');
%     [nrows,ncols] = size(Tc);
%     format = [repmat('%s ',1,ncols),' \n'];
%     for row = 1:nrows
%         fprintf(fid,format,Tc{row,:});
%     end
%     fclose(fid);
% end