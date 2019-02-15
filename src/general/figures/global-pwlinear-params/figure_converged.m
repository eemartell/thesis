
clear
clc
close all

npers = 120;
nsimsall = 100;
nsims = 50;
nZLBpers = 6;
nqs = [0,6,12,18,24,30];

%Load in V from misspecified model
load('../../../Results-nocap/options.mat','V','P','F')


filename = ['../../../global-pwlinear-params/me5/states/converged.txt'];
if exist(filename, 'file') == 2
    converged = textread(filename,'%f');
    converged = reshape(converged,[nsims,nZLBpers]);
end


50 - sum(converged,1)

% Loop through ZLB durations
for inq = 1:numel(nqs)
    % ZLB duration
    nq = nqs(inq);
    disp(['ZLB duration = ' num2str(nq)])
    for isim = 1:nsims
        %Load posterior distributions they exist, otherwise keep as
        %NaNs
        filename = ['../../../Estimation-pwlinear/','zlb',num2str(inq),'/mean/',num2str(isim),'-postmean.txt'];
        if exist(filename, 'file') == 2
            %Calculate mean across posterior distribution
            means(:,isim,inq) = textread(filename,'%f');
        end
    end
end

T = table(converged(:),'VariableNames',{'converged'});
for ivar = 1:F.nparam
    temp = means(ivar,:,:);
    T= [T,table(temp(:),'VariableNames',{F.priors{ivar,1}})];
end

%T.sigsum = T.sigs + T.sigg + T.sigmp;
tree = fitctree(T,'converged','MaxNumSplits',6);
view(tree,'Mode','graph')
y_hat = predict(tree, T);
cm = confusionmat(T.converged,y_hat)