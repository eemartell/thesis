clear
clc

% Set random number seed
mtstream = RandStream('mt19937ar','seed',123456);
RandStream.setGlobalStream(mtstream);

% Load options
load('options')
V.nstate = V.nvar-V.nfore;

% Load solution
load('solution_test.mat')

%% Stochastic steady state
T = 20000;
epsg = zeros(T,1);
epss = zeros(T,1);
epsmp = zeros(T,1);
sims = simulation(pf,P,S,G,V,epsg,epss,epsmp);
stochss = sims(end,:);
format longg
disp(stochss')
% Ergodic distribution
nergdraws = 10000;
nergpers = 500;
epsg = randn(nergpers,nergdraws);
epss = randn(nergpers,nergdraws);
epsmp = randn(nergpers,nergdraws);
ergsims = simulation(pf,P,S,G,V,epsg,epss,epsmp);
ergdraws = squeeze(ergsims(end,:,:));
isnumidx = ones(nergdraws,1);
for idraw = 1:nergdraws
    for istate = 1:V.nstate
        if isnan(ergdraws(istate,idraw))
            isnumidx(idraw) = 0;
            break
        end
    end
end
nergdrawsisnum = sum(isnumidx);
ergdrawsisnum = zeros(V.nstate,nergdrawsisnum);
idrawisum = 0;
for idraw = 1:nergdraws
    if isnumidx(idraw) == 1
        idrawisum = idrawisum + 1;
        ergdrawsisnum(:,idrawisum) = ergdraws(:,idraw);
    end
end

format shortg
% in, c, k, x
disp('Ergodic dist. credible sets of endog. vars')

disp('  notional rate 1/99 sims')
disp(quantile(ergdrawsisnum(V.in,:),[.001,.999]))
disp('  notional rate bounds')
disp([G.inmin,G.inmax])

disp('  consumption 1/99 sims')
disp(quantile(ergdrawsisnum(V.c,:),[.001,.999]))
disp('  consumption bounds')
disp([G.cmin,G.cmax])

%% Find simulations with singular ZLB event with specified duration
T = 120;
nsimstemp = 100;
nsims = 10;
% Duration of ZLB event(s) as percentage of simulation length
ZLBpers = [0,.05];
nZLBpers = numel(ZLBpers);
sims = NaN(T,V.nvar-V.nfore,nsims,nZLBpers);
for iZLBper = 1:nZLBpers
    ZLBnq = ZLBpers(iZLBper)*T;
    if ZLBnq == 0
        inbound = 1.00125;
    elseif ZLBnq > 0
        inbound = 1.00;
    end
    tic
    isim = 0;
    while isim < nsims
        % Simulate
        simstemp1 = ergdraws(:,ceil(rand(nsimstemp,1)*nergdraws));
        epsg = randn(T,nsimstemp);
        epss = randn(T,nsimstemp);
        epsi = randn(T,nsimstemp);
        simstemp = simulation(pf,P,S,G,V,epsg,epss,epsi,simstemp1);

        % Keep simulations that meet criteria
        binds = squeeze(simstemp(:,V.in,:) <= inbound); % use where in Fortran
        bindsdiff = binds(2:T,:)-binds(1:T-1,:);
        for isimtemp = 1:nsimstemp
            temp = simstemp(:,:,isimtemp);
            isnantemp = any(isnan(temp(:)));
            if ~isnantemp
                nZLBqtemp = sum(binds(:,isimtemp));
                savesim = 0;
                if ZLBnq == 0
                    savesim = nZLBqtemp == 0 && isim < nsims;
                elseif ZLBnq > 0
                    enters = bindsdiff(:,isimtemp) == 1; % use where in Fortran
                    nenters = sum(enters);
                    exits = bindsdiff(:,isimtemp) == -1; % use where in Fortran
                    nexits = sum(exits);
                    savesim = ((nenters == 1 && nexits == 1) || ...
                            (nenters == 1 && nexits == 0) || ...
                            (nenters == 0 && nexits == 1)) && ...
                        nZLBqtemp == ZLBnq && ...
                        isim < nsims;
                end
                % Save simulation
                if savesim
                    isim = isim + 1;
                    sims(:,:,isim,iZLBper) = simstemp(:,:,isimtemp);
                    disp(isim)
                end
            end
        end
    end
end
toc
binds = squeeze(sims(:,V.in,:,:) <= 1);
bindsdiff = diff(binds,1,1);

% Test Fotran sims
simstemp = textread('Fortran/sims.txt','%f');
sims = reshape(simstemp,[120,V.nvar-V.nfore,nsims,nZLBpers]);
binds = squeeze(sims(:,V.in,:,:) <= 1);
bindsdiff = diff(binds,1,1);