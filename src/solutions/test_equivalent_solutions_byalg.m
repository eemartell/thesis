% Tests whether fixed point and time iteration solutions are numerically equivalent
load('solutionfpART.mat')
perbind_ALG = perbind;
pf_ALG.c = pf.hh;
pf_ALG.pigap = pf.firm;
load('solutionfpGust.mat')
perbind_Gust = perbind;
pf_Gust.c = 1/pf.hh;
pf_Gust.pigap = (1+sqrt((P.varphi+4*pf.firm)/P.varphi))/2;
if all(pf_Gust.c(:) - pf_ALG.c(:) < P.tol) && all(pf_Gust.pigap(:) - pf_ALG.pigap(:) < P.tol)
    disp(['The solutions are equivalent up to a tolerance of: ' num2str(P.tol)]);
end
if perbind_Gust == perbind_ALG
    disp('ZLB percent bindings are equivalent between solutions')
end
disp(['Max difference in consumption policy: ',num2str(max(pf_Gust.c(:) - pf_ALG.c(:)))])
disp(['Max difference in inflation gap policy: ',num2str(max(pf_Gust.pigap(:) - pf_ALG.pigap(:)))])