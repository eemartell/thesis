% Tests whether fixed point and time iteration solutions are numerically equivalent
load('solutionti1.mat')
perbind_ti = perbind;
pf_ti = pf;
load('solutionfp1.mat')
perbind_fp = perbind;
pf_fp = pf;
if all(pf_fp.c(:) - pf_ti.c(:) < P.tol) && all(pf_fp.pigap(:) - pf_ti.pigap(:) < P.tol)
    disp(['The solutions are equivalent up to a tolerance of: ' num2str(P.tol)]);
end
if perbind_fp == perbind_ti
    disp('ZLB percent bindings are equivalent between solutions')
end
disp(['Max difference in consumption policy: ',num2str(max(pf_fp.c(:) - pf_ti.c(:)))])
disp(['Max difference in inflation gap policy: ',num2str(max(pf_fp.pigap(:) - pf_ti.pigap(:)))])