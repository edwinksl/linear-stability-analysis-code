function [func,timeperexpmv] = expcomp(h,A,B,vb,expmvtol,acclejas,EXTREIGS,vn)
%
% expcomp compares the four algorithms for approximating exp(h*inv(B)*A)*vb
% to the prescribed relative change stopping criterion, and returns the
% fastest method and the time needed for each exp(h*inv(B)*A)*vb evaluation
%
% Authors of the MATLAB code:  Minghao W. Rostami and Fei Xue
% Last date of chage: May 8th, 2018
%
% Input:
% for most input parameters, see expm_LR_eigs.m for details
%
% acclejas is an array containing 2 values >= 1; default [1 1];
%
% The 2 values are the speedup factor of P-Leja and R-Leja. By using a
% stepfactor > 1, the corresponding method adopts a substep size that is 
% stepfactor times as large as the automatically found substep size
%
% Our experience is that the speedup factor for R-Leja could be > 1,
% sometimes much larger, and the factor for P-Leja could be slightly > 1;
% however, for Arnoldi and RD, this factor probably has to be 1 for them 
% to approximate the action of matrix exponentials accurately
%
% Output:
% func          the function handle of the fastest method approximating 
%               exp(h*inv(B)*A)*vb to the prescribed tolerance
% timeperexpmv  the time cost for approximating exp(h*inv(B)*A)*vb

if ~exist('acclejas','var')
    acclejas = [1 1];
elseif length(acclejas) ~= 2 || ~all(acclejas >= 1)
    error('An array of 2 scalars >= 1 is needed for input @acclejas.\n');
end
if ~exist('vn','var')
    vn = zeros(length(A),1);
end

%%% factorization of B for polynomial methods
if norm(B-speye(size(B)),'fro') ~= 0
    if norm(B-B','fro')/norm(B,'fro') > 2*eps
        [pL,pU,pP,pQ] = lu(B,0.25);
        lufactors.pL = pL;  lufactors.pU = pU;  lufactors.pP = pP;
        lufactors.pQ = pQ;  lufactors.eyeB = false; lufactors.symB = false;
    else
        [L,D,P,S] = ldl(B,0.1);
        lufactors.pL = L;    lufactors.pLt = L';
        lufactors.pinvD = blkdiaginv(D);
        lufactors.pSP = S*P; lufactors.pPtS = P'*S;
        lufactors.eyeB = false;     lufactors.symB = true;
    end
else
    lufactors.eyeB = true;
end


% RD rational (shift-invert Arnoldi)
steppow_a = -5;     steppow_b = 5;
% default maximum size of the RD space 
ell_rd = 100;
% default shift s = h/ell_rd (a heuristic shift for RD to converge most
% rapidly to exp(h*A)*v in ell_rd steps without restart)
% 1e-4 for tol4000, 6.5e-3 for aerofoil_atf, and 1e-2 for others (expmvtol = 1e-9)
s = 1e-2;   %6.5e-3;   
fprintf('\nAutomatic tuning of parameters for RD-rational with shift = %.2e ...\n',s);
fprintf('Feel free to tune the shift for RD by trial and error to maximize its speed.\n\n');
[rL,rU,rP,rQ,rR] = lu(B-s*A,0.25);
lufactors.rdL = rL;   lufactors.rdU = rU;    lufactors.rdP = rP;
lufactors.rdQ = rQ;   lufactors.rdR = rR;
[~,flag_a] = BRDexpmv(2^steppow_a,s,A,B,vb,vn,1,ell_rd,expmvtol,lufactors);
[~,flag_b] = BRDexpmv(2^steppow_b,s,A,B,vb,vn,1,ell_rd,expmvtol,lufactors);

fprintf('Initial substep size is %.2e.\n',2^steppow_a);
% first locate an interval on which the maximum feasible tau for tau*A can
% be computed within the required tolerance for exp(tau*A)*b
for ii = 1 : 3
    if ~flag_a
        steppow_b = steppow_a;
        steppow_a = steppow_a - 10;
        [~,flag_a] = BRDexpmv(2^steppow_a,s,A,B,vb,vn,1,ell_rd,expmvtol,lufactors);
        fprintf('Initial substep size updated to %.2e. A smaller shift for RD might work better.\n',2^steppow_a);
    elseif flag_b
        steppow_a = steppow_b;
        steppow_b = steppow_b + 10;
        [~,flag_b] = BRDexpmv(2^steppow_b,s,A,B,vb,vn,1,ell_rd,expmvtol,lufactors);
        fprintf('Initial substep size updated to %.2e.\n',2^steppow_a);
    end
end

% bisection
while steppow_b - steppow_a > 1e-2
    steppow_c = (steppow_a+steppow_b)/2;
    [~,flag_c] = BRDexpmv(2^steppow_c,s,A,B,vb,vn,1,ell_rd,expmvtol,lufactors);
    if flag_c
        steppow_a = steppow_c;
        fprintf('substep size updated to %.2e.\n',2^steppow_a);
        if 2^(steppow_a-1e-2) >= h
            break;
        end
    else
        steppow_b = steppow_c;
    end
end
h_powbnd_rd = steppow_a;

substeps_rd = max([1 ceil(h/(2^(h_powbnd_rd-1e-2)))]);
fprintf('RD-rational(%d) takes %d substeps to compute exp(h*invB*A)*v (h = %.2e) to tol %.2e.\n',...
    ell_rd,substeps_rd,h,expmvtol);
afun4 = @(v)BRDexpmv(h,s,A,B,v,vn,substeps_rd,ell_rd,expmvtol,lufactors);



% rational Leja
steppow_a = -5;     steppow_b = 5;
ell_rl = 50;
fprintf('\nAutomatic tuning of parameters for rational Leja ...\n');
[~,flag_a] = raexplejaV2(2^steppow_a,A,B,vb,vn,1,ell_rl,expmvtol);
[~,flag_b] = raexplejaV2(2^steppow_b,A,B,vb,vn,1,ell_rl,expmvtol);
fprintf('Initial substep size is %.2e.\n',2^steppow_a);
for ii = 1 : 3
    if ~flag_a
        steppow_b = steppow_a;
        steppow_a = steppow_a - 10;
        [~,flag_a] = raexplejaV2(2^steppow_a,A,B,vb,vn,1,ell_rl,expmvtol);
        fprintf('Initial substep size updated to %.2e.\n',2^steppow_a);
    elseif flag_b
        steppow_a = steppow_b;
        steppow_b = steppow_b + 10;
        [~,flag_b] = raexplejaV2(2^steppow_b,A,B,vb,vn,1,ell_rl,expmvtol);
        fprintf('Initial substep size updated to %.2e.\n',2^steppow_a);
    end
end
while steppow_b - steppow_a > 1e-2
    if steppow_b - steppow_a < 1 && steppow_b < h_powbnd_rd - 3
        fprintf('R-Leja(%d) unlikely much faster than RD. No more tuning of R-Leja.\n',ell_rl);
        break;
    end
    steppow_c = (steppow_a+steppow_b)/2;
    [~,flag_c] = raexplejaV2(2^steppow_c,A,B,vb,vn,1,ell_rl,expmvtol);
    if flag_c
        steppow_a = steppow_c;
        fprintf('substep size updated to %.2e.\n',2^steppow_a);
        if 2^(steppow_a-1e-2) >= h
            break;
        end
    else
        steppow_b = steppow_c;
    end
end
h_powbnd_ra = steppow_a;

substeps_rl = max([1 ceil(h/(2^(h_powbnd_ra-1e-2)))]);
if acclejas(2) > 1
    fprintf('R-Leja(%d) adopts a stepfactor %d.\n',ell_rl,acclejas(2));
    substeps_rl = max([1 floor(substeps_rl/acclejas(2))]);
end
fprintf('R-Leja(%d) takes %d substeps to compute exp(h*invB*A)*v (h = %.2e) to tol %.2e.\n',...
    ell_rl,substeps_rl,h,expmvtol);

%%% IMPORTANT: make sure a0 and c0 are the same in raexplejaV2.m and here
a0 = 50;    c0 = 0;     tau = h/substeps_rl;
[rL,rU,rP,rQ,rR] = lu(tau*A-(a0+c0)*B,0.25);
lufactors.rL = rL;  lufactors.rU = rU;  lufactors.rP = rP;
lufactors.rQ = rQ;  lufactors.rR = rR;
afun3 = @(v)raexplejaV2(h,A,B,v,vn,substeps_rl,ell_rl,expmvtol,lufactors);



%%% polynomial Leja %%%
ell_pl = 100;
fprintf('\nAutomatic tuning of parameters for P-Leja ...\n');
if ~exist('EXTREIGS','var')
    [~,~,~,EXTREIGS] = explejaV2(2^steppow_a,A,B,1,ell_pl,expmvtol,vb,lufactors);
end
[~,flag_a] = explejaV2(2^steppow_a,A,B,1,ell_pl,expmvtol,vb,lufactors,EXTREIGS,vn);
[~,flag_b] = explejaV2(2^steppow_b,A,B,1,ell_pl,expmvtol,vb,lufactors,EXTREIGS,vn);
fprintf('Initial substep size is %.2e.\n',2^steppow_a);
for ii = 1 : 3
    if ~flag_a
        steppow_b = steppow_a;
        steppow_a = steppow_a - 10;
        [~,flag_a] = explejaV2(2^steppow_a,A,B,1,ell_pl,expmvtol,vb,lufactors,EXTREIGS,vn);
        fprintf('Initial substep size updated to %.2e.\n',2^steppow_a);
    elseif flag_b
        steppow_a = steppow_b;
        steppow_b = steppow_b + 10;
        [~,flag_b] = explejaV2(2^steppow_b,A,B,1,ell_pl,expmvtol,vb,lufactors,EXTREIGS,vn);
        fprintf('Initial substep size updated to %.2e.\n',2^steppow_a);
    end
end
while steppow_b - steppow_a > 1e-2
    if steppow_b - steppow_a < 1 && steppow_b < h_powbnd_rd - 3
        fprintf('P-Leja(%d) unlikely much faster than RD. No more tuning of P-Leja.\n',ell_pl);
        break;
    end
    steppow_c = (steppow_a+steppow_b)/2;
    [~,flag_c] = explejaV2(2^steppow_c,A,B,1,ell_pl,expmvtol,vb,lufactors,EXTREIGS,vn);
    if flag_c
        steppow_a = steppow_c;
        fprintf('substep size updated to %.2e.\n',2^steppow_a);
        if 2^(steppow_a-1e-2) >= h
            break;
        end
    else
        steppow_b = steppow_c;
    end
end
h_powbnd_pl = steppow_a;

substeps_pl = max([1 ceil(h/(2^(h_powbnd_pl-1e-2)))]);
if acclejas(1) > 1
    fprintf('P-Leja(%d) adopts a stepfactor %d.\n',ell_pl,acclejas(1));
    substeps_pl = max([1 floor(substeps_pl/acclejas(1))]);
end
fprintf('P-Leja(%d) takes %d substeps to compute exp(h*invB*A)*v (h = %.2e) to tol %.2e.\n',...
    ell_pl,substeps_pl,h,expmvtol);
afun1 = @(v)explejaV2(h,A,B,substeps_pl,ell_pl,expmvtol,v,lufactors,EXTREIGS,vn);


% regular Arnoldi
steppow_a = -5;     steppow_b = 5;
ell_ar = 100;
fprintf('\nAutomatic tuning of parameters for Arnoldi ...\n');
[~,flag_a] = arnoldiexpv(2^steppow_a,A,B,vb,vn,1,ell_ar,expmvtol,lufactors);
[~,flag_b] = arnoldiexpv(2^steppow_b,A,B,vb,vn,1,ell_ar,expmvtol,lufactors);

fprintf('Initial substep size is %.2e.\n',2^steppow_a);

for ii = 1 : 3
    if ~flag_a
        steppow_b = steppow_a;
        steppow_a = steppow_a - 10;
        [~,flag_a] = arnoldiexpv(2^steppow_a,A,B,vb,vn,1,ell_ar,expmvtol,lufactors);
        fprintf('Initial substep size updated to %.2e.\n',2^steppow_a);
    elseif flag_b
        steppow_a = steppow_b;
        steppow_b = steppow_b + 10;
        [~,flag_b] = arnoldiexpv(2^steppow_b,A,B,vb,vn,1,ell_ar,expmvtol,lufactors);
        fprintf('Initial substep size updated to %.2e.\n',2^steppow_a);
    end
end

% bisection
while steppow_b - steppow_a > 1e-2
    if steppow_b - steppow_a < 2 && steppow_b < h_powbnd_rd - 1
        fprintf('Arnoldi(%d) unlikely much faster than RD. No more tuning of Arnoldi.\n',ell_ar);
        break;
    end
    steppow_c = (steppow_a+steppow_b)/2;
    [~,flag_c] = arnoldiexpv(2^steppow_c,A,B,vb,vn,1,ell_ar,expmvtol,lufactors);
    if flag_c
        steppow_a = steppow_c;
        fprintf('substep size updated to %.2e.\n',2^steppow_a);
        if 2^(steppow_a-1e-2) >= h
            break;
        end
    else
        steppow_b = steppow_c;
    end
end
h_powbnd_ar = steppow_a;

substeps_ar = max([1 ceil(h/(2^(h_powbnd_ar-1e-2)))]);
fprintf('Arnoldi(%d) takes %d substeps to compute exp(h*invB*A)*v (h = %.2e) to tol %.2e.\n',...
    ell_ar,substeps_ar,h,expmvtol);
afun2 = @(v)arnoldiexpv(h,A,B,v,vn,substeps_ar,ell_ar,expmvtol,lufactors);


fprintf('\nComparing the algorithms ... \n\n');

costpermvp1 = nnz(A);
if ~lufactors.eyeB
    if ~lufactors.symB
        costpermvp1 = costpermvp1 + nnz(lufactors.pU) + nnz(lufactors.pL);
    else
        costpermvp1 = costpermvp1 + 2*nnz(lufactors.pL) + 2*nnz(lufactors.pSP) + nnz(lufactors.pinvD);
    end
end
cost1 = costpermvp1*substeps_pl*ell_pl;

cost2 = costpermvp1*substeps_ar*ell_ar;
% the additional factor 2 here refers to double MGS in Arnoldi
cost2 = cost2 + 2*substeps_ar*(ell_ar-1)*ell_ar*length(A);

costpermvp3 = nnz(lufactors.rL) + nnz(lufactors.rU) + nnz(A) + nnz(B);
cost3 = costpermvp3*substeps_rl*45; %*ell_rl;

costpermvp4 = nnz(lufactors.rdL) + nnz(lufactors.rdU) + nnz(B);
cost4 = costpermvp4*substeps_rd*ell_rd;
% the additional factor 2 here refers to double MGS in RD's Arnoldi
cost4 = cost4 + 2*substeps_rd*(ell_rd-1)*ell_rd*length(A);

fprintf('Estimated flops per expmv:\n');
fprintf('RD-rational(%d)\t %.2e\n',ell_rd,cost4);
fprintf('R-Leja(%d)\t\t %.2e\n',ell_rl,cost3);
fprintf('P-Leja(%d)\t\t %.2e\n',ell_pl,cost1);
fprintf('Arnoldi(%d)\t\t %.2e\n\n',ell_ar,cost2);

[sortedcost,idx] = sort([cost1 cost2 cost3 cost4]);
cellfunc = cell(4,2);
cellfunc{1,1} = afun1;  cellfunc{1,2} = sprintf('P-Leja(%d)',ell_pl);
cellfunc{2,1} = afun2;  cellfunc{2,2} = sprintf('Arnoldi(%d)',ell_ar);
cellfunc{3,1} = afun3;  cellfunc{3,2} = sprintf('R-Leja(%d)',ell_rl);
cellfunc{4,1} = afun4;  cellfunc{4,2} = sprintf('RD-rational(%d)',ell_rd);
num_competitors = 1;
for ii = 2 : 4
    if sortedcost(ii) < 2*sortedcost(1)
        num_competitors = ii;
    else
        break;
    end
end

tstart = uint64(zeros(num_competitors,1));
telapsed = zeros(num_competitors,1);
%if num_competitors > 1
for ii = 1 : num_competitors
    func = cellfunc{idx(ii),1};
    name = cellfunc{idx(ii),2};
    tstart(ii) = tic;
    [~,flag,mvps] = func(vb);
    telapsed(ii) = toc(tstart(ii));
    str = sprintf(' returned flag %d in %.3d seconds (%d mvps).',flag,telapsed(ii),mvps);
    fprintf([name,str,'\n']);
    if flag ~= 1
        fprintf('CAUTION: this expmv algorithm did not converge to tolerance %.2e.\n',expmvtol);
        %telapsed(ii) = realmax;
    end
end
%end

[timeperexpmv,fastest] = min(telapsed);
func = cellfunc{idx(fastest),1};
name = cellfunc{idx(fastest),2};
fprintf(['\nChoosing the fastest method ',name,'\n\n']);

return
