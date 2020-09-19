function [expmv,flag,mvps] = BRDexpmv(h,s,A,B,F,vn,substeps,m,reltol,lufs)

% BRDexpmv is the RD (restricted denominator, i.e., shift-invert Arnoldi)
% method approximating exp(h*inv(B)*A)*F by the Krylov subspace K_m(L,F),
% where the linear operator L = inv(I-s*inv(B)*A) = inv(B-s*A)*B
%
% Authors of the MATLAB code:  Minghao W. Rostami and Fei Xue
% Last date of chage: May 8th, 2018
%
% Reference: I. Moret & P. Novati, RD-Rational Approximations of the Matrix
%            Exponential, BIT Numer. Math., Vol. 44 (2004), pp 595?615
%
% Input:
% h         positive scaling factor in exp(h*inv(B)*A)*F
% A,B       matrix A  or matrix pencil (A,B) of interest
% F         the right-hand side in exp(h*inv(B)*A)*F
% s         positive shift for the RD method; crucial for the efficiency
% vn        the shared null space basis of A and B
% substeps  the number of substeps for approximating exp(h*inv(B)*A)*F
%           in each substep, we approximate exp(tau*inv(B)*A)*K, where 
%           tau = h/substeps and K is the result from the previous substep
% m         the maximum Krylov subspace dimension (typically <= 100)
% reltol    the relative change stopping criterion for all expmv algorithms
% lufs      precomputed LU factors for B-s*A, if provided
% 
% Output:
% expmv     approximation to exp(h*inv(B)*A)*F
% flag      whether RD has converged numerically; true or false
% mvps      the number of linear solves involving B-s*A

if h < 0
    error('Input h must be nonnegative.');
elseif h == 0
    expmv = F;  flag = 1;   mvps = 0;   return;
end
warning('off','MATLAB:nearlySingularMatrix');
% we assume that A is nonsymmetric, and so is B-s*A. Therefore no need to
% perform ldl factorization of B-s*A
if ~exist('lufs','var')
    [L,U,P,Q,R] = lu(B-s*A,0.25);
    lufs.rdL = L;	lufs.rdU = U;	lufs.rdP = P;	lufs.rdQ = Q;	lufs.rdR = R;
end
% lufs must have the lu factors of B-s*A;
p = size(F,2);
V = zeros(length(A),(m+1)*p);
H = zeros((m+1)*p,m*p);
[expmv,R] = qr(F,0);
scaling = R;
flags = zeros(substeps,1);
regularvn = (rank(vn,sqrt(eps)) == size(vn,2));
tau = h/substeps;
mvps = 0;
relchange = zeros(substeps*m,1);
normhis = zeros(substeps*m+1,1);
normtype = inf;
normhis(1) = norm(F,normtype);
if nnz(B-speye(size(B))) == 0,  idB = true;
else,    idB = false;
end
erratic = 0;
iterref = (1 < 1);  % iterative refinement needed? (false by default)
if iterref, fprintf('Iterative refinement enabled.\n'); end
lastokvsize = 0;    lastokexpm = zeros(0,p);

for iter = 1 : substeps
    H(:,:) = 0;     V(:,:) = 0;
    oldprjexpm = zeros(0,p);
    [V(:,1:p),R] = qr(expmv,0);
    scaling = R*scaling;
    for ii = 1 : m
        iiblkidx = ((ii-1)*p+1):ii*p;
        
        if idB, rhs = V(:,iiblkidx);
        else,   rhs = B*V(:,iiblkidx);
        end
        v = lufs.rdQ*(lufs.rdU\(lufs.rdL\(lufs.rdP*(lufs.rdR\rhs))));
        mvps = mvps + p;
        % iterative refinement, if necessary
        if iterref
            rhs = rhs - (B*v-s*(A*v));
            deltav = lufs.rdQ*(lufs.rdU\(lufs.rdL\(lufs.rdP*(lufs.rdR\rhs))));
            v = v + deltav;     mvps = mvps + p;
        end
        
        if regularvn
            v = v - vn*((vn'*vn)\(vn'*v));  v = v - vn*((vn'*vn)\(vn'*v));
        end
        for jj = 1 : ii
            jjblkidx = ((jj-1)*p+1):jj*p;
            H(jjblkidx,iiblkidx) = V(:,jjblkidx)'*v;
            v = v-V(:,jjblkidx)*H(jjblkidx,iiblkidx);
        end
        % reorthogonalization seems necessary for this algorithm
        for jj = 1 : ii
            jjblkidx = ((jj-1)*p+1):jj*p;
            dh = V(:,jjblkidx)'*v;
            H(jjblkidx,iiblkidx) = H(jjblkidx,iiblkidx)+dh;
            v = v-V(:,jjblkidx)*dh;
        end
        if regularvn
            v = v - vn*((vn'*vn)\(vn'*v));  v = v - vn*((vn'*vn)\(vn'*v));
        end
        [V(:,iiblkidx+p),R] = qr(v,0);
        H(iiblkidx+p,iiblkidx) = R;
        
        prjexpm = expm((tau/s)*(eye(ii*p)-H(1:ii*p,1:ii*p)\eye(ii*p)));
        normhis(nnz(normhis)+1) = norm(prjexpm(:,1:p),normtype);
        tmp_relchange = norm([oldprjexpm; zeros(p,p)]-prjexpm(:,1:p),normtype)/norm(prjexpm(:,1:p),normtype);
        relchange(nnz(relchange)+1) = tmp_relchange;
        
        oldprjexpm = prjexpm(:,1:p);
        if norm(prjexpm(:,1:p),normtype) <= realmax
            lastokexpm = prjexpm(:,1:p);
            lastokvsize = ii*p;
            erratic = 0;
            if tmp_relchange <= reltol
                flags(iter) = 1;
                break;
            end
        else
            erratic = erratic + 1;
            if erratic >= 100
                fprintf('NaN or Inf found in 100 consecutive iters at step %d.\n',ii);
                fprintf('Returning the result at previous good step.\n');
                break;
            end
        end
    end
    expmv = V(:,1:lastokvsize)*lastokexpm;
end
flag = (nnz(flags) >= substeps*1/2);

if regularvn
    expmv = expmv - vn*((vn'*vn)\(vn'*expmv));
    expmv = expmv - vn*((vn'*vn)\(vn'*expmv));
end

expmv = expmv*scaling;
warning('on','MATLAB:nearlySingularMatrix');
