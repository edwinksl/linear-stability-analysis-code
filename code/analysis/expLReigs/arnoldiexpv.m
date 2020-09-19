function [expmv,flag,mvps] = arnoldiexpv(h,A,B,b,vn,substeps,m,reltol,lufactors)

% arnoldiexpv implements the standard Arnoldi's method for approximating
% exp(h*inv(B)*A)*vb to the relative change tolerance reltol
%
warning('off','MATLAB:nearlySingularMatrix');
if ~exist('lufactors','var') && norm(B-speye(size(B)),'fro') ~= 0
    if norm(B-B','fro')/norm(B,'fro') > 2*eps
        [L,U,P,Q] = lu(B,0.25);
        lufactors.pL = L;    lufactors.pU = U;    lufactors.pP = P;
        lufactors.pQ = Q;    lufactors.eyeB = false;     lufactors.symB = false;
    else
        [L,D,P,S] = ldl(B,0.1);
        lufactors.pL = L;    lufactors.pLt = L';
        lufactors.pinvD = blkdiaginv(D);
        lufactors.pSP = S*P; lufactors.pPtS = P'*S;
        lufactors.eyeB = false;     lufactors.symB = true;
    end
end
V = zeros(length(A),m+1);
H = zeros(m+1,m);
expmv = b/norm(b);
scaling = norm(b);
%oldexpmv = zeros(length(A),1);
flags = zeros(substeps,1);
regularvn = (rank(vn,sqrt(eps)) == size(vn,2));
tau = h/substeps;
mvps = 0;

for iter = 1 : substeps
    H(:,:) = 0;
    V(:,:) = 0;
    oldprjexpm = zeros(0,1);
    scaling = scaling*norm(expmv);
    V(:,1) = expmv/norm(expmv);
    for ii = 1 : m
        if lufactors.eyeB
            v = A*V(:,ii);
        elseif ~lufactors.symB
            v = lufactors.pQ*(lufactors.pU\(lufactors.pL\(lufactors.pP*(A*V(:,ii)))));
        else
            v = lufactors.pSP*(lufactors.pLt\(lufactors.pinvD*(lufactors.pL\(lufactors.pPtS*(A*V(:,ii))))));
        end
        mvps = mvps+1;
        if regularvn
            v = v - vn*((vn'*vn)\(vn'*v));
            v = v - vn*((vn'*vn)\(vn'*v));
        end
        for jj = 1 : ii
            H(jj,ii) = V(:,jj)'*v;
            v = v-V(:,jj)*H(jj,ii);
        end
        for jj = 1 : ii
            dh = V(:,jj)'*v;
            H(jj,ii) = H(jj,ii)+dh;
            v = v-V(:,jj)*dh;
        end
        if regularvn
            v = v - vn*((vn'*vn)\(vn'*v));
            v = v - vn*((vn'*vn)\(vn'*v));
        end
        H(ii+1,ii) = norm(v);
        V(:,ii+1) = v/norm(v);
        
        prjexpm = expm(tau*H(1:ii,1:ii));
%         [HU,HD] = eig(tau*H(1:ii,1:ii));
%         prjexpm = HU*diag(exp(diag(HD)))/HU;
%         if realprob
%             prjexpm = real(prjexpm);
%         end
        
        %expmv = V(:,1:ii)*prjexpm(:,1);
        if norm([oldprjexpm; 0]-prjexpm(:,1))/norm(prjexpm(:,1)) < reltol
            flags(iter) = 1;
            break;
        else
            oldprjexpm = prjexpm(:,1);
        end
    end
    expmv = V(:,1:ii)*prjexpm(:,1);
end
flag = (nnz(flags) > substeps*3/4);

if regularvn
    expmv = expmv - vn*((vn'*vn)\(vn'*expmv));
    expmv = expmv - vn*((vn'*vn)\(vn'*expmv));
end

expmv = expmv*scaling;
warning('on','MATLAB:nearlySingularMatrix');
