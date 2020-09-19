function [V,D,failflag,iter] = arnoldi(afun,v0,m,p,eigtol)

% Standard Arnoldi's method for eigenvalue computation
%
% Input:
% afun      linear operator (v -> exp(h*inv(B)*A)*v in our exp eig setting)
% v0        starting vector of Arnoldi
% m         maximum Arnoldi steps
% p         number of desired dominant eigenpairs
% eigtol    relative tolerance of the dominant invariant subspace
% 
% Output:
% V,D       disired dominant eigenvector/eigenvalue approximations  
% failflag  whether Arnoldi's method converged in m steps to eigtol
% iter      number of Arnoldi steps taken upon exit

failflag = true;
W = zeros(length(v0),m+1);
H = zeros(m+1,m);
W(:,1) = v0/norm(v0);
eigres = realmax;
% teststep = 25;
% testtol = 2.5e-2;
for ii = 1 : m
    w = afun(W(:,ii));
    for jj = 1 : ii
        H(jj,ii) = W(:,jj)'*w;
        w = w - W(:,jj)*H(jj,ii);
    end
    for jj = 1 : ii
        dh = W(:,jj)'*w;
        w = w - W(:,jj)*dh;
        H(jj,ii) = H(jj,ii) + dh;
    end
    H(ii+1,ii) = norm(w);
    w = w/H(ii+1,ii);
    nmw = norm(w);
    H(ii+1,ii) = H(ii+1,ii)*nmw;
    W(:,ii+1) = w/nmw;
    
    [U,D] = eig(H(1:ii,1:ii));
    dg = diag(D);
    [~,idx] = sort(abs(dg),'descend');
    U = U(:,idx);   dg = dg(idx);
    D = diag(dg);
    q = min([ii p]);
    U = U(:,1:q);   D = D(1:q,1:q);
    V = W(:,1:ii)*U;
    eigres = H(ii+1,ii)*norm(U(end,:))/norm(U*D);
    if ii == p || (ii >= p && mod(ii,5) == 0) || eigres <= eigtol
        fprintf('Arnoldi step %d: relative eig-res %.2e.\n',ii,eigres);
    end
    if ii >= p && eigres <= eigtol
        failflag = false;
        fprintf('Eigen tol %.2e satisfied at step %d. Arnoldi terminated.\n',eigtol,ii);
        break;
%     elseif ii >= teststep && eigres >= testtol
%         failflag = true;
%         fprintf('Note: Arnoldi did not decrease eigen residual to %.2e in %d steps.\n',...
%             testtol,ii);
        %return;
    end
end
%V = W(:,1:ii)*U;
iter = ii;
if eigres > eigtol
    fprintf('Arnoldi failed to converge to tolerance %.2e in %d steps.\n',eigtol,ii);
end