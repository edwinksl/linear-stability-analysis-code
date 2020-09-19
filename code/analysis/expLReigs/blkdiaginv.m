function invD = blkdiaginv(D)

% returns the inverse of the block diagonal matrix D, with 1x1 or 2x2
% blocks on the diagonal; used for more efficient linear solves based on
% the LDL factorization of real symmetric matrices (by polynomial methods
% Leja and Arnoldi for approximating exmpv)
invD = D;
n = length(D);
ii = 1;
while ii < n
    if D(ii+1,ii) ~= 0
        invD(ii:ii+1,ii:ii+1) = D(ii:ii+1,ii:ii+1)\speye(2);
        ii = ii + 2;
    else
        invD(ii,ii) = D(ii,ii)\1;
        ii = ii + 1;
    end
end
if ii == n
    invD(n,n) = D(n,n)\1;
end