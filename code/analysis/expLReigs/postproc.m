function [V,D,relres] = postproc(A,B,V,dscp)

% postproc works on an approximate invariant subspace basis V, and outputs
% the desired eigenvectors described by dscp

    if isreal(A) && isreal(B)
        kk = 1;
        while kk <= size(V,2)
            if isreal(V(:,kk))
                kk = kk+1;
            else
                V(:,kk:kk+1) = [real(V(:,kk)) imag(V(:,kk))];
                V(:,kk) = V(:,kk)/norm(V(:,kk));    V(:,kk+1) = V(:,kk+1)/norm(V(:,kk+1));
                kk = kk+2;
            end
        end
    end
    RQ = (B*V)\(A*V); 
    [U,S] = eig(RQ);
    ritzval = diag(S);
    if strcmpi(dscp,'lr')
        [~,idx] = sort(real(ritzval),'descend');
    elseif strcmpi(dscp,'sr')
        [~,idx] = sort(real(ritzval));
    elseif strcmpi(dscp,'li2')
        [~,idx] = sort(abs(imag(ritzval)),'descend');
    elseif strcmpi(dscp,'li')
        [~,idx] = sort(imag(ritzval),'descend');
    elseif strcmpi(dscp,'si')
        [~,idx] = sort(imag(ritzval));
    else
        error('\nDescription not supported!');
    end
    ritzval = ritzval(idx);
    U = U(:,idx);
    V = V*U;
    for kk = 1 : size(V,2)
        V(:,kk) = V(:,kk)/norm(V(:,kk));
    end
    D = diag(ritzval);
    RES = A*V-B*(V*D); AVNM = sqrt(sum(conj(A*V).*(A*V))); RESNM = sqrt(sum(conj(RES).*(RES)));
    relres = (RESNM./AVNM)';
    fprintf('\nThe %d ''%s'' eigenvalues of (A,B) and relative residual norms :\n',length(ritzval),dscp);
    disp_evalres(ritzval,relres);
    fprintf('\n');
end