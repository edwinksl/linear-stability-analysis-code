function [phikAv,flag,totalmvps,extreigs] = explejaV2(h,A,B,nsteps,maxm,NewTol,varargin)

% Original code by M. Calari & P. Kandolf, modified by Fei Xue to have the 
% same input parameters as other algorithms approximating exp(h*inv(B)*A)*v

%%
%EXPLEJA   Matrix exponential times vector or matrix.
%   [PHIKAV, ERREST, INFO, EXTREIGS] = EXPLEJA (H, A, V, TOL, EXTREIGS)
%
%   Compute EXPM(H * A) * V without forming EXP(H*A). A is a matrix 
%   *preferably* sparse and with eigenvalues with negative real part. 
%   V is a vector or a matrix.
%       If length (TOL) == 1, then TOL is the absolute tolerance. If
%   length (TOL) == 2, then TOL(1) is the absolute tolerance and TOL(2)
%   the tolerance relative to the current approximation of EXPM(H * A)*V. 
%   If length (TOL) == 3, then TOL(3) specifies the norm, by default this
%   is INF. If nothing is provided TOL = [0, 2^(-24), INF].
%       If EXTREIGS is not given, the spectrum of A is estimated by
%   Gershgorin's disks theorem. Otherwise, it is estimated by EXTREIGS.SR
%   (smaller real part of the eigenvalues), EXTREIGS.LR (largest real
%   part) and EXTREIGS.LI2 (squared largest imaginary part).
%       On output, sum (INFO) is the total number of iterations
%   (== matrix-vector products) and sum(ERREST) the estimated error (in
%   the specified norm). On output, EXTREIGS is like above.
%
%   The code is based on PHILEJA provided by Marco Caliari.

%   Reference: M. Caliari, P. Kandolf, A. Ostermann and S. Rainer, 
%   Comparison of software for computing the action of the matrix 
%   exponential. BIT, 2013, DOI: 10.1007/s10543-013-0446-0

%   Peter Kandolf, June 26, 2014, cooperation with Marco Caliari,
%   Alexander Ostermann and Stefan Rainer

%   Additional notes:
%   - The error is estimated during the newton interpolation, which is
%   performed until 
%        ||errorestimate||_TOL(3) > max(TOL(1),TOL(2)*||Y||_TOL(3))
%   is satisfied. 
%   - The number of inner steps is estimated like described in the above
%   paper, for sampling tolerances 10.^[-4,-6,-8,-10,-12].
%   - The necessary divided differences are computed by the algorithm
%   of M. Caliari, see nested functions below (line 323).
% 

%   Minimal example:
%       A=diag([-10:0])+diag(ones(1,10),1);
%       v=ones(11,1); h=1;
%       y=expleja(h,A,v);
%%

warning('off','MATLAB:nearlySingularMatrix');
if length(varargin) < 3 && norm(B-speye(size(B)),'fro') ~= 0
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
elseif norm(B-B(1,1)*speye(size(B)),'fro') == 0
    lufactors.pL = speye(size(A));  lufactors.pU = B;
    lufactors.pQ = speye(size(A));  lufactors.pQ = speye(size(A));
    if B(1,1) == 1, lufactors.eyeB = true;
    else,           lufactors.eyeB = false;
    end
else
    lufactors = varargin{2};
end
%% Spectrum estimate and parameter check
if (isfloat (A))
    % check consistency of matrix and vector sizes
    v = varargin{1};
    if (size (v,1) ~= size (A, 2)), size(v), size(A)
        error ('Inconsistent matrix and vector sizes')
    end
    %Check if tolerance is given and set to default otherwise
    
    % get spectral estimate
    if (nargin >= 8) && ~isempty(varargin{3})
        extreigs = varargin{3}; 
    else
        fprintf('No spectrum estimate. Gershgorin disc used to estimate the spectrum (may not be accurate).\n');
        if norm(B-B(1,1)*speye(size(B)),'fro') > 0
            fprintf('Matrix B is not a multiple of identity. Please DO provide spectrum estimate.\n');
            fprintf('Returned spectrum estimate unreliable!\n');
        end
        extreigs = gersh(A);
    end
else
    error('This version is only supporting matrices');
end
SR = extreigs.SR * max(h); LR = extreigs.LR * max(h); 
LI2 = extreigs.LI2 * max(h) ^ 2;

% h and/or A are zero - computation is finished
if (abs(SR)+abs(LR)+abs(LI2) == 0)
    phikAv = v;         return
end

%% Point selection and divided differences
% compute the properties of the ellipse with smallest capacity 
% circumscribing the box given by extreigs
d = (SR + LR) / 2;
C2 = ((LR - SR) / 2) ^ 2;
c2 = (nthroot (C2, 3) + nthroot (LI2, 3)) * (C2 ^ (2/3) - LI2 ^ (2/3));
gamma = (nthroot (C2, 3) * sqrt (nthroot (C2, 3) + nthroot (LI2, 3)) + ...
    sqrt (nthroot (C2 * LI2 ^ 2, 3) + LI2)) / 2;
gammafoc = sqrt (abs (c2)) / 2;

if (c2 >= 0)
    % real case
    % the focal interval should be contained in the left half plane
    gammafoc = min(gammafoc,-d/2);
    % gammafoc cannot be zero (constant diagonal case, ellipse=circle, ...)
    if (gammafoc < 1e-4), gammafoc = max(gammafoc,-d/100); end
else
    % Complex case
    % gammafoc cannot be zero (constant diagonal case, ellipse=circle, ...)
    if (gammafoc < 1e-4), gammafoc = max(gammafoc,gamma/2); end
end

% Compute the maximal and minimal amount of interpolation points and the
% number of inner steps necessary for the computation
%[mt,nsteps]=getsubsteps(tol,SR,sqrt(LI2),c2);
%maxm=max(mt);
%m=min(mt);
%if isnan(m), error('stepsize is to big'); end

% Selection of Leja points
if (c2 >= 0)
    % real Leja points
    xi = leja(maxm) + d / gammafoc;
    newt = @newton_fix;
else
    % symmetriezed complex Leja points
    xi = 1i * lejas(maxm) + d / gammafoc;
    newt = @newtons_fix;
end

% Compute the divided differences depending on inner steps and 
% interpolation points
dd = exptaylordd(xi,gammafoc/nsteps);

%% Inner steps loop for the newton interpolation
totalmvps = 0;
phikAv = v;
flags = zeros(nsteps,1);
for j = 1:nsteps
    [phikAv,flags(j),mvps] = newt(h,A,lufactors,phikAv,xi,dd,gammafoc,NewTol,maxm);
    totalmvps = totalmvps + mvps;
%    errest(j) = perrest; info(j) = pinfo;
    if length(varargin) > 3
        vn = varargin{4};
        if norm(vn) > 0
            phikAv = phikAv - vn*((vn'*vn)\(vn'*phikAv));
            phikAv = phikAv - vn*((vn'*vn)\(vn'*phikAv));
        end
    end
end
flag = (nnz(flags) >= nsteps*1/2);
warning('on','MATLAB:nearlySingularMatrix');
end

%% Nested functions 

%% Spectral estimate
%
%   Computaiton of Gerschgorin disks for the given matrix to get upper 
%   bound for the eigenvalues of A
%
function extreigs = gersh (A,varargin)
if issparse(A) && nargin == 1
    %Original version - works well for sparse matrices
    AS = (A + A') / 2;
    radius = sum (abs (AS), 2);
    extreigs.SR = full (min (diag (AS) - (radius - abs (diag (AS)))));
    extreigs.LR = min(full(max(diag (AS) + (radius - abs(diag(AS))))), 0);
    AS = (A - AS);
    extreigs.LI2 = full (max(sum (abs (AS), 2)))^2;
else
    radius=zeros(size(A,2),2);
    for i=1:size(A,2)
        radius(i,1)=sum(abs(A(i,:)'+A(:,i)))/2;
        radius(i,2)=sum(abs(A(i,:)'-A(:,i)))/2;
    end
    d=real(diag(A));
    extreigs.SR=(min(d - (radius(:,1) - abs(d))));
    extreigs.LR=min((max(d + (radius(:,1) - abs(d)))), 0);
    extreigs.LI2=(max(radius(:,2)))^2;
end
end

%% Get substeps
%
%   Computes the number of substeps needed to solve the corresponding
%   matrix interpolation by interpolation according to the eigenvalues of
%   of the problem.
%
% function [stepperit,nsteps] = getsubsteps(tol, sr, li,c2,varargin)
% %Printing is disabled by default
% print=false; if nargin==5, print=varargin{1}; end
% 
% %Select tolerance and if non is found modify the value of the
% %multiplicative factor
% maxnsteps=200;
% factor=[1.2,1];
% t=min(tol(tol(1:2)~=0)); if isempty(t), t=1e-13; end
% Tol=[1e-4,1e-6,1e-8,1e-10,1e-12];
% tmp=find(Tol<=t,1,'first');
% if isempty(tmp)
%     factor=factor+0.5;
%     if print
%     fprintf(1,'Nothing specified for such a tolerance,\n');
%     fprintf(1,'taking %1.1e with multiplication factors of %1.2e %1.2e\n',...
%         Tol(end),factor);
%     end
%     tmp=length(Tol);
% end
% 
% %Make sure SR and LI are in the the correct range
% sr=min(sr,-1/100);
% li=max(li,1/100);
% if c2>=0
%     [F,boarder]=getReal(tmp);
%     maxstepperit=[151,151,151,151,131];
%     nsteps=max(1,max(ceil([sr;li]./boarder)));
% else
%     [F,boarder]=getComplex(tmp);
%     maxstepperit=[151,151,151,151,151];
%     nsteps=max(1,max(ceil([li;sr]./boarder)));
% end
% 
% for nsteps=nsteps:maxnsteps
%     stepperit=ceil(F(-sr/nsteps,li/nsteps).*factor);
%     if max(stepperit)<maxstepperit(tmp)
%         break;
%     end
% end
% 
% if isempty(nsteps) || (nsteps==maxnsteps && max(stepperit)>maxstepperit(tmp))
%     error('reached substep limit of %d', maxnsteps);
% end
% 
% return
% 
% function [a,boarder]=getComplex(i)
%     trash=[15,15,15,17,19];
%     Boarder=[ 131, 126, 120, 113, 105;
%              -140,-140,-130,-120,-100];
%     b=[ 5.68434188608080e-14  4.52545986795450e-16 ...
%         4.43036619711279e-03  4.46228318841690e-03  4.20670442924376e-16;
%         6.75712808019284e-01  6.11845812932171e-01 ...
%         5.05042700959962e-01  4.75754824714841e-01  5.26935625152687e-01;
%         1.81008852641003e+00  2.14221265568663e+00 ...
%         2.65326854924768e+00  2.81090758079188e+00 2.64959160998887e+00];
%     c=[-2.34707286779767e-02 -2.47730137683950e-02 ...
%        -2.46857045250150e-02 -2.38912597716627e-02 -2.27766004357710e-02;
%         7.23785030253681e-01  6.63155993922102e-01 ...
%         6.10937560490686e-01  5.88477884219582e-01  5.67201960858340e-01;
%         1.46849475831187e+00  1.77424777523952e+00 ...
%         2.02988678277987e+00  2.16784432688276e+00 2.29418423275041e+00];
%     a=@(sr,li)max(trash(i),exp([log(sr),log(li),1]*[b(:,i),c(:,i)]));
%     boarder=Boarder(i);
% end
% function [a,boarder]=getReal(i)
%     trash=[14,14,16,17,19];
%     Boarder=[-951,-754,-492,-459,-347;
%               100,  95,  90,  85,  80];
%     b=[ 4.74135221075085e+00  4.46897082443911e+00  ...
%         4.36059821239678e+00  4.26390052999110e+00  4.37962073223809e+00;
%        -1.16291831007402e-02 -1.76554580857044e-02  ...
%        -1.43808682344115e-03  1.87156853583432e-02  2.06467834519798e-02;
%        -1.57331341950896e+00  2.66150937596763e-01  ...
%         1.34891825526006e+00  2.89005385937215e+00  3.36918671323585e+00];
%     c=[ 4.05472340854117e+00  3.98701095841291e+00  ...
%         4.33524202253646e+00  4.14639299189636e+00  3.63972462754257e+00;
%        -3.14595617603228e-01 -2.82468750420033e-01  ...
%        -4.50642222746083e-01 -4.94693570373209e-01 -4.15753762304312e-01;
%        -2.37322347225680e+00 -1.21325609751378e+00  ...
%        -1.47339078957474e+00  1.85656738380089e-01  2.44349225765701e+00];
%     a=@(sr,li)max(trash(i),exp(sqrt([log(sr),log(li),1]*[b(:,i),c(:,i)])));
%     boarder=Boarder(:,i);
% end
% end


%% Newton
%
% Compute the Newton interpolation polynomial in real Leja points for the 
% matrix function specified with the divided differences DD applied to the
% right hand side V of the operator A*H*V as Y=P_m(H*A)V. The nodes are in 
% [-2,2] + d/gammafoc and the interpolation stops when it reached the 
% maximal steps iter.
%
% The result is stored in Y the estimated error in NORMERREST and the
% number of steps in INFO. if the maximal number of iterations is reached
% but the desired error is not reached INFO contains -MAX_NUMBER_OF_STEPS.
%
function [y, flag, mvps] = newton_fix(h, A,lufactors,v, xi, dd, gammafoc,NewTol,iter)

flag = 0;
mvps = 0;
w = v;
y = w * dd(1);
if any(isnan(dd)) || any(isinf(dd))
    fprintf('Polynomial Leja failed: Newton table has NaN or Inf due to large tau.\n');
    return;
end
for m = 2:iter
    if lufactors.eyeB
        w = (A*w*h)/gammafoc - xi(m-1)*w;
    elseif ~lufactors.symB
        %w = (Qm*(Um\(Lm\(Pm*(A*w))))*h) / gammafoc  - xi(m-1) * w;
        w = lufactors.pQ*(lufactors.pU\(lufactors.pL\(lufactors.pP*(A*w*h)))) / gammafoc  - xi(m-1) * w;
    else
        w = lufactors.pSP*(lufactors.pLt\(lufactors.pinvD*(lufactors.pL\(lufactors.pPtS*(A*w*h)))))/gammafoc - xi(m-1)*w;
    end
    mvps = mvps + 1;
    updatevec = w * dd(m);
    y = y + updatevec;
    normwdd = norm(updatevec)/norm(y);
    if normwdd < NewTol && normwdd > eps
        %fprintf('converged at inner step %d\n',m)
        flag = 1;
        break
    end
end
% if normwdd > NewTol
%     fprintf('Polynomial Leja did not converge at inner step %d\n',m)
% end
   
end
%
% Compute the Newton interpolation polynomial in conjugate pairs of  Leja
% points for the matrix function specified with the divided differences DD
% applied to the right hand side V of the operator A*V*H as Y=P_m(A*H)V.
% The nodes are in 1i*[-2,2] + d/gammafoc and the interpolation stops when 
% it reached the maximal steps iter.
%
% The result is stored in Y the estimated error in NORMERREST and the
% number of steps in INFO. if the maximal number of iterations is reached
% but the desired error is not reached INFO contains -MAX_NUMBER_OF_STEPS.
%
function [y, flag, mvps] = newtons_fix (h, A,lufactors, v, xi, dd, gammafoc, NewTol,iter)

w = v;
y = w * real (dd(1)); % dd(1,:) should be real
%w = (Qm*(Um\(Lm\(Pm*(A*w))))*h) / gammafoc - xi(1) * w;

if lufactors.eyeB
    w = (A*w*h)/gammafoc - xi(1)*w;
elseif ~lufactors.symB
    w = lufactors.pQ*(lufactors.pU\(lufactors.pL\(lufactors.pP*(A*w*h)))) / gammafoc  - xi(1) * w;
else
    w = lufactors.pSP*(lufactors.pLt\(lufactors.pinvD*(lufactors.pL\(lufactors.pPtS*(A*w*h)))))/gammafoc - xi(1)*w;
end

flag = 0;
mvps = 1;
if any(isnan(dd)) || any(isinf(dd))
    fprintf('Polynomial Leja failed: Newton table has NaN or Inf due to large tau.\n');
    return;
end
for m=3:2:iter
    %wtilde = (Qm*(Um\(Lm\(Pm*(A*w))))*h) / gammafoc  - real (xi(m-1)) * w;
    if lufactors.eyeB
        wtilde = (A*w*h)/gammafoc - real(xi(m-1))*w;
    elseif ~lufactors.symB
        %w = (Qm*(Um\(Lm\(Pm*(A*w))))*h) / gammafoc  - xi(m-1) * w;
        wtilde = lufactors.pQ*(lufactors.pU\(lufactors.pL\(lufactors.pP*(A*w*h)))) / gammafoc  - real(xi(m-1))*w;
    else
        wtilde = lufactors.pSP*(lufactors.pLt\(lufactors.pinvD*(lufactors.pL\(lufactors.pPtS*(A*w*h)))))/gammafoc - real(xi(m-1))*w;
    end
    mvps = mvps + 1;
    updatewdd = w*real(dd(m-1))+wtilde*real(dd(m));
    y = y + updatewdd;
    normwdd = norm(updatewdd)/norm(y);
      if normwdd < NewTol && normwdd > eps
          %fprintf('converged at inner step %d\n',m)
         flag = 1;
         break
      end
    w = (A*wtilde*h) / gammafoc - real (xi(m-1)) * wtilde + ...
        imag (xi(m-1))^2 * w;    
end
% if normwdd > NewTol
%     fprintf('Polynomial Leja did not converge at inner step %d\n',m)
% end

end

%% Divided differences
%
%   Compute accurate divided differences for \exp at nodes XI(:). If XI(:)
%   are in  [-2,2]+D/GAMMA then they are in the divided differences for
%   EXP(C+GAMMA*XI) at nodes XI(:).
%   S specifies the amount of substeps, if non is provided it is computed
%   automatically.
%
%   The computation is carried out with scaling and squaring and by Taylor
%   expansion.
%   There is a second version included that is faster in Octave and my be
%   beneficial in a different MATLAB version. 
%
%   Reference: M. Caliari, Accurate evaluation of divided differences for 
%   polynomial interpolation of exponential propagators, Computing 2007, 
%   vol. 80(2): pp. 189â€?01.
%
function d = exptaylordd (xi, gamma, varargin)
l = max (gamma * (abs (xi(1:3)) + 1));
s = max(floor (2 * l),1);
s = min([s 200]);
% Edited by F. Xue to ensure that this function goes through quickly, so
% that polynomial Leja claims failure rapidly for computing exp(h*A)*b for
% a large h with a limited number of substeps (the original code would try
% hard to set up a Newton table of divided difference to high accuracy, but
% that costs a lot of time if h is large, and does not lead to convergence
% if the number of substeps is limited). In the current version, we specify 
% the number of substeps manualy.
m = length (xi);
gamma = gamma / s;
F = toeplitz(gamma .^ (0:m-1)./[1,1,cumprod(2:m-1)]);
% pointers to main diagonal
iddiag = (1:m + 1:m ^ 2);
for l = 1:16 % 16 is enough to approximate by Taylor's series with xi
    for j = 1:m-1
        F(j, j) = gamma * xi(j) * F(j, j) / l;
        for i = j+1:m
            F(j, i) = gamma * (xi(i) * F(j, i) + F(j, i - 1)) / (l + i - j);
            F(i, j) = F(i, j) + F(j, i);
        end
    end
end
F(iddiag) = exp (gamma * xi);
Fk = tril(F);
d = Fk(:, 1);
for i = 1:s-1
    d = Fk * d;
end

% The original code above is mathematically equivalent to the following
% direct evaluation of Newton's finite difference table
%
% nt_table = zeros(length(xi),length(xi));
% nt_table(:,1) = exp (gamma * xi);     % original gamma (not divided by s)
% for k = 2 : length(xi)
%     nt_table(k:end,k) = (nt_table(k:end,k-1)-nt_table(k-1:end-1,k-1))...
%         ./(xi(k:length(xi))-xi(1:length(xi)-k+1));
% end
% d = diag(nt_table);

end
% % compute amount of substeps
% if (isempty(varargin))
%   l = max (gamma * (abs (xi(1:3)) + 1));
%   s = max(floor (2 * l), 1);
% else s=varargin{1};
% end
% m = length (xi);
% gamma = gamma / s;
% Fk = zeros(m);
% 
% F = toeplitz(gamma .^ (0:m-1)./[1,1,cumprod(2:m-1)]); 
% 
% d = F(:,1);
% error = 1;
% l = 1;
% xi = xi(:).';
% 
% iddiag = (1:m + 1:m ^ 2);
% 
% while (error > 0 || F(1,1) > eps) 
%   F(iddiag) = (gamma / l) * F(iddiag) .* xi;
%   idx1 = iddiag;
%   for i = 1:m - 1
% % pointers to sup diagonals
%     idx2 = idx1(1:m - i);
%     idx1 = idx1(2:m - i + 1) - 1;
%     F(idx1) = (gamma / (l + i)) * (F(idx1) .* xi(i + 1:m) + F(idx2));
%   end  
%   F = F + triu (F).' - diag (F(iddiag));
%   error = max (abs (F(2:m, 1) - d(2:m)) ./ abs (F(2:m, 1)));
%   d = F(:, 1);
%   l = l + 1;
% end
% F(iddiag) = exp (gamma * xi);
% Fk(1:m, 1:m) = tril (F);
% d = Fk(:, 1);
% for i = 1:s-1
%   d = Fk * d;
% end
% d = d(1:m,:);
% end


%% Leja points
%
% Returns the first M (real) Leja points in [-2,2]. Only the first 200 Leja
% points are availabel --> M<=200.
%
function xi = leja(m)
xi=[0.20000000000000000E+01;
 -0.20000000000000000E+01;
  0.00000000000000000E+00;
  0.11547005389072722E+01;
 -0.13174131921908798E+01;
  0.16785083483460106E+01;
 -0.17400142991341268E+01;
 -0.61122666051585473E+00;
  0.64341522576591714E+00;
  0.18859583645359135E+01;
 -0.19053465427082013E+01;
 -0.95882465993895383E+00;
  0.14252772819409354E+01;
  0.31191873036563822E+00;
 -0.15497446836687399E+01;
  0.19589552375933266E+01;
 -0.32233053684798985E+00;
 -0.19666526192556142E+01;
  0.92274120557475403E+00;
  0.17837856425002223E+01;
 -0.11437941710547146E+01;
 -0.18251194862987536E+01;
  0.12985070496776170E+01;
 -0.15963594650138674E+00;
  0.48461307997069958E+00;
 -0.14445887774072832E+01;
  0.19853372477170090E+01;
 -0.78157903185550548E+00;
  0.15642613771121558E+01;
 -0.19881134949997448E+01;
  0.79515378837945172E+00;
 -0.16527693558289136E+01;
  0.18400967618043405E+01;
 -0.47031502482404541E+00;
  0.16237742431025759E+00;
 -0.19362973057205193E+01;
  0.19283923034097246E+01;
  0.10488606867030650E+01;
 -0.10567169290471541E+01;
  0.14988429403852255E+01;
 -0.17846584511097803E+01;
 -0.12391389737003360E+01;
  0.17305853530680211E+01;
  0.39781032979286135E+00;
 -0.69722681619218396E+00;
  0.19949158007765879E+01;
 -0.18707564721468248E+01;
  0.12274622398609041E+01;
 -0.23862075208070349E+00;
 -0.16002204184230693E+01;
  0.72036969721353816E+00;
 -0.19958685914305654E+01;
  0.16242140212548106E+01;
 -0.87412702016542365E+00;
  0.19724261744339364E+01;
  0.80404158966704414E-01;
 -0.13876623788610769E+01;
  0.13631392018454864E+01;
 -0.19529498748553720E+01;
  0.98318328054576365E+00;
 -0.40340085143054333E+00;
  0.19071556630242674E+01;
 -0.16975573834583724E+01;
  0.56064317542589448E+00;
 -0.14958571038756896E+01;
  0.18124992798591713E+01;
 -0.78233177400729326E-01;
 -0.19792196201336756E+01;
  0.11048530744290983E+01;
 -0.10091469774874193E+01;
  0.19981611763672014E+01;
 -0.54661155815383133E+00;
 -0.18490785255385729E+01;
  0.15324942269853394E+01;
  0.24047293749664256E+00;
 -0.11945330474062112E+01;
  0.18638789913676468E+01;
  0.85523120854811341E+00;
 -0.19205417848388144E+01;
  0.17047215720365889E+01;
 -0.82681093074815104E+00;
 -0.16269396953652424E+01;
  0.19454358501538969E+01;
  0.44049818260554907E+00;
 -0.19985529667308248E+01;
  0.12646226928741116E+01;
 -0.13516689386557659E+01;
 -0.27959567945593367E+00;
  0.14595518634318552E+01;
 -0.17632647864745934E+01;
  0.12100991164387512E+00;
  0.19796700401009231E+01;
 -0.65578679811245655E+00;
  0.68259629139607525E+00;
 -0.18882965699713949E+01;
  0.17587184762659824E+01;
 -0.11021736407306699E+01;
  0.15977510949799085E+01;
 -0.15224243565360993E+01;
 -0.11713917398659571E+00;
  0.10156251201074136E+01;
 -0.19920609183244737E+01;
  0.19909754389942678E+01;
 -0.91656729720428121E+00;
  0.35125666439067899E+00;
  0.13342493738720118E+01;
 -0.18057573015625830E+01;
  0.19178282223700041E+01;
 -0.12783352957127079E+01;
 -0.50636888139925262E+00;
  0.88711907622909425E+00;
 -0.19728690859025828E+01;
  0.16529163699606504E+01;
  0.36324223875185128E-01;
 -0.16772535961126018E+01;
  0.18269295872608657E+01;
  0.59828071959314566E+00;
 -0.14185314526121953E+01;
  0.11892745028984089E+01;
 -0.74004480758623958E+00;
  0.19655366634970863E+01;
 -0.19448214118826375E+01;
 -0.36473230837564752E+00;
  0.13971063858841664E+01;
 -0.17201596479296610E+01;
  0.20638660864548197E+00;
  0.19993615581122488E+01;
 -0.11694111992233767E+01;
  0.76118357247152990E+00;
 -0.19840843366776040E+01;
  0.18751188672466657E+01;
 -0.15743849287764438E+01;
 -0.19837146670101830E+00;
  0.10794685548932486E+01;
 -0.18600756919982304E+01;
  0.17448849163392661E+01;
 -0.98435277450936931E+00;
  0.52115164438684380E+00;
  0.14804723469509264E+01;
 -0.58000573293431446E+00;
 -0.19994907653772023E+01;
  0.19378502518760174E+01;
 -0.14698147331485774E+01;
 -0.38665298263104081E-01;
  0.15813962348335746E+01;
 -0.19130149913849928E+01;
  0.95174241520235947E+00;
 -0.85012679882534403E+00;
  0.17978673155225322E+01;
  0.27790115008970273E+00;
 -0.12976347570633995E+01;
  0.19884036833991652E+01;
 -0.18356760596259141E+01;
 -0.43591547045755408E+00;
  0.12458760307759087E+01;
 -0.19600030758784839E+01;
  0.82352085825767629E+00;
  0.18970546328941940E+01;
 -0.10796895476808857E+01;
 -0.16398694765802584E+01;
  0.16657857586328304E+01;
  0.37556569558272346E+00;
 -0.71851159293948874E+00;
 -0.17522623869784195E+01;
  0.19528019576324605E+01;
  0.11318200384276298E+01;
 -0.13701244963593509E+01;
 -0.13843992971968599E+00;
 -0.19941800376844832E+01;
  0.13806290736436702E+01;
  0.19967087118551494E+01;
  0.62076650375226761E+00;
 -0.12179461714907480E+01;
 -0.18961260995395672E+01;
  0.59321622300415357E-01;
  0.18515979143014407E+01;
 -0.34333214126906741E+00;
  0.15174104510155204E+01;
 -0.15363309140757493E+01;
 -0.93608754429389618E+00;
  0.17169049281230784E+01;
 -0.19295036355898909E+01;
  0.74040169612773110E+00;
 -0.63351253860727552E+00;
  0.19761830519314971E+01;
 -0.17951877483432130E+01;
  0.46286329034297946E+00;
  0.13156820343408633E+01;
 -0.19974323573365131E+01;
  0.18374074686749642E+00;
  0.17720726595862875E+01;
 -0.11234112999456882E+01;
 -0.15879222569701819E+01;
  0.10318521087330463E+01;
 -0.25886694816712452E+00;
  0.19232144716501198E+01;
 -0.17087331818172968E+01;
  0.14421895057088641E+01;
 -0.80240512714215129E+00;
 -0.19760626278667983E+01;
  0.16368140642612501E+01;
  0.54134829760904868E+00;
 -0.14046511952885661E+01;
  0.19931734850658200E+01;
 -0.52540321013695324E+00;
  0.90471021449558708E+00;
 -0.18791077893082284E+01;
 -0.57496557396576657E-01;
  0.12062363247183083E+01;
 -0.10328450798122886E+01;
  0.18914741137953366E+01;
 -0.18156421615632654E+01;
  0.25994109504469465E+00;
 -0.12592408239980593E+01;
  0.15495374423760373E+01;
  0.19997771269354820E+01;
 -0.19565437514122395E+01;
  0.66510360596327245E+00;
 -0.14828859488893442E+01;
 -0.45255027071268694E+00;
  0.18197318521595538E+01;
 -0.16661085174995318E+01;
  0.11705194181332477E+01;
  0.10248896482064303E+00;
 -0.19900489631934479E+01;
  0.19687310193550804E+01;
 -0.89439995684080076E+00;
  0.16918063576731750E+01;
 -0.21660687463806830E+00;
  0.83931153119878499E+00;
 -0.17737863422802569E+01;
 -0.13342974127878295E+01;
  0.13487982074272646E+01;
 -0.67663150865144928E+00;
  0.19417324619830225E+01;
 -0.19406612086140811E+01;
  0.41871445631709769E+00;
  0.16104042052877992E+01;
 -0.16133660507757301E+01;
  0.96821857939159217E+00;
 -0.19998216287681223E+01;
  0.18577978199406866E+01;
 -0.38481592732544789E+00;
 -0.11816558948490286E+01;
  0.33020376949941965E+00;
  0.19827672639234795E+01;
 -0.76159920845964324E+00;
 -0.18429308141537173E+01;
  0.12812943279381548E+01;
 -0.17421475485852884E-01;
 -0.14561587444279338E+01;
  0.17907669681095686E+01;
  0.77832848803746946E+00;
 -0.19696673266172788E+01;
  0.14125372999016872E+01;
 -0.99696993264214240E+00;
  0.19987753701654216E+01;
 -0.17301164448261215E+01;
  0.57968332715471649E+00;
 -0.30085287689251250E+00;
  0.10660513668629648E+01;
 -0.19247788397316214E+01;
  0.19121492714876558E+01;
 -0.56471201638588608E+00;
 -0.15101018918392919E+01;
  0.17377548010153787E+01;
  0.14268736140062088E+00;
 -0.19859808012663684E+01;
  0.15079242541034690E+01;
 -0.10909202482605411E+01;
  0.19619255586965365E+01;
 -0.98568340014749664E-01;
 -0.16873779171436061E+01;
  0.11184110878802367E+01;
  0.70159825402612641E+00;
 -0.18656326928837734E+01;
 -0.12492464317924714E+01;
  0.16449731476462048E+01;
 -0.83846533209976970E+00;
  0.18342414869889936E+01;
  0.50176886337164861E+00;
 -0.15623753554246487E+01;
 -0.19967152509394768E+01;
  0.19897363639201788E+01;
 -0.48850557480153239E+00;
  0.99903538669353953E+00;
  0.22326043560717895E+00;
 -0.19007727069998448E+01;
  0.14698034001810949E+01;
 -0.14301917806537630E+01;
  0.19331056769471811E+01;
 -0.17875649221010009E+00;
 -0.94762700509045317E+00;
  0.15730360944594053E+01;
 -0.19490628039745230E+01;
  0.87169272306901247E+00;
 -0.59757074156589451E+00;
  0.18800175148995142E+01;
 -0.17899392769268649E+01;
  0.12173334522913510E+01;
 -0.13076123189439273E+01];
xi = xi(1:m);
end
%
% Returns the first M (symmetric complex conjugate) Leja points in 
% 1Ã®*[-2,2]. Only the first 200 Leja points are availabel --> M<=200.
%
function xi = lejas(m)
% symmetric
xi = [0.00000000000000000E+00;
 -0.20000000000000000E+01;
  0.20000000000000000E+01;
 -0.11547005389072722E+01;
  0.11547005389072722E+01;
 -0.16798869581724092E+01;
  0.16798869581724092E+01;
 -0.55212517153625673E+00;
  0.55212517153625673E+00;
 -0.18860377814061269E+01;
  0.18860377814061269E+01;
 -0.86926567048810433E+00;
  0.86926567048810433E+00;
 -0.14648353835608652E+01;
  0.14648353835608652E+01;
 -0.25358605654518529E+00;
  0.25358605654518529E+00;
 -0.19608080183933749E+01;
  0.19608080183933749E+01;
 -0.17855360051319420E+01;
  0.17855360051319420E+01;
 -0.13076036318707878E+01;
  0.13076036318707878E+01;
 -0.70927408113420620E+00;
  0.70927408113420620E+00;
 -0.19857559840157573E+01;
  0.19857559840157573E+01;
 -0.15746119124852025E+01;
  0.15746119124852025E+01;
 -0.38641404790579925E+00;
  0.38641404790579925E+00;
 -0.10252007303405923E+01;
  0.10252007303405923E+01;
 -0.19238218904644926E+01;
  0.19238218904644926E+01;
 -0.11356505078523751E+00;
  0.11356505078523751E+00;
 -0.18328648391574527E+01;
  0.18328648391574527E+01;
 -0.13869875745457783E+01;
  0.13869875745457783E+01;
 -0.17299132735089220E+01;
  0.17299132735089220E+01;
 -0.79110322966313817E+00;
  0.79110322966313817E+00;
 -0.19950094737219035E+01;
  0.19950094737219035E+01;
 -0.12255720758043205E+01;
  0.12255720758043205E+01;
 -0.46869711976191569E+00;
  0.46869711976191569E+00;
 -0.16249179463353562E+01;
  0.16249179463353562E+01;
 -0.19439167904762229E+01;
  0.19439167904762229E+01;
 -0.95416908607531581E+00;
  0.95416908607531581E+00;
 -0.18233461816615804E+00;
  0.18233461816615804E+00;
 -0.18604221037194564E+01;
  0.18604221037194564E+01;
 -0.15165949191908132E+01;
  0.15165949191908132E+01;
 -0.63208883307573149E+00;
  0.63208883307573149E+00;
 -0.19757001452981686E+01;
  0.19757001452981686E+01;
 -0.10945725314996602E+01;
  0.10945725314996602E+01;
 -0.17583097219194987E+01;
  0.17583097219194987E+01;
 -0.32085293148252758E+00;
  0.32085293148252758E+00;
 -0.13487547756021807E+01;
  0.13487547756021807E+01;
 -0.19982428230270117E+01;
  0.19982428230270117E+01;
 -0.51257930118453976E-01;
  0.51257930118453976E-01;
 -0.19054515139305148E+01;
  0.19054515139305148E+01;
 -0.14300340970675731E+01;
  0.14300340970675731E+01;
 -0.91121106035049493E+00;
  0.91121106035049493E+00;
 -0.16533046750271887E+01;
  0.16533046750271887E+01;
 -0.59277547325880930E+00;
  0.59277547325880930E+00;
 -0.18114456322046653E+01;
  0.18114456322046653E+01;
 -0.12609004281591163E+01;
  0.12609004281591163E+01;
 -0.19904989930840369E+01;
  0.19904989930840369E+01;
 -0.75108525679913374E+00;
  0.75108525679913374E+00;
 -0.15471104012629255E+01;
  0.15471104012629255E+01;
 -0.19526142499017065E+01;
  0.19526142499017065E+01;
 -0.42646412321259564E+00;
  0.42646412321259564E+00;
 -0.10611747892375527E+01;
  0.10611747892375527E+01;
 -0.17068490728341663E+01;
  0.17068490728341663E+01;
 -0.21711234979315699E+00;
  0.21711234979315699E+00;
 -0.18737209418635996E+01;
  0.18737209418635996E+01;
 -0.11885139040381980E+01;
  0.11885139040381980E+01;
 -0.19695392132630480E+01;
  0.19695392132630480E+01;
 -0.83065749929225863E+00;
  0.83065749929225863E+00;
 -0.15998430128521397E+01;
  0.15998430128521397E+01;
 -0.50969282235341318E+00;
  0.50969282235341318E+00;
 -0.19329922003399136E+01;
  0.19329922003399136E+01;
 -0.14887639818661527E+01;
  0.14887639818661527E+01;
 -0.82500260421577265E-01;
  0.82500260421577265E-01;
 -0.19993753417027305E+01;
  0.19993753417027305E+01;
 -0.99031285197543495E+00;
  0.99031285197543495E+00;
 -0.17984694537070336E+01;
  0.17984694537070336E+01;
 -0.12856068973614110E+01;
  0.12856068973614110E+01;
 -0.67083107410985932E+00;
  0.67083107410985932E+00;
 -0.18467090979098812E+01;
  0.18467090979098812E+01;
 -0.29086225313969671E+00;
  0.29086225313969671E+00;
 -0.19813115154694114E+01;
  0.19813115154694114E+01;
 -0.14082128486138015E+01;
  0.14082128486138015E+01;
 -0.11248472454724956E+01;
  0.11248472454724956E+01;
 -0.17439345363407788E+01;
  0.17439345363407788E+01;
 -0.35608248193825248E+00;
  0.35608248193825248E+00;
 -0.19141794488578103E+01;
  0.19141794488578103E+01;
 -0.16393881747858858E+01;
  0.16393881747858858E+01;
 -0.89037695394600458E+00;
  0.89037695394600458E+00;
 -0.19930460273640178E+01;
  0.19930460273640178E+01;
 -0.14722764520692921E+00;
  0.14722764520692921E+00;
 -0.13303270238321043E+01;
  0.13303270238321043E+01;
 -0.18954957715477345E+01;
  0.18954957715477345E+01;
 -0.73002695721571742E+00;
  0.73002695721571742E+00;
 -0.15322007488124774E+01;
  0.15322007488124774E+01;
 -0.53093030417538989E+00;
  0.53093030417538989E+00;
 -0.19652740490686988E+01;
  0.19652740490686988E+01;
 -0.16937068017621355E+01;
  0.16937068017621355E+01;
 -0.10432560583884489E+01;
  0.10432560583884489E+01;
 -0.23104523884455920E-01;
  0.23104523884455920E-01;
 -0.19969351840405918E+01;
  0.19969351840405918E+01;
 -0.12072030165460474E+01;
  0.12072030165460474E+01;
 -0.17725812834896795E+01;
  0.17725812834896795E+01;
 -0.14479963855736195E+01;
  0.14479963855736195E+01;
 -0.44735129673309959E+00;
  0.44735129673309959E+00;
 -0.18233260155856839E+01;
  0.18233260155856839E+01;
 -0.81116207966067522E+00;
  0.81116207966067522E+00;
 -0.19387483888387484E+01;
  0.19387483888387484E+01;
 -0.97174274699549334E+00;
  0.97174274699549334E+00;
 -0.15872920982131618E+01;
  0.15872920982131618E+01;
 -0.61287752274714813E+00;
  0.61287752274714813E+00;
 -0.19880574665310511E+01;
  0.19880574665310511E+01;
 -0.13681257452441069E+01;
  0.13681257452441069E+01;
 -0.23539334410327134E+00;
  0.23539334410327134E+00;
 -0.18671540015844768E+01;
  0.18671540015844768E+01;
 -0.11699972963064553E+01;
  0.11699972963064553E+01;
 -0.16671457859178513E+01;
  0.16671457859178513E+01;
 -0.19997819819489917E+01;
  0.19997819819489917E+01;
 -0.33879692067760181E+00;
  0.33879692067760181E+00;
 -0.93193510956905001E+00;
  0.93193510956905001E+00;
 -0.17193556727601473E+01;
  0.17193556727601473E+01;
 -0.19564712782337348E+01;
  0.19564712782337348E+01;
 -0.68810730108543350E+00;
  0.68810730108543350E+00;
 -0.15017426944300611E+01;
  0.15017426944300611E+01;
 -0.13044591694027682E+00;
  0.13044591694027682E+00;
 -0.19189501051074362E+01;
  0.19189501051074362E+01;
 -0.12449462828393627E+01;
  0.12449462828393627E+01;
 -0.18400500193660232E+01;
  0.18400500193660232E+01;
 -0.48936152483709394E+00;
  0.48936152483709394E+00;
 -0.11091391592482718E+01;
  0.11091391592482718E+01;
 -0.19785744492297805E+01;
  0.19785744492297805E+01;
 -0.15609602234431856E+01;
  0.15609602234431856E+01;
 -0.77167718033496546E+00;
  0.77167718033496546E+00;
 -0.17656073238286094E+01;
  0.17656073238286094E+01;
 -0.66682609227377843E-01;
  0.66682609227377843E-01;
 -0.14186948906067163E+01;
  0.14186948906067163E+01;
 -0.19960372860004174E+01;
  0.19960372860004174E+01;
 -0.57293757345009144E+00;
  0.57293757345009144E+00;
 -0.18807232956436515E+01;
  0.18807232956436515E+01;
 -0.10084464422900927E+01;
  0.10084464422900927E+01;
 -0.16133381879660718E+01;
  0.16133381879660718E+01;
 -0.40436849775440142E+00;
  0.40436849775440142E+00;
 -0.19481910481270022E+01;
  0.19481910481270022E+01;
 -0.12967412655742256E+01;
  0.12967412655742256E+01;
 -0.85067147549021660E+00;
  0.85067147549021660E+00;
 -0.18048858729861337E+01;
  0.18048858729861337E+01;
 -0.27280542126864110E+00;
  0.27280542126864110E+00;
 -0.19727206400250756E+01;
  0.19727206400250756E+01;
 -0.14763121584831804E+01;
  0.14763121584831804E+01;
 -0.10782667121173539E+01;
  0.10782667121173539E+01;
 -0.19006344714612742E+01;
  0.19006344714612742E+01;
 -0.16697056908903357E+00;
  0.16697056908903357E+00;
 -0.16869141167379311E+01;
  0.16869141167379311E+01;
 -0.65170298750894351E+00;
  0.65170298750894351E+00;
 -0.19988653296263976E+01;
  0.19988653296263976E+01;
 -0.13583663084152575E+01;
  0.13583663084152575E+01;
 -0.11407374157356993E+01;
  0.11407374157356993E+01;
 -0.17373227019455160E+01;
  0.17373227019455160E+01;
 -0.19837088666911726E+01;
  0.19837088666911726E+01;
 -0.37147733137072370E+00;
  0.37147733137072370E+00;
 -0.12721419615906764E+01;
  0.12721419615906764E+01];
xi = xi(1:m); 
end
