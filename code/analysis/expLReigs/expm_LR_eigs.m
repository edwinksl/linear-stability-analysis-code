function [V,D,eigres] = expm_LR_eigs(h,A,B,p,expmvtol,acclejas,EXTREIGS,vn)

% expm_LR_eigs is the main file for computing a few eigenvalues of the
% largest real part (rightmost) of a large sparse matrix A or matrix pencil
% (A,B), by mapping these desired eigenvalues to the dominant eigenvalues
% of exp(h*inv(B)*A) (h > 0), called the exponential transformation.

% Authors of the MATLAB code:  Minghao W. Rostami and Fei Xue
% Last date of chage: May 8th, 2018

% Input:
% h         positive scaling factor in exp(h*inv(B)*A); a larger h leads to
%           faster convergence rate of the eigensolver (outer) iteration,  
%           but it is more expensive to approximate exp(h*inv(B)*A)*v;
%           *** if the rightmost eigenvalues are close in real part, use **
%           *** a larger h (a rule of thumb here: if these eigenvalues'  **
%           *** real part differ at level 1e-2, h = 1 is sufficient; if  **
%           *** at 1e-3, h = 10 may be appropriate, and so on)           **
% A,B       the matrix A or matrix pencil (A,B) of interest; to explore
%           a matrix A alone, let B = speye(size(A))
%
%           *** The input matrix (pencil) should have vast majority of  ***
%           *** eigenvalues of negative real part, and may have a small ***
%           *** portion of eigenvalues of small positive real part      ***
%
% p         the number of rightmost eigenvalues wanted; p <= 10 recommended
% expmvtol  the relative change stopping criterion for all algorithms that
%           approximate exp(h*inv(B)*A)*v; by default, the relative eigen
%           tolerance is 10 times expmvtol
% acclejas  the factors that enforce a larger substep size for polynomial 
%           Leja and rational Leja; [1 1] by default
% EXTREIGS  the spectral estimate of A or (A,B); needed for poly Leja alone
%           EXTREIGS.SR is a lower bound on the real part of the leftmost
%           eigenvalue, EXTREIGS.LR is an upper bound on the real part of
%           the rightmost eigenvalue (typically set as 0), and EXTREIGS.LI2
%           is the _square_ of an upper bound on the largest imaginary part
%           of all eigenvalues of A or (A,B)
% vn        the shared null space of A and B
%
% Output:
% V and D   the desired rightmost eigenvectors and associated eigenvalues
% eigres    relative residual norm of the desired eigenpairs
%
% Comments:
% Four algorithms for approximating exp(h*inv(B)*A)*v are compared and the
% fastest one is chosen to work with eigs. These algorithms are: polynomial
% Leja, Arnoldi, rational Leja, and the RD-rational method. Usually but not
% always, RD with an appropriate shift s > 0 is the most efficient method,
% and Arnoldi is almost always the slowest one. 
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESSED 
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABI-
% LITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT 
% SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES 
% OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
% ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
% OTHER DEALINGS IN THE SOFTWARE.

% choose whether to use unrestarted Arnoldi (0) or eigs (1); 1 by default
eigsflag = (1 > 0);

n = length(A);
vb = ones(n,1);     vb = vb/norm(vb);
if nargin < 8
    vn = zeros(n,1);
    if nargin < 7
        [~,~,~,EXTREIGS] = explejaV2(2^(-5),A,B,1,10,1e-3,vb);
        fprintf('\nATTENTION: spectrum estimate needed by not provided for P-Leja.\n');
        fprintf('\t   This method may not be as efficient as it could be.\n');
        if nargin < 6
            acclejas = [1 1];
            if nargin < 5
                expmvtol = 1e-9;
            end
        elseif length(acclejas) ~= 2 || ~all(acclejas >= 1)
            error('An array of 2 scalars >= 1 is needed for input ''acclejas''.\n');
        end
    end
end

if p > 10
    fprintf('p = %d rightmost eigenpairs sought. This may take a longer while.\n',p);
    fprintf('We usually let p <= 5 to get the computation done relatively quickly.\n');
end

[afun,timeperexpmv] = expcomp(h,A,B,vb,expmvtol,acclejas,EXTREIGS,vn);
opts.disp = 0;	opts.testtol = 2.5e-2;    opts.v0 = vb;
opts.p = max([25 3*p]);     arnsteps = opts.p;      opts.maxit = 1;
maxhtest = 6;
for ii = 1 : maxhtest
    fprintf('Testing if h = %.2e leads to acceptable rate of convergence for the eigensolver.\n',h);
    [V,~,failflag,arn_iter] = arnoldi(afun,vb,arnsteps,p,opts.testtol);
    
    if failflag
        fprintf('\nEigensolver (outer) iteration converges too slow with h = %.2e.\n',h);
        h = h*2;
        fprintf('Doubling h to %.2e to speed up eigensolver (outer) iteration.\n',h);
        fprintf('Comparing algorithms for the action of matrix exponentials again.\n');
        fprintf('If doubling h happens a few times, consider a slight increase of the shift for RD.\n');
        [afun,timeperexpmv] = expcomp(h,A,B,vb,expmvtol,acclejas,EXTREIGS,vn);
    else
        fprintf('h = %.2e leads to eigensolver''s satisfactory rate of convergence -- moving on.\n\n',h);
        opts.disp = 0;	opts.tol = max([min([10*expmvtol 1e-7]) 1e-14]);
        opts.p = max([25 3*p]);   opts.maxit = 25;    opts.v0 = vb;
        arnsteps = opts.p + (opts.p-p)*(opts.maxit-1);
        restart_estimate = log(opts.tol)/log(opts.testtol^(opts.p/arn_iter));
        fprintf('\nESTIMATE: the eigensolver could take about %.2f secs but might run more quickly.\n\n',...
            restart_estimate*opts.p*timeperexpmv);
        if eigsflag
            fprintf('eigs started for computing the rightmost %d eigenvalues to tol %.2e.\n',p,opts.tol);
            fprintf('NOTE: if eigs complains about the starting vector (MATLAB R2017+), increase\n');
            fprintf('      the corresponding tolerance in the KrylovSchur function in eigs.\n');
            tic;
            [V,~,~] = eigs(afun,length(A),p,'lm',opts);
            telapsed = toc;
            fprintf('eigs took %.2f secs.\n',telapsed);
        else
            fprintf('Arnoldi started for computing the rightmost %d eigenvalues to tol %.2e.\n',p,opts.tol);
            tic;
            [V,~,~] = arnoldi(afun,vb,arnsteps,p,opts.tol);
            telapsed = toc;
            fprintf('Arnoldi took %.2f secs.\n',telapsed);
        end
        break;
    end
end
if ii == maxhtest && failflag
    fprintf('h has been doubled %d times, but eigensolver still exhibits slow convergence.\n',maxhtest);
    fprintf('Choose a large h > %d and try again (or this problem is too hard to solve).\n',h);
end
[V,D,eigres] = postproc(A,B,V,'LR');

end