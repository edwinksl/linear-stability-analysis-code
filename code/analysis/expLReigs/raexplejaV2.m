function [expmv,flag,mvps] = raexplejaV2(h,A,B,vb,vn,substeps,ell,reltol,lufactors)

% raexplejaV2 is the single-pole rational Leja method for approximating
% exp(h*inv(B)*A)*vb to the prescribed relative change tolerance
%
% Authors of the MATLAB code:  Minghao W. Rostami and Fei Xue
% Last date of chage: May 8th, 2018
%
% Input:
% h         positive scaling factor in exp(h*inv(B)*A)
% A,B       the matrix or matrix pencil of interest
% vb        the right-hand side of exp(h*inv(B)*A)*vb
% vn        the shared null space basis of (A,B); set to zero for regulars
% substeps  the number of substeps for approximating exp(h*inv(B)*A)*F
%           in each substep, we approximate exp(tau*inv(B)*A)*K, where 
%           tau = h/substeps and K is the result from the previous substep
% ell       the maximum # of Leja points used in each substep (50 is ideal)
% reltol    the relative change stopping criterion for all expmv algorithms
% lufactors precomputed LU factors of tau*A-(a0+c0)*B
%
% Output:
% expmv     approximation to exp(h*inv(B)*A)*vb
% flag      whether the algorithm has converged upon exit; true or false
% mvps      the number of linear solves involving tau*A-(a0+c0)*B
% 
% Comments:
% Let z = a0*(w-b1)/(w+b2)+c0, where b1 and b2 are sufficiently close to 2.
% Since w \in [-b2,b1], we have z \in (-infty,c0]. Therefore exp(z) = g(w),
% where w = [(b2*c0-b1*a0)-b2*z]/[z-(a0+c0)], is a function of w defined on
% [-b2,b1] that can be well approximated by a polynomial p_m(w), i.e., an
% (m,m)-type rational function of z. We use Leja points {w_k} on [-b2,b1]
% to construct the polynomial interpolation p_m(w) of g(w)(hence the name
% "Rational Leja"). To approximate exp(h*inv(B)*A)*vb, replacing z with
% tau*inv(B)*A gives w = (tau*A-(a0+c0)*B)\((b2*c0-b1*a0)*B-b2*tau*A).
% We need to perform the action of the Leja polynomial p_m(w) on vb. 
% This is a single-pole rational method. Newton polynomial based on Leja 
% points is the way we do so.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESSED 
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABI-
% LITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT 
% SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES 
% OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
% ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
% OTHER DEALINGS IN THE SOFTWARE.

warning('off','MATLAB:nearlySingularMatrix');
nt_table = zeros(ell,ell);

flags = zeros(substeps,1);
regularvn = (rank(vn,sqrt(eps)) == size(vn,2));

pmnorms = zeros(ell-1,1);
tau = h/substeps;

%%% IMPORTANT: a0 and c0 MUST BE identical in expcomp.m and here
a0 = 50;	c0 = 0;
    
if ~exist('lufactors','var')
    [L,U,P,Q,R] = lu(tau*A-(a0+c0)*B,0.25);
    lufactors.rL = L;   lufactors.rU = U;   lufactors.rP = P;
    lufactors.rQ = Q;   lufactors.rR = R;
end

b1 = (2+sqrt(eps)/10);  b2 = (2+sqrt(eps)/10);
sh_1 = b2*c0-b1*a0;     sh_2 = b2;

fun = @(s)exp(a0*(s-b1)./(s+b2)+c0);
xi = zeros(ell,1);
load('xi.mat','xi');
xi = xi(1:ell);
nt_table(:,1) = fun(xi);
for k = 2 : ell
    nt_table(k:end,k) = (nt_table(k:end,k-1)-nt_table(k-1:end-1,k-1))...
        ./(xi(k:length(xi))-xi(1:length(xi)-k+1));
end

expmv = vb;
mvps = 0;
for iter = 1 : substeps
    rm = expmv;     expmv = nt_table(1,1)*expmv;
    for m = 1 : ell-1
        rm = lufactors.rQ*(lufactors.rU\(lufactors.rL\(lufactors.rP*(lufactors.rR\(B*(rm*sh_1)-A*(rm*(tau*sh_2)))))))-xi(m)*rm;
        if regularvn
            rm = rm - vn*((vn'*vn)\(vn'*rm));
            rm = rm - vn*((vn'*vn)\(vn'*rm));
        end
        mvps = mvps + 1;
        delta_pm = nt_table(m+1,m+1)*rm;
        expmv = expmv + delta_pm;
        pmnorms(m) = norm(expmv);
        rel_change = norm(delta_pm)/pmnorms(m);
        if rel_change < reltol && rel_change > eps
            flags(iter) = 1;
            break;
        end
    end
end
flag = (nnz(flags) >= substeps*1/2);
if regularvn
    expmv = expmv - vn*((vn'*vn)\(vn'*expmv));
    expmv = expmv - vn*((vn'*vn)\(vn'*expmv));
end

warning('on','MATLAB:nearlySingularMatrix');
return