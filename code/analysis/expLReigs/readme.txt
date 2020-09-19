INTRODUCTION

expm_LR_eigs is a MATLAB software package for reliable computation of a few eigenvalues of largest real part (rightmost) of a large sparse matrix A or a matrix pencil (A,B). It is based on the exponential transformation A --> exp(h*inv(B)*A) (h > 0), which maps the rightmost eigenvalues of A to the dominant (largest in modulus) eigenvalues of exp(h*inv(B)*A). The matrix exponential itself is not formed explicitly; instead, only its action on vectors is needed to invoke MATLAB's eigs.

Please read the documentation of expm_LR_eigs.m for more details about this package.

Software Author: M. W. Rostami and F. Xue

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESSED 
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABI-
% LITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT 
% SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES 
% OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
% ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
% OTHER DEALINGS IN THE SOFTWARE.

We also provide a few examples for users to run. Once you download this package and the example problems, run the following commands

load tols4000.mat;
[V,D,eigres] = expm_LR_eigs(0.2,A,B,5,1e-9,[1 1],EXTREIGS);

load tubular_10000.mat; % rightmost eigenvalues have positive real part
[V,D,eigres] = expm_LR_eigs(1,A,B,5,1e-9,[1 1],EXTREIGS);

load cryg10000.mat; % rightmost eigenvalues have positive real part
[V,D,eigres] = expm_LR_eigs(1,A,B,5,1e-9,[1 1],EXTREIGS);

load aerofoil.mat; % we can speed up the rational Leja method by a factor of 100
[V,D,eigres] = expm_LR_eigs(8,A,B,5,1e-9,[1 100],EXTREIGS);

load obstacle_350.mat;
[V,D,eigres] = expm_LR_eigs(2,A,B,5,1e-9,[1 1],EXTREIGS);

load cavity_7800.mat; % A and B share a one-dimensional null space 
[V,D,eigres] = expm_LR_eigs(5,A,B,5,1e-9,[1 1],EXTREIGS,vn);

NOTE: for new problems, it might not be easy to obtain a spectrum estimate for polynomial Leja (see EXTREIGS for existing problems); you may simply run a command like

[V,D,eigres] = expm_LR_eigs(2,A,B,5,1e-9,[1 1]); 

If spectrum estimate is not provided or is not relatively accurate, polynomial Leja is most probably not as efficient as it could be. But for problems of wide spectrum, this method (and standard Arnoldi) is not competitive anyway.

There are also a few artificial problems (names ending with "_atf") that impose challenges to traditional methods because their rightmost eigenvalues have large imaginary part. The RD method is typically most efficient for computing the matrix exponentials for these problems. Run the artificial problems to see the reliability of this package based on the exponential transformation; for example:

load tubular_10000_atf.mat;
[V,D,eigres] = expm_LR_eigs(1,A,B,5,1e-9,[1 1],EXTREIGS);

If you observe any bugs or have questions, please contact M. W. Rostami at mwrostam@syr.edu or F. Xue at fxue@clemson.edu.

