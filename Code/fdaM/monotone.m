function [Wfd, beta, Fstr, iternum, iterhist] = ...
          monotone(x, y, wt, Wfd, zmat, Lfdobj, lambda, ...
                   conv, iterlim, active, dbglev)
%  Smooths the relationship of Y to X by fitting a monotone function
%     of the form
%                   f(x) = b_0 + b_1 D^{-1} exp W(x)
%     where  W  is a function defined over the same range as X,
%  W + ln b_1 = log Df and w = D W = D^2f/Df.
%  The constant term b_0 in turn can be a linear combinations of covariates:
%     b_0 = zmat * c
%  The fitting criterion is penalized mean squared error:
%    PENSSE(lambda) = \sum [y_i - f(x_i)]^2 +
%                     \lambda * \int [L W(x)]^2 dx
%  The function W(x) is expanded by the basis in functional data object
%    Wfd.   The coefficients of this expansion are called "coefficients"
%    in the comments, while the b's are called "regression coefficients"
%
%  Arguments are ...
%  X       ...  argument value array
%  Y       ...  function value array (the values to be fit)
%  WT      ...  a vector of weights
%  WFD     ...  functional data object for W.  It's coefficient array
%               has a single column, and these are the starting values
%               for the iterative minimization of mean squared error.
%  ZMAT    ...  a matrix of covariate values for the constant term.
%               It defaults to a column of one's;
%  LFDOBJ  ...  linear differential opr defining roughness penalty to
%               be applied to W.  This may be either a functional data
%               object defining a linear differential operator, or a
%               nonnegative integer.  If the latter, it specifies the
%               order of derivative to be penalized.
%               LFDOBJ = 1 by default, corresponding to L = D.
%  LAMBDA  ...  smoothing parameter determining the amount of penalty,
%               0 by default.
%  CONV    ...  convergence criterion, 0.0001 by default
%  ITERLIM ...  maximum number of iterations, 20 by default
%  ACTIVE  ... indices among 1:NBASIS of parameters to optimize 
%  DBGLEV  ...  Controls the level of output on each iteration.  If 0,
%               no output, if 1, output at each iteration, if higher, output
%               at each line search iteration. 1 by default.
%
%  Returns are:
%  WFD     ...  Functional data object for W.  It's coefficient vector
%               contains the optimized coefficients.
%  BETA    ...  The regression coefficients b_0 and b_1.
%  FNEW    ...  The minimized penalized mean squared error.
%  MSG     ...  The mean of the squares of the components of the gradient
%  ITERNUM ...  The total number of iterations taken.
%  ITERHIST...  An array with number of rows equal to number of iterations,
%               and columns corresponding to it. number, function value,
%               mean squared gradient, and beta values.

%  Last modified 30 January 2003

%  make sure that both X and Y are column vectors

x = x(:);
y = y(:);

%  check WFD

if ~isa_fd(Wfd)
    error('WFD is not a functional data object.');
end

nobs     = length(x);              %  number of observations
oneobs   = ones(nobs,1);
cvec     = getcoef(Wfd);           %  initial coefficients
basisobj = getbasis(Wfd);          %  basis for Wfd
nbasis   = getnbasis(basisobj);    %  no. basis functions
type     = getbasistype(basisobj); %  type of basis

%  initialize arguments that are not included

if nargin < 11
    dbglev = 1;
end
if nargin < 10
    active=1:nbasis;
end
if nargin < 9
    iterlim = 20;
end
if nargin < 8
    conv = 0.0001;
end
if nargin < 7
    lambda  = 0;
end
if nargin < 6
    Lfdobj = int2Lfd(1);
end
if nargin < 5 | zmat == 0 | zmat == NaN
    zmat = oneobs;
end

%  check LFDOBJ

Lfdobj = int2Lfd(Lfdobj);

%  check some arguments

if any(wt) < 0
    error('One or more weights are negative.');
end
zdim = size(zmat);
if zdim(1) ~= nobs
    error('First dimension of ZMAT not correct.')
end

%  set up some variables

ncov    = zdim(2);                %  number of covariates
ncovp1  = ncov + 1;
wtroot  = sqrt(wt);
wtrtmt  = wtroot*ones(1,ncovp1);
yroot   = y.*wtroot;
climit  = 100.0.*([-1; 1]*ones(1,nbasis));
inactive  = ones(1,nbasis);
inactive(active) = 0;
inactive  = find(inactive);
ninactive = length(inactive);

%  set up cell for storing basis function values

JMAX = 15;
basiscell = cell(1,15);

%  initialize matrix Kmat defining penalty term

if lambda > 0
    Kmat = lambda*eval_penalty(basisobj, Lfdobj);
else
    Kmat  = zeros(nbasis,nbasis);
end

%  Compute initial function and gradient values

[Fstr, beta, Dyhat] = fngrad(y, x, zmat, wt, Wfd, ...
    lambda, Kmat, basiscell, inactive);

%  compute the initial expected Hessian

hessmat = hesscal(beta, Dyhat, wtroot, lambda, Kmat, inactive);

%  evaluate the initial line search direction vector
[deltac, cosangle] = linesearch(Fstr, hessmat, dbglev);

%  initialize iteration status arrays

iternum = 0;
status = [iternum, Fstr.f, Fstr.norm, beta'];
if dbglev >= 1
    fprintf('\nIter.   PENSSE   Grad Length Intercept   Slope\n')
    fprintf('%3.f %10.4f %10.4f %10.4f %10.4f\n', ...
        [status(1:4),beta(ncovp1)]);
end
if dbglev > 2
    for i = 1:nbasis, fprintf('%10.4f%', cvec(i)); end
    fprintf('\n');
    for i = 1:3,      fprintf('%10.4f%', beta(i)); end
    fprintf('\n');
end
iterhist = zeros(iterlim+1,length(status));
iterhist(1,:)  = status;
if iterlim == 0, return;  end

%  -------  Begin main iterations  -----------

MAXSTEPITER = 10;
MAXSTEP = 100;
trial   = 1;
reset   = 0;
linemat = zeros(3,5);
betaold = beta;
cvecold = cvec;
Foldstr = Fstr;
dbgwrd  = dbglev >= 2;
%  ---------------  beginning of optimization loop  -----------
for iter = 1:iterlim
    iternum = iternum + 1;
    %  initialize logical variables controlling line search
    dblwrd = [0,0];  limwrd = [0,0];  stpwrd = 0;  ind = 0; ips = 0;
    %  compute slope at 0 for line search
    linemat(2,1) = sum(deltac.*Fstr.grad);
    %  normalize search direction vector
    sdg          = sqrt(sum(deltac.^2));
    deltac       = deltac./sdg;
    linemat(2,1) = linemat(2,1)/sdg;
    % initialize line search vectors
    linemat(:,1:4) = [0; linemat(2,1); Fstr.f]*ones(1,4);
    stepiter  = 0;
    if dbglev >= 2
        fprintf('                 %3.f %10.4f %12.6f %12.6f\n', ...
            [stepiter, linemat(:,1)']);
    end
    %  return with error condition if initial slope is nonnegative
    if linemat(2,1) >= 0
        if dbglev >= 2, disp('Initial slope nonnegative.'); end
        ind = 3;
        break;
    end
    %  return successfully if initial slope is very small
    if linemat(2,1) >= -1e-7;
        if dbglev >= 2, disp('Initial slope too small'); end
        ind = 0;
        break;
    end
    %  first step set to trial
    linemat(1,5)  = trial;
    %  -------  Begin line search iterations  -----------
    cvecnew = cvec;
    Wfdnew  = Wfd;
    for stepiter = 1:MAXSTEPITER
        limflg = 0;
        %  check that step size does not go beyond limits on parameters
        [linemat(1,5), ind, limwrd] = ...
            stepchk(linemat(1,5), cvec, deltac, limwrd, ind, ...
            climit, active, dbgwrd);
        if ind == 1, break; end  % break of limit hit twice in a row
        if linemat(1,5) <= 1e-7
            %  Current step size too small ... terminate
            if dbglev >= 2
                fprintf('Stepsize too small: %15.7f\n', linemat(1,5));
            end
            if limflg, ind = 1; else ind = 4; end
            break;
        end
        %  compute new function value and gradient
        cvecnew = cvec + linemat(1,5).*deltac;  %  update coefficients
        Wfdnew  = putcoef(Wfd, cvecnew);   %  update function W
        [Fstr, beta, Dyhat] = fngrad(y, x, zmat, wt, Wfdnew, ...
            lambda, Kmat, basiscell, inactive);
        linemat(3,5) = Fstr.f;
        %  compute new directional derivative
        linemat(2,5) = sum(deltac.*Fstr.grad);
        if dbglev >= 2
            fprintf('                 %3.f %10.4f %12.6f %12.6f\n', ...
                [stepiter, linemat(:,5)']);
        end
        %  compute next line search step, also testing for convergence
        [linemat, ips, ind, dblwrd] = ...
            stepit(linemat, ips, ind, dblwrd, MAXSTEP, dbglev);
        trial  = linemat(1,5);
        if trial == MAXSTEP, break; end
        %  ind == 0 means convergence
        if ind == 0 | ind == 5, break; end
    end
    %  -------  End of line search iterations  -----------
    cvec = cvecnew;
    Wfd  = Wfdnew;
    %  check that function value has not increased
    if Fstr.f > Foldstr.f
        %  if it has, terminate iterations with a warning
        if dbglev >= 2
            fprintf('Criterion increased:');
            fprintf('%10.4f %10.4f\n',[Foldstr.f, Fstr.f]);
        end
        %  reset parameters and fit
        beta   = betaold;
        cvec   = cvecold;
        Wfd    = putcoef(Wfd, cvecold);
        Fstr   = Foldstr;
        deltac = -Fstr.grad;
        if dbglev > 2
            for i = 1:nbasis, fprintf('%10.4f%', cvec(i)); end
            fprintf('\n');
            for i = 1:3,      fprintf('%10.4f%', beta(i)); end
            fprintf('\n');
        end
        if reset == 1
            %  This is the second time in a row that this
            %     has happened ...  quit
            if dbglev >= 2
                fprintf('Reset twice, terminating.\n');
            end
            return;
        else
            reset = 1;
        end
    else
        %  function value has not increased,  check for convergence
        if abs(Foldstr.f-Fstr.f) < conv
            status = [iternum, Fstr.f, Fstr.norm, beta'];
            iterhist(iter+1,:) = status;
            if dbglev >= 1
                fprintf('%3.f %10.4f %10.4f %10.4f %10.4f\n', ...
                    [status(1:4),beta(ncovp1)]);
            end
            break;
        end
        %  update old parameter vectors and fit structure
        cvecold = cvec;
        betaold = beta;
        Foldstr = Fstr;
        %  update the expected Hessian
        hessmat = hesscal(beta, Dyhat, wtroot, lambda, Kmat, inactive);
        %  update the line search direction vector
        [deltac, cosangle] = linesearch(Fstr, hessmat, dbglev);
        reset = 0;
    end
    %  store iteration status
    status = [iternum, Fstr.f, Fstr.norm, beta'];
    iterhist(iter+1,:) = status;
    if dbglev >= 1
        fprintf('%3.f %10.4f %10.4f %10.4f %10.4f\n', ...
            [status(1:4),beta(ncovp1)]);
    end
end

%  ----------------------------------------------------------------

function [deltac, cosangle] = linesearch(Fstr, hessmat, dbglev)
deltac = -symsolve(hessmat,Fstr.grad);
cosangle  = -sum(Fstr.grad.*deltac)./sqrt(sum(Fstr.grad.^2)*sum(deltac.^2));
if dbglev >= 2
    fprintf('Cos(angle) = %8.4f\n', cosangle);
end
if cosangle < 1e-7
    if dbglev >=2, fprintf('\nCosine of angle too small\n'); end
    deltac = -Fstr.grad;
end

%  ----------------------------------------------------------------

function [Fstr, beta, Dyhat] = fngrad(y, x, zmat, wt, Wfd, ...
    lambda, Kmat, basiscell, inactive)
ncov   = size(zmat,2) + 1;
nrep   = size(y,2);
nobs   = length(x);
cvec   = getcoef(Wfd);
nbasis = size(cvec,1);
f      =    monfn(x, Wfd, basiscell);
Dyhat  = mongrad1(x, Wfd, basiscell);
xmat   = [zmat,f];
Dxmat  = zeros([nobs,ncov,nbasis]);
Dxmat(:,ncov,:) = Dyhat;
wtroot = sqrt(wt);
wtrtmt = wtroot*ones(1,ncov);
yroot  = y.*wtroot;
xroot  = xmat.*wtrtmt;
ninactive = length(inactive);
%  compute regression coefs.
beta = LSfit(yroot, zmat, f, wtrtmt);
%  update fitted values
yhat   = xmat*beta;
%  update residuals and function values
res    = y - yhat;
Fstr.f = mean(res.^2.*wt);
grad   = zeros(nbasis,1);
for j=1:nbasis
    Dxroot  = squeeze(Dxmat(:,:,j)).*wtrtmt;
    yDx     = yroot'*Dxroot*beta;
    xDx     = xroot'*Dxroot;
    grad(j) = beta'*(xDx+xDx')*beta - 2*yDx;
end
Fstr.grad = grad/nobs;
if lambda > 0
    Fstr.grad = Fstr.grad + 2 .* Kmat * cvec;
    Fstr.f    = Fstr.f + cvec' * Kmat * cvec;
end
if ninactive > 0, Fstr.grad(inactive) = 0; end
Fstr.norm = sqrt(sum(Fstr.grad.^2));   %  gradient norm

%  ----------------------------------------------------------------

function hessmat = hesscal(beta, Dyhat, wtroot, lambda, Kmat, inactive)
nbet = length(beta);
[nobs, nbasis] = size(Dyhat);
temp = beta(nbet).*Dyhat;
temp = temp.*(wtroot*ones(1,nbasis));
hessmat = 2.*temp'*temp./nobs;
%  adjust for penalty
if lambda > 0, hessmat = hessmat + 2.*Kmat; end
%  adjust for inactive coefficients
ninactive = length(inactive);
if ninactive > 0
    hessmat(inactive,:    ) = 0;
    hessmat(:    ,inactive) = 0;
    hessmat(inactive,inactive) = eye(ninactive);
end

%  ----------------------------------------------------------------

function beta = LSfit(yroot, zmat, f, wtrtmt)
xmat  = [zmat,f];
ncol  = size(xmat,2);
xroot = xmat.*wtrtmt;
[Q, R, E] = qr(xroot);
tol       = size(xmat,1)*3e-16*abs(R(1,1));
y         = Q'*yroot;
for j=1:ncol
    if abs(R(j,j)) < tol
        if R(j,j) < 0, R(j,j) = -tol;  else R(j,j) = tol; end
    end
end
beta      = R\y;
beta      = E*beta;


