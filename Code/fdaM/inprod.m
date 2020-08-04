function  ss = inprod(fdstr1, fdstr2, Lfdobj1, Lfdobj2, rng, wtfd)
%  INPROD   Computes matrix of inner products of functions.
%    If both functions have the same B-spline basis, both Lfdobj's are
%       numeric, there is no wtfd argument, and the ranges are the same,
%       the inner products are exact, and computed by inprod_bspline.
%    Otherwise, by numerical integration using Romberg integration 
%       with the trapezoidal rule.

%  Arguments:
%  FD1STR and FDSTR2 ...  these may be either functional data or basis function
%                    objects.  In the latter case, a functional data object
%                    is created from a basis function object by using the
%                    identity matrix as the coefficient matrix.
%                    Both functional data objects must be univariate.
%                    If  inner products for multivariate objects are needed,
%                    use a loop and call inprod(fdstr1(i),fdstr2(i)).
%  LFDOBJ1 and LFDOBJ2 ...  linear differential operators for inner product for
%                    FD1 and FD2, respectively.  Default to 0.  
%  RNG  ...  Limits of integration
%  WTFD ...  A functional data object defining a weight
%  Return:
%  A NREP1 by NREP2 matrix SS of inner products for each possible pair
%  of functions.

%  last modified 16 January 2003

%  set up default values of arguments

if nargin < 6, wtfd = 0;               end
if nargin < 4, Lfdobj2 = int2Lfd(0);   end
if nargin < 3, Lfdobj1 = int2Lfd(0);   end

%  check LFDOBJ1 and LFDOBJ2

Lfdobj1 = int2Lfd(Lfdobj1);
Lfdobj2 = int2Lfd(Lfdobj2);
  
%  set constants determining Richardson extrapolation
  
EPS  = 1e-4;  %  convergence criterion
JMAX = 15;    %  maximum number of iterations
JMIN =  5;    %  minimum number of iterations

%  check    WTFD

if  isa_fd(wtfd)
    coefw = getcoef(wtfd);
    coefd = size(coefw);
    if coefd(2) > 1
        error('Argument WTFD is not a single function');
    end
end
  
%  check FDSTR1

fdclass = 1;
if  isa_fd(fdstr1)
    coef1 = getcoef(fdstr1);
elseif isa_basis(fdstr1)
    coef1  = eye(getnbasis(fdstr1));
    temp1  = fd(coef1, fdstr1);
    fdstr1 = temp1;
else
    fdclass = 0;
end

%  check FDSTR2

if  isa_fd(fdstr2)
    coef2 = getcoef(fdstr2);
elseif isa_basis(fdstr2)
    coef2  = eye(getnbasis(fdstr2));
    temp2  = fd(coef2, fdstr2);
    fdstr2 = temp2;
else
    fdclass = 0;
end

if ~fdclass
    error (['The two first arguments are not', ...
            ' functional data or basis objects.']);
end

%  determine NREP1 and NREP2, and check for common range

coefd1 = size(coef1);
coefd2 = size(coef2);
if length(coefd1) > 2 | length(coefd2) > 2
    error('Functional data objects must be univariate');
end
nrep1 = coefd1(2);
nrep2 = coefd2(2);
basisobj1 = getbasis(fdstr1);
basisobj2 = getbasis(fdstr2);
type1  = getbasistype(basisobj1);
type2  = getbasistype(basisobj1);
range1 = getbasisrange(basisobj1);
range2 = getbasisrange(basisobj2);
if nargin < 5, rng = range1; end
if rng(1) < range1(1) | rng(2) > range1(2)
    error('Limits of integration are inadmissible.');
end
  
%  call B-spline version if both bases are the same B-spline basis
%  else proceed with the use of the Romberg integration
  
if strcmp(type1,'bspline')       & strcmp(type2,'bspline')      & ...
    basisobj1 == basisobj2       & ...
    strcmp(class(Lfdobj1),'double') & strcmp(class(Lfdobj1),'double') & ...
    nargin < 6                   & ...
    all(rng == range1)
 
    ss = inprod_bspline(fdstr1, fdstr2, Lfdobj1, Lfdobj2);
    
else
 
  %  set up first iteration

    iter  = 1;
    width = rng(2) - rng(1);
    JMAXP = JMAX + 1;
    h     = ones(JMAXP,1);
    h(2)  = 0.25;
    s = reshape(zeros(JMAXP*nrep1*nrep2,1),[JMAXP,nrep1,nrep2]);
    %  the first iteration uses just the endpoints
    fx1 = eval_fd(rng, fdstr1, Lfdobj1);
    fx2 = eval_fd(rng, fdstr2, Lfdobj2);
    if ~isnumeric(wtfd)
        wtd = eval_fd(wtfd, rng);
        fx2 = (wtd * ones(1,nrep2)) .* fx2;
    end
    s(1,:,:)  = width.*(fx1' * fx2)./2;
    tnm = 0.5;

    %  now iterate to convergence

    for iter = 2:JMAX
        tnm = tnm.*2;
        del = width./tnm;
        x   = rng(1)+del/2:del:rng(2);
        fx1 = eval_fd(x, fdstr1, Lfdobj1);
        fx2 = eval_fd(x, fdstr2, Lfdobj2);
        if ~isnumeric(wtfd)
            wtd = eval_fd(wtfd, x);
            fx2 = (wtd * ones(1,nrep2)) .* fx2;
        end
        chs = reshape(width.*(fx1' * fx2)./tnm,[1,nrep1,nrep2]);
        s(iter,:,:) = (s(iter-1,:,:) + chs)./2;
        if iter >= 5
            ind = (iter-4):iter;
            ya = s(ind,:,:);
            xa = h(ind);
            absxa = abs(xa);
            [absxamin, ns] = min(absxa);
            cs = ya;
            ds = ya;
            y  = squeeze(ya(ns,:,:));
            ns = ns - 1;
            for m = 1:4
                for i = 1:(5-m)
                    ho      = xa(i);
                    hp      = xa(i+m);
                    w       = (cs(i+1,:,:) - ds(i,:,:))./(ho - hp);
                    ds(i,:,:) = hp.*w;
                    cs(i,:,:) = ho.*w;
                end
                if 2*ns < 5-m
                    dy = squeeze(cs(ns+1,:,:));
                else
                    dy = squeeze(ds(ns,:,:));
                    ns = ns - 1;
                end
                y = y + dy;
            end
            ss = reshape(y, nrep1, nrep2);
            errval = max(max(abs(dy)));
            ssqval = max(max(abs(ss)));
            if all(ssqval > 0)
                crit = errval./ssqval;
            else
                crit = errval;
            end
            if crit < EPS & iter >= JMIN
                return
            end
        end
        s(iter+1,:,:) = s(iter,:,:);
        h(iter+1)   = 0.25.*h(iter);
    end
    disp(['No convergence after',num2str(JMAX),' steps in INPROD']);
  
end

