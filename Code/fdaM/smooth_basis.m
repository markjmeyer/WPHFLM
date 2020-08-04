function smoothlist = smooth_basis(y, argvals, basis, wtvec, Lfdobj, ...
                                   lambda, fdnames)
%SMOOTH_BASIS  Smooths discrete curve values using penalized basis functions
%  Arguments for this function:
%
%  Y        ... an array containing values of curves
%               If the array is a matrix, rows must correspond to argument
%               values and columns to replications, and it will be assumed
%               that there is only one variable per observation.
%               If Y is a three-dimensional array, the first dimension
%               corresponds to argument values, the second to replications,
%               and the third to variables within replications.
%               If Y is a vector, only one replicate and variable are assumed.
%  ARGVALS  ... A set of argument values, set by default to equally spaced on
%               the unit interval (0,1).
%  BASIS    ... A basis.fd object created by function create_basis.fd.
%  WTVEC    ... A vector of N weights, set to one by default, that can
%               be used to differentially weight observations.
%  LFDOBJ   ... The order of derivative or a linear differential
%               operator to be penalized.
%  LAMBDA   ... The smoothing parameter determining the weight to be
%               placed on the roughness penalty.  
%               This is 0 by default.
%  FDNAMES  ... A cell of length 3 with names for
%               1. argument domain, such as 'Time'
%               2. replications or cases
%               3. the function.
%  Returns a cell object containing:
%    FDOBJ ...  an object of class fd containing coefficients
%    DF    ...  a degrees of freedom measure
%    GCV   ...  a measure of lack of fit discounted for df.

%  last modified 27 January 2003

if nargin < 3
    error('There is not at least three arguments.');
end

n = length(argvals);

%  set default argument values

if nargin < 7
    fdnames{1} = 'time';
    fdnames{2} = 'reps';
    fdnames{3} = 'values';
end

if nargin < 6, lambda = 0;           end
if nargin < 5, Lfdobj = int2Lfd(2);  end;
if nargin < 4, wtvec  = ones(n,1);   end

%  check LFD

Lfdobj = int2Lfd(Lfdobj);
nderiv = getnderiv(Lfdobj);

%  check BASIS

if ~isa_basis(basis)
    error('BASIS is not a basis object.');
end

nbasis   = getnbasis(basis);
onebasis = ones(1,nbasis);

%  check WTVEC

sizew = size(wtvec);
if (length(sizew) > 1 & sizew(1) > 1 & sizew(2) > 1) | ...
      length(sizew) > 2
    error ('WTVEC must be a vector.');
end
if length(sizew) == 2 & sizew(1) == 1
    wtvec = wtvec';
end
if length(wtvec) ~= n
    error('WTVEC of wrong length');
end
if min(wtvec) <= 0
    error('All values of WTVEC must be positive.');
end

if lambda < 0
    warning ('Value of LAMBDA was negative, and 0 used instead.');
    lambda = 0;
end

sizec = size(y);
ndim  = length(sizec);

if sizec(1) ~= n
    error('Number of arguments differs from first dimension of matrix of values');
end

%  set number of curves and number of variables

switch ndim
    case 1
        ncurves = 1;
        nvar    = 1;
    case 2
        ncurves = sizec(2);
        nvar    = 1;
    case 3
        ncurves = sizec(2);
        nvar    = sizec(3);
    otherwise
        error('Second argument must not have more than 3 dimensions');
end

basismat  = getbasismatrix(argvals, basis);

if n >= nbasis | lambda > 0
    
    %  The following code is for the coefficients completely determined
    
    basisw = basismat .* (wtvec * ones(1,nbasis));
    Bmat = basisw' * basismat;
    
    if ndim < 3
        Dmat = basisw' * y;
    else
        Dmat = zeros(nbasis,ncurves,nvar);
        for ivar = 1:nvar
            Dmat(:,:,ivar) = basisw' * y(:,:,ivar);
        end
    end
    
    df   = nbasis;
    
    if lambda > 0
        %  smoothing required, set up coefficient matrix for normal equations
        afdcell = getafd(Lfdobj);  % multiplier(s) of forcing function(s)
        ufdcell = getufd(Lfdobj);  % forcing function(s)
        if ~isnumeric(Lfdobj) & ~isempty(afdcell) & ~isempty(ufdcell)
            %  here the linear differential operator is not an integer
            %  and is not homogeneous.  
            %  first set up the homogeneous counterpart
            Lfdhom = Lfd(nderiv, getwfd(Lfdobj));
            %  evaluate the penalty matrix for the homogeneous operator
            penmat = eval_penalty(basis, Lfdhom);
            %  set up the part of the roughness penalty affected by the
            %  presence of forcing function(s)
            penvec = zeros(nbasis,1);
            for k=1:nderiv
                afdk = afdcell{k};
                ufdk = ufdcell{k};
                ffdk = times(afdk,ufdk);
                fvec = inprod(basis, ffdk, Lfdhom, int2Lfd(0));
                penvec = penvec - fvec;
            end
        else
            %  here the linear differential operator is homogeneous
            %  only the penalty matrix is needed.
            penmat = eval_penalty(basis, Lfdobj);
            penvec = zeros(nbasis,1);
        end
        Bnorm   = sqrt(sum(sum(Bmat.^2)));
        pennorm = sqrt(sum(sum(penmat.^2)));
        condno  = pennorm/Bnorm;
        if lambda*condno > 1e12
            lambda = 1e12/condno;
            disp(['lambda reduced to ',num2str(lambda),...
                    ' to prevent overflow']);
        end
        Bmat = Bmat + lambda .* penmat;
        if ~all(penvec == 0)
            %  if the linear differential operator is nonhomogeneous
            %  use PENVEC to alter the right side of the equation.
            if ndim < 3
                Dmat = Dmat + lambda.*(penvec*ones(1,ncurves));
            else
                for icurve=1,ncurves
                    for ivar=1:nvar
                        Dmat(:,icurve,ivar) = Dmat(:,icurve,ivar) + lambda .* penvec;
                    end
                end
            end
        end
    end
    
    %  compute inverse of Bmat
    
    if is_diag(Bmat)
        Bmatinv = diag(1/diag(Bmat));
    else
        Lmat    = chol(Bmat);
        Lmatinv = inv(Lmat);
        Bmatinv = Lmatinv * Lmatinv';
    end
    
    %  compute degrees of freedom of smooth
    
    df = sum(diag(Bmatinv * Bmat));
    
    %  solve normal equations for each observation
    
    if ndim < 3
        coef = Bmatinv * Dmat;
    else
        for ivar = 1:nvar
            coef(:,:,ivar) = Bmatinv * Dmat(:,ivar);
        end
    end
    
else
    
    %  The following code is for the underdetermined coefficients:
    %     the number of basis functions exceeds the number of argument values.
    %  No smoothing is used.  
    
    [Qmat,Rmat] = qr(basismat');
    Q1mat  = Qmat(:,1:n);
    Q2mat  = Qmat(:,((n+1):nbasis));
    Hmat   = eval_penalty(basis);
    Q2tHmat   = Q2mat' * Hmat;
    Q2tHQ2mat = Q2tHmat * Q2mat;
    Q2tHQ1mat = Q2tHmat * Q1mat;
    if ndim < 3
        z1mat = symsolve(Rmat, y);
        z2mat = symsolve(Q2tHQ2mat, Q2tHQ1matz1mat);
        coef = Q1mat * z1mat + Q2mat * z2mat;
    else
        for ivar = 1:nvar
            z1mat = symsolve(Rmat, y(:,:,ivar));
            z2mat = symsolve(Q2tHQ2mat, Q2tHQ1mat*z1mat);
            coef(:,:,ivar) = Q1mat * z1mat + Q2mat * z2mat;
        end
    end
    df = n;
end

  %  compute  GCV index

if df < n
    if ndim < 3
        yhat = basismat * coef;
        SSE = sum((y - yhat).^2);
    else
        SSE = 0;
        for ivar = 1:nvar
            yhat = basismat * coef(:,:,ivar);
            SSE = SSE + sum((y(:,:,ivar) - yhat).^2);
        end
    end
    gcv = (SSE/n)/(nvar*(n - df)/n)^2;
else
    gcv = NaN;
end

fdobj = fd(coef, basis, fdnames);

smoothlist.fdobj = fdobj;
smoothlist.df    = df;
smoothlist.gcv   = gcv;


