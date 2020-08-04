function rdividefd = rdivide(fd1, fd2)
%  RDIVIDE  Pointwise quotient of two functional data objects, or
%    the quotient of a scalar and a functional data object.
%  Two functional data objects need not have the same basis, but must
%    have the same number of replicates and the same number of functions.
%  Functional data object RDIVIDEFD inherits the fdnames of the
%    first functional data object and the basis of the object with the
%    most basis functions.

%  last modified 28 February 2002

if ~(isa_fd(fd1) | isa_fd(fd2))
    error('Neither argument for * is a functional data object.');
end

if ~(isnumeric(fd1) | isnumeric(fd2))
    %  both arguments are functional data objects
    if (~(isa_fd(fd1) & isa_fd(fd2)))
        error('Both arguments are not functional data objects.');
    end
    coef1  = getcoef(fd1);
    coef2  = getcoef(fd2);
    coefd1 = size(coef1);
    coefd2 = size(coef2);
    ndim1  = length(coefd1);
    ndim2  = length(coefd2);
    if coefd1(2) ~= coefd2(2) 
        error('Number of replications are not equal.');
    end
    if ndim1 > 2 & ndim2 > 2 & ndim1 ~= ndim2
        error(['Both arguments multivariate, but involve different numbers'...
               ' of functions.']);
    end
    basisobj1 = getbasis(fd1);
    basisobj2 = getbasis(fd2);
    nbasis1   = getnbasis(basisobj1);
    nbasis2   = getnbasis(basisobj2);
    rangeval1 = getbasisrange(basisobj1);
    rangeval2 = getbasisrange(basisobj2);
    if (any(rangeval1 ~= rangeval2))
        error('The ranges of the arguments are not equal.');
    end
    neval = max(10*max(nbasis1,nbasis2) + 1, 101);
    neval = min(neval,201);
    evalarg   = linspace(rangeval1(1),rangeval2(2), neval);
    fdarray1  = eval_fd(fd1, evalarg);
    fdarray2  = eval_fd(fd2, evalarg);
    if (ndim1 <= 2 & ndim2 <= 2) | (ndim1 > 2 & ndim2 > 2)
        fdarray = fdarray1./fdarray2;
    end
    if ndim1 == 2 & ndim2 > 2
        fdarray = zeros(coefd2);
        for ivar = 1:coefd2(3)
            fdarray(:,:,ivar) = fdarray1 ./ fdarray2(:,:,ivar);
        end
    end
    if ndim1 > 2 & ndim2 == 2
        fdarray = zeros(coefd1);
        for ivar = 1:coefd1(3)
          fdarray(:,:,ivar) = fdarray1(:,:,ivar) ./ fdarray2;
        end
    end
    if nbasis1 > nbasis2
        basisobj = basisobj1;
        coefquot = project_basis(fdarray, evalarg, basisobj);
    else
        basisobj = basisobj2;
        coefquot = project_basis(fdarray, evalarg, basisobj);
    end
    
    fdnames1 = fd1.fdnames;
    fdnames2 = fd2.fdnames;
    fdnames  = fdnames1;
    fdnames{3} = [fdnames1{3},'*',fdnames2{3}];
    
 else
    %  one argument is numeric and the other is functional
    if ~(isnumeric(fd1) | isnumeric(fd2))
        error('Neither argument for * is numeric.');
    end
    if isnumeric(fd1)
        fac = fd1;
        fd  = fd2;
    else
        fac = fd2;
        fd  = fd1;
        if fac == 0
          error('Divide by zero.');
        end
    end
    coef     = getcoef(fd);
    coefd    = size(coef);
    basisobj = getbasis(fd);
    nbasis   = getnbasis(basisobj);
    rangeval = getbasisrange(basisobj);
    neval    = max(10*nbasis + 1,101);
    neval    = min(neval,201);
    evalarg  = linspace(rangeval(1),rangeval(2), neval);
    if isnumeric(fd1)
        denom   = eval(fd, evalarg);
        fdarray = fac ./ denom;
    else
        numer   = eval(fd, evalarg);
        fdarray = numer ./ fac;
    end
    coefquot = project_basis(fdarray, evalarg, basisobj);
    fdnames  = fd.fdnames;
    if isnumeric(fd1)
        fdnames{3} = [num2str(fac),'/',fdnames{3}];
    else
        fdnames{3} = [fdnames{3},'/',num2str(fac)];
    end
end

rdividefd.coef     = coefquot;
rdividefd.basisobj = basisobj;
rdividefd.fdnames  = fdnames;

rdividefd = class(rdividefd, 'fd');


