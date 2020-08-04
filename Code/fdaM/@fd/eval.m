function evalarray = eval(fd, x, Lfd)
%  EVAL  Evaluates the value of the derivative or the linear
%  differential operator defined in Lfd for a functional data
%  object FD at argument values in array XS.
%  The default for Lfd is 0, meaning that the function itself
%  is evaluated.

%  last modified 1 July 1998

  if nargin < 3
    Lfd = 0;
  end

  sizex = size(x);
  ndim  = length(sizex);
  switch ndim
    case 2
      if sizex(1) > 1 & sizex(2) > 1
        error('Second argument must be a vector');
      else
        n = length(x);
      end
    case 1
      n = length(x);
      x = x';
    otherwise
      error('Second argument must be a vector');
  end

  if ~isa_fd(fd)
    error ('Argument FD is not a functional data object.');
  end

  if ~isa_Lfd(Lfd)
    error (['Argument Lfd is neither a functional data object', ...
             ' nor an integer.']);
  end

  basisobj = getbasis(fd);
  nbasis   = getnbasis(basisobj);

  if strcmp(class(Lfd),'double')
    if length(Lfd) == 1
      nderiv = round(Lfd);
      if nderiv ~= Lfd
        error('Order of derivative must be an integer');
      end
      if nderiv < 0
        error('Order of derivative must be 0 or positive');
      end
    else
      error('Order of derivative must be a single number');
    end
  else
    Lcoef  = getcoef(Lfd);
    Lcoefd = size(Lcoef);
    nderiv = Lcoefd(2);
  end

  basismat = getbasismatrix(x, basisobj, nderiv);

  %  If a differential operator has been defined in DERIVWTFD, compute
  %  the weighted combination of derivatives

  if nderiv > 0 & strcmp(class(Lfd), 'fd')
    Lbasisobj = Lfd.basisobj;
    Lbasismat = getbasismatrix(x, Lbasisobj);
    Lfdmat = Lbasismat * Lcoef;
    onerow = ones(1,nbasis);
    for j = 1:nderiv
      if any(abs(Lfdmat(:,j)) > 1e-7)
        basismat = basismat + (Lfdmat(:,j)*onerow) .* ...
          getbasismatrix(x, basisobj, j-1);
      end
    end
  end

  %  evaluate the functions at arguments in X

  coef  = getcoef(fd);
  coefd = size(coef);
  ndim  = length(coefd);

  if ndim <= 2
    evalarray = basismat * coef;
  else
    ncurves   = coefd(2);
    nvar      = coefd(3);
    evalarray = zeros(n,ncurves,nvar);
    for ivar = 1:nvar
      evalarray(:,:,ivar) = basismat * coef(:,:,ivar);
    end
  end
