function evalarray = eval(basisobj, x, Lfd)
%  EVAL  Evaluates the value of the derivative or the linear
%  differential operator defined in Lfd for a basisobj
%  object BASISOBJ at argument values in array XS.
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

  if ~isa_basis(basisobj)
    error ('Argument BASISOBJ is not a functional data object.');
  end

  if ~isa_Lfd(Lfd)
    error (['Argument Lfd is neither a functional data object', ...
             ' nor an integer.']);
  end

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

  evalarray = getbasismatrix(x, basisobj, nderiv);

  %  If a differential operator has been defined in DERIVWTFD, compute
  %  the weighted combination of derivatives

  if nderiv > 0 & strcmp(class(Lfd), 'fd')
    Lbasisobj = Lfd.basisobj;
    Lbasismat = getbasismatrix(x, Lbasisobj);
    Lfdmat = Lbasismat * Lcoef;
    onerow = ones(1,nbasis);
    for j = 1:nderiv
      if any(abs(Lfdmat(:,j)) > 1e-7)
        evalarray = evalarray + (Lfdmat(:,j)*onerow) .* ...
          getbasismatrix(x, basisobj, j-1);
      end
    end
  end

