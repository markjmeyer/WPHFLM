function derivfd = deriv(fdobj, Lfd)
%  DERIV  Applies LDO Lfd to functional data object FDOBJ.
%  Lfd is either a positive integer or a
%    fd object defining a linear differential operator.

%  last modified 27 May 1999

  if nargin < 2
    Lfd = 1;
  end
  if ~isa_fd(fdobj)
    error('Argument  FD not a functional data object.');
  end

  basisobj = getbasis(fdobj);
  nbasis   = getnbasis(basisobj);
  rangeval = getbasisrange(basisobj);

  evalarg  = linspace(rangeval(1), rangeval(2), 10*nbasis+1);
  Lfdmat   = eval_fd(fdobj, evalarg, Lfd);

  Lfdcoef  = project_basis(Lfdmat, evalarg, basisobj);

  Dfdnames    = getnames(fdobj);
  Dfdnames{3} = ['D',Dfdnames{3}];

  derivfd = fd(Lfdcoef, basisobj, Dfdnames);


