function minusfd = uminus(fd)
% Unary minus or negative of functional data object.

%  last modified 29 June 1998

  if (~(isa_fd(fd1) & isa_fd(fd2))
    error('Both arguments are not functional data objects.');
  end
  coef     = getcoef(fd);
  basisobj = getbasis(fd);
  fdnames  = getnames(fd);
  minusfd  = fd(-coef, basisobj, fdnames);

