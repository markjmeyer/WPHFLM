function plusfd = plus(fd1, fd2)
% PLUS Sum of two functional data objects having the same basis.
% Sum functional data object PLUSFD inherits the fdnames field of the
%  first argument.

%  last modified 1 July 1998

  if ~(isa_fd(fd1) & isa_fd(fd2))
    error('Both arguments are not functional data objects.');
  end
  basisobj1 = getbasis(fd1);
  basisobj2 = getbasis(fd2);
  type1   = getbasistype(basisobj1);
  type2   = getbasistype(basisobj2);
  nbasis1 = getnbasis(basisobj1);
  nbasis2 = getnbasis(basisobj2);
  range1  = getbasisrange(basisobj1);
  range2  = getbasisrange(basisobj2);
  params1 = getbasispar(basisobj1);
  params2 = getbasispar(basisobj2);
  if ~strcmp(type1, type2)  | ...
     any(range1  ~= range2) | ...
         nbasis1 ~= nbasis2 | ...
     any(params1 ~= params2)
    error('Basis structures are not equal.');
  end
  coef1 = getcoef(fd1);
  coef2 = getcoef(fd2);
  if (any(size(coef1) ~= size(coef2)))
    error('Coefficient arrays are not of same dimensions.');
  end

  fdnames = getnames(fd1);
  plusfd.coef = coef1+coef2;
  plusfd.basisobj = basisobj1;
  plusfd.fdnames  = fdnames;
  plusfd = class(plusfd, 'fd');

