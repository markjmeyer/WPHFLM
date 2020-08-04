function  smthfd = smooth(fdobj, lambda, Lfd, rebase)
%SMOOTH  Smooths a functional data object.
%  Arguments for this function:
%
%  FDOBJ   ... A functional data object.
%
%  LAMBDA  ... The smoothing parameter determining the weight to be
%              placed on the size of the derivative = smoothing.  This
%              is 0 by default, but this will produce a warning
%              message that no smoothing has been carried out.
%
%  Lfd     ... The order of derivative or a linear differential
%              operator to be penalized = the smoothing phase.
%              By default Lfd is set = function GETBASISPEnaLTY
%
%  If rebase=1 and the basis type is 'polyg' then the basis
%    is changed to a cubic bspline  basis and before smoothing
%
%  Returns a functional data object containing a smoothed version
%    of the input functional data object
%
%  Last modified:  30 December 2000
%
% Rebase to default B spline basis if rebase is T and basistype is
%    polygonal.  Then test to see if any smoothing is actually required.
%

  if nargin < 4
    rebase = 1;
  end

  if nargin < 3
    Lfd = 0;
  end

  if nargin < 2
    lambda = 0;
  end

  basis = getbasis(fdobj);
  type  = getbasistype(basis);
  if rebase == 1 & strcmp(type,'polyg')
    fdnames = getnames(fdobj);
    params = getbasispar(basis);
    fdobj = data2fd(getcoef(fdobj), params, fdnames);
    basis = getbasis(fdobj);
  end

  if lambda <= 0
    warning('LAMBDA was not positive. No smoothing carried out.');
    smthfd = fdobj;
    return;
  end
%
%  Main smoothing step
%
  coef  = getcoef(fdobj);
  coefd = size(coef);
  ndim  = length(coefd);
  if ndim == 3
     nvar = coefd(3);
  else
     nvar = 1;
  end
  Bmat  = inprod(basis, basis);
%
%  set up coefficient matrix for normal equations
%
  penmat = getbasispenalty(basis, Lfd);
  Cmat   = Bmat + lambda .* penmat;
%
%  solve normal equations for each observation
%
  if(ndim < 3)
     Dmat = inprod(basis, fdobj);
     coef = symsolve(Cmat, Dmat);
  else
    for ivar = 1:nvar
       Dmat = inprod(basis, fdobj(:,ivar));
       coef(:,:,ivar) = symsolve(Cmat, Dmat);
    end
  end
%
%  replace coefficient matrix = fd, leaving other properties alone
%
  smthfd = fd(coef, basis, getnames(fdobj));


