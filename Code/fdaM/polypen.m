function penaltymatrix = polypen(basisobj, Lfd)
%  POLYPEN  Computes the monomial penalty matrix.
%  Arguments:
%  BASISFD  ... a basis.fd object
%  Lfd     ... either the order of derivative or a
%               linear differential operator to be penalized.
%  Returns a list the first element of which is the basis matrix
%   and the second element of which is the diagonal of the penalty matrix.

%  Last modified:  15 March 99

  if ~strcmp(class(basisobj), 'basis')
    error('First argument is not a basis.fd object.');
  end

  if nargin < 2
    Lfd = 2;
  end

  type = getbasistype(basisobj);
  if ~strcmp(type, 'poly')
    error('BASISOBJ not of type poly');
  end

  penaltymatrix = inprod(basisobj, basisobj, Lfd, Lfd);

