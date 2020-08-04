function basisobj = create_polygon_basis(argvals)
%  CREATE_POLYGON_BASIS Creates a polygonal basis
%  Argument:
%  ARGVALS  ... strictly increasing argument values
%  Return:
%  BASIS_FD ... a functional data basis object of type 'polygon'

%  last modified 30 January 2003

%  check that argument values are strictly increasing

if min(diff(argvals)) <= 0
    error('ARGVALS are not strictly increasing.');
end

type     = 'polyg';
nbasis   = length(argvals);
rangeval = [min(argvals), max(argvals)];
params   = argvals;

basisobj = basis(type, rangeval, nbasis, params);

