function basisobj = create_polynomial_basis(rangeval, nbasis, ctr)
%  CREATE_POLY_BASIS  Creates a monomial basis:, 1, x, ..., x^{nbasis-1}
%  Arguments:
%  RANGEVAL ... an array of length 2 containing the lower and upper
%               boundaries for the rangeval of argument values.  If a 
%               single value is input, it must be positive and the lower
%               limit of the range is set to 0.
%  NBASIS   ... the number of basis functions
%  CTR      ... If 1, center the range before evaluating.
%  Return:
%  BASIS.FD  ... a functional data basis object of type 'constant'

%  last modified 30 January 2003

%  check RANGEVAL

if length(rangeval) == 1
    if rangeval <= 0
        error('RANGEVAL a single value that is not positive.');
    end
    rangeval = [0,rangeval];
end

if rangechk(rangeval) ~= 1
    error('RANGEVAL is not a legitimate range.');
end

if nargin < 3, ctr = 0;          end
if nargin < 2, nbasis = 2;       end
if nargin < 1, rangeval = [0,1]; end
type    = 'poly';
params  = ctr;

basisobj = basis(type, rangeval, nbasis, params);

