function basisobj = create_exponential_basis(rangeval, nbasis, ratevec)
%  CREATE_EXPONENTIAL_BASIS  Creates a exponential basis: 
%            exp[RATEVEC(1)*x], exp[RATEVEC(2)*x], ...
%  Argument:
%  RANGEVAL ... an array of length 2 containing the lower and upper
%               boundaries for the rangeval of argument values.  If a 
%               single value is input, it must be positive and the lower
%               limit of the range is set to 0.
%  NBASIS    ... number of basis functions
%  RATEVEC   ... an array of NBASIS rate values
%  Return:
%  BASIS     ... a functional data basis object of type 'expon'

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

if nargin < 2, nbasis = 1;    end
if nargin < 3, ratevec = 0:(nbasis-1); end

% check if there are duplicate ratevec

if min(diff(sort(ratevec))) <= 0
    error('There are duplicate ratevec.');
end

type   = 'expon';
params = ratevec;

basisobj = basis(type, rangeval, nbasis, params);

