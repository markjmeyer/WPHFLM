function basismat = polyg(x, argvals, nderiv)
%  POLYG Evaluates the basis for a linear interpolant or its first derivative.
%  It calls function spcol.
%  Arguments are as follows:
%  EVALARG ... A vector of values at which the spline functions are to
%              evaluated
%  ARGVAL  ... a STRICTLY INCREASING sequence of argument values.
%  NDERIV  ... Either 0 or 1.  0 means only function values
%  Return is a matrix with length(EVALARG) rows and number of columns equal to
%             number of argument values

%  last modified 30 January 2003

if nargin < 3
    nderiv = 0;
end

argvals = argvals(:);
nargvals = length(argvals);

evalarg = evalarg(:);
n  = length(evalarg);

if (max(evalarg) > max(argvals)) | (min(evalarg) < min(argvals)) 
    error('ARGVALS do not span the values of EVALARG.');
end

if (min(diff(argvals)) <= 0 )
    error('Break-points are not strictly increasing');
end

if (~(nderiv == 0 | nderiv == 1))
    error('NDERIV is neither 0 nor 1.');
end

derivs = nderiv .* ones(n,1);
nbasis = length(argvals);

knots    = [argvals(1), argvals, argvals(nbasis)];

if nderiv == 0
    tau    = reshape(evalarg,[1,n]);
else
    onevec = ones(nderivp1,1);
    evalargmat   = reshape(evalarg,[1,n]);
    tau    = reshape(onevec * evalargmat,[1,n*nderivp1]);
end

%  call function SPCOL

norder = 2;
basismat = spcol(knots,norder,tau);

if nderiv > 0
    index = 2:2:(n*2);
    basismat = basismat(index,:);
end


