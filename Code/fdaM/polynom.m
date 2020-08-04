function polynommat = polynom(x, norder, nderiv, ctr)
%  POLYNOM  Values of monomials, or their derivatives.
%  The powers of X are the NBASIS nonnegative integers in EXPONENTS.
%  The default is 1, meaning X itself.
%  Arguments are as follows:
%  X      ... array of values at which the polynomials are to
%                evaluated
%  NORDER ... highest degree plus one
%  NDERIV ... order of derivative to be returned.
%  CTR    ... a constant shift that helps to keep the polynomials from
%                getting too ill-conditioned.  A good choice is the mid-range.
%  Return is:
%  A matrix with length(X) rows and NBASIS columns containing
%    the values of the monomials or their derivatives

%  last modified 15 June 1999

  n = length(x);
  xdim = size(x);
  if length(xdim) == 2
    if xdim(2) == n
      x = x';
      xdim = size(x);
    end
    if xdim(1) ~= n | xdim(2) ~= 1
      error('X is not a column vector.');
    end
  end

  % set default arguments
  if nargin < 4, ctr = (min(x) + max(x))/2; end
  if nargin < 3, nderiv = 0;                end
  if nargin < 2, exponents = 1;             end

  ndegree = norder - 1;
  if nderiv > ndegree
    error('NDERIV exceeds highest degree of polynomial.');
  end

  nbasis = norder;
  polynommat = zeros(n,nbasis);
  exponents  = 0:ndegree;
  if nderiv == 0
    for ibasis=1:nbasis
      polynommat(:,ibasis) = (x-ctr).^exponents(ibasis);
    end
  else
    for ibasis=1:nbasis
      degree = exponents(ibasis);
      if nderiv <= degree
        fac = degree;
        for ideriv=2:nderiv
          fac = fac*(degree-ideriv+1);
        end
        polynommat(:,ibasis) = fac.*x.^(degree-nderiv);
      end
    end
  end

