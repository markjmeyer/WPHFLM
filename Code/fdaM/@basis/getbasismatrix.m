function basismat = getbasismatrix(evalarg, basisobj, nderiv)
%  GETBASISMATRIX   Computes the basis matrix evaluated at arguments in
%    EVALARG associated with basis.fd object BASISOBJ.
%    The returned basis matrix BASISMAT contains the basis
%    derivatives of order NDERIV (0 by default).

%  last modified 11 March 2003

if nargin < 3,  nderiv = 0;  end

if ~isa_basis(basisobj)
    error('Argument BASISOBJ is not a functional basis object');
end

type   = getbasistype(basisobj);
nbasis = getnbasis(basisobj);

switch type
    case 'fourier'
        period   = basisobj.params(1);
        nbasis   = basisobj.nbasis;
        basismat = fourier(evalarg, nbasis, period, nderiv);
    case 'bspline'
        rangex   = basisobj.rangeval;
        breaks   = [rangex(1), basisobj.params, rangex(2)];
        norder   = basisobj.nbasis - length(breaks) + 2;
        basismat = bsplineFDA(evalarg, breaks, norder, nderiv);
    case 'polyg'
        basismat = polyg(evalarg, basisobj.params);
    case 'power'
        exponents = getbasispar(basisobj);
        norder    = getnbasis(basisobj);
        basismat  = powerbasis(evalarg, exponents, nderiv);
    case 'expon'
        exponents = basisobj.params;
        basismat = expon(evalarg, exponents, nderiv);
    case 'const'
        basismat = ones(length(evalarg),1);
    otherwise
        error('Basis type not recognizable')
end

