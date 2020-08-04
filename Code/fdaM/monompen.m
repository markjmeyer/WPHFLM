function penaltymat = monompen(basisobj, Lfdobj)
%  MONOMPEN  Computes the monomial penalty matrix.
%  Arguments:
%  BASISFD  ... a monomial basis object
%  LFDOBJ   ... either the order of derivative or a
%               linear differential operator to be penalized.
%  Returns a list the first element of which is the basis matrix
%   and the second element of which is the diagonal of the penalty matrix.

%  Last modified:  16 January 2003

%  check BASISOBJ

if ~strcmp(class(basisobj), 'basis')
    error('First argument is not a basis.fd object.');
end
type = getbasistype(basisobj);
if ~strcmp(type, 'monom')
    error('BASISOBJ not of type monom');
end

%  set default value of LFDOBJ

if nargin < 2
    Lfdobj = int2Lfd(2);
end

%  check LFDOBJ

Lfdobj = int2Lfd(Lfdobj);

%  get basis information

nbasis     = getnbasis(basisobj);
xrange     = getbasisrange(basisobj);
exponents  = getbasispar(basisobj);

%  compute the penalty matrix

penaltymat = zeros(nbasis);
for ibasis=1:nbasis
    ideg = exponents(ibasis);
    if nderiv == 0
        ifac = 1;
    else
        ifac = ideg;
        for k=2:nderiv
            ifac = ifac*(ideg - k + 1);
        end
    end
    for jbasis=1:ibasis
        jdeg = exponents(jbasis);
        if nderiv == 0
            jfac = 1;
        else
            jfac = jdeg;
            for k=2:nderiv
                jfac = jfac*(jdeg - k + 1);
            end
        end
        if ideg >= nderiv & jdeg >= nderiv
            penaltymat(ibasis,jbasis) = ifac*jfac* ...
                (xrange(2)^(ideg+jdeg-2*nderiv+1) -  ...
                xrange(1)^(ideg+jdeg-2*nderiv+1));
            penaltymat(jbasis,ibasis) = penaltymat(ibasis,jbasis);
        end
    end
end
else
    penaltymat = inprod(basisobj, basisobj, Lfdobj, Lfdobj);
end
