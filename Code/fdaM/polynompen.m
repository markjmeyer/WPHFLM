function penaltymat = polynompen(basisobj, Lfd)
%  POLYNOMPEN  Computes the monomial penalty matrix.
%  Arguments:
%  BASISFD  ... a monomial basis object
%  Lfd     ... either the order of derivative or a
%               linear differential operator to be penalized.
%  Returns a list the first element of which is the basis matrix
%   and the second element of which is the diagonal of the penalty matrix.

%  Last modified:  15 June 99

  if ~strcmp(class(basisobj), 'basis')
    error('First argument is not a basis.fd object.');
  end

  if nargin < 2
    Lfd = 2;
  end

  type = getbasistype(basisobj);
  if ~strcmp(type, 'poly')
    error('BASISOBJ not of type POLY');
  end

  ctr    = getbasispar(basisobj);
  nbasis = getnbasis(basisobj);
  ndegre = nbasis - 1;
  exponents = 0:ndegre;

  if ~isa_Lfd(Lfd)
    error (['Argument Lfd is neither a functional data object', ...
             ' nor an integer.']);
  end

  if strcmp(class(Lfd),'double')
    if length(Lfd) == 1
      nderiv = round(Lfd);
      if nderiv ~= Lfd
        error('Order of derivative must be an integer');
      end
      if nderiv < 0
        error('Order of derivative must be 0 or positive');
      end
    else
      error('Order of derivative must be a single number');
    end
    nbasis = getnbasis(basisobj);
    penaltymat = zeros(nbasis);
    xrange = getbasisrange(basisobj);
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
            ((xrange(2)-ctr)^(ideg+jdeg-2*nderiv+1) -  ...
             (xrange(1)-ctr)^(ideg+jdeg-2*nderiv+1));
          penaltymat(jbasis,ibasis) = penaltymat(ibasis,jbasis);
        end
      end
    end
  else
    penaltymat = inprod(basisobj, basisobj, Lfd, Lfd);
  end
