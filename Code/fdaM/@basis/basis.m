function basisobj = basis(fdtype, rangeval, nbasis, params)
%  BASIS  Creates a functional data basis.
%  Arguments
%  FDTYPE   ... a string indicating the type of basis.  This may be one of
%               'Fourier', 'fourier', 'Fou', 'fou',
%               'Bspline', 'bspline', 'Bsp', 'bsp',
%               'pol', 'poly', 'polynomial',
%               'con', 'const', 'constant'
%               'exp', 'exponen', 'exponential'
%               'polyg' 'polygon', 'polygonal'
%  RANGEVAL ... an array of length 2 containing the lower and upper
%               boundaries for the rangeval of argument values
%  NBASIS   ... the number of basis functions
%  PARAMS   ... If the basis is 'fourier', this is a single number indicating
%                 the period.  That is, the basis functions are periodic on
%                 the interval (0,PARAMS) or any translation of it.
%               If the basis is 'bspline', the values are interior points at
%                 which the piecewise polynomials join.
%                 Note that the number of basis functions NBASIS is equal
%                 to the order of the Bspline functions plus the number of
%                 interior knots, that is the length of PARAMS.
%               This means that NBASIS must be at least 1 larger than the
%                 length of PARAMS.
%  Returns
%  BASIS_fd  ... a functional data basis object
%  An alternative name for this function is CREATE_BASIS, but PARAMS argument
%     must be supplied.
%  Specific types of bases may be set up more conveniently using functions
%  CREATE_BSPLINE_basis  ...  creates a b-spline basis
%  CREATE_FOURIER_basis  ...  creates a fourier basis
%  CREATE_POLYGON_basis  ...  creates a polygonal basis
%  CREATE_MONOM_BASIS    ...  creates a monomial basis
%  CREATE_CONSTANT_BASIS ...  creates a constant basis

%  last modified 22 January 2002

  if nargin==0
    basisobj.type     = 'bspline';
    basisobj.rangeval = [0,1];
    basisobj.nbasis   = 2;
    basisobj.params   = [];
    basisobj = class(basisobj, 'basis');
    return;
  end

  if isa(fdtype, 'basis')
    basisobj = fdtype;
    return;
  end

  fdtype = use_proper_basis(fdtype);
  if strcmp(fdtype,'unknown')
    error ('TYPE unrecognizable.');
  end

  switch fdtype
   case 'fourier'
     paramvec   = rangeval(2) - rangeval(1);
     period     = params(1);
     if (period <= 0)
       error ('Period must be positive for a Fourier basis');
     end
     params = period;
     if (2*floor(nbasis/2) == nbasis)
       nbasis = nbasis + 1;
     end
   case 'bspline'
     allbrks  = linspace(rangeval(1), rangeval(2), nbasis-2);
     intbrks  = allbrks(2:(nbasis-3));
     paramvec = intbrks;
     breaks   = params;
     nbreaks  = length(breaks);
     if (nbreaks < 1)
       error ('No values=BREAKS');
     end
     if (nbasis <= nbreaks)
       error ('NBASIS must exceed number of internal knots');
     end
     if (breaks(1) <= rangeval(1))
       error('Smallest value=BREAKS not within RANGEVAL');
     end
     if (breaks(nbreaks) >= rangeval(2))
       error('Largest value=BREAKS not within RANGEVAL');
     end
     if min(diff(breaks)) <= 0
       error('Values=BREAKS not strictly increasing');
     end
     params = breaks;
   case 'expon'
     if (length(params) ~= nbasis)
       error(['No. of parameters not equal to no. of basis fns ',  ...
              'for exponential basis.']);
     end
   case 'polyg'
     if (length(params) ~= nbasis)
       error(...
       'No. of parameters not equal to no. of basis fns for polygonal basis.');
     end
   case 'power'
    if length(params) ~= nbasis
       error(...
       'No. of parameters not equal to no. of basis fns for power basis.');
     end
   case 'const'
     params = 0;
   otherwise
     error('Unrecognizable basis');
  end

  basisobj.type     = fdtype;
  basisobj.rangeval = rangeval;
  basisobj.nbasis   = nbasis;
  basisobj.params   = params;

  basisobj = class(basisobj, 'basis');

