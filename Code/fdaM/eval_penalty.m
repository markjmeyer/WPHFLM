function penaltymat = eval_penalty(basisobj, Lfdobj)
%  EVAL_PENALTY evaluates a the inner products of a linear 
%  differential operator L defined by LFDOBJ applied to a set of 
%  basis functions defined by BASISOBJ.
%
%  LFDOBJ is a functional data object defining the order m 
%  NONHOMOGENEOUS linear differential operator of the form
%  Lx(t) = w_0(t) x(t) + w_1(t) Dx(t) + ... 
%          w_{m-1}(t) D^{m-1}x(t) + D^m x(t) + ...
%          a_1(t) u_1(t)  + ... + a_k(t) u_k(t).
%  This is a change from previous usage where LFDOBJ was assumed to 
%  define a HOMOGONEOUS differential operator.  See function
%  @Lfd/Lfd() for details.
%
%  Arguments:
%  BASISOBJ ... A basis object
%  LFDOBJ   ... A linear differential operator object
%              applied to the functions that are evaluated.
%
%  Note that the first two arguments may be interchanged.
%
%  Returns:  An array of function values corresponding to the evaluation
%              arguments in EVALARG

%  last modified 12 January 2003

if ~isa_basis(basisobj)
    error('Argument BASISOBJ is not a functional basis object.');
end

%  set up default value for Lfdobj

if nargin < 2 
    Lfdobj = int2Lfd(2);  
end

%  deal with the case where LFDOBJ is an integer

if isnumeric(Lfdobj)
    nderiv = Lfdobj;
    if nderiv ~= round(Lfdobj)
        error('LFDOBJ numeric but not an integer.');
    end
    if nderiv < 0
        error('LFDOBJ an integer but negative.');
    end
    Lfdobj = int2Lfd(nderiv);
end
    
%  check LFDOBJ

if ~isa_Lfd(Lfdobj)
    error (['Argument LFDOBJ is not a ', ...
            'linear differential operator object.']);
end

%  determine basis type

type = getbasistype(basisobj);

%  choose appropriate penalty matrix function
switch type
    case 'fourier'
        penaltymat = fourierpen(basisobj, Lfdobj);
    case 'bspline'
        penaltymat = bsplinepen(basisobj, Lfdobj);
    case 'power'
        penaltymat = powerpen(basisobj, Lfdobj);
    case 'polyg'
        penaltymat = polygpen(basisobj, Lfdobj);
    case 'expon'
        penaltymat = exponpen(basisobj, Lfdobj);
    case 'const'
        penaltymat = basisobj.rangeval(2) - basisobj.rangeval(1);
    otherwise
        error('Basis type not recognizable');
end

