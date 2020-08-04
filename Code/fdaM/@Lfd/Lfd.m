function Lfdobj = Lfd(m, wfdcell, afdcell, ufdcell)
%  LFD creates a linear differential operator object of the form
%
%  Lx(t) = w_0(t) x(t) + w_1(t) Dx(t) + ... 
%          w_{m-1}(t) D^{m-1}x(t) + D^m x(t) + ...
%          a_1(t) u_1(t)  + ... + a_k(t) u_k(t).
%  
%  Function x(t) is operated on by this operator L, and the operator
%  computes a linear combination of the function and its first m
%  derivatives.  The function x(t) must be scalar.  This part
%  part of the operator (the first two lines in the above equation)
%  is called the HOMOGENEOUS part of the operator.
%
%  The linear combination of derivatives is defined by the weight 
%  or coefficient functions w_j(t), and these are assumed to vary  
%  over t, although of course they may also be constant as a 
%  special case.  
%
%  The operator $L$ is potentially defined by one or more known
%  scalar forcing functions u_1(t), ..., u_k(t), each multiplied by  
%  a weight or coefficient function a(t).  If the forcing functions
%  are defined, then the operator is called NONHOMOGENEOUS, and the
%  nonhomogenous part of the operator is the third lined in the 
%  above equation.
%
%  It may be required that within any of these three groups the 
%  functions will vary in complexity.  Consequently, each group
%  is input to constructor function Lfd() as a cell object, and
%  each individual function within the group is defined within
%  a cell member.
%  
%  The evaluation of the linear differential operator L is applied to 
%  basis functions takes place in function EVAL_BASIS().  Here only 
%  the homogeneous part of the operator is used, defined in WFDCELL.

%  The evaluation of the linear differential operator L applied to
%  functional data objects takes placed in function EVAL_FD(). Here
%  the entire operator is applied, including any forcing functions
%  defined in UFDCELL and weighted by AFDCELL.
%
%  The inner products of the linear differential operator L 
%  applied to basis functions is evaluated in the functions
%  called in function EVAL_PENALTY().  Here only the homogeneous part 
%  of the operator is used, defined in WFDCELL.
%
%  Some important functions also have the capability of allowing
%  the argument that is an LFD object be an integer. They convert 
%  the integer internally to an LFD object by INT2LFD().  These are:
%     EVAL_FD()
%     EVAL_MON()
%     EVAL_POS()
%     EVAL_BASIS()
%     EVAL_PENALTY()
%
%  Arguments:
%
%  M       ... the order of the operator, NDERIV, that is,
%          the highest order of derivative.
%  WFDCELL ... A cell vector object with m cells.
%  AFDCELL ... A cell vector object with k cells,
%          where k is the number of forcing functions.  
%          k may be zero, in which case AFDCELL will be an empty cell,
%          and this is the default if AFDCELL is not supplied.
%  UFDCELL ... A cell vector containing the forcing functions.
%          If UFDCELL is not supplied but AFDCELL is, then the 
%          default is the unit function using the constant basis.
%
%  Simple cases: 
%
%  All this generality may not be needed, and, for example, 
%  often the linear differential operator will be 
%  simply L = D^m, defining Lx(t) = D^mx(t).  Or the weights and 
%  forcing functions may all have the same bases, in which case 
%  it is simpler to use a functional data objects to define them.  
%  These situations cannot be accommodated within Lfd(), but
%  there is function int2Lfd(m) that converts a nonnegative 
%  integer m into an Lfd object equivalent to D^m. 
%  There is also fd2cell(fdobj) and that converts a functional 
%  data object into cell object, which can then be used as
%  an argument of Lfd().
%
%  Returns:
%
%  LFDOBJ ... a functional data object

%  last modified 14 January 2003

% superiorto('double', 'sparse', 'struct', 'cell', 'char', ...
%    'inline', 'basis');

superiorto('double', 'struct', 'cell', 'char', ...
    'inline', 'basis');


%  check m

if ~isnumeric(m)
    error('Order of operator is not numeric.');
end
if m ~= round(m)
    error('Order of operator is not an integer.');
end
if m < 0
    error('Order of operator is negative.');
end

%  check that WFDCELL is a cell object

if ~iscell(wfdcell)
    error('WFDCELL not a cell object.');
end

wfdsize = size(wfdcell);

%  WFDCELL is one-dimensional.  Only possibilities are M = 0 or 1;
if length(wfdsize) == 1
    if m > 1
        error('Dimension of WFDCELL not compatible with M.');
    end
end

%  WFDCELL two-dimensional.  
%  Only possibilities are (1) N > 1, M = 1, and 
%                         (2) N = 1, M >= 1;
if length(wfdsize) == 2
    if wfdsize(1) > 1 &  wfdsize(2) > 1
        error('WFDCELL is not a vector.');
    else
        if m > 0
            if wfdsize(1) ~= m & wfdsize(2) ~= m
                error('Dimension of WFDCELL not compatible with M.');
            end
        end
    end
end

%  WFDCELL has more than two dimensions.  
if length(wfdsize) > 2
    error('WFDCELL has more than two dimensions.');
end

%  find the range, which must be common to all functions

wrange = getbasisrange(getbasis(wfdcell{1}));
wnames = getnames(wfdcell{1});

%  check AFDCELL

if nargin >= 3
    if ~iscell(afdcell)
        error('AFDCELL is not a cell object.');
    end
    if ~isempty(afdcell)
        afdsize = size(afdcell);
        if length(afdsize) > 3
            error('AFDCELL has more than two dimensions');
        end
        if length(afdsize) == 2
            if afdsize(1) > 1 & afdsize(2) > 1
                error('AFDCELL has both dimensions greater than one.');
            end
            if afdsize(1) == 1
                k = afdsize(2);
            else
                k = afdsize(1);
            end
            for j=1:k
                arange = getbasisrange(getbasis(afdcell{j}));
                if any(arange ~= wrange)
                    error('WRANGE and ARANGE do not match.');
                end
                acoef = getcoef(afdcell{j});
                if ~size(acoef,2) == 1
                    error('AFDCELL is not a single function');
                end
            end
        end
    end
else
    afdcell = {};
end

%  set up the default unit functions for UFDCELL

if nargin == 3      
    ubasis = create_constant_basis(wrange);
    unames = wnames;
    unames{2} = 'forcing fn.';
    unames{3} = [wnames{3}, 'forcing fn.'];
    ufd = fd(1, ubasis, unames);
    for i=1:k
        ufdcell{i} = ufd;
    end
end

%  check UFDCELL

if nargin >= 3
    if ~iscell(ufdcell)
        error('UFDCELL is not a cell object.');
    end
    %  check each cell of UFDCELL
    if ~isempty(ufdcell)
        for i=1:k
            ubasis = getbasis(ufdcell{i});
            urange = getbasisrange(ubasis);
            if any(urange ~= wrange)
                error('WRANGE and URANGE do not match for a UFDCELL.');
            end
            ucoef = getcoef(ufdcell{i});
            if ~size(ucoef,2) == 1
                error('A UFDCELL is not a single function');
            end
        end
    end
else
    ufdcell = {};
end

%  set up the Lfd object

Lfdobj.nderiv  = m;
Lfdobj.wfdcell = wfdcell;
Lfdobj.afdcell = afdcell;
Lfdobj.ufdcell = ufdcell;

Lfdobj = class(Lfdobj, 'Lfd');

