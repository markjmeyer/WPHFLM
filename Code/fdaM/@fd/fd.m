function fdobj = fd(coef, basisobj, fdnames)
%  FD   Creates a functional data object.
%    A functional data object consists of a basis for expanding a functional
%    observation and a set of coefficients defining this expansion.
%    The basis is contained=a 'basis' object; that is, a realization
%    of the 'basis' class.
%  Arguments
%  COEF ... An array containing coefficient values for the expansion of each
%             set of function values=terms of a set of basis functions.
%           If COEF is a three-way array, then the first dimension
%             corresponds to basis functions, the second to replications,
%             and the third to variables.
%           If COEF is a matrix, it is assumed that there is only
%             one variable per replication, and then
%                 rows    correspond to basis functions
%                 columns correspond to replications
%           If COEF is a vector, it is assumed that there is only one
%             replication and one variable.
%  BASISOBJ ... a functional data basis object
%  FDNAMES  ... A cell of length 3 with members containing
%               1. a name for the argument domain, such as 'Time'
%               2. a name for the replications or cases
%               3. a name for the function
%  Returns:
%  FD ... a functional data object

%  last modified 29 March 99

%  superiorto('double', 'sparse', 'struct', 'cell', 'char', 'inline', ...
%             'basis');

  superiorto('double', 'struct', 'cell', 'char', 'inline', ...
             'basis');

  if nargin == 0
    fdobj.coef     = zeros(1,1);
    fdobj.basisobj = basis('bspline',[0,1],3,0.5);
    fdnames{1} = 'time';
    fdnames{2} = 'reps';
    fdnames{3} = 'values';
    fdobj.fdnames  = fdnames;
    fdobj = class(fdobj, 'fd');
    return;
  end

  if isa(coef, 'fd')
    fdobj = coef;
    return; 
  end

  coefd = size(coef);
  ndim  = length(coefd);

  if (ndim > 3)
    error('First argument not of dimension 1, 2 or 3');
  end

  if ~strcmp(class(basisobj), 'basis')
    error('Argument BASISOBJ must be of basis class');
  end

  nbasis = getnbasis(basisobj);
  if coefd(1) ~= nbasis
     error( ...
     'Number of coefficients does not match number of basis functions.');
  end

  if (ndim > 1)
    ncurve = coefd(2);
  else
    ncurve = 1;
  end
  if (ndim > 2)
    nvar = coefd(3);
  else
    nvar = 1;
  end

  %  set up default fdnames

  if nargin < 3
    fdnames{1} = 'time';
    fdnames{2} = 'reps';
    fdnames{3} = 'values';
  end

  fdobj.coef     = coef;
  fdobj.basisobj = basisobj;
  fdobj.fdnames  = fdnames;

  fdobj = class(fdobj, 'fd');

