
function fdtype = use_proper_basis(fdtype)
%  USE_PROPER_BASIS recognizes type of basis by use of several variant spellings

  %  Last modified 22 January 2002

  switch fdtype
  case 'Fourier'
     fdtype = 'fourier';
  case 'fourier'
     fdtype = 'fourier';
  case 'Fou'
     fdtype = 'fourier';
  case 'fou'
     fdtype = 'fourier';

  case 'bspline'
    fdtype = 'bspline';
  case 'Bspline'
    fdtype = 'bspline';
  case 'Bsp'
    fdtype = 'bspline';
  case 'bsp'
    fdtype = 'bspline';

  case 'power'
    fdtype = 'power';
  case 'pow'
    fdtype = 'power';

  case 'polyg'
    fdtype = 'polyg';
  case 'polygon'
    fdtype = 'polyg';
  case 'polygonal'
    fdtype = 'polyg';

  case 'exp'
    fdtype = 'expon';
  case 'expon'
    fdtype = 'expon';
  case 'exponential'
    fdtype = 'expon';

  case 'con'
    fdtype = 'const';
  case 'const'
    fdtype = 'const';
  case 'const'
    fdtype = 'const';

  otherwise
    fdtype = 'unknown';

  end
