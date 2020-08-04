function display(bifd)
%  DISPLAY  Display a functional data object.

%  Last modified 6 July 1998

  if strcmp(class(bifd), 'bifd')
    fprintf('Dimensions of data:\n');
    fprintf(bifd.fdnames);
  else
    error('Argument not a functional data object');
  end

  fbdo = bifd.basisobj;
  fprintf('\nBasis:\n');
  fprintf(['  Type:', fbdo.type,'\n']);
  fprintf(['  Range:',fbdo.rangeval[1],'to',fbdo.rangeval[2],'\n']);
  fprintf(['  Number of basis functions:',  fbdo.nbasis,     '\n']);
  type = getbasistype(fbdo);
  if  strcmp(type, 'fourier')
    fprintf(['  Period:',fbdo.params,'\n']);
  end
  if strcmp(type, 'bspline')
    fprintf('  Interior knots\n');
    disp(fbdo.params);
  end
  if strcmp(type, 'poly')
    fprintf('  Polynomial coefficients\n');
    disp(fbdo.params);
  end
  if strcmp(type, 'polyg')
    fprintf('  Argument values\n');
    disp(fbdo.params);
  end
  if strcmp(type, 'expon')
    fprintf('  Rate coefficients\n');
    disp(fbdo.params);
  end

