function display(basis)
%  DISPLAY  Display a functional data basis object.

%  Last modified 22 January 2002

  if ~strcmp(class(basis), 'basis')
    error('Argument not a functional data object');
  end

  fprintf('\nBasis:\n');
  type = getbasistype(basis);
  fprintf(['  Type: ', type,'\n']);
  fprintf(['  Range: ', num2str(basis.rangeval(1)), ...
           ' to ',      num2str(basis.rangeval(2)),'\n']);
  fprintf(['  Number of basis functions: ', ...
                        num2str(basis.nbasis),     '\n']);
  if  strcmp(type, 'fourier')
    fprintf(['  Period: ',num2str(basis.params),'\n']);
  end
  if strcmp(type, 'bspline')
    fprintf('  Interior knots\n');
    disp(basis.params);
  end
  if strcmp(type, 'polyg')
    fprintf('  Argument values\n');
    disp(basis.params);
  end
  if strcmp(type, 'expon')
    fprintf('  Rate coefficients\n');
    disp(basis.params);
  end
  if strcmp(type, 'power')
    fprintf('  Exponents\n');
    disp(basis.params);
  end
  if strcmp(type, 'const')
    fprintf('  No parameters\n');
  end
