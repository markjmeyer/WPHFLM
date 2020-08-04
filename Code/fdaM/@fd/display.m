function display(fd)
%  DISPLAY  Display a functional data object.

%  Last modified 27 May 1999

  if strcmp(class(fd), 'fd')
    fprintf('Dimensions of data:\n');
    fdnames = getnames(fd);
    fprintf(['   ',fdnames{1},'\n']);
    fprintf(['   ',fdnames{2},'\n']);
    fprintf(['   ',fdnames{3},'\n']);
  else
    error('Argument not a functional data object');
  end

  fbdo = getbasis(fd);
  fprintf('\nBasis:\n');
  fprintf(['  Type: ', getbasistype(fbdo),'\n']);
  rangeval = getbasisrange(fbdo);
  fprintf(['  Range: ',num2str(rangeval(1)), ...
           ' to ',     num2str(rangeval(2)),'\n']);
  fprintf(['  Number of basis functions: ', ...
           num2str(getnbasis(fbdo)),     '\n']);
  type = getbasistype(fbdo);
  if  strcmp(type, 'fourier')
    fprintf(['  Period: ',  ...
             num2str(getbasispar(fbdo)),'\n']);
  end
  if  strcmp(type, 'cosine')
    fprintf(['  Period: ',  ...
             num2str(getbasispar(fbdo)),'\n']);
  end
  if strcmp(type, 'bspline')
    fprintf('  Interior knots\n');
    disp(getbasispar(fbdo));
  end
  if strcmp(type, 'poly')
    fprintf('  Polynomial coefficients\n');
    disp(getbasispar(fbdo));
  end
  if strcmp(type, 'polyg')
    fprintf('  Argument values\n');
    disp(getbasispar(fbdo));
  end
  if strcmp(type, 'expon')
    fprintf('  Rate coefficients\n');
    disp(getbasispar(fbdo));
  end
  if strcmp(type, 'monom')
    fprintf('  Exponents\n');
    disp(getbasispar(fbdo));
  end
  if strcmp(type, 'const')
    fprintf('  No parameters\n');
  end

