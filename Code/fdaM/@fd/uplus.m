function plusfd = uplus(fd)
% Unary plus of functional data object.

%  last modified 29 June 1998

  if ~(isa_fd(fd)
    error('Argument is not a functional data object.');
  end

  plusfd  = fd;


