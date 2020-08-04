function isafd = isa_fd(fdstr)
%  ISA_FD  checks a struct object for fields 'coef' and 'basisstr'

%  last modified 23 May 1998

  isafd = 1;
  if ~strcmp(class(fdstr), 'fd')
    isafd = 0;
    return;
  end

