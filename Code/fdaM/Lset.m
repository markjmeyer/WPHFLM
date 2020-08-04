function [nderiv, derivcoef, Lconstwrd] = Lset(Lfd)
%  LSET  Determine order and coefficients of linear differential operator

%  last modified 1 July 1998

  derivcoef = 0;
  Lconstwrd = 0;
  switch class(Lfd)
    case 'double'
      if length(Lfd) == 1
        nderiv = Lfd;
        if nderiv ~= round(nderiv)
          error('Order of derivative must be an integer');
        end
        if nderiv < 0
          error('Order of derivative must be 0 or positive');
        end
      else
        error('Order of derivative must be a single number');
      end
    case 'fd'
      if isa_fd(Lfd)
        derivcoef  = getcoef(Lfd);
        derivcoefd = size(derivcoef);
        nderiv     = derivcoefd(2);
        basisobj   = getbasis(Lfd);
        type       = getbasistype(basisobj);
        if strcmp(type,'const')
          Lconstwrd = 1;
        end
      else
        error(['Second argument a struct object', ...
               ' but not a functional object.']);
      end
    otherwise
      error('Second argument neither an integer nor a functional data object');
  end


