function wfd = pdan(fdobj, difeorder, wbasisfd, estimate, ...
                    lambda, wfd0, n)
%  PDAN  Compute the basis function expansion of the
%    estimate of the coefficient functions
%    for a linear differential operator of degree NORDER that comes as
%    close as possible in a least squares sense to annihilating the NREP
%    functions in multivariate functional data object FDOBJ.
%    This version differs from function PDA in that the operator for each
%    includes contributions from the derivatives of all other functions,
%    as would be appropriate for 3D curves where the coordinate system is
%    only known to within an linear transformation.  
%  Arguments:
%  FDOBJ     ...  functional data object, or a matrix of function values
%  DIFEORDER ...  order of the linear differential operator, that is, the order
%                 of the highest derivative.
%  WBASISFD  ...  basis object for weight functions
%  ESTIMATE  ...  logical array of length DIFEORDER, if a value is T, the
%                 corresponding coefficient function is estimated, otherwise
%                 the target value is used.
%  LAMBDA    ...  penalty parameter for penalizing the departure of the
%                 estimated weight functions from those defined = WFD0
%  WFD0      ...  A specification of a functional data object that is used for
%                 those weight functions not estimated, or as target functions
%                 toward which the estimated weight functions are smoothed. WFD0
%                 can either be a vector of DIFEORDER constants, or a functional
%                 data object with the same structure as WFN that is returned
%                 by this function.
%  N         ...  number of sampling points for numerical integration

%  Returns:
%  WFN       ...  estimated weight functional data object.  It has DIFEORDER
%                 replicates, and the corresponding linear differential operator
%                 is  L = w_1 I + w_2 D + ... + w_DIFEORDER D^DIFEORDER-1 
%                       + D^DIFEORDER

%  last modified 8 March 2001

  if isa(fdobj, 'fd')
      basisfd  = getbasis(fdobj);
      nbasis   = getnbasis(basisfd);
      rangeval = getbasisrange(basisfd);
  else
      nbasis = 0;
  end

  if nargin < 7
      if isa(fdobj, 'fd')
          n = 5*nbasis; 
      else
          n = size(fdobj,1);
      end
  end
  if nargin < 6, wfd0 =    zeros(difeorder, 1); end
  if nargin < 5, lambda =  zeros(difeorder, 1); end
  if nargin < 4, estimate = ones(difeorder, 1); end
  
  ncoef = sum(estimate);
  coefindx = 1:ncoef;
  if ncoef < difeorder
    m = 0;
    for j=1:difeorder
      if estimate(j)
        m = m + 1;
        coefindx(m) = j;
      end
    end
  end

  if isa(fdobj, 'fd')
      coef  = getcoef(fdobj);
  else
      coef  = fdobj;
  end
  coefd = size(coef);
  ndim  = length(coefd);
  ncurve = coefd(2);
  if ndim >= 3
    nvar = coefd(3);
  else
    nvar = 1;
  end
  nwgt = nvar.*ncoef;

  typew   = getbasistype(wbasisfd);
  nbasisw = getnbasis(wbasisfd);
  rangew  = getbasisrange(wbasisfd);
  indexw  = 1:nbasisw;
  if isa(fdobj, 'fd')
      if any(rangew ~= rangeval)
          error('Weight function range not equal to FD range');
      end
      delta = (rangew(2)-rangew(1))/(n-1);
      x     = rangew(1):delta:rangew(2);
  else
      x = linspace(rangew(1),rangew(2),n)';
      wbasismat = getbasismatrix(x, wbasisfd);
      onebas = ones(1,nbasisw);
  end


  %  compute expected products of functions for range of derivatives
  
  rarray = zeros(n, nwgt.*(nwgt+1)./2);
  m = 0;
  for j1=coefindx
    for l1=1:nvar
        if isa(fdobj, 'fd')
            ymat1 = eval(fdobj(:,l1), x, j1-1);
        else
            ymat1 = squeeze(fdobj(:,:,l1,j1));
        end
        for j2=coefindx;
            if j2 <= j1
                if j1 == j2, l2lim = l1; else, l2lim = nvar;  end
                for l2=1:l2lim
                    if isa(fdobj, 'fd')
                        ymat2 = squeeze(eval(fdobj(:,l2), x, j2-1));
                    else
                        ymat2 = squeeze(fdobj(:,:,l2,j2));
                    end
                    m = m + 1;
                    if ncurve > 1
                        rarray(:,m) = mean((ymat1.*ymat2)')';
                    else
                        rarray(:,m) = ymat1.*ymat2;
                    end
                end
            end
        end
    end
  end
  
  if isa(fdobj, 'fd')
      rfd = data2fd(rarray, x, basisfd);
  end

  sarray = zeros(n, nwgt, nvar);
  for l2=1:nvar
      if isa(fdobj, 'fd')
          ymat2 = squeeze(eval(fdobj(:,l2), x, difeorder));
      else
          ymat2 = squeeze(fdobj(:,:,l2,difeorder+1));
      end
      m = 0;
      for j1=coefindx
          for l1=1:nvar
              m = m + 1;
              if isa(fdobj, 'fd')
                  ymat1 = squeeze(eval(fdobj(:,l1), x, j1-1));
              else
                  ymat1 = squeeze(fdobj(:,:,l1,j1));
              end
              if ncurve > 1
                  sarray(:,m,l2) = mean((ymat1.*ymat2)')';
              else
                  sarray(:,m,l2) = ymat1.*ymat2;
              end
          end
      end
  end
  
  if isa(fdobj, 'fd')
      sfd = data2fd(sarray, x, basisfd);
  end
             
  npars = ncoef*nbasisw; 
  onefd = fd(1,create_constant_basis(rangew));
  Cmat = zeros(npars, npars);
  Dmat = zeros(npars, nvar);
  index1 = indexw;
  m1 = 0;
  m = 0;
  for j1=1:ncoef
    for l1=1:nvar
      m1 = m1 + 1;
      for l2=1:nvar
          if isa(fdobj, 'fd')
              Dmat(index1,l2) = inprod(onefd,wbasisfd,0,0,rangew,sfd(m1,l2))';
          else
              Dmat(index1,l2) = mean(wbasismat.*(sarray(:,m1,l2)*onebas))';
          end
      end
      index2 = indexw;
      for j2 = 1:j1
        if j2 == j1, l2lim = l1; else, l2lim = nvar; end
        for l2=1:l2lim
            m = m + 1;
            if isa(fdobj, 'fd')
                Cmat(index1,index2) = inprod(wbasisfd,wbasisfd,0,0,rangew,rfd(m));
            else
                Cmat(index1,index2) = ((wbasismat.*(rarray(:,m)*onebas))'*wbasismat)./n;
            end
            Cmat(index2,index1) = Cmat(index1,index2);
          index2 = index2 + nbasisw;
        end
      end
      index1 = index1 + nbasisw;
    end
  end

  if any(lambda > 0)
    if ~strcmp(class(wfn0), 'fd') & isnumeric(wfn0)
      if length(wfn0) ~= difeorder 
        error('WFN0 is a vector of =correct length');
      end
      wbasis0 = create_constant_basis(rangew);
      wfn0 = fd(wfn0, wbasisfd0);
    else
      error('WFN0 is neither a vector nor a FD object');
    end
    Gmat   = getpenaltymatrix(basisfd);
    index = indexw;
    for j = 1:ncoef
      Hmat = inprod(basisfd,wfn0(:,j));
      for l=1:nvar
        if lambda(i) > 0
          Cmat(index,index) = Cmat(index,index) - lambda(i).*Gmat;
          for l2=1:nvar;
            Dmat(index,l2) = Dmat(index,l2) + lambda(i).*Hmat;
          end
        end
        index = index + nbasisw;
      end
    end
  end
  
  dvec = -Cmat\Dmat;
  
  coefmat = zeros(nbasisw,nvar.*difeorder,nvar);
  index = indexw;
  for j=coefindx
    for l1=1:nvar;
      m = (j-1).*nvar + l1;
      for l2=1:nvar;
        coefmat(:,m,l2) = dvec(index,l2);
      end
      index = index + nbasisw;
    end
  end

  wfdnames{1} = 'Time';
  wfdnames{2} = 'Weight functions';
  wfdnames{3} = 'Weight value';
  wfd = fd(coefmat, wbasisfd, wfdnames);

