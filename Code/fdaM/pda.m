function wfd = pda(fdobj, difeorder, wbasisfd, ...
                   lambda, wfd0, estimate, n)
%  PDA  Compute the basis function expansion of the
%    estimate of the coefficient functions
%    for a linear differential operator of degree NORDER that comes as
%    close as possible in a least squares sense to annihilating the NREP
%    functions in array Y sampled at the argument values in array X
%  Arguments:
%  FDOBJ     ...  functional data object
%  DIFEORDER ...  order of the linear differential operator, that is, the order
%                 of the highest derivative.
%  WBASISFD  ...  basis object for weight functions
%  LAMBDA    ...  penalty parameter for penalizing the departure of the
%                 estimated weight functions from those defined = WFD0
%  WFD0      ...  A specification of a functional data object that is used for
%                 those weight functions not estimated, or as target functions
%                 toward which the estimated weight functions are smoothed. WFD0
%                 can either be a vector of DIFEORDER constants, or a functional
%                 data object with the same structure as WFN that is returned
%                 by this function.
%  ESTIMATE  ...  logical array of length DIFEORDER, if a value is T, the
%                 corresponding coefficient function is estimated, otherwise
%                 the target value is used.
%  N         ...  number of sampling points for numerical integration

%  Returns:
%  WFN       ...  estimated weight functional data object.  It has DIFEORDER
%                 replicates, and the corresponding linear differential operator
%                 is  L = w_1 I + w_2 D + ... + w_DIFEORDER D^DIFEORDER-1end + D^DIFEORDER

%  last modified 18 April 2001

  basisfd  = getbasis(fdobj);
  nbasis   = getnbasis(basisfd);
  rangeval = getbasisrange(basisfd);

  if nargin < 7
    n = 5*nbasis;
  end

  if nargin < 6
    estimate = ones(difeorder, 1);
  end

  if nargin < 5
    wfd0 = zeros(1,difeorder);
  end

  if ~strcmp(class(wfd0), 'fd') & ~isnumeric(wfd0)
      error('WFD0 is neither a vector nor a FD object');
  end
  if ~strcmp(class(wfd0), 'fd') 
      if length(wfd0) ~= difeorder 
          error('WFN0 is a vector of incorrect length');
      end
      wbasis0 = create_constant_basis(rangew);
      wfd0    = fd(wfd0, wbasis0);
  end

  if nargin < 4
    lambda = zeros(difeorder,1);
  end

  ncoef = sum(estimate);

  coef  = getcoef(fdobj);
  coefd = size(coef);
  ndim  = length(coefd);
  ncurve = coefd(2);
  if ndim == 3
    nvar = coefd(3);
  else
    nvar = 1;
  end

  typew   = getbasistype(wbasisfd);
  nbasisw = getnbasis(wbasisfd);
  rangew  = getbasisrange(wbasisfd);

  if any(rangew ~= rangeval)
    error('Weight function range not equal to FD range');
  end

  if strcmp(typew, 'bspline')
    params   = getbasispar(wbasisfd);
    nbreaksw = length(params);
    difeorderw  = nbasisw - nbreaksw;
  end

  delta = (rangew(2)-rangew(1))/(n-1);
  x     = rangew(1):delta:rangew(2);

  nordp1   = difeorder + 1;
  if ndim <= 2
    yarray = zeros([n,ncurve,nordp1]);
    for j = 0:difeorder
      yarray(:,:,j+1) = eval(fdobj, x, j);
    end
  else
    yarray = zeros([n,ncurve,nvar,nordp1]);
    for j = 0:difeorder
      yarray(:,:,:,j+1) = eval(fdobj, x, j);
    end
  end

  basismat = getbasismatrix(x, wbasisfd);

  if ndim <= 2

  %  --------------  univariate case  -------------------------

    dmat = zeros(nbasisw,difeorder);
    if ncurve == 1
      DV = delta*yarray(:,:,nordp1);
      IV = zeros(n,ncoef*nbasisw);
      onenb = ones(1, nbasisw);
      mi = 0;
      for i = 1:difeorder
        if estimate(i)
          mi = mi + 1;
          index = (1 + (mi-1)*nbasisw):(mi*nbasisw);
          IV(:,index) = delta.*(squeeze(yarray(:,1,i))*onenb).*basismat;
        end
      end
      dvec   = -IV\DV;
    else
      mi   = 0;
      mij  = 0;
      Swgt = zeros(n,ncoef);
      Rwgt = zeros(n,ncoef*(ncoef+1)/2);
      for i = 1:difeorder
        if estimate(i)
          mi = mi + 1;
          index = (1 + (mi-1)*nbasisw):(mi*nbasisw);
          Swgt(:,mi) = mean((yarray(:,:,i).*yarray(:,:,nordp1))')';
          mj = 0;
          for j = 1:i
            if estimate(j)
              mij = mij + 1;
              Rwgt(:,mij) = mean((yarray(:,:,i).*yarray(:,:,j))')';
            end
          end
        end
      end

      [Dmat, Cmat] = SRsetup(ncoef, nbasisw, Swgt, Rwgt, basismat);
      if any(lambda > 0)
        Hmat   = getbasispenalty(wbasisfd);
        for i = 1:ncoef
          index = (1 + (i-1)*nbasisw):(i*nbasisw);
          if lambda(i) > 0
            Cmat(index,index) = Cmat(index,index) - lambda(i).*Hmat;
            if any(getcoef(wfd0))
              Dmat(index,1) = Dmat(index,1) + ...
                lambda(i).*inprod(basisfd,wfd0(i));
            end
          end
        end
      end
      dvec = -Cmat\Dmat;
    end

    mi  = 0;
    for i = 1:difeorder
      if estimate(i)
        mi = mi + 1;
        index = (1 + (mi-1)*nbasisw):(mi*nbasisw);
        dmat(:,i) = dvec(index);
      end
    end

  %  --------------  multivariate case  -------------------------

  else
    dmat = zeros(nbasisw,difeorder,nvar);
    for ivar=1:nvar
      if ncurve == 1
        DV = delta*yarray(:,:,ivar,nordp1);
        IV = zeros(n,ncoef*nbasisw);
        onenb = ones(1, nbasisw);
        mi = 0;
        for i = 1:difeorder
          if estimate(i)
            mi = mi + 1;
            index = (1 + (mi-1)*nbasisw):(mi*nbasisw);
            IV(:,index) = delta.*(squeeze(yarray(:,1,ivar,i))*onenb).*basismat;
          end
        end
        dvec   = -IV\DV;
      else
        mi   = 0;
        mij  = 0;
        Swgt = zeros(n,ncoef);
        Rwgt = zeros(n,ncoef*(ncoef+1)/2);
        for i = 1:difeorder
          if estimate(i)
            mi = mi + 1;
            index = (1 + (mi-1)*nbasisw):(mi*nbasisw);
            Swgt(:,mi) = mean((yarray(:,:,ivar,i).*yarray(:,:,ivar,nordp1))')';
            mj = 0;
            for j = 1:i
              if estimate(j)
                mij = mij + 1;
                Rwgt(:,mij) = mean((yarray(:,:,ivar,i).*yarray(:,:,ivar,j))')';
              end
            end
          end
        end

        [Dmat, Cmat] = SRsetup(ncoef, nbasisw, Swgt, Rwgt, basismat);
        if any(lambda > 0)
          Hmat = getbasispenalty(wbasisfd);
          for i = 1:ncoef
            index = (1 + (i-1)*nbasisw):(i*nbasisw);
            if lambda(i) > 0
              Cmat(index,index) = Cmat(index,index) - lambda(i).*Hmat;
              if any(getcoef(wfd0))
                Dmat(index,1) = Dmat(index,1) + ...
                   lambda(i).*inprod(basisfd,wfd0(:,i));
              end
            end
          end
        end
        dvec = -Cmat\Dmat;
      end
      mi  = 0;
      for i = 1:difeorder
        if estimate(i)
          mi = mi + 1;
          index = (1 + (mi-1)*nbasisw):(mi*nbasisw);
          dmat(:,i,ivar) = dvec(index);
        end
      end
    end
  end

  wfdnames = getnames(fdobj);
  wfdnames{2} = 'Weight functions';
  wfdnames{3} = 'Weight value';
  wfd = fd(dmat, wbasisfd, wfdnames);

%  ----------------------------------------------------------------------------

function [Dmat, Cmat] = SRsetup(ncoef, nbasisw, Swgt, Rwgt, basismat)
% SRSETUP sets up coefficient matrices for basis expansion of weight functions

  Dmat = zeros(ncoef*nbasisw, 1);
  Cmat = zeros(ncoef*nbasisw, ncoef*nbasisw);
  [n, m1] = size(basismat);
  m  = 0;
  onen  = ones(n,1);
  onem1 = ones(1,m1);
  ind1n = [1,n];
  index = 1:nbasisw;
  for i = 1:ncoef
    indexi = index + (i-1)*nbasisw;
    temp     = basismat .* (Swgt(:,i) * onem1);
    temp(ind1n,:) = temp(ind1n,:)./2;
    Dmat(indexi) = temp' * onen;
    for j = 1:i
      m = m + 1;
      indexj = index + (j-1)*nbasisw;
      temp     = basismat .* (Rwgt(:,m) * onem1);
      temp(ind1n,:) = temp(ind1n,:)./2;
      Cmat(indexi,indexj) = temp' * basismat;
      if i ~= j
        Cmat(indexj,indexi) = Cmat(indexi,indexj);
      end
    end
  end

