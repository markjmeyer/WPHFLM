function [gval, tval, fval] = mongrad1(x, Wfd, basiscell, EPS, JMIN, JMAX)
% MONGRAD1 evaluates the gradient of a single monotone function of the form
  %             h(x) = (D^{-1} exp Wfd)(x)
  %  where  D^{-1} means taking the indefinite integral.
  %  The interval over which the integration takes places is defined in
  %  the basis object in Wfd.
  %  If WFD has multiple curves, use MONGRAD instead.
  %  Arguments:
  %  X         ...  Vector of argument values at which gradient is evaluated
  %  WFD       ...  Functional data object defining monotone function
  %  BASISCELL ...  A cell object containing basis function values
  %  EPS       ...  Relative error needed for convergence
  %  JMIN      ...  Minimum number of step halving steps
  %  JMAX      ...  Maximum number of step halving steps
  %
  %  Returns:
  %  GVAL  ... values of derivatives in NBASIS cols
  %  TVAL  ... Arguments used for trapezoidal approximation to integral
  %  FVAL  ... Values of exp Wfd corresponding to TVAL

  %  Last modified 16 February 2001

  if nargin < 2,  error('There are less than two arguments');  end
  
  %  set some constants
  
  if nargin < 6,  EPS  = 1e-5;  end
  if nargin < 5,  JMIN = 11;    end
  if nargin < 4,  JMAX = 15;    end
  if JMIN   < 5,  JMIN =  5;    end
  
  %  get coefficient matrix and check it
  
  coef  = getcoef(Wfd);
  coefd = size(coef);
  ndim  = length(coefd);
  if ndim > 1 & coefd(2) ~= 1 
    error('WFD is not a single function');
  end

  basis    = getbasis(Wfd);
  rangeval = getbasisrange(basis);
  nbasis   = getnbasis(basis);
  onebas   = ones(1,nbasis);
  width    = rangeval(2) - rangeval(1);

  %  set up first iteration

  JMAXP = JMAX + 1;
  h     = ones(JMAXP,1);
  h(2)  = 0.25;
  %  matrix SMAT contains the history of discrete approximations to the
  %    integral
  smat = zeros(JMAXP,nbasis);
  %  array TVAL contains the argument values used = the approximation
  %  array FVAL contains the integral values at these argument values,
  %     rows corresponding to argument values
  %  the first iteration uses just the endpoints
  j    = 1;
  tj   = rangeval;
  tval = tj;
  if nargin == 3
     if isempty(basiscell{j})
        bmat = eval(basis, tj);
        basiscell{j} = bmat;
     else
        bmat = basiscell{j};
     end
  else
     bmat = eval(basis, tj);
  end
  fx   = exp(bmat*coef);
  grad = (fx * onebas).*bmat;
  fval = grad;
  smat(1,:)  = width*sum(grad)./2;
  tnm = 0.5;
  %  now iterate to convergence
  for j = 2:JMAX
    tnm  = tnm*2;
    del  = width/tnm;
    hdel = del/2;
    tj   = (rangeval(1)+hdel):del:(rangeval(2)-hdel);
    tval = [tval, tj];
    if nargin == 3
       if isempty(basiscell{j})
          bmat = eval(basis, tj);
          basiscell{j} = bmat;
       else
          bmat = basiscell{j};
       end
    else
       bmat = eval(basis, tj);
    end
    fx   = exp(bmat*coef);
    grad = (fx*onebas).*bmat;
    fval = [fval; grad];
    smat(j,:) = (smat(j-1,:) + del.*sum(grad))./2;
    if j >= max([JMIN,5])
      ind = (j-4):j;
      [ss,dss] = polintmat1(h(ind),smat(ind,:),0);
      if all(abs(dss) < EPS*max(abs(ss))) | j == JMAX
        % successful convergence
        % sort argument values and corresponding function values
        [tval,ordind] = sort(tval);
        fval  = fval(ordind,:);
        % set up partial integral values
        integfval = (tval(2) - tval(1)).*cumtrapz(fval);
        gval = zeros(length(x),nbasis);
        gval = interp1(tval, integfval, x, 'cubic');
        return;
      end
    end
    h(j+1) = 0.25*h(j);
  end
 %error(['No convergence after ',num2str(JMAX),' steps in MONGRAD1']);
