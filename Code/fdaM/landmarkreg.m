function regstr = landmarkreg(fdobj, ximarks, x0marks, wbasis, Lfdobj, ...
                               wLfdobj, lambda, monwrd)
%  LANDMARKREG ... Register curves using landmarks.
%  Arguments:
%  FDOBJ        ...  a functional data object for curves to be registered
%  XIMARKS   ...  N by NL array of times of interior landmarks for
%                 each observed curve
%  XOMARKS   ...  vector of length NL of times of interior landmarks for
%                 target curve
%  LFDOBJ    ...  integer or functional data object defining derivative
%                 or linear differential operator value to be registered.
%                 Default:  0
%  WBASIS    ... an optional basis object used for estimating warp
%                 functions.  If not supplied the basis for FDOBJ is used.
%                 Default:  Order 4 B-spline with knots at landmarks.
%  WLFDOBJ   ...  integer or functional data object defining derivative
%                 or linear differential operator value to be penalized
%                 in estimating the warping function
%                 Default:  2
%  LAMBDA    ...  smoothing parameter used in estimating the warping function
%                 Default:  1e-10
%  MONWRD    ...  If 1, warping functions are estimated by monotone smoothing,
%                 otherwise by regular smoothing.  The latter is faster, but
%                 not guaranteed to produce a strictly monotone warping 
%                 function.  If MONWRD is 0 and an error message results 
%                 indicating nonmonotonicity, rerun with MONWRD = 1.
%                 Default:  1
%  Returns:
%  REGSTR  ...  A struct object with fields
%    REGSTR.REGFD  ...  A functional data object for the registered curves
%    REGSTR.WARPFD ...  A functional data object for the warping functions
%    REGSTR.WFD    ...  A Functional data object for function W defining 
%                         warping fns

%  last modified 16 January 2003

%  check first two arguments for being functional data objects

if ~isa_fd(fdobj)
    error ('Argument FDOBJ is not a functional data object.');
end

%  extract information from curve functional data object and its basis

coef       = getcoef(fdobj);
coefd      = size(coef);
ndim       = length(coefd);
ncurve     = coefd(2);
fdbasisobj = getbasis(fdobj);
fdnames    = getnames(fdobj);
fdnbasis   = getnbasis(fdbasisobj);
rangeval   = getbasisrange(fdbasisobj);

%  check landmarks

ximarksd = size(ximarks);
if ximarksd(1) ~= ncurve
    error('Number of rows of second argument wrong.');
end
nlandm = ximarksd(2);

if nargin < 3
    x0marks = mean(ximarks);
end
if (length(x0marks) ~= nlandm)
    error(['Number of target landmarks not equal to ', ...
            'number of curve landmarks.']);
end
if size(x0marks,2) == 1, x0marks = x0marks';  end

%  set default argument values

if nargin < 8, monwrd  = 1;           end
if nargin < 7, lambda  = 1e-10;       end
if nargin < 6, wLfdobj = int2Lfd(2);  end
if nargin < 5, Lfdobj  = int2Lfd(0);  end
if nargin < 4  
    %  default basis is order 4 with knots at landmarks
    breaks = [rangeval(1), x0marks, rangeval(2)];
    wbasis = create_bspline_basis(rangeval, nlandm+4, 4, breaks); 
end

%  check linear differential operators

Lfdobj  = int2Lfd(Lfdobj);
wLfdobj = int2Lfd(wLfdobj);

%  check landmark values

if any(ximarks <= rangeval(1)) | any(ximarks >= rangeval(2))
    error('Some landmark values are not within the range.');
end

n   = min([101,10*fdnbasis]);
x   = linspace(rangeval(1),rangeval(2),n)';
wtn = ones(n,1);
y   = eval_fd(fdobj, x, Lfdobj);
yregmat = y;
hfunmat = zeros(n,ncurve);
lambda  = max([lambda,1e-10]);

xval = [rangeval(1),x0marks,rangeval(2)]';
nval = length(xval);
wval = ones(nval,1);
nwbasis = getnbasis(wbasis);
Wfd0    = fd(zeros(nwbasis,1),wbasis);
Wcoef   = zeros(nwbasis,ncurve);

for icurve = 1:ncurve
    
    %  set up landmark times for this curve
    disp(['Curve ',num2str(icurve)]);
    yval = [rangeval(1),ximarks(icurve,:),rangeval(2)]';
    
    %  smooth relation between this curve's values and target's values
    %  warpfd is the functional data object for the warping functions
    %  h is a vector of warping function values corresponding to arguments x
    if monwrd
        %  use monotone smoother
        Wfd = warpsmth(xval, yval, wval, Wfd0, wLfdobj, lambda);
        h   = monfn(x, Wfd);
        %  normalize h
        h = h.*(rangeval(2)-rangeval(1))./(h(n)-h(1));
        h = h - h(1) + rangeval(1);
        h(1) = rangeval(1);
        h(n) = rangeval(2);
        warpfd = data2fd(h, x, wbasis);
        wcoefi = getcoef(Wfd);
        Wcoef(:,icurve) = wcoefi;
    else
        %  use regular smoother
        smoothlist = smooth_basis(yval, xval, wbasis, wval, wLfdobj, lambda);
        warpfd     = smoothlist.fdobj;
        Wfd        = warpfd;
        %  set up warping function sampling values
        h          = eval_fd(warpfd, x);
        %  normalize h
        h = h.*(rangeval(2)-rangeval(1))./(h(n)-h(1));
        h = h - h(1) + rangeval(1);
        h(1) = rangeval(1);
        h(n) = rangeval(2);
        %  check for monotonicity because regular smooth may not be monotone
        deltah = diff(h);
        if any(deltah <= 0) ...
                error(['Non-increasing warping function estimated for curve',...
                    num2str(icurve),'.\n',...
                    'Try setting MONWRD to 1.']);
        end
    end
    hfunmat(:,icurve) = h;
    
    %  compute h-inverse in order to register curves
    
    if monwrd
        Wfdinv = fd(-wcoefi,wbasis);
        Wfdinv = warpsmth(h, x, wtn, Wfdinv, wLfdobj, lambda);
        hinv   = monfn(x, Wfdinv);
        hinv   = hinv.*(rangeval(2)-rangeval(1))./(hinv(n)-hinv(1));
        hinv   = hinv - hinv(1) + rangeval(1);
        hinv(n) = rangeval(2);
        hinv(1) = rangeval(1);
    else
        smoothlist = smooth_basis(x, h, wbasis, wtn, wLfdobj, lambda);
        hinvfd = smoothlist.fdobj;
        hinv   = eval_fd(hinvfd, x);
        hinv   = hinv.*(rangeval(2)-rangeval(1))./(hinv(n)-hinv(1));
        hinv   = hinv - hinv(1) + rangeval(1);
        hinv(n) = rangeval(2);
        hinv(1) = rangeval(1);
        %  check for monotonicity because regular smooth may not be monotone
        deltahinv = diff(hinv);
        if any(deltahinv <= 0) ...
                error(['Non-increasing inverse warping function estimated for curve',...
                    num2str(icurve)]);
        end
    end
    
    %  compute registered curves
    
    if ndim == 2
        %  single variable case
        smoothlist = smooth_basis(y(:,icurve), hinv, fdbasisobj, wtn, ...
            int2Lfd(2), 1e-10);
        yregfd = smoothlist.fdobj;
        yregmat(:,icurve) = eval_fd(yregfd, x);
    end
    if ndim == 3
        %  multiple variable case
        for ivar = 1:nvar
            % evaluate curve as a function of h at sampling points
            smoothlist = smooth_basis(y(:,icurve,ivar), hinv, fdbasisobj, wtn, ...
                int2Lfd(2), 1e-10);
            yregfd = smoothlist.fdobj;
            yregmat(:,icurve,ivar) = eval_fd(yregfd, x);
        end
    end
end

%  create functional data objects for the registered curves

yregcoef      = project_basis(yregmat, x, fdbasisobj);
regfdnames    = getnames(fdobj);
regfdnames{3} = ['Registered ',regfdnames{3}];
regfd         = fd(yregcoef, fdbasisobj, regfdnames);

%  create functional data objects for the warping functions

warpcoef       = project_basis(hfunmat, x, wbasis);
warpfdnames    = getnames(fdobj);
warpfdnames{3} = ['Warped ',warpfdnames{3}];
warpfd         = fd(warpcoef, wbasis, warpfdnames);

%  create functional data object for Wfd if monwrd is true

Wfd = fd(Wcoef,wbasis);

%  set up REGSTR  ...  A struct object with fields
%    REGSTR.REGFD  ...  A functional data object for the registered curves
%    REGSTR.WARPFD ...  A functional data object for the warping functions
%    REGSTR.WFD    ...  A Functional data object for function W defining 
%                         warping fns

regstr.regfd  = regfd;
regstr.warpfd = warpfd;
regstr.Wfd    = Wfd;


