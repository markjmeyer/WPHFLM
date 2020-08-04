Windows:

addpath ('c:\matlab\fdaM')
addpath ('c:\matlab\fdaM\examples\growth')

Unix:

addpath('examples/growth')

%  Last modified 15 January 2003

%  -----------------------------------------------------------------------
%                    Berkeley Growth Data
%  -----------------------------------------------------------------------

%  ------------------------  input the data  -----------------------

ncasem = 39;
ncasef = 54;
nage   = 31;

fid = fopen('hgtm.dat','rt');
hgtmmat = reshape(fscanf(fid,'%f'),[nage,ncasem]);

fid = fopen('hgtf.dat','rt');
hgtfmat = reshape(fscanf(fid,'%f'),[nage,ncasef]);

age = [ 1:0.25:2, 3:8, 8.5:0.5:18 ]';
rng = [1,18];

%  --------------  Smooth the data nonmonotonically  --------------
%  This smooth uses the usual smoothing methods to smooth the data,
%  but is not guaranteed to produce a monotone fit.  This may not
%  matter much for the estimate of the height function, but it can
%  have much more serious consequences for the velocity and
%  accelerations.  See the monotone smoothing method below for a
%  better solution, but one with a much heavier calculation overhead.

%  -----------  Create fd objects   ----------------------------
%  A B-spline basis with knots at age values and order 6 is used

knots  = age;
norder = 6;
nbasis = length(knots) + norder - 2;
hgtbasis = create_bspline_basis(rng, nbasis, norder, knots);

hgtmfd = data2fd(hgtmmat, age, hgtbasis);

hgtffd = data2fd(hgtfmat, age, hgtbasis);

%  --- Smooth these objects, penalizing the 4th derivative  --
%  This gives a smoother estimate of the acceleration functions

Lfd    = int2Lfd(4);
lambda = 1e-2;

hgtmfd = smooth_fd(hgtmfd, lambda, Lfd);
hgtffd = smooth_fd(hgtffd, lambda, Lfd);

%  plot data and smooth, residuals, velocity, and acceleration

agefine = linspace(1,18,101)';

%  Males:

hgtmfit = eval_fd(hgtmfd, age);
hgtmhat = eval_fd(hgtmfd, agefine);
velmhat = eval_fd(hgtmfd, agefine, 1);
accmhat = eval_fd(hgtmfd, agefine, 2);

for i = 1:ncasem
  subplot(2,2,1)
  plot(age, hgtmmat(:,i), 'o', agefine, hgtmhat(:,i), '-')
  axis([1,18,60,200]);
  xlabel('Years');  title(['Height for male ',num2str(i)])
  subplot(2,2,2)
  resi = hgtmmat(:,i) - hgtmfit(:,i);
  ind  = find(resi >= -.7 & resi <= .7);
  plot(age(ind), resi(ind), 'o-', rng, [0,0], '--')
  axis([1,18,-.7,.7])
  xlabel('Years');  title('Residuals')
  subplot(2,2,3)
  ind = find(velmhat(:,i) >= 0 & velmhat(:,i) <= 20);
  plot(agefine(ind), velmhat(ind,i), '-')
  axis([1,18,0,20])
  xlabel('Years');  title('Velocity')
  subplot(2,2,4)
  ind = find(accmhat(:,i) >= -6 & accmhat(:,i) <= 6);
  plot(agefine(ind), accmhat(ind,i), '-', rng, [0,0], 'r--')
  axis([1,18,-6,6]),
  xlabel('Years');  title('Acceleration')
  pause
end

% Females:

hgtffit = eval_fd(hgtffd, age);
hgtfhat = eval_fd(hgtffd, agefine);
velfhat = eval_fd(hgtffd, agefine, 1);
accfhat = eval_fd(hgtffd, agefine, 2);

for i = 1:ncasef
  subplot(2,2,1)
  plot(age, hgtfmat(:,i), 'o', agefine, hgtfhat(:,i), '-')
  axis([1,18,60,200]);
  xlabel('Years');  title(['Height for female ',num2str(i)])
  subplot(2,2,2)
  resi = hgtfmat(:,i) - hgtffit(:,i);
  ind  = find(resi >= -.7 & resi <= .7);
  plot(age(ind), resi(ind), 'o-', rng, [0,0], '--')
  axis([1,18,-.7,.7])
  xlabel('Years');  title('Residuals')
  subplot(2,2,3)
  ind = find(velfhat(:,i) >= 0 & velfhat(:,i) <= 20);
  plot(agefine(ind), velfhat(ind,i), '-')
  axis([1,18,0,20])
  xlabel('Years');  title('Velocity')
  subplot(2,2,4)
  ind = find(accfhat(:,i) >= -6 & accfhat(:,i) <= 6);
  plot(agefine(ind), accfhat(ind,i), '-', rng, [0,0], 'r--')
  axis([1,18,-6,6]),
  xlabel('Years');  title('Acceleration')
  pause
end


%  -------  Compute monotone smooths of the data  -----------

%  These analyses use a function written entirely in S-PLUS called
%  smooth.monotone that fits the data with a function of the form
%                   f(x) = b_0 + b_1 D^{-1} exp W(x)
%     where  W  is a function defined over the same range as X,
%                 W + ln b_1 = log Df and w = D W = D^2f/Df.
%  The constant term b_0 in turn can be a linear combinations of covariates:
%                         b_0 = zmat * c.
%  The fitting criterion is penalized mean squared error:
%    PENSSE(lambda) = \sum [y_i - f(x_i)]^2 +
%                     \lambda * \int [L W(x)]^2 dx
%  where L is a linear differential operator defined in argument LFD.
%  The function W(x) is expanded by the basis in functional data object
%  Because the fit must be calculated iteratively, and because S-PLUS
%  is so slow with loopy calculations, these fits are VERY slow.  But
%  they are best quality fits that I and my colleagues, notably
%  R. D. Bock, have been able to achieve to date.
%  The Matlab version of this function is much faster.

%  ------  First set up a basis for monotone smooth   --------
%  We use b-spline basis functions of order 6
%  Knots are positioned at the ages of observation.

norder = 6;
nbasis = nage + norder - 2;
wbasis = create_bspline_basis(rng, nbasis, norder, knots);

%  starting values for coefficient

cvec0 = zeros(nbasis,1);
Wfd0  = fd(cvec0, wbasis);

%  starting values for coefficient

zmat  = ones(nage,1);
wgt   = zmat;

Lfd    = int2Lfd(3);  %  penalize curvature of velocity
lambda = 10^(-0.5);   %  smoothing parameter

% -----------------  Male data  --------------------

cvecm = zeros(nbasis, ncasem);
betam = zeros(2,      ncasem);
RMSEm = zeros(1,      ncasem);

index = 1:ncasem;

for icase=index
  hgt = hgtmmat(:,icase);
  [Wfd, beta, Fstr, iternum, iterhist] = ...
             smooth_monotone(age, hgt, wgt, Wfd0, zmat, Lfd, lambda);
  cvecm(:,icase) = getcoef(Wfd);
  betam(:,icase) = beta;
  hgthat = beta(1) + beta(2).*monfn(age, Wfd);
  RMSEm(icase) = sqrt(mean((hgt - hgthat).^2));
  fprintf('%5.f %5.f %12.5f %10.4f\n', [icase, iternum, Fstr.f, RMSEm(icase)])
end

% -----------------  Female data  --------------------

cvecf = zeros(nbasis, ncasef);
betaf = zeros(2,      ncasef);
RMSEf = zeros(1,      ncasef);

index = 1:ncasef;

for icase=index
  hgt = hgtfmat(:,icase);
  [Wfd, beta, Fstr, iternum, iterhist] = ...
             smooth_monotone(age, hgt, wgt, Wfd0, zmat, Lfd, lambda);
  cvecf(:,icase) = getcoef(Wfd);
  betaf(:,icase) = beta;
  hgthat = beta(1) + beta(2).*monfn(age, Wfd);
  RMSEf(icase) = sqrt(mean((hgt - hgthat).^2));
  fprintf('%5.f %5.f %12.5f %10.4f\n', [icase, iternum, Fstr.f, RMSEf(icase)])
end

%  plot data and smooth, residuals, velocity, and acceleration

%  Males:

index = 1:ncasem;
for i = index
  Wfd  = fd(cvecm(:,i),wbasis);
  beta = betam(:,i);
  hgtmhat   = beta(1) + beta(2).*monfn(age, Wfd);
  Dhgtmhat  = beta(2).*eval_mon(age, Wfd, 1);
  D2hgtmhat = beta(2).*eval_mon(age, Wfd, 2);
  subplot(2,2,1)
  plot(age, hgtmmat(:,i), 'go', age, hgtmhat, '-')
  axis([1, 18, 60, 200]);
  xlabel('Years');  title(['Height for male ',num2str(i)])
  resi = hgtmmat(:,i) - hgtmhat;
  subplot(2,2,2)
  plot(age, resi, '-o',     [1,18], [0,0], 'r--')
  axis([1,18,-1,1]);
  xlabel('Years');  title('Residuals')
  subplot(2,2,3)
  plot(age, Dhgtmhat, '-',  [1,18], [0,0], 'r--')
  axis([1,18,0,15]);
  xlabel('Years');  title('Velocity')
  subplot(2,2,4)
  plot(age, D2hgtmhat, '-')
  axis([1,18,-6,6]);
  xlabel('Years') ;  title('Acceleration')
  pause;
end

% Females:

index = 1:ncasef;
for i = index
  Wfd  = fd(cvecf(:,i),wbasis);
  beta = betaf(:,i);
  hgtfhat   = beta(1) + beta(2).*monfn(age, Wfd);
  Dhgtfhat  = beta(2).*eval_mon(age, Wfd, 1);
  D2hgtfhat = beta(2).*eval_mon(age, Wfd, 2);
  subplot(2,2,1)
  plot(age, hgtfmat(:,i), 'go', age, hgtfhat, '-')
  axis([1, 18, 60, 200]);
  xlabel('Years');  title(['Height for female ',num2str(i)])
  resi = hgtfmat(:,i) - hgtfhat;
  subplot(2,2,2)
  plot(age, resi, '-o',     [1,18], [0,0], 'r--')
  axis([1,18,-1,1]);
  xlabel('Years');  title('Residuals')
  subplot(2,2,3)
  plot(age, Dhgtfhat, '-',  [1,18], [0,0], 'r--')
  axis([1,18,0,15]);
  xlabel('Years');  title('Velocity')
  subplot(2,2,4)
  plot(age, D2hgtfhat, '-')
  axis([1,18,-6,6]);
  xlabel('Years') ;  title('Acceleration')
  pause;
end

%  ---------------------------------------------------------------------
%            Register the velocity curves for the girls
%  ---------------------------------------------------------------------

nbasisw = 15;
norder  = 5;
basisw  = create_bspline_basis([1,18], nbasisw, norder);

index = 1:ncasef;

agefine = linspace(1, 18, 101)';
agemat  = agefine * ones(1,length(index));

y0fd = deriv(mean(hgtffd), 1);
yfd  = deriv(hgtffd(index), 1);

y0vec = eval_fd(y0fd, agefine);
yvec  = eval_fd(yfd,  agefine);

coef0 = zeros(nbasisw,length(index));
Wfd0  = fd(coef0, basisw);

Lfd    = int2Lfd(2);
lambda = 1;

regstr = registerfd(y0fd, yfd, Wfd0, Lfd, lambda);

yregfd  = regstr.regfd;
yregmat = eval_fd(yregfd,agefine);
Wfd     = regstr.Wfd;
warpmat = monfn(agefine, Wfd);
warpmat = 1 + 17.*warpmat./(ones(101,1)*warpmat(101,:));

for i = 1:length(index)
   subplot(1,2,1)
   plot(agefine, yvec(:,i), '-', agefine, y0vec, '--', agefine, yregmat(:,i), '-');
   axis('square')
   title(['Case ',num2str(i)])
   subplot(1,2,2)
   plot(agefine, warpmat(:,i), '-', agefine, agefine, '--')
   axis('square')
   pause
end

%  ---------------------------------------------------------------------
%        Monotone smooth of short term height measurements
%  ---------------------------------------------------------------------

%  ---------------- input the data  ----------------------------------

clear;
fid  = fopen('onechild.dat','rt');
temp = fscanf(fid,'%f');
n    = 83;
data = reshape(temp, [n, 2]);
day  = data(:,1);
hgt  = data(:,2);
rng  = [day(1), day(n)];
wgt  = ones(n,1);
zmat = wgt;

%  set up the basis

nbasis   = 43;
norder   = 4;
hgtbasis = create_bspline_basis(rng, nbasis, norder);

%  set up the functional data object for W = log Dh

cvec0 = zeros(nbasis,1);
Wfd   = fd(cvec0, hgtbasis);

%  set parameters for the monotone smooth

Lfd    = int2Lfd(2);
lambda = 1;

%  carry out the monotone smooth

[Wfd, beta, Fstr, iternum, iterhist] = ...
    smooth_monotone(day, hgt, wgt, Wfd, zmat, Lfd, lambda);

%  plot the function W = log Dh

subplot(1,1,1)
plot(Wfd);

%  plot the data plus smooth

dayfine  = linspace(day(1),day(n),151)';
yhat     = beta(1) + beta(2).*eval_mon(day, Wfd);
yhatfine = beta(1) + beta(2).*eval_mon(dayfine, Wfd);
plot(day, hgt, 'o', dayfine, yhatfine, '-')
xlabel('\fontsize{12} Day')
ylabel('\fontsize{12} Centimeters')

%  plot growth velocity

Dhgt = beta(2).*eval_mon(dayfine, Wfd, int2Lfd(1));
plot(dayfine, Dhgt)
xlabel('\fontsize{12} Days')
ylabel('\fontsize{12} Centimeters/day')

