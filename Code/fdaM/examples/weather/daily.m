Windows:

addpath ('c:\matlab\fdaM')
addpath ('c:\matlab\fdaM\examples\weather')

Unix:

addpath('/export/home/steve/ramsay/u1/fdaM')
addpath('/export/home/steve/ramsay/u1/fdaM/examples/weather')

%  Last modified 15 January 2003

%  -----------------------------------------------------------------------
%                     Daily Weather Data
%  -----------------------------------------------------------------------

%  ------------------------  input the data  -----------------------

fid = fopen('dailtemp.dat','rt');
tempav = fscanf(fid,'%f');
tempav = reshape(tempav, [365,35]);

fid = fopen('dailprec.dat','rt');
precav = fscanf(fid,'%f');
precav = reshape(precav, [365,35]);

daytime = (1:365)-0.5;

place = [ ...
'arvida  '; 'bagottvi'; 'calgary '; 'charlott'; 'churchil'; 'dawson  '; ...
'edmonton'; 'frederic'; 'halifax '; 'inuvik  '; 'iqaluit '; 'kamloops'; ...
'london  '; 'montreal'; 'ottawa  '; 'princeal'; 'princege'; 'princeru'; ...
'quebec  '; 'regina  '; 'resolute'; 'scheffer'; 'sherbroo'; 'stjohns '; ...
'sydney  '; 'thepas  '; 'thunderb'; 'toronto '; 'uraniumc'; 'vancouvr'; ...
'victoria'; 'whitehor'; 'winnipeg'; 'yarmouth'; 'yellowkn'];

%  -------------  set up fourier basis  --------------------------

nbasis = 65;
daybasis = create_fourier_basis([0,365], nbasis);

%  ---------  create fd objects for temp. and prec. ---------------

daytempfd = data2fd(tempav, daytime, daybasis);
daytempfd_fdnames{1} = 'Day';
daytempfd_fdnames{2} = 'Station';
daytempfd_fdnames{3} = 'Deg C';
daytempfd = putnames(daytempfd, daytempfd_fdnames);

dayprecfd = data2fd(precav, daytime, daybasis);
dayprecfd_fdnames{1} = 'Day';
dayprecfd_fdnames{2} = 'Station';
dayprecfd_fdnames{3} = 'Deg C';
dayprecfd = putnames(dayprecfd, dayprecfd_fdnames);

%  Plot temperature curves and values

plotfit_fd(tempav, daytime, daytempfd, place)

%  Plot residuals for three best fits and three worst fits

casenames = place;
varnames  = 'Temperature';
rng       = [0,365];
index     = [1,2,3,33,34,35];
residual  = 1;
sortwrd   = 1;

plotfit_fd(tempav, daytime, daytempfd, casenames, varnames, ...
           residual, sortwrd, rng, index)

%  ---------------  interactive plots  ---------------------------------

global tempstr precstr nbasis place

tempstr.y = tempav;
tempstr.c = getcoef(daytempfd);
tempstr.mean = mean(daytempfd);

precstr.y = precav;
precstr.c = getcoef(dayprecfd);
precstr.mean = mean(dayprecfd);

dailymenu

%  these smoothing parameter values probably undersmooth the data,
%  but we can impose further smoothness on the results of analyses

%  set up the harmonic acceleration operator

Lbasis  = create_constant_basis([0,365]);  %  create a constant basis
Lcoef   = [0,(pi/6)^2,0];                 %  set up three coefficients
wfd     = fd(Lcoef,Lbasis);   % define an FD object for weight functions
wfdcell = fd2cell(wfd);       % convert the FD object to a cell object
harmaccelLfd = Lfd(3, wfdcell);  %  define the operator object

templambda = 1e1;
preclambda = 1e5;

daytempfd = smooth_fd(daytempfd, templambda, harmaccelLfd);
dayprecfd = smooth_fd(dayprecfd, preclambda, harmaccelLfd);

%  plot each pair of functions along with raw data

tempmat = eval_fd(daytempfd, daytime);
precmat = eval_fd(dayprecfd, daytime);

index = 1:35
for i = index
  subplot(2,1,1)
  plot(daytime,tempav(:,i),'go',daytime,tempmat(:,i),'-')
  axis([0 365 -35 25])
  xlabel('Day')
  ylabel('Temperature (deg. C)')
  title(place(i,:))
  subplot(2,1,2)
  plot(daytime,precav(:,i),'go',daytime,precmat(:,i),'-')
  axis([0 365 0 13])
  xlabel('Day')
  ylabel('Precipitation (mm)')
  pause
end

%  plot all the functions

subplot(1,1,1)
plot(daytempfd);
  axis([0 365 -35 25])

plot(dayprecfd);
  axis([0 365 0 13])

%  -------------------  do a PCA of temperature  -------------------

nharm  = 4;
Lfd    = harmaccelLfd;
lambda = 1e4;

daytemppcastr = pca(daytempfd, nharm, lambda, Lfd);
daytemppcastr = varmx_pca(daytemppcastr);

%  plot harmonics

subplot(1,1,1)
plot_pca(daytemppcastr)

%  plot log eigenvalues

daytempharmeigval = daytemppcastr.eigvals;
x = ones(16,2);
x(:,2) = reshape((5:20),[16,1]);
y = log10(daytempharmeigval(5:20));
c = x\y;
subplot(1,1,1)
plot(1:20,log10(daytempharmeigval(1:20)),'-o', ...
     1:20, c(1)+ c(2).*(1:20), ':')
xlabel('Eigenvalue Number')
ylabel('Log10 Eigenvalue')

%  plot factor scores

harmscr = daytemppcastr.harmscr;

plot(harmscr(:,1), harmscr(:,2), 'o')
xlabel('Harmonic I')
ylabel('Harmonic II')
text(harmscr(:,1), harmscr(:,2), place)

%  -----  model log10 total precipitation using temperature  -------------

annualprec = log10(sum(precav))';

placewt = ones(35,1);
xLfd    = int2Lfd(2);
xlambda = 1e6;
yLfd    = int2Lfd(0);
ylambda = 0;

linmodstr = linmod(daytempfd, annualprec, placewt, ...
                       xLfd, yLfd, xlambda, ylambda);

alpha = linmodstr.alpha;
regfd = linmodstr.reg;
yhat  = linmodstr.yhat;

plot(regfd)

plot(yhat, annualprec, 'o', yhat, yhat, 'r-')
text(yhat, annualprec, place)

SSE = sum((annualprec - yhat).^2);
SSY = sum((annualprec - mean(annualprec)).^2);
Rsq = (SSY - SSE)/SSY;
disp(Rsq);
