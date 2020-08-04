addpath ('c:\matlab6p5\fdaM')
addpath ('c:\matlab6p5\fdaM\examples\lip')

%  Last modified 15 January 2003


%  -----------------------------------------------------------------------
%                       Lip Movement Data
%  -----------------------------------------------------------------------

%  ----------------  input the data  ------------------------

fid = fopen('lip.dat','rt');
lipmat = reshape(fscanf(fid,'%f'), [51, 20]);
nobs = size(lipmat,2);  %  number of replications

liptime  = (0:0.02:1)';  %  sampling points for each curve

%  ----------  set up the b-spline basis object  ------------
%       use order 6 splines so we can look at acceleration

lipbasis = create_bspline_basis([0,1], 31, 6);

%  -----  create the fd object  ---------

lipfd = data2fd(lipmat, liptime, lipbasis);
lipfd_fdnames{1} = 'Normalized time';
lipfd_fdnames{2} = 'Replications';
lipfd_fdnames{3} = 'mm';
lipfd = putnames(lipfd, lipfd_fdnames);

%  -----------  apply light smoothing --------------

lipfd = smooth_fd(lipfd, 1e-12, int2Lfd(4));
lipmeanfd = mean(lipfd);  % set mean lip position

%  ------------  plot data and fit  ----------------

plotfit_fd(lipmat, liptime, lipfd)

%  ------------- plot residuals  ------------

plotfit_fd(lipmat, liptime, lipfd, [], [], 1, 1)

%  ---------  summarize and plot the functions and their accelerations  -----

subplot(2,1,1)
plot(lipfd);
title('Lip position')

subplot(2,1,2)
plot(lipfd, 1, 1, 2);
ylabel('mm/t/t')
title('Lip acceleration')

%  ---- register the data using the minimum and the elbow as landmarks  ------
%               manually identify these points in each curve

%  there are two landmarks, in addition to the beginning and end

%  set up the matrix of acceleration values for identifying landmarks

D2lipmat = eval_fd(lipfd, liptime, int2Lfd(2));

%  plot each acceleration curve, and click on the two
%    maxima, near t = .4, the other near t = .75

nmarks   = 2;
lipmarks = zeros(nobs,nmarks);
index    = zeros(nmarks,1);
subplot(1,1,1)
for i = 1:nobs
  plot(liptime, D2lipmat(:,i), 'o', [0,1], [0,0], ':')
  title(['Curve ',num2str(i)])
  for j = 1:nmarks
    [x y] = ginput(1);  %  input two clicks on points here
    index(j) = round(x*51);
  end
  lipmarks(i,:) = liptime(index)';
end

save lipmarks lipmarks  % save lip marks in case you need them for future work
load lipmarks

lipmeanmarks = mean(lipmarks);

%  -----------   register the curves  ---------------------

%  create a basis object for the warping function
%  it has order 4 (piecewise cubic) and two interior knots
%  positioned at the mean landmark values since
%  NBASIS = NORDER + # interior knots

nbasis = 6;
norder = 4;
breaks = [0,lipmeanmarks,1];
warpbasis = create_bspline_basis([0,1], nbasis, norder, breaks);
%  plot the basis
plot(warpbasis)  %  plot of B-spline basis functions
line([lipmeanmarks(1),lipmeanmarks(1)],[0,1])  % first knot
line([lipmeanmarks(2),lipmeanmarks(2)],[0,1])  % second knot

%  call landmark registration function to set up struct LMRKSTR
lmrkstr = landmarkreg(lipfd, lipmarks, lipmeanmarks, warpbasis);

lipregfd  = lmrkstr.regfd;  %  registered curves
lipwarpfd = lmrkstr.warpfd; %  warping functions

%  plot unregistered and registered curves

subplot(1,2,1)
plot(lipfd)
title('Unregistered')

subplot(1,2,2)
plot(lipregfd)
title('Registered')

%  plot unregistered and registered accelerations

subplot(1,2,1)
plot(lipfd,int2Lfd(2),1,1)
title('Unregistered')

subplot(1,2,2)
plot(lipregfd,int2Lfd(2),1,1)
title('Registered')


%  plot warping functions

subplot(1,1,1)
plot(lipwarpfd)
title('Warping Functions')

%  plot deformation functions: warp(t) - t

lipwarpmat = eval_fd(lipwarpfd,liptime);
lipdefmat  = lipwarpmat - liptime*ones(1,nobs);
plot(liptime, lipdefmat, '-', [0,1], [0,0], ':')
title('Deformation Functions')

%  ---------- principal component analysis --------------------------

nharm = 4;
lippcastr = pca(lipfd, nharm);

%  plot unrotated harmonics

subplot(1,1,1)
plot_pca(lippcastr)

%  rotate harmonics

lippcastr = varmx_pca(lippcastr);

%  plot rotated harmonics

plot_pca(lippcastr)
  
%  plot log eigenvalues

lipeigvals = lippcastr.eigvals;
plot(1:19,log10(lipeigvals(1:19)),'-o')
xlabel('Eigenvalue Number')
ylabel('Log10 Eigenvalue')

%  -------------  principal differential analysis  -------------------

%  compute the weight functions by the basis expansion method

difeorder  = 2;  %  order of equation

awtcell = {};
ufdcell = {};

nwbasis = 21;
wbasis = create_bspline_basis([0,1], nwbasis);
wcoef0 = zeros(nwbasis,1);
bwtstruct.fd       = fd(wcoef0, wbasis); 
bwtstruct.estimate = ones(difeorder,1); 
bwtcell{1,1,1} = bwtstruct;
bwtcell{1,1,2} = bwtstruct;

xfdcell{1} = lipregfd;

%  carry out principal differential analysis

[bfdcell, resfdcell, afdcell] = ...
    pdacell(xfdcell, ufdcell, awtcell, bwtcell, difeorder);

%  plot the weight functions

for j=1:2
    subplot(2,1,j)
    plot(bfdcell{1,1,j});
    ylabel(['Weight function ',num2str(j-1)]);
end

%  set up a linear differential operator 

wcoef = [getcoef(bfdcell{1}), getcoef(bfdcell{2})];
wfd = fd(wcoef,wbasis);
lipLfd = Lfd(difeorder, fd2cell(wfd));

%  compute forcing functions

force = eval_fd(lipfd, liptime, lipLfd);
%  plot the forcing functions for each curve
plot(liptime, force);
axis([0,1,-1e3,1e3]);

%  plot the mean forcing function along with second deriv.

forcemean = mean(force')';
D2mean = eval_fd(mean(lipregfd),liptime,int2Lfd(2));

plot(liptime, forcemean, '-', liptime, D2mean, ':')
axis([0,1,-700,700])

%  solve equation

global wfd  %  this is necessary for function derivs

ystart = eye(2);
[tp1, yp1] = ode45('derivs', liptime, ystart(:,1));
[tp2, yp2] = ode45('derivs', liptime, ystart(:,2));

%  plot the two solutions

umat = [yp1(:,1),yp2(:,1)];
subplot(2,1,1)
plot(liptime, umat, liptime, zeros(51,1), ':'), title('Function');
Dumat = [yp1(:,2),yp2(:,2)];
subplot(2,1,2)
plot(liptime, Dumat, liptime, zeros(51,1), ':'), title('Derivative');

%  plot fit to each curve ...  hit any key after each plot

index = 1:nobs;
lipmat   = eval_fd(lipfd, liptime);
D2lipmat = eval_fd(lipfd, liptime, int2Lfd(2));
for i = index
   subplot(2,1,1)
   %  solid line is forcing function, dashed 2nd deriv.
   plot(liptime, force(:,i),    '-',  ...
        liptime, D2lipmat(:,i), '--', ...
        liptime, zeros(51,1),   ':')
   axis([0,1,-1000,1000])
   title(['Record ',num2str(i),' Forcing Fn.'])
   xhat = umat * (umat\lipmat(:,i));
   subplot(2,1,2)
   %  solid is fit, dashed is actual
   plot(liptime, xhat, '-', liptime, lipmat(:,i), '--')
   title('Function')
   pause;
end

