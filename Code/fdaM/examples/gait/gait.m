addpath ('c:\matlab\fdaM')
addpath ('c:\matlab\fdaM\examples\gait')

%  Last modified 15 January 2003

%  -----------------------------------------------------------------------
%                        Gait data
%  -----------------------------------------------------------------------

%  -------------  input the data for the two measures  ---------------

fid  = fopen('hip.dat','rt');
hip  = reshape(fscanf(fid,'%f'), [20,39]);
fid  = fopen('knee.dat','rt');
knee = reshape(fscanf(fid,'%f'), [20,39]);

gaittime = (1:20)/21;

gaitarray = zeros(20, 39, 2);
gaitarray(:,:,1) = hip;
gaitarray(:,:,2) = knee;

%  ---------------  set up the fourier basis  ------------------------

gaitbasis = create_fourier_basis([0,1], 21);

%  -----------  create the fd object (no smoothing)  -----------------

gaitfd = data2fd(gaitarray, gaittime,  gaitbasis);
gaitfd_fdnames{1} = 'Normalized time';
gaitfd_fdnames{2} = 'Boys';
gaitfd_fdnames{3} = 'Angle (deg.)';
gaitfd = putnames(gaitfd, gaitfd_fdnames);

% -----------  set up the harmonic acceleration operator  ----------

Lbasisobj = create_constant_basis([0,1]);
Lcoef     = [0, (2*pi)^2, 0];
wfd       = fd(Lcoef, Lbasisobj);
wfdcell   = fd2cell(wfd);       % convert the FD object to a cell object
harmaccelLfd = Lfd(3, wfdcell);  %  define the operator object

%   smooth the data a bit penalizing harmonic acceleration

lambda = 1e-9;
gaitfd = smooth_fd(gaitfd, lambda, harmaccelLfd);

%  plot for FDA lecture

subplot(2,1,1)
plot(gaitfd(:,1))
xlabel('')
title('\fontsize{12} Knee Angle')
subplot(2,1,2)
plot(gaitfd(:,2))
title('\fontsize{12} Hip Angle')
xlabel('')

print -dpsc2 'c:/MyFiles/talks/fdacourse/figs/gaitangles.ps'

%  --------  plot curves and their first derivatives  ----------------

%  plot each pair of curves interactively

casenames = [];
varnames  = ['Knee angle';'Hip angle ';];
plotfit_fd(gaitarray, gaittime, gaitfd, casenames, varnames)

%  plot the residuals, sorting cases by residual sum of squares

residual = 1;
sortwrd  = 1;
plotfit_fd(gaitarray, gaittime, gaitfd, casenames, varnames, residual, sortwrd)

%  plot first derivative of all curves

plot(gaitfd, 1)

%  -----  plot curves as cycles  --------

subplot(1,1,1);
cycleplot(gaitfd, 0);

%  -----  plot the mean functions and their first two derivatives

gaitmeanfd = mean(gaitfd);

plot(gaitmeanfd, 0)

plot(gaitmeanfd, 1)

plot(gaitmeanfd, 2)

%  plot of gait cycle for FDA lecture

gaitvec     = squeeze(eval_fd(gaitfd(1,:), gaittime));
gaitmeanvec = squeeze(eval_fd(gaitmeanfd,  gaittime));
gaitlet = ['A', 'B', 'C', 'D', 'E'];
gaitind = [1,4,7,12,16];

subplot(1,1,1)
plot(gaitvec(:,1), gaitvec(:,2), '.-', ...
     gaitmeanvec(:,1), gaitmeanvec(:,2), '.--')
xlabel('\fontsize{12} Knee Angle')
xlabel('\fontsize{12} Hip Angle')
axis([0,50,0,80])
hold on
for i=1:5
    text(gaitvec(gaitind(i),1),     gaitvec(gaitind(i),2),     gaitlet(i))
    text(gaitmeanvec(gaitind(i),1), gaitmeanvec(gaitind(i),2), gaitlet(i))
end
hold off

print -dpsc2 'c:/MyFiles/talks/fdacourse/figs/gaitloop.ps'



% ---------------  do a PCA of gait data  -------------------------------

%  do the PCA with varimax rotation

nharm   = 4;
lambda  = 1e-9;
gaitpca = pca(gaitfd, nharm, lambda, harmaccelLfd);
gaitpca = varmx_pca(gaitpca);

%  plot harmonics using cycle plots

subplot(1,1,1)
plot_pca(gaitpca, 101, 1, 0, 0, 1);

%  plot eigenvalues

gaiteigvals = gaitpca.eigvals;
x = ones(16,2);
x(:,2) = reshape((5:20),[16,1]);
y = log10(gaiteigvals(5:20));
c = x\y;
subplot(1,1,1)
plot(1:20,log10(gaiteigvals(1:20)),'-o', ...
     1:20, c(1)+ c(2).*(1:20), ':')
xlabel('Eigenvalue Number')
ylabel('Log10 Eigenvalue')

gaitharmfd = gaitpca.harmfd;
gaitharmmat = eval_fd(gaitharmfd,gaittime);
gaitvarprop = gaitpca.varprop;
gaitmeanvec = squeeze(eval_fd(gaitmeanfd,gaittime));

con = 5.*ones(1,4);
for j=1:4
    subplot(2,2,j)
    yplus = gaitmeanvec + con(j).*squeeze(gaitharmmat(:,j,:));
    plot(gaitmeanvec(:,1),gaitmeanvec(:,2),'g.')
    hold on
    for i=1:20
        plot([gaitmeanvec(i,1),yplus(i,1)],...
             [gaitmeanvec(i,2),yplus(i,2)],'b-')
    end
    hold off
    xlabel('Hip Angle')
    ylabel('Knee Angle')
    title(['PC ',num2str(j),' (',num2str(round(gaitvarprop(j)*1000)/10),'%)'])
    axis([-20,60,0,80])
end

print -dpsc2 'c:/Myfiles/talks/fdacourse/figs/gaitpca.ps'

%  ------  do a canonical correlation analysis of knee-hip curves  ------

%  first penalize the second derivative to get the results in the book

ncan    = 3;
lambda  = 7e-4;
gaitcca = cca(gaitfd, ncan, lambda);

plot_cca(gaitcca,1)

gaitcca.corr

%  now penalize the harmonic acceleration

lambda  = 1e-6;
gaitcca = cca(gaitfd, ncan, lambda, harmaccelLfd);

plot_cca(gaitcca, 1)

gaitcca.corr

%  ----------  compute the variance and covariance functions  -------

gaitvarbifd = var(gaitfd);

gaitvararray = eval_bifd(gaitvarbifd, gaittime, gaittime);

subplot(2,3,1)
contour(gaitvararray(:,:,1,1))
title('Knee - Knee')

subplot(2,3,2)
contour(gaitvararray(:,:,1,2))
title('Knee - Hip')

subplot(2,3,3)
contour(gaitvararray(:,:,1,3))
title('Hip - Hip')

subplot(2,3,4)
surf(gaitvararray(:,:,1,1))
title('Knee - Knee')

subplot(2,3,5)
surf(gaitvararray(:,:,1,2))
title('Knee - Hip')

subplot(2,3,6)
surf(gaitvararray(:,:,1,3))
title('Hip - Hip')

print -dpsc2 'c:/MyFiles/talks/fdacourse/figs/gaitcorr.ps'

%  ----  register the first derivative of the gait data  

%  set up basis for warping function

nbasis = 7;
wbasis = create_fourier_basis([0,1],nbasis);

%  set parameters for registerfd

Lfd      = 3;
lambda   = 1e-3;
conv     = 1;
periodic = 1;
iterlim  = 20;
dbglev   = 1;
crit     = 2;

index = 1:39;

Dgaitfd = deriv(gaitfd(index,:),1);
y0fd = mean(Dgaitfd);
index = 1;  %   select a subset of curves to register here
yfd  = Dgaitfd(index);
xfine = linspace(0,1,101)';
ofine = ones(101,1);
y0vec = squeeze(eval_fd(y0fd, xfine));
yvec  = eval_fd(yfd, xfine);

cvec0  = zeros(nbasis,length(index));
Wfd0   = fd(cvec0, wbasis);

regstr = registerfd(y0fd, yfd, Wfd0, Lfd, lambda, periodic, ...
   iterlim, dbglev, conv, crit);

yregfd  = regstr.regfd;
yregmat = eval_fd(yregfd,xfine);
Wfd     = regstr.Wfd;
shift   = regstr.shift;
warpmat = monfn(xfine, Wfd);
warpmat = ofine*shift' + warpmat./(ofine*warpmat(101,:));

%  plot the registered gait functions

plot(yregfd)

%  plot the knee and hip angles for each case

for i = index
   subplot(1,2,1)
   plot(xfine, yvec(:,i,1), '-', xfine, y0vec(:,1), '--', xfine, yregmat(:,i,1), '-');
   axis('square')
   title(['Knee derivative ',num2str(index(i))])
   subplot(1,2,2)
   plot(xfine, yvec(:,i,2), '-', xfine, y0vec(:,2), '--', xfine, yregmat(:,i,2), '-');
   axis('square')
   title(['Hip derivative ',num2str(index(i))])
   pause
end

%  Plot the warping functions and display shifts
%  Note case 4, for which knee remains way out of
%    of phase with target, and case 31 with a large shift.

subplot(1,1,1)
for i = index
   plot(xfine, warpmat(:,i), '-', xfine, xfine+shift(i), '--')
   axis('square')
   title(['Case ',num2str(index(i)),' shift = ',num2str(shift(i))])
   pause
end

