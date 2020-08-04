addpath ('c:\matlab\fdaM')
addpath ('c:\matlab\fdaM\examples\handwrit')

%  Last modified 15 January 2003


%  -----------------------------------------------------------------------
%                 Registered Handwriting Data
%  -----------------------------------------------------------------------

fid = fopen('fdareg.dat','rt');
fdarray = reshape(fscanf(fid,'%f'), [20,2,1401]);
fdarray = permute(fdarray,[3,1,2]);

fdatime   = linspace(0, 2.3, 1401);

fdarange  = [0, 2.3];
fdatype   = 'bspline';

%  After some experimentation, 205 order 6 splines without smoothing
%    seemed to do an adequate job of representing the original data.
%  Order 6 was used to get a reasonable estimate of the third deriv.

fdabasisobj = create_bspline_basis([0, 2.3], 205, 6);

%  set up the functional data structure

fdafd = data2fd(fdarray, fdatime, fdabasisobj);
fdafd_fdnames{1} = 'Seconds';
fdafd_fdnames{2} = 'Replications';
fdafd_fdnames{3} = 'mm';
fdafd = putnames(fdafd, fdafd_fdnames);

%  plot all curves

plot(fdafd)

%  plot individual curves, including both sampling points and fit

fdaeval = eval_fd(fdafd, fdatime);
fdamean = squeeze(eval_fd(mean(fdafd), fdatime));

subplot(1,1,1);
for i = 1:20
  plot(fdarray(:,i,1), fdarray(:,i,2), 'go', ...
       fdaeval(:,i,1), fdaeval(:,i,2), '-', ...
       fdamean(:,1), fdamean(:,2), 'r--');
  axis([-40, 40,-40, 40]);
  title(['Record ', num2str(i)]);
  pause;
end

%  plot the acceleration records

i = 1;
fdasmthfd = smooth_fd(fdafd, 1e-12, int2Lfd(4));
subplot(1,1,1);
accel = squeeze(eval_fd(fdasmthfd(i,:), fdatime, int2Lfd(2)));
plot(fdatime, accel);
axis([0, 2.3, -8000, 8000]);
xlabel('\fontsize{12} Seconds')
ylabel('\fontsize{12} mm/sec/sec')
legend('X', 'Y')

print -dpsc2 'c:/MyFiles/talks/fdacourse/figs/handwritaccel.ps'

i=1;
subplot(2,1,1);
for j = 1:2
  subplot(2,1,j);
  plot(fdafd(i,j), int2Lfd(2));
  axis([0, 2.3, -20000, 20000]);
  ylabel('mm/sec/sec')
  title(['Acceleration for Variable ',num2str(j)]);
end

%  plot the acceleration magnitudes

D2fdamat = eval_fd(fdafd, fdatime, int2Lfd(2));
D2mag = sqrt(D2fdamat(:,:,1).^2 + D2fdamat(:,:,2).^2);

subplot(1,1,1);
plot(fdatime, D2mag)
axis([0,2.3,0,15000]);
xlabel('Seconds')
ylabel('mm/sec/sec')
title('Acceleration Magnitude');
