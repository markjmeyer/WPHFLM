function Hline = line(fd, Lfd, nx)
%  LINE adds lines to a plot of a functional data object.

%  last modified 20 October 2000

if nargin < 3, nx = 101;  end
if nargin < 2, Lfd = int2Lfd(0);   end

if ~isa_Lfd(Lfd)
    error ('LFD is not a linear differential operator object.');
end

coefd = size(getcoef(fd));
ndim  = length(coefd);
if ndim < 3
   rangex = getbasisrange(getbasis(fd));
   x      = linspace(rangex(1),rangex(2),nx);
   fdmat  = eval(fd, x, Lfd);
   Hline = line(x, fdmat);
else
   error('Too many dimensions');
end

